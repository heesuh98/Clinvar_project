#!/bin/python

import pandas as pd
import json
import time
import logging
import argparse
from pathlib import Path
from tqdm import tqdm
from openai import OpenAI
from datetime import datetime

# 로깅 설정
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# OpenAI 설정
OPENAI_API_KEY = "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"
OPENAI_MODEL = "gpt-5.2"


class PhenotypePatientExtractor:
    """변이 정보에 phenotype과 patient 정보를 추가하는 클래스"""

    def __init__(self, api_key, model):
        self.client = OpenAI(api_key=api_key)
        self.model = model
        self.stats = {
            'total_variants': 0,
            'phenotype_found': 0,
            'patient_found': 0,
            'both_found': 0,
            'none_found': 0,
            'api_errors': 0
        }



    def create_extraction_prompt(self, fulltext, variant_c, variant_p, gene):
        """Phenotype과 patient 정보 추출 프롬프트"""

        variant_info = f"cDNA: {variant_c}, Protein: {variant_p}" if variant_c != "Not specified" or variant_p != "Not specified" else "variant information mentioned in paper"

        prompt = f"""
    You are an expert clinical geneticist specializing in genotype-phenotype correlations.

    Your task: Extract **clinical symptoms (phenotypes)** and **patient information** for a specific variant in the **{gene}** gene.

    ---

    ## TARGET VARIANT
    Gene: {gene}
    Variant: {variant_info}

    ---

    ## TASK 1 — PHENOTYPE EXTRACTION

    Extract ONLY **clinical symptoms observed in patients**.

    ### Allowed phenotype examples
    - seizures
    - epilepsy
    - developmental delay
    - intellectual disability
    - hypotonia
    - movement disorder
    - microcephaly
    - ataxia
    - speech delay

    ### IMPORTANT
    Phenotypes must be **clinical symptoms only**.

    DO NOT include:
    - variant descriptions
    - functional study results
    - protein effects
    - pathogenicity classification
    - inheritance (de novo, familial)
    - population frequency
    - EEG / MRI findings
    - treatment information

    If a symptom is explicitly **absent**, record it like:
    - "No intellectual disability"
    - "No developmental delay"

    ---

    ## TASK 2 — PATIENT EXTRACTION

    Extract individual patient data ONLY if explicitly described.

    For each patient extract:
    - patient_id
    - sex
    - age
    - seizure_onset
    - seizure_status
    - phenotypes (clinical symptoms only)

    Merge duplicate patients mentioned in different sections.

    ---

    ## CRITICAL RULES

    1. Extract ONLY patients with variants in **{gene}**.
    2. Phenotypes must be **clinical symptoms only**.
    3. Do NOT include variant strings in phenotype text.
    4. Merge duplicate patients across sections.

    ---

    ## OUTPUT FORMAT

    Return ONLY a valid JSON object.

    {{
    "phenotype": {{
        "found": true/false,
        "data": "clinical symptoms only (comma separated) OR 'Not found'",
        "location": "where in paper",
        "reasoning": "short explanation"
    }},
    "patient": {{
        "found": true/false,
        "count": number,
        "data": [
        {{
            "patient_id": "string",
            "sex": "string",
            "age": "string",
            "seizure_onset": "string",
            "seizure_status": "string",
            "phenotypes": "clinical symptoms only"
        }}
        ],
        "location": "where in paper",
        "reasoning": "short explanation"
    }}
    }}

    ---

    ## PAPER TEXT

    {fulltext}

    ---

    Extract phenotype and patient information for the specified {gene} variant.
    """
        return prompt
    
    def extract_json_safely(self, text):
        """응답에서 JSON을 안전하게 추출"""
        try:
            clean_text = text.replace("```json", "").replace("```", "").strip()
            start_idx = clean_text.find('{')
            end_idx = clean_text.rfind('}')

            if start_idx != -1 and end_idx != -1:
                json_str = clean_text[start_idx:end_idx+1]
                return json.loads(json_str)
            else:
                return None

        except json.JSONDecodeError as e:
            logging.warning(f"JSON 파싱 실패: {e}")
            return None
        except Exception as e:
            logging.error(f"JSON 추출 중 오류: {e}")
            return None

    def extract_phenotype_patient(self, paper_id, fulltext, variant_c, variant_p, gene):
        """한 변이에 대해 phenotype과 patient 추출"""

        # Fulltext 체크
        if pd.isna(fulltext) or not fulltext or len(str(fulltext)) < 500:
            logging.debug(f"Paper {paper_id}: Full text 없음")
            return {
                'phenotype_data': 'Not available (no fulltext)',
                'phenotype_location': 'N/A',
                'phenotype_reasoning': 'No fulltext available',
                'patient_count': 0,
                'patient_data': '',
                'patient_location': 'N/A',
                'patient_reasoning': 'No fulltext available'
            }

        # 프롬프트 생성
        prompt = self.create_extraction_prompt(fulltext, variant_c, variant_p, gene)

        # API 호출 (최대 3회 재시도)
        for attempt in range(3):
            try:
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {
                            "role": "system",
                            "content": "You are an expert clinical geneticist specializing in genotype-phenotype correlations and patient data extraction."
                        },
                        {
                            "role": "user",
                            "content": prompt
                        }
                    ],
                    max_completion_tokens=3000
                )

                raw_response = response.choices[0].message.content
                logging.debug(f"Paper {paper_id}: Response preview: {raw_response[:300]}...")

                data = self.extract_json_safely(raw_response)

                if not data:
                    raise ValueError("Failed to parse JSON response")

                # Phenotype 정보 추출
                phenotype = data.get('phenotype', {})
                phenotype_found = phenotype.get('found', False)
                phenotype_data = phenotype.get('data', 'Not found')
                phenotype_location = phenotype.get('location', 'N/A')
                phenotype_reasoning = phenotype.get('reasoning', 'N/A')

                # Patient 정보 추출
                patient = data.get('patient', {})
                patient_found = patient.get('found', False)
                patient_count = patient.get('count', 0)
                patient_list = patient.get('data', [])
                patient_location = patient.get('location', 'N/A')
                patient_reasoning = patient.get('reasoning', 'N/A')

                # Patient 데이터를 JSON string으로 변환
                patient_data_str = json.dumps(patient_list, ensure_ascii=False) if patient_list else ''

                # 통계 업데이트
                if phenotype_found and patient_found:
                    self.stats['both_found'] += 1
                elif phenotype_found:
                    self.stats['phenotype_found'] += 1
                elif patient_found:
                    self.stats['patient_found'] += 1
                else:
                    self.stats['none_found'] += 1

                logging.info(f"Paper {paper_id}: Phenotype={phenotype_found}, Patient={patient_found} (count={patient_count})")

                return {
                    'phenotype_data': phenotype_data,
                    'phenotype_location': phenotype_location,
                    'phenotype_reasoning': phenotype_reasoning,
                    'patient_count': patient_count,
                    'patient_data': patient_data_str,
                    'patient_location': patient_location,
                    'patient_reasoning': patient_reasoning
                }

            except Exception as e:
                logging.warning(f"Paper {paper_id} API 호출 실패 (시도 {attempt+1}/3): {e}")
                if attempt < 2:
                    time.sleep(2 ** attempt)
                else:
                    self.stats['api_errors'] += 1
                    return {
                        'phenotype_data': f'Error: {str(e)}',
                        'phenotype_location': 'N/A',
                        'phenotype_reasoning': f'API error: {str(e)}',
                        'patient_count': 0,
                        'patient_data': '',
                        'patient_location': 'N/A',
                        'patient_reasoning': f'API error: {str(e)}'
                    }

        # 여기 도달하면 안됨
        return None


def main():
    parser = argparse.ArgumentParser(
        description='변이 CSV에 phenotype과 patient 정보를 추가합니다.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('variants_csv', type=str, help='변이 정보 CSV (paper_id, variant_c, variant_p 포함)')
    parser.add_argument('fulltext_csv', type=str, help='Fulltext CSV (paper_id, fulltext 포함)')
    parser.add_argument('--sample_num', type=int, default=5, help='처리할 변이 개수 (기본값: 5)')
    parser.add_argument('--all', action='store_true', help='전체 변이 처리')
    parser.add_argument('--start_row', type=int, default=0, help='시작 row (0-indexed)')

    args = parser.parse_args()

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # 결과 디렉토리 생성
    script_dir = Path(__file__).parent
    workspace_dir = script_dir.parent
    output_dir = workspace_dir / 'results' / f'variants_phenotype_patient_{timestamp}'
    output_dir.mkdir(parents=True, exist_ok=True)

    output_csv = output_dir / 'variants_with_phenotype_patient.csv'
    metadata_json = output_dir / 'metadata.json'

    logging.info("="*60)
    logging.info("Phenotype & Patient 정보 추가 시작")
    logging.info(f"변이 CSV: {args.variants_csv}")
    logging.info(f"Fulltext CSV: {args.fulltext_csv}")
    logging.info(f"출력 디렉토리: {output_dir}")
    logging.info(f"출력 파일: {output_csv}")
    logging.info("="*60)

    # CSV 로드
    variants_df = pd.read_csv(args.variants_csv)
    fulltext_df = pd.read_csv(args.fulltext_csv)

    logging.info(f"변이 데이터: {len(variants_df)}개")
    logging.info(f"Fulltext 데이터: {len(fulltext_df)}개")

    # paper_id 기준으로 fulltext 매칭
    # variants_df의 paper_id 컬럼과 fulltext_df의 pmcid 컬럼 매칭
    fulltext_dict = dict(zip(fulltext_df['pmcid'], fulltext_df['fulltext']))

    # 처리할 범위 결정
    if args.all:
        process_df = variants_df.iloc[args.start_row:].copy()
        logging.info(f"전체 모드: Row {args.start_row+1}부터 끝까지 ({len(process_df)}개)")
    else:
        end_row = min(args.start_row + args.sample_num, len(variants_df))
        process_df = variants_df.iloc[args.start_row:end_row].copy()
        logging.info(f"샘플 모드: Row {args.start_row+1}부터 {args.sample_num}개 ({len(process_df)}개)")

    # Extractor 초기화
    extractor = PhenotypePatientExtractor(OPENAI_API_KEY, OPENAI_MODEL)

    # 새 컬럼 초기화
    process_df['phenotype_data'] = ''
    process_df['phenotype_location'] = ''
    process_df['phenotype_reasoning'] = ''
    process_df['patient_count'] = 0
    process_df['patient_data'] = ''
    process_df['patient_location'] = ''
    process_df['patient_reasoning'] = ''

    # 추출 실행
    logging.info(f"\nPhenotype & Patient 추출 시작... ({len(process_df)}개 변이)")

    for idx, row in tqdm(process_df.iterrows(), total=len(process_df), desc="변이 처리"):
        extractor.stats['total_variants'] += 1

        paper_id = row.get('paper_id', '')
        variant_c = row.get('variant_c', 'Not specified')
        variant_p = row.get('variant_p', 'Not specified')
        gene = row.get('gene', 'EEF1A2')

        # Fulltext 가져오기
        fulltext = fulltext_dict.get(paper_id, '')

        # 추출
        result = extractor.extract_phenotype_patient(paper_id, fulltext, variant_c, variant_p, gene)

        if result:
            process_df.at[idx, 'phenotype_data'] = result['phenotype_data']
            process_df.at[idx, 'phenotype_location'] = result['phenotype_location']
            process_df.at[idx, 'phenotype_reasoning'] = result['phenotype_reasoning']
            process_df.at[idx, 'patient_count'] = result['patient_count']
            process_df.at[idx, 'patient_data'] = result['patient_data']
            process_df.at[idx, 'patient_location'] = result['patient_location']
            process_df.at[idx, 'patient_reasoning'] = result['patient_reasoning']

        # Rate limiting
        time.sleep(1)

    # 결과 저장
    process_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    logging.info(f"\n결과 저장 완료: {output_csv}")

    # 메타데이터 저장
    metadata = {
        'timestamp': timestamp,
        'input_variants_csv': str(args.variants_csv),
        'input_fulltext_csv': str(args.fulltext_csv),
        'output_dir': str(output_dir),
        'output_csv': str(output_csv),
        'model': OPENAI_MODEL,
        'stats': extractor.stats,
        'total_variants_processed': len(process_df),
        'total_patients_found': int(process_df['patient_count'].sum()) if 'patient_count' in process_df.columns else 0
    }

    with open(metadata_json, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    logging.info(f"메타데이터 저장 완료: {metadata_json}")

    # 통계 출력
    logging.info("\n" + "="*60)
    logging.info("추출 완료!")
    logging.info("="*60)
    logging.info(f"처리한 변이: {extractor.stats['total_variants']}개")
    logging.info(f"Phenotype & Patient 모두 발견: {extractor.stats['both_found']}개")
    logging.info(f"Phenotype만 발견: {extractor.stats['phenotype_found']}개")
    logging.info(f"Patient만 발견: {extractor.stats['patient_found']}개")
    logging.info(f"둘 다 없음: {extractor.stats['none_found']}개")
    logging.info(f"API 오류: {extractor.stats['api_errors']}개")
    logging.info(f"총 Patient 수: {metadata['total_patients_found']}명")
    logging.info(f"\n출력 파일:")
    logging.info(f"  - CSV: {output_csv}")
    logging.info(f"  - Metadata: {metadata_json}")
    logging.info("="*60)


if __name__ == "__main__":
    main()