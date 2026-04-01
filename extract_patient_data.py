#!/usr/bin/env python3
"""
Full text가 포함된 CSV에서 OpenAI API를 사용하여 환자 임상 데이터를 추출하는 스크립트

사용법:
    # 처음 5개 논문만 처리 (기본값)
    python extract_patient_data.py

    # 특정 개수의 논문만 처리
    python extract_patient_data.py --sample_num 10

    # 전체 논문 처리 (494개)
    python extract_patient_data.py --all

    # Row 62부터 시작 (EEF1A2가 실제로 언급되는 논문들)
    python extract_patient_data.py --sample_num 5 --start_row 62

출력 파일 (모두 extraction_partials_{timestamp}/ 폴더 내):
    - final_merged.csv : 최종 추출된 전체 데이터 (환자 없음 포함)
    - patients_only.csv : 실제 환자 데이터만 (N/A 제외)
    - metadata.json : 전체 통계 + 환자 전용 통계 통합
    - partial_0001.csv, partial_0002.csv, ... : 중간 저장 파일들 (50개마다)

중간 저장 기능:
    - 50개 논문마다 자동으로 중간 결과 저장
    - 오류 발생 시 partial_*.csv 파일에서 복구 가능
    - Ctrl+C로 중단해도 중간 결과 보존됨

주의사항:
    - OpenAI API 키가 필요합니다 (OPENAI_API_KEY 변수 설정)
    - 전체 처리 시 약 2시간 소요, 비용 약 $0.52 예상
    - Rate limit을 고려하여 0.5초 간격으로 API 호출합니다
"""

import pandas as pd
import json
import time
import logging
import argparse
import os
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
OPENAI_MODEL = "gpt-5-mini"  # 현재는 저렴한 버전(테스트)
TARGET_GENE = "EEF1A2"


class PatientDataExtractor:
    """OpenAI API를 사용해 논문에서 환자 데이터를 추출하는 클래스"""

    def __init__(self, api_key, model, gene_name):
        self.client = OpenAI(api_key=api_key)
        self.model = model
        self.gene_name = gene_name
        self.stats = {
            'total_papers': 0,
            'successful_extractions': 0,
            'no_patients_found': 0,
            'api_errors': 0,
            'total_patients_extracted': 0
        }

    def create_extraction_prompt(self, fulltext, target_gene):
        """환자 데이터 추출을 위한 프롬프트 생성"""

        prompt = f"""You are an expert clinical geneticist and data extractor specializing in epilepsy and neurodevelopmental disorders.

Your task is to analyze the provided research paper and extract **ALL individual patient clinical data** for patients with variants in the **{target_gene}** gene.

## TARGET GENE TO FIND:
- Gene: {target_gene}

**IMPORTANT**: Extract ALL patients who have variants in the {target_gene} gene mentioned in this paper.

---

## CRITICAL EXTRACTION RULES

### 1. Gene-Specific Extraction
- Extract ALL patients who have variants in the {target_gene} gene
- For EACH patient, record:
  - Which variant they have (both cDNA and protein notation if available)
  - Their clinical phenotypes
  - All demographic information
- If you find NO patients with {target_gene} variants, return empty list

### 2. Patient Identification
- Extract individual patients, NOT aggregate statistics
- ❌ DO NOT extract: "50% of patients had seizures"
- ✅ DO extract: "Patient 1 had seizures at 6 months"
- ❌ DO NOT extract: Data from cited studies (e.g., "Smith et al. reported...")
- ✅ DO extract: New patient data from THIS paper

### 3. Phenotype Extraction (Clean Data)
- Extract ONLY clinical symptoms and observations
- ❌ FORBIDDEN: Do NOT include variant notation in phenotypes field
- Example: If text says "Seizures (c.123A>G)", extract:
  - phenotypes: "Seizures"
  - variant_c: "c.123A>G"

### 4. Variant Notation Extraction

**Case A: Standard SNVs/Indels**
- If variant is already in HGVS format, use as-is:
  - "c.123A>G" → variant_c: "c.123A>G"
  - "p.Arg41His" → variant_p: "p.Arg41His"

- If table splits data (Position, Ref, Alt columns), CONSTRUCT HGVS:
  - Position: 123, Ref: A, Alt: G → variant_c: "c.123A>G"
  - Position: 41, Ref: Arg, Alt: His → variant_p: "p.Arg41His"

**Case B: Structural Variants (CNVs, Deletions)**
- DO NOT force into 'c.' format
- Keep description exactly as written:
  - "Exon 1-3 deletion"
  - "Copy number gain"
  - "Genomic deletion chr20:62,120,000-62,125,000"

### 5. Data Source Tracking (REASONING)
**This is CRITICAL - you MUST provide detailed reasoning for EVERY patient**

Include in reasoning field:
a) **Where you found the data**:
   - "Table 1, row 3"
   - "Main text Results section, page 4, paragraph 2"
   - "Supplementary Table S2"
   - "Case report section describing Patient 1"

b) **What evidence supports your extraction**:
   - Quote exact phrases: "Paper states: 'Patient 1 (3y/M) with c.123A>G mutation...'"
   - Describe table structure: "Patient ID in first column, genotype in column 3"

c) **Confidence level**:
   - High: All fields clearly stated with no ambiguity
   - Medium: Some fields require minor inference
   - Low: Significant ambiguity or missing information

d) **Any concerns or caveats**:
   - "Age not explicitly stated, inferred from 'toddler'"
   - "Two different variants mentioned - selected c.123A>G as primary"
   - "Gender not specified"

**Example reasoning:**
"Extracted from Table 2, row 5. Paper explicitly states 'Patient 5 (4y/F)' with variant 'c.1061G>A (p.Ser354Asn)' in the Genotype column. Phenotype 'developmental delay and seizures' described in Clinical Features column. Seizure onset '8 months' found in main text page 6 paragraph 3 referring to Patient 5. High confidence - all data clearly documented."

---

## OUTPUT FORMAT

Return ONLY valid JSON in this EXACT structure:

{{
  "patients": [
    {{
      "patient_id": "string (e.g., 'Patient 1', 'P1', 'Case 1', 'Subject III-2')",
      "gene": "string (gene symbol, e.g., 'SCN2A', 'EEF1A2', 'STXBP1')",
      "sex": "string (Male/Female/M/F/Unknown)",
      "age": "string (e.g., '3 years', '6 months', '2y5m', 'Unknown')",
      "seizure_onset": "string (age at first seizure, e.g., '6 months', 'birth', 'Unknown')",
      "seizure_status": "string (e.g., 'Controlled', 'Refractory', 'Seizure-free', 'Unknown')",
      "phenotypes": "string (comma-separated clinical features, NO variant notation)",
      "variant_c": "string (cDNA notation, e.g., 'c.123A>G' or 'Exon 3 deletion')",
      "variant_p": "string (protein notation, e.g., 'p.Arg41His' or 'Not specified')",
      "summary": "string (ONE-LINE summary: e.g., 'Patient found in Table 2' or 'No patients - review article' or 'No patients - mouse model study')",
      "reasoning": "string (detailed explanation of WHERE and WHY you extracted this data, with confidence level)"
    }}
  ]
}}

**If NO patients found**: Return a special explanation record:

{{"patients": [
  {{
    "patient_id": "N/A",
    "gene": "{target_gene}",
    "sex": "N/A",
    "age": "N/A",
    "seizure_onset": "N/A",
    "seizure_status": "N/A",
    "phenotypes": "N/A",
    "variant_c": "N/A",
    "variant_p": "N/A",
    "summary": "One-line reason why no patients found (e.g., 'Review article', 'Mouse model study', 'No {target_gene} mention', 'Only aggregate statistics')",
    "reasoning": "Explain WHY no patients were found. Examples: 'This is a review article with no original patient data', 'This is a mouse model study, no human patients', 'Paper mentions {target_gene} but only in references/discussion, no individual patient cases', 'Only aggregate statistics provided, no individual patient data'"
  }}
]}}

---

## PAPER TO ANALYZE:

{fulltext}

---

Now extract all patient data following the rules above. Remember:
1. ONE-LINE summary is mandatory for EVERY record (patient found or not)
2. Detailed reasoning for EACH patient is mandatory
3. If NO patients found, provide both a summary and detailed explanation
"""
        return prompt

    def extract_json_safely(self, text):
        """응답에서 JSON을 안전하게 추출"""
        try:
            # Remove markdown code blocks if present
            clean_text = text.replace("```json", "").replace("```", "").strip()

            # Find JSON object
            start_idx = clean_text.find('{')
            end_idx = clean_text.rfind('}')

            if start_idx != -1 and end_idx != -1:
                json_str = clean_text[start_idx:end_idx+1]
                return json.loads(json_str)
            else:
                return {"patients": []}

        except json.JSONDecodeError as e:
            logging.warning(f"JSON 파싱 실패: {e}")
            return {"patients": []}
        except Exception as e:
            logging.error(f"JSON 추출 중 오류: {e}")
            return {"patients": []}

    def extract_from_paper(self, paper_id, fulltext, fulltext_source, target_gene, ground_truth_hgvsc, ground_truth_hgvsp):
        """한 논문에서 환자 데이터 추출 (paper_id는 PMCID)"""

        # Full text가 없거나 실패한 경우
        if not fulltext or fulltext == 'Not available' or fulltext_source == 'Failed':
            logging.debug(f"Paper {paper_id}: Full text 없음 (source: {fulltext_source})")
            self.stats['no_patients_found'] += 1
            return []

        # 너무 짧은 텍스트 (에러 메시지일 가능성)
        if len(fulltext) < 500:
            logging.debug(f"Paper {paper_id}: 텍스트가 너무 짧음 ({len(fulltext)} chars)")
            self.stats['no_patients_found'] += 1
            return []

        # 프롬프트 생성 (gene만 전달, variant는 LLM이 찾도록)
        prompt = self.create_extraction_prompt(fulltext, target_gene)

        # API 호출 (최대 3회 재시도)
        for attempt in range(3):
            try:
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {
                            "role": "system",
                            "content": "You are an expert clinical geneticist and medical data extractor. You extract structured patient data from research papers with high precision."
                        },
                        {
                            "role": "user",
                            "content": prompt
                        }
                    ],
                    # temperature 파라미터 제거 (gpt-5-mini는 기본값만 지원)
                    max_completion_tokens=4000
                )

                # 응답 파싱
                raw_response = response.choices[0].message.content

                # [디버깅] Raw response 로깅
                logging.info(f"Paper {paper_id}: LLM Raw Response (처음 500자):\n{raw_response[:500]}...")

                data = self.extract_json_safely(raw_response)

                patients = data.get("patients", [])

                # PMCID와 metadata 태깅 (ground truth 포함)
                for patient in patients:
                    patient['pmcid'] = paper_id
                    patient['fulltext_source'] = fulltext_source
                    patient['target_gene'] = target_gene
                    patient['ground_truth_hgvsc'] = ground_truth_hgvsc
                    patient['ground_truth_hgvsp'] = ground_truth_hgvsp

                # 환자 데이터 확인
                if patients:
                    # N/A 레코드인지 실제 환자인지 확인
                    first_patient = patients[0]
                    is_na_record = (first_patient.get('patient_id') == 'N/A')

                    if is_na_record:
                        logging.info(f"Paper {paper_id}: 환자 데이터 없음")
                        logging.info(f"  LLM Reasoning: {first_patient.get('reasoning', 'N/A')[:200]}...")
                        self.stats['no_patients_found'] += 1
                    else:
                        logging.info(f"Paper {paper_id}: {len(patients)}명 환자 추출 성공")
                        self.stats['successful_extractions'] += 1
                        self.stats['total_patients_extracted'] += len(patients)
                else:
                    # LLM이 빈 리스트 반환 시 fallback (이상적으로는 발생하지 않아야 함)
                    logging.warning(f"Paper {paper_id}: LLM이 빈 리스트 반환 (프롬프트 오류 가능성)")
                    self.stats['no_patients_found'] += 1

                    patients.append({
                        'pmcid': paper_id,
                        'patient_id': 'N/A',
                        'gene': target_gene,
                        'sex': 'N/A',
                        'age': 'N/A',
                        'seizure_onset': 'N/A',
                        'seizure_status': 'N/A',
                        'phenotypes': 'N/A',
                        'variant_c': 'N/A',
                        'variant_p': 'N/A',
                        'summary': 'ERROR: LLM returned empty list',
                        'reasoning': f"ERROR: LLM returned empty list. Raw response: {raw_response[:500]}",
                        'fulltext_source': fulltext_source,
                        'target_gene': target_gene,
                        'ground_truth_hgvsc': ground_truth_hgvsc,
                        'ground_truth_hgvsp': ground_truth_hgvsp
                    })

                return patients

            except Exception as e:
                logging.warning(f"Paper {paper_id} API 호출 실패 (시도 {attempt+1}/3): {e}")
                if attempt < 2:
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    self.stats['api_errors'] += 1
                    return []

        return []


def main():
    """메인 실행 함수"""

    # 커맨드라인 인자 파싱
    parser = argparse.ArgumentParser(
        description='OpenAI API를 사용하여 논문에서 환자 임상 데이터를 추출합니다.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  python extract_patient_data.py                          # 처음 5개 논문만 처리 (기본값)
  python extract_patient_data.py --sample_num 10          # 10개 논문 처리
  python extract_patient_data.py --all                    # 전체 논문 처리
  python extract_patient_data.py --start_row 62 --sample_num 5  # Row 62부터 5개
        """
    )
    parser.add_argument('--sample_num', type=int, default=5,
                        help='처리할 논문 개수 (기본값: 5)')
    parser.add_argument('--all', action='store_true',
                        help='전체 논문 처리 (--sample_num 무시)')
    parser.add_argument('--start_row', type=int, default=0,
                        help='시작 row 번호 (0-indexed, 기본값: 0)')

    args = parser.parse_args()

    # Timestamp 생성
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # 입력/출력 파일 경로
    input_csv = '/Users/huios/workspace/litvar_df_with_fulltext_added_0227.csv'

    # 중간 결과 저장용 폴더 생성
    partial_dir = Path(f'/Users/huios/workspace/extraction_partials_{timestamp}')
    partial_dir.mkdir(exist_ok=True)

    logging.info("="*60)
    logging.info("환자 데이터 추출 시작")
    logging.info(f"입력 파일: {input_csv}")
    logging.info(f"출력 폴더: {partial_dir}")
    logging.info(f"타겟 유전자: {TARGET_GENE}")
    logging.info(f"OpenAI 모델: {OPENAI_MODEL}")
    logging.info("="*60)

    # CSV 로드
    try:
        df = pd.read_csv(input_csv)
        logging.info(f"총 {len(df)}개 논문 로드 완료")
    except FileNotFoundError:
        logging.error(f"파일을 찾을 수 없습니다: {input_csv}")
        logging.error("먼저 fetch_fulltext.py를 실행하여 full text를 다운로드하세요.")
        return

    # Extractor 초기화
    extractor = PatientDataExtractor(
        api_key=OPENAI_API_KEY,
        model=OPENAI_MODEL,
        gene_name=TARGET_GENE
    )

    # 모든 추출된 환자 데이터 저장
    all_patients = []

    # 처리할 논문 범위 결정
    if args.all:
        # 전체 논문 처리
        process_df = df.iloc[args.start_row:]
        logging.info(f"\n환자 데이터 추출 중... (전체 모드: Row {args.start_row+1}부터 {len(df)}까지)")
    else:
        # 샘플 개수만큼 처리
        end_row = min(args.start_row + args.sample_num, len(df))
        process_df = df.iloc[args.start_row:end_row]
        logging.info(f"\n환자 데이터 추출 중... (샘플 모드: Row {args.start_row+1}부터 {args.sample_num}개)")

    # 중간 저장 설정
    SAVE_INTERVAL = 50  # 50개마다 저장
    partial_counter = 0
    batch_counter = 0

    for idx, row in tqdm(process_df.iterrows(), total=len(process_df), desc="논문 처리"):
        pmcid = row.get('pmcid', 'Unknown')
        fulltext = row.get('fulltext', '')
        fulltext_source = row.get('fulltext_source', 'PMC')  # 기본값 PMC (컬럼 없을 경우)

        # Target gene과 ground truth variant 정보 추출
        target_gene = row.get('gene', 'Unknown')
        ground_truth_hgvsc = row.get('norm_hgvsc', 'Unknown')
        ground_truth_hgvsp = row.get('norm_hgvsp', 'Unknown')

        extractor.stats['total_papers'] += 1

        # 환자 데이터 추출 (gene만 전달, variant는 LLM이 찾도록)
        patients = extractor.extract_from_paper(pmcid, fulltext, fulltext_source,
                                                target_gene, ground_truth_hgvsc, ground_truth_hgvsp)
        all_patients.extend(patients)

        batch_counter += 1

        # 중간 저장 (50개마다)
        if batch_counter >= SAVE_INTERVAL:
            partial_counter += 1
            partial_file = partial_dir / f'partial_{partial_counter:04d}.csv'

            # 현재까지 누적된 데이터 저장
            if all_patients:
                temp_df = pd.DataFrame(all_patients)
                temp_df.to_csv(partial_file, index=False, encoding='utf-8-sig')
                logging.info(f"\n💾 중간 저장: {partial_file.name} ({len(all_patients)}개 레코드)")

            batch_counter = 0  # 카운터 리셋

        # API Rate limit (OpenAI는 보통 분당 요청 수 제한)
        time.sleep(0.5)

    # 최종 결과 저장 (환자 없어도 항상 저장)
    if all_patients:
        result_df = pd.DataFrame(all_patients)

        # 컬럼 순서 정리 (ground truth와 LLM 추출 결과 비교용)
        column_order = [
            'pmcid', 'target_gene', 'ground_truth_hgvsc', 'ground_truth_hgvsp',
            'patient_id', 'gene', 'sex', 'age',
            'seizure_onset', 'seizure_status', 'phenotypes',
            'variant_c', 'variant_p',
            'summary', 'reasoning', 'fulltext_source'
        ]

        # 존재하는 컬럼만 선택
        available_columns = [col for col in column_order if col in result_df.columns]
        result_df = result_df[available_columns]

        # 최종 CSV 저장 (partial_dir 내에만)
        final_merged = partial_dir / 'final_merged.csv'
        result_df.to_csv(final_merged, index=False, encoding='utf-8-sig')
        logging.info(f"\n✅ 최종 결과 저장 완료: {final_merged}")
        logging.info(f"   총 {len(result_df)}개 레코드 저장 (환자 없음 포함)")
    else:
        logging.warning("\n⚠️ 처리된 논문이 없습니다.")

    # Merge 및 통계 생성 (partial_dir 내에서)
    if all_patients:
        # 실제 환자만 필터링
        patients_only_df = result_df[
            (result_df['patient_id'].notna()) &
            (result_df['patient_id'] != '') &
            (result_df['patient_id'] != 'N/A')
        ].copy()

        # 통계 계산
        total_records = len(result_df)
        na_records = len(result_df) - len(patients_only_df)
        patient_records = len(patients_only_df)
        unique_pmcids = result_df['pmcid'].nunique()

        if patient_records > 0:
            unique_patient_ids = patients_only_df['patient_id'].nunique()
            gene_counts = patients_only_df['gene'].value_counts().to_dict()

            # 실제 환자만 저장 (partial_dir 내)
            patients_only_csv = partial_dir / 'patients_only.csv'
            patients_only_df.to_csv(patients_only_csv, index=False, encoding='utf-8-sig')
            logging.info(f"\n✅ 실제 환자 데이터만 저장: {patients_only_csv}")
            logging.info(f"   환자 레코드: {patient_records}개 (Unique 환자 ID: {unique_patient_ids}개)")
        else:
            unique_patient_ids = 0
            gene_counts = {}

        success_rate = (patient_records / total_records * 100) if total_records > 0 else 0.0

        # 메타데이터 저장 (partial_dir 내) - patients_only 정보도 통합
        metadata_json_final = partial_dir / 'metadata.json'

        # PMCID 리스트 (환자 발견된 논문만)
        pmcid_list = []
        if patient_records > 0:
            pmcid_list = sorted(patients_only_df['pmcid'].unique().tolist())

        metadata = {
            'extraction_config': {
                'input_csv': input_csv,
                'output_csv': str(partial_dir / 'final_merged.csv'),
                'patients_only_csv': str(partial_dir / 'patients_only.csv') if patient_records > 0 else None,
                'metadata_json': str(metadata_json_final),
                'target_gene': TARGET_GENE,
                'openai_model': OPENAI_MODEL,
                'partial_dir': str(partial_dir)
            },
            'all_records_statistics': {
                'total_records': total_records,
                'unique_pmcids': unique_pmcids,
                'patient_records': patient_records,
                'na_records': na_records,
                'success_rate_percent': round(success_rate, 1)
            },
            'patients_only_statistics': {
                'total_patient_records': patient_records,
                'unique_pmcids': len(pmcid_list),
                'unique_patient_ids': unique_patient_ids,
                'gene_distribution': gene_counts,
                'pmcid_list': pmcid_list
            },
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }

        with open(metadata_json_final, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

        logging.info(f"\n✅ 메타데이터 저장 완료: {metadata_json_final}")

        # 최종 통계 출력
        logging.info("\n" + "="*60)
        logging.info("최종 통계:")
        logging.info(f"  총 레코드: {total_records:,}개")
        logging.info(f"  Unique 논문: {unique_pmcids}개")
        logging.info(f"  환자 발견: {patient_records}개 ({success_rate:.1f}%)")
        logging.info(f"  환자 없음: {na_records}개")
        if patient_records > 0:
            logging.info(f"  Unique 환자 ID: {unique_patient_ids}개")
        logging.info("="*60)
        logging.info("완료!")


if __name__ == "__main__":
    main()
