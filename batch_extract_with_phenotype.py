#!/usr/bin/env python3
"""
환자별 변이 + Phenotype 정보 추출 (Responses API)

기능:
    - 환자별로 변이와 phenotype 정보를 함께 추출
    - Seizure types, age of onset, response to ASM 등
"""

import json
import argparse
import logging
import time
from pathlib import Path
from openai import OpenAI
import os
from datetime import datetime
import csv

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "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")
TARGET_GENE = "EEF1A2"

# Responses API가 지원하는 파일 확장자
SUPPORTED_EXTENSIONS = {
    '.pdf', '.docx', '.doc', '.pptx', '.ppt',
    '.xlsx', '.xls', '.csv', '.tsv',
    '.txt', '.md'
}


class PhenotypeExtractor:
    """환자별 변이 + Phenotype 추출 클래스"""

    def __init__(self, api_key, gene_name):
        self.client = OpenAI(api_key=api_key)
        self.gene_name = gene_name
        self.stats = {
            'total_pmcs': 0,
            'successful_pmcs': 0,
            'failed_pmcs': 0,
            'total_patients': 0,
            'api_errors': 0,
        }
        self.results = []

    def create_extraction_prompt(self, target_gene):
        """환자별 변이 + phenotype 추출 프롬프트"""

        return f"""You are an expert clinical geneticist. Analyze ALL attached files and extract patient-level information for the **{target_gene}** gene.

## TASK:
Extract information for EACH PATIENT with a {target_gene} variant. Each patient should have:
1. Genetic variant information (variant_c, variant_p)
2. Clinical phenotype information (seizures, age of onset, response to treatment, etc.)

## OUTPUT FORMAT (JSON only):

{{
  "patients": [
    {{
      "patient_id": "Patient 1 / Case 1 / ID from paper",
      "variant_c": "c.123A>G or 'Not specified'",
      "variant_p": "p.Arg41His or 'Not specified'",
      "phenotype": {{
        "seizure_presence": "Yes / No / Not specified",
        "age_of_onset": "6 months / 2 years / Not specified",
        "seizure_types": "Focal motor / Epileptic spasms / Generalized tonic-clonic / Not specified",
        "response_to_asm": "Controlled / Refractory / Unknown / Not specified",
        "asm_used": "VPA, LEV / Not specified",
        "responsive_asm": "VPA / Not specified",
        "epilepsy_type": "DEE / Focal epilepsy / Not specified",
        "other_phenotypes": "Developmental delay, intellectual disability, etc."
      }},
      "location_in_paper": "Table 1, row 3 / Case report section / Supplementary Table S1",
      "file_source": "Which file this information was found in",
      "confidence": "High / Medium / Low"
    }}
  ]
}}

## IMPORTANT NOTES:
- Extract ONLY patients from THIS paper's own data (not cited literature)
- Each patient should have BOTH variant AND phenotype information
- If phenotype information is not available for a specific field, use "Not specified"
- Look for information in:
  * Main text (case reports, patient descriptions)
  * Tables (patient summary tables, variant tables)
  * Supplementary materials (detailed patient data)
  * Figures (pedigrees, clinical timelines)

## Common phenotype fields to look for:
- Seizure types: focal, generalized, epileptic spasms, absence, myoclonic, etc.
- Age of onset: when seizures first appeared
- Response to ASM (Anti-Seizure Medications): controlled, refractory, partial response
- ASM used: VPA, LEV, LTG, OXC, TPM, etc.
- Epilepsy syndrome: DEE, West syndrome, Dravet syndrome, etc.
- Developmental outcomes: normal, delayed, intellectual disability
- Other features: autism, movement disorders, etc.

If NO patients found: return {{"patients": []}}

Return ONLY JSON, no additional text.
"""

    def get_pmc_files(self, pmc_dir):
        """PMC 디렉토리에서 지원되는 파일들 가져오기"""
        pmc_path = Path(pmc_dir)

        if not pmc_path.exists():
            logging.warning(f"디렉토리가 존재하지 않습니다: {pmc_dir}")
            return []

        files = []
        for file in pmc_path.iterdir():
            if file.is_file() and file.suffix.lower() in SUPPORTED_EXTENSIONS:
                files.append(file)

        return sorted(files)

    def extract_from_pmc(self, pmc_dir, pmcid):
        """단일 PMC에서 환자별 정보 추출"""
        self.stats['total_pmcs'] += 1

        logging.info(f"\n{'='*60}")
        logging.info(f"Processing {pmcid}")
        logging.info(f"Directory: {pmc_dir}")
        logging.info(f"{'='*60}")

        # 1. 파일 목록 가져오기
        files = self.get_pmc_files(pmc_dir)

        if not files:
            logging.warning(f"{pmcid}: 처리 가능한 파일이 없습니다")
            return {
                'pmcid': pmcid,
                'status': 'no_files',
                'patients': [],
                'error': 'No processable files found'
            }

        logging.info(f"📁 Found {len(files)} file(s):")
        for f in files:
            logging.info(f"   - {f.name} ({f.suffix})")

        # 2. 파일 업로드
        uploaded_files = []

        for file_path in files:
            try:
                logging.info(f"📤 Uploading {file_path.name}...")
                with open(file_path, "rb") as f:
                    file_obj = self.client.files.create(
                        file=f,
                        purpose="user_data"
                    )
                uploaded_files.append({
                    'file_id': file_obj.id,
                    'filename': file_path.name
                })
                logging.info(f"   ✅ {file_obj.id}")

            except Exception as e:
                logging.error(f"   ❌ 업로드 실패: {e}")
                continue

        if not uploaded_files:
            logging.error(f"{pmcid}: 파일 업로드 실패")
            return {
                'pmcid': pmcid,
                'status': 'upload_failed',
                'patients': [],
                'error': 'All file uploads failed'
            }

        # 3. Responses API 호출
        try:
            prompt_text = self.create_extraction_prompt(self.gene_name)

            # content 배열 구성
            content = []
            for file_info in uploaded_files:
                content.append({
                    "type": "input_file",
                    "file_id": file_info['file_id']
                })

            content.append({
                "type": "input_text",
                "text": prompt_text
            })

            logging.info(f"🤖 Analyzing {len(uploaded_files)} file(s) with Responses API...")

            response = self.client.responses.create(
                model="gpt-5.4",
                input=[{
                    "role": "user",
                    "content": content
                }]
            )

            response_text = response.output[0].content[0].text

            # JSON 파싱
            result = json.loads(response_text)

            patients = result.get('patients', [])
            self.stats['total_patients'] += len(patients)

            logging.info(f"\n📄 Response preview:")
            logging.info(json.dumps(result, indent=2, ensure_ascii=False)[:500] + "...")

            logging.info(f"✅ Found {len(patients)} patient(s)")

            self.stats['successful_pmcs'] += 1

            # 결과 구성
            extraction_result = {
                'pmcid': pmcid,
                'status': 'success',
                'files_processed': [f['filename'] for f in uploaded_files],
                'patients': patients,
                'timestamp': datetime.now().isoformat()
            }

        except json.JSONDecodeError as e:
            logging.error(f"❌ JSON 파싱 실패: {e}")
            logging.error(f"Response: {response_text[:500]}")
            self.stats['failed_pmcs'] += 1
            self.stats['api_errors'] += 1

            extraction_result = {
                'pmcid': pmcid,
                'status': 'json_error',
                'patients': [],
                'error': str(e),
                'raw_response': response_text[:1000]
            }

        except Exception as e:
            logging.error(f"❌ API 호출 실패: {e}")
            self.stats['failed_pmcs'] += 1
            self.stats['api_errors'] += 1

            extraction_result = {
                'pmcid': pmcid,
                'status': 'api_error',
                'patients': [],
                'error': str(e)
            }

        finally:
            # 4. 업로드한 파일 삭제
            for file_info in uploaded_files:
                try:
                    self.client.files.delete(file_info['file_id'])
                except Exception as e:
                    logging.warning(f"파일 삭제 실패: {e}")

        return extraction_result

    def save_results(self, output_dir):
        """결과 저장 (JSON + CSV)"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # 1. 전체 결과 JSON
        all_results_path = output_path / 'all_patients.json'
        with open(all_results_path, 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)

        # 2. 환자 테이블 CSV
        csv_path = output_path / 'patients_phenotype.csv'
        csv_records = []

        for result in self.results:
            if result['status'] != 'success':
                continue

            pmcid = result['pmcid']
            for patient in result.get('patients', []):
                phenotype = patient.get('phenotype', {})

                csv_records.append({
                    'pmcid': pmcid,
                    'patient_id': patient.get('patient_id', ''),
                    'variant_c': patient.get('variant_c', ''),
                    'variant_p': patient.get('variant_p', ''),
                    'seizure_presence': phenotype.get('seizure_presence', ''),
                    'age_of_onset': phenotype.get('age_of_onset', ''),
                    'seizure_types': phenotype.get('seizure_types', ''),
                    'response_to_asm': phenotype.get('response_to_asm', ''),
                    'asm_used': phenotype.get('asm_used', ''),
                    'responsive_asm': phenotype.get('responsive_asm', ''),
                    'epilepsy_type': phenotype.get('epilepsy_type', ''),
                    'other_phenotypes': phenotype.get('other_phenotypes', ''),
                    'location_in_paper': patient.get('location_in_paper', ''),
                    'file_source': patient.get('file_source', ''),
                    'confidence': patient.get('confidence', '')
                })

        if csv_records:
            with open(csv_path, 'w', encoding='utf-8-sig', newline='') as f:
                fieldnames = [
                    'pmcid', 'patient_id', 'variant_c', 'variant_p',
                    'seizure_presence', 'age_of_onset', 'seizure_types',
                    'response_to_asm', 'asm_used', 'responsive_asm',
                    'epilepsy_type', 'other_phenotypes',
                    'location_in_paper', 'file_source', 'confidence'
                ]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_records)

        # 3. 개별 PMC JSON 파일
        for result in self.results:
            pmcid = result['pmcid']
            pmc_file = output_path / f"{pmcid}.json"
            with open(pmc_file, 'w', encoding='utf-8') as f:
                json.dump(result, f, indent=2, ensure_ascii=False)

        logging.info(f"\n결과 저장: {output_dir}")
        logging.info(f"  - all_patients.json (전체 데이터)")
        logging.info(f"  - patients_phenotype.csv (환자 테이블)")
        logging.info(f"  - {{pmcid}}.json (개별 결과)")

    def print_stats(self):
        """통계 출력"""
        logging.info(f"\n{'='*60}")
        logging.info("FINAL STATISTICS")
        logging.info(f"{'='*60}")
        logging.info(f"Total PMCs processed: {self.stats['total_pmcs']}")
        logging.info(f"Successful: {self.stats['successful_pmcs']}")
        logging.info(f"Failed: {self.stats['failed_pmcs']}")
        logging.info(f"Total patients found: {self.stats['total_patients']}")
        logging.info(f"API errors: {self.stats['api_errors']}")
        logging.info(f"")


def main():
    parser = argparse.ArgumentParser(description='환자별 변이 + Phenotype 추출')
    parser.add_argument('--pmcid', help='단일 PMCID 처리')
    parser.add_argument('--gene', default=TARGET_GENE, help=f'유전자명 (기본: {TARGET_GENE})')
    parser.add_argument('--input_dir', default='pdfs/pmc_pubmed',
                        help='PMC 폴더가 있는 디렉토리')

    args = parser.parse_args()

    # Extractor 초기화
    extractor = PhenotypeExtractor(
        api_key=OPENAI_API_KEY,
        gene_name=args.gene
    )

    logging.info("="*60)
    logging.info("Patient-level Phenotype Extraction Starting")
    logging.info("="*60)

    # 단일 PMCID 처리
    if args.pmcid:
        pmc_dir = Path(args.input_dir) / args.pmcid
        result = extractor.extract_from_pmc(pmc_dir, args.pmcid)
        extractor.results.append(result)

        # 결과 저장
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_dir = f"results/{timestamp}_phenotype_{args.pmcid}"

    extractor.save_results(output_dir)
    extractor.print_stats()

    logging.info(f"\n결과 저장: {output_dir}")
    logging.info(f"{'='*60}")


if __name__ == '__main__':
    main()
