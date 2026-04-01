import os
import argparse
import logging
import pandas as pd
import glob
import re
import fitz  # PyMuPDF
import tiktoken
import datetime
from tqdm import tqdm
from typing import Optional, List, Dict

# 로컬 설정 임포트 (config.py가 없을 경우 기본값 사용)
try:
    import config
except ImportError:
    class Config:
        OPENAI_API_KEY = "sk-..." 
        ENTREZ_EMAIL = "email@example.com"
        DOWNLOAD_ROOT = "downloads"
        GENE_NAME = "SCN2A"
        OPENAI_MODEL_NAME = "gpt-5.1" 
    config = Config()

from data_processor import AdvancedDataProcessor
from analyzer import PatientDataAnalyzer

# === [Configuration] ===
# GPT-5.1 (400k Context) 기준 입력 제한 설정
# 250k 토큰을 Input으로 사용하고, 나머지 150k는 Output 및 시스템 여유분으로 확보
MAX_INPUT_TOKENS = 250000 

# === [Helper Functions] 로그 및 유틸리티 ===

def setup_logging(gene_name: str, mode: str) -> str:
    """
    로그 디렉토리를 생성하고 파일 핸들러를 설정합니다.
    파일명 형식: logs/GeneName_mode_timestamp.log
    """
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = f"{gene_name}_{mode}_{timestamp}.log"
    log_path = os.path.join(log_dir, log_filename)

    # 기존 핸들러 제거 (중복 로깅 방지)
    root_logger = logging.getLogger()
    if root_logger.handlers:
        root_logger.handlers = []

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path, encoding='utf-8'),
            logging.StreamHandler()
        ]
    )
    return log_path

def get_token_count(text: str) -> int:
    """
    텍스트의 토큰 수를 계산합니다 (GPT-4o/5.1 기준).
    오류 발생 시 글자 수 기반으로 근사치를 반환합니다.
    """
    try:
        encoding = tiktoken.encoding_for_model("gpt-5.1")
        return len(encoding.encode(text))
    except:
        # Fallback: 영어 약 4자=1토큰, 한글 포함 시 보수적 계산
        return max(1, len(text) // 3)


class SmartBatcher:
    """
    1. 대용량 처리 속도 최적화 (토큰 계산 최소화)
    2. 배치에 포함된 파일명 로깅 기능 추가
    """
    def __init__(self, max_tokens=MAX_INPUT_TOKENS): 
        self.max_tokens = max_tokens
        self.current_batch = []
        self.current_tokens = 0
        self.batches = []
        # [NEW] 현재 배치에 들어간 파일명들을 추적하기 위한 리스트
        self.current_batch_files = []

    def add_content(self, filename: str, content: str, header_info: str = ""):
        remaining_content = content
        remaining_header = header_info
        is_continuation = False

        # 헤더의 대략적 토큰 수 미리 계산
        base_header_tokens = get_token_count(f"\n\n--- FILE: {filename} ---\n{header_info}\n")

        while remaining_content:
            space_left = self.max_tokens - self.current_tokens
            
            # 공간이 너무 적으면 배치를 마감하고 새로 시작
            if space_left < 2000:
                self._finalize_batch()
                space_left = self.max_tokens
                space_left -= (base_header_tokens + 100) 

            # 남은 공간에 들어갈 수 있는 최대 글자 수를 추산
            estimated_chars = int(space_left * 3.5)
            
            # [Case 1] 남은 내용이 추산치보다 작으면 -> 통째로 추가
            if len(remaining_content) <= estimated_chars:
                candidate_text = remaining_content
                candidate_tokens = get_token_count(candidate_text)
                
                if self.current_tokens + candidate_tokens + base_header_tokens > self.max_tokens:
                    pass # 예측 실패 시 아래 슬라이싱 로직으로 이동
                else:
                    final_text = f"\n\n--- FILE: {filename}{' (Continuation)' if is_continuation else ''} ---\n{remaining_header}\n{candidate_text}"
                    final_tokens = get_token_count(final_text)
                    
                    self.current_batch.append({"type": "text", "content": final_text})
                    self.current_tokens += final_tokens
                    # [NEW] 파일명 기록 (중복 방지 없이 순서대로 기록하거나, set으로 관리 가능. 여기선 순서대로 기록)
                    self.current_batch_files.append(filename if not is_continuation else f"{filename}(Part)")
                    break

            # [Case 2] 잘라야 하는 경우 (Slice)
            slice_idx = estimated_chars
            chunk = remaining_content[:slice_idx]
            
            # 정밀 조정
            while True:
                temp_text = f"\n\n--- FILE: ... ---\n{remaining_header}\n{chunk}"
                chunk_tokens = get_token_count(temp_text)
                
                if self.current_tokens + chunk_tokens <= self.max_tokens:
                    break
                
                slice_idx = int(slice_idx * 0.9)
                if slice_idx == 0: break
                chunk = remaining_content[:slice_idx]

            batch_text = f"\n\n--- FILE: {filename}{' (Continuation)' if is_continuation else ''} ---\n{remaining_header}\n{chunk}\n...(Truncated)..."
            
            self.current_batch.append({"type": "text", "content": batch_text})
            self.current_tokens += get_token_count(batch_text)
            # [NEW] 파일명 기록
            self.current_batch_files.append(filename if not is_continuation else f"{filename}(Part)")

            remaining_content = remaining_content[slice_idx:]
            is_continuation = True
            
            self._finalize_batch()

            if "Continuation" not in remaining_header:
                remaining_header = f"[Continuation of {filename}]\n{header_info}"

    def _finalize_batch(self):
        if self.current_batch:
            # [NEW] 파일명 리스트를 문자열로 변환
            files_list_str = ", ".join(self.current_batch_files)
            
            logging.info(f"   [Batch Info] Finalized (Tokens: {self.current_tokens}/{self.max_tokens}) | Files: [{files_list_str}]")
            
            self.batches.append(self.current_batch)
            self.current_batch = []
            self.current_tokens = 0
            # [NEW] 파일명 리스트 초기화
            self.current_batch_files = []
            
    def get_batches(self):
        self._finalize_batch()
        return self.batches

# === [File Processing] 파일 읽기 ===

def read_file_content(file_path: str, target_gene: str = "") -> tuple:
    """
    파일 확장자에 따라 적절한 방식으로 내용을 읽어 텍스트로 반환합니다.
    
    Returns:
        tuple: (header_info, content_text)
        - header_info: 파일의 메타데이터(컬럼명 등)로, 파일이 분할될 때 문맥 유지를 위해 사용됩니다.
    """
    filename = os.path.basename(file_path).lower()
    header_info = ""
    content_text = ""

    try:
        # 1. PDF 처리 (텍스트 추출)
        if filename.endswith('.pdf'):
            doc = fitz.open(file_path)
            full_text = []
            for page in doc:
                full_text.append(page.get_text())
            content_text = "\n".join(full_text)
            doc.close()
            header_info = "[Type: PDF Document]"

        # 2. Excel/CSV 처리
        elif filename.endswith(('.xlsx', '.xls', '.csv')):
            if filename.endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                dfs = pd.read_excel(file_path, sheet_name=None)
                df = pd.concat(dfs.values(), ignore_index=True)
            
            # 헤더 정보 추출
            cols = ", ".join([str(c) for c in df.columns])
            header_info = f"[Columns: {cols}]"
            
            # 내용 텍스트화
            content_text = df.to_string(index=False)
            
            # (옵션) 초대용량 파일 필터링
            if len(content_text) > 1000000 and target_gene:
                 logging.info(f"   [Filter] Large file detected ({len(content_text)} chars). Applying keyword filter for '{target_gene}'.")
                 filtered_rows = []
                 for row in df.astype(str).values:
                     if target_gene.lower() in " ".join(row).lower():
                         filtered_rows.append(row)
                 if filtered_rows:
                     new_df = pd.DataFrame(filtered_rows, columns=df.columns)
                     content_text = new_df.to_string(index=False)
                     header_info += " [Filtered by Target Gene]"

        # 3. 일반 텍스트 파일
        elif filename.endswith('.txt'):
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content_text = f.read()

    except Exception as e:
        logging.warning(f"[Read Error] Failed to read {filename}: {e}")
        return "", ""

    return header_info, content_text

# === [Main Logic] ===

def deduplicate_patients_advanced(patients_list: list) -> list:
    """
    추출된 환자 데이터 리스트에서 중복을 제거합니다.
    Variant, Patient ID, Sex 정보를 종합적으로 비교하여 동일 환자로 추정되면 정보를 병합합니다.
    """
    merged_list = []
    
    def normalize(s): 
        return re.sub(r'[^a-z0-9]', '', str(s).lower()) if s and str(s).lower() not in ['nan','not specified'] else ""
    
    def extract_nums(s): 
        return "".join(re.findall(r'\d+', str(s)))

    for new_p in patients_list:
        is_merged = False
        for exist_p in merged_list:
            nc_new, nc_exist = normalize(new_p.get('variant_c')), normalize(exist_p.get('variant_c'))
            id_new, id_exist = extract_nums(new_p.get('patient_id')), extract_nums(exist_p.get('patient_id'))
            
            variant_match = (nc_new and nc_exist and nc_new == nc_exist)
            id_match = (id_new and id_exist and id_new == id_exist)
            
            should_merge = False
            
            # 병합 조건 판단
            if variant_match:
                if id_match or (not id_new or not id_exist): 
                    should_merge = True
            elif id_match:
                sex_new, sex_exist = normalize(new_p.get('sex')), normalize(exist_p.get('sex'))
                if not sex_new or not sex_exist or sex_new == sex_exist: 
                    should_merge = True

            if should_merge:
                # 정보 병합 (더 상세한 정보 우선)
                for key in new_p:
                    val_new = str(new_p.get(key, "")).strip()
                    val_exist = str(exist_p.get(key, "")).strip()
                    
                    if val_new and val_new.lower() not in ["not specified", "n/a", ""]:
                        if not val_exist or val_exist.lower() in ["not specified", "n/a"] or len(val_new) > len(val_exist):
                            exist_p[key] = val_new
                is_merged = True
                break
        
        if not is_merged: merged_list.append(new_p)
    return merged_list

def main():
    parser = argparse.ArgumentParser(description="Phase 2: Download & Analyze")
    parser.add_argument("gene_name", type=str, help="Target gene name (e.g., SCN2A)")
    parser.add_argument("--mode", type=str, choices=['download', 'analyze'], required=True, help="Execution mode")
    args = parser.parse_args()
    
    # 로깅 초기화
    log_file = setup_logging(args.gene_name, args.mode)
    logging.info(f"=== Start Process: {args.gene_name} ({args.mode}) ===")
    logging.info(f"Log File: {log_file}")
    logging.info(f"Max Input Tokens: {MAX_INPUT_TOKENS}")

    # 설정 객체 생성
    conf_dict = {k: v for k, v in config.__dict__.items() if not k.startswith('__')}
    conf_dict['GENE_NAME'] = args.gene_name
    
    processor = AdvancedDataProcessor(conf_dict)
    analyzer = PatientDataAnalyzer(conf_dict)

    # 타겟 PMID 로드 (Phase 1 결과물 CSV 활용)
    target_pmids = []
    csv_path = f"./results/{args.gene_name}_final_analysis.csv"
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        df.columns = [c.lower().strip() for c in df.columns]
        if 'pmid' in df.columns:
            target_pmids = [str(p).split('.')[0] for p in df['pmid'].dropna().unique()]
    
    logging.info(f"Target PMIDs Count: {len(target_pmids)}")

    # --- Mode 1: Download ---
    if args.mode == 'download':
        logging.info("=== Mode: Download ===")
        for pmid in tqdm(target_pmids):
            try:
                logging.info(f"[Processing] PMID: {pmid}")
                processor.download_raw_data(pmid)
            except Exception as e:
                logging.error(f"[Error] Failed to download PMID {pmid}: {e}")

    # --- Mode 2: Analyze ---
    elif args.mode == 'analyze':
        logging.info("=== Mode: Analyze ===")
        all_results = []
        raw_root = processor.raw_data_root

        available_pmids = [p for p in target_pmids if os.path.exists(os.path.join(raw_root, str(p)))]
        
        for i, pmid in enumerate(tqdm(available_pmids, desc="Analyzing Papers")):
            logging.info(f"\n[PMID: {pmid}] Processing ({i+1}/{len(available_pmids)})")
            
            base_path = os.path.join(raw_root, str(pmid))
            all_files = glob.glob(os.path.join(base_path, "*"))
            
            # 유효 파일 선별
            valid_files = [f for f in all_files if f.lower().endswith(('.pdf', '.xlsx', '.xls', '.csv', '.docx', '.txt'))]
            if not valid_files:
                logging.warning(f"[Warning] No valid files found for PMID {pmid}")
                continue

            # 1. 파일 내용 로드 및 스마트 배칭
            batcher = SmartBatcher(max_tokens=MAX_INPUT_TOKENS)
            
            # Main PDF 우선 처리
            pdfs = [f for f in valid_files if f.lower().endswith('.pdf')]
            main_pdf = next((f for f in pdfs if '_main.pdf' in f.lower()), max(pdfs, key=os.path.getsize) if pdfs else None)
            
            sorted_files = []
            if main_pdf: sorted_files.append(main_pdf)
            sorted_files.extend([f for f in valid_files if f != main_pdf])

            for f_path in sorted_files:
                header, content = read_file_content(f_path, args.gene_name)
                if content:
                    batcher.add_content(os.path.basename(f_path), content, header)
            
            # 2. 배치별 LLM 분석
            batches = batcher.get_batches()
            pmid_extracted = []
            
            logging.info(f"   [Batch] Generated {len(batches)} batches.")
            
            for b_idx, batch_content in enumerate(batches):
                logging.info(f"   [LLM Request] Batch {b_idx+1}/{len(batches)}")
                res = analyzer.analyze_with_text_content(batch_content)
                
                if res.get("patients"):
                    logging.info(f"     [Result] Extracted {len(res['patients'])} patients")
                    pmid_extracted.extend(res["patients"])
                else:
                    logging.info("     [Result] No data found.")

            # 3. 중복 제거 및 결과 병합
            if pmid_extracted:
                unique = deduplicate_patients_advanced(pmid_extracted)
                for p in unique: 
                    p['pmid'] = pmid 
                all_results.extend(unique)
                logging.info(f"   [Save] Merged {len(unique)} unique patients for PMID {pmid}")
            
            # 중간 결과 저장 (10건 단위)
            if len(all_results) > 0 and i % 10 == 0:
                temp_df = pd.DataFrame(all_results)
                temp_df.to_csv(f"./results/{args.gene_name}_temp_extracted.csv", index=False, encoding='utf-8-sig')

        # 최종 저장
        if all_results:
            df_final = pd.DataFrame(all_results)
            df_final.fillna("Not specified", inplace=True)
            output_csv = f"./results/{args.gene_name}_hybrid_extracted_data.csv"
            df_final.to_csv(output_csv, index=False, encoding='utf-8-sig')
            logging.info(f"[Complete] Data saved to {output_csv}")
        else:
            logging.warning("[Complete] No data extracted.")

if __name__ == "__main__":
    main()