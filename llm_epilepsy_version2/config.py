import os
# ==============================
# [기본 설정]
# ==============================
GENE_NAME = "EEF1A2"   # 기본 타겟 유전자
DOWNLOAD_ROOT = "downloads" # 다운로드 폴더 자동 생성
ENTREZ_EMAIL = "email@example.com" # 굳이 입력 안하시고 이대로 사용해도 문제는 없습니다.
FINAL_ANALYSIS_CSV_TPL = "./results/{gene_name}_final_analysis.csv" #main.py 실행 결과 저장 파일명
# ==============================
# [API 키 설정]
# ==============================
OPENAI_API_KEY = "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"

# ==============================
# [LLM 모델 설정]
# ==============================
OPENAI_MODEL_NAME = "gpt-5.1"