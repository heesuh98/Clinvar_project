import argparse
import logging
import os
import datetime
from pipeline import AnalysisPipeline

def setup_logging(gene_name: str):
    """
    로그 설정을 초기화합니다.
    - 포맷: Clinvar_{gene_name}_{timestamp}.log
    - 위치: logs 폴더 (없으면 생성)
    """
    # 로그 디렉토리 생성
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)

    # 타임스탬프 생성
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = f"Clinvar_{gene_name}_{timestamp}.log"
    log_path = os.path.join(log_dir, log_filename)

    # 기존 핸들러 제거 (중복 방지)
    root_logger = logging.getLogger()
    if root_logger.handlers:
        root_logger.handlers = []

    # 로깅 설정 (파일 + 콘솔 동시 출력)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path, encoding='utf-8'), # 파일 저장
            logging.StreamHandler()                          # 화면 출력
        ],
        force=True # 파이썬 3.8+ : 다른 모듈의 설정을 덮어씀
    )
    
    logging.info(f"=== 로그 기록 시작: {log_path} ===")
    return log_path

def main():
    """
    커맨드 라인에서 유전자 이름을 입력받아 파이프라인을 실행합니다.
    """
    parser = argparse.ArgumentParser(description="유전자 관련 논문 분석 파이프라인")
    parser.add_argument("gene_name", type=str, help="분석할 유전자의 이름 (예: SCN2A)")
    
    args = parser.parse_args()
    
    # [중요] 로깅 설정 먼저 수행
    setup_logging(args.gene_name)
    
    try:
        logging.info(f"분석 시작: Gene Name = {args.gene_name}")
        pipeline = AnalysisPipeline(gene_name=args.gene_name)
        pipeline.run()
    except Exception as e:
        logging.error(f"파이프라인 실행 중 치명적 오류 발생: {e}", exc_info=True)

if __name__ == "__main__":
    main()