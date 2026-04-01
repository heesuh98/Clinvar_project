import pandas as pd
import logging
import config
from data_fetcher import DataFetcher
from publication_retriever import PublicationRetriever

class AnalysisPipeline:
    def __init__(self, gene_name: str):
        self.gene_name = gene_name
        self.data_fetcher = DataFetcher(self.gene_name) 
        self.pub_retriever = PublicationRetriever(config.ENTREZ_EMAIL)
        # print -> logging.info 변경
        logging.info(f"분석 파이프라인 초기화: '{self.gene_name}' (ClinVar Only / Intro 파일 제외)")

    def run(self):
        # Step 1. ClinVar 수집
        all_pmids, clinvar_df, _, _, *_ = self.data_fetcher.fetch_and_save_data()
        
        # 데이터가 없는 경우 처리
        if not all_pmids and clinvar_df.empty:
            logging.warning("수집된 데이터가 없어 파이프라인을 종료합니다.")
            return

        # Step 2. 상세 정보 수집 (Abstract, Intro 확인)
        papers_df = self.pub_retriever.fetch_details(list(all_pmids))
        
        # papers_df가 비어있거나 abstract가 없는 경우 처리
        if papers_df.empty:
            logging.warning("논문 상세 정보를 가져오지 못했습니다.")
            valid_papers_df = pd.DataFrame()
        else:
            valid_papers_df = papers_df[papers_df["abstract"].notna()].copy()
        
        if valid_papers_df.empty:
            logging.warning("유효한 초록(Abstract)이 있는 논문이 없습니다.")
            return

        # Step 3. 전체 사용
        final_df = valid_papers_df.copy()
        final_df['pmid'] = final_df['pmid'].astype(str)

        logging.info("최종 결과 파일 병합 중...")

        # Step 4. ClinVar 병합
        clinvar_df_final = pd.DataFrame()
        if not clinvar_df.empty:
            clinvar_df['pmid'] = clinvar_df['pmid'].astype(str)
            clinvar_df_renamed = clinvar_df.rename(columns={
                'hgvsc': 'clinvar_cDNA',
                'hgvsp': 'clinvar_protein',
                'rsid': 'clinvar_rsid'
            })
            clinvar_df_final = clinvar_df_renamed[[
                'pmid', 'clinical_significance', 'clinvar_cDNA', 'clinvar_protein', 'clinvar_rsid'
            ]]

        # Step 5. Merge
        if not clinvar_df_final.empty:
            final_df = pd.merge(final_df, clinvar_df_final, on='pmid', how='left')
            
        # Step 6. 출처 마킹
        clinvar_pmids = set(clinvar_df['pmid']) if not clinvar_df.empty else set()
        final_df['In_Clinvar'] = final_df['pmid'].apply(lambda x: 'O' if x in clinvar_pmids else 'X')
        
        # Step 7. 컬럼 정리 및 저장
        new_cols = ['clinvar_rsid', 'clinvar_cDNA', 'clinvar_protein', 'clinical_significance']
        for col in new_cols:
            if col not in final_df.columns: final_df[col] = ''
            else: final_df[col] = final_df[col].fillna('') 

        # introduction 컬럼 제외
        cols_order = [
            'pmid', 
            'title', 
            'journal', 
            'pubdate', 
            'abstract',      
            'clinical_significance',
            'clinvar_rsid', 'clinvar_cDNA', 'clinvar_protein',
            'In_Clinvar'
        ]
        
        existing_cols = [c for c in cols_order if c in final_df.columns]
        final_df = final_df.reindex(columns=existing_cols)
        
        output_path = config.FINAL_ANALYSIS_CSV_TPL.format(gene_name=self.gene_name)
        final_df.to_csv(output_path, index=False, encoding='utf-8-sig')
        logging.info(f"==> 완료! '{output_path}' 저장됨.")

        return final_df