#!/bin/python

import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time
from tqdm import tqdm
from Bio import Entrez
import logging

litvar_result="C:\python\Litvar\output\litvar_vv_gnomad_hgvs_result.csv"
litvar_df=pd.read_csv(litvar_result)


litvar_df['pmcid'] = litvar_df['pmcid'].fillna("")
litvar_df['pmcid'] = litvar_df['pmcid'].str.split(",")
litvar_df = litvar_df.explode('pmcid')
litvar_df['pmcid'] = litvar_df['pmcid'].str.strip()
print(litvar_df)

path= r"C:\python\Litvar\output"
pmcid_output=path+"\litvar_pmcid_match.csv"
litvar_df.to_csv(pmcid_output, index=False)



Entrez.email = "heesuh98@gmail.com"
Entrez.api_key = "17b7a8642acb4874f49ddced130220e1eb09" 
Entrez.tool = "eef1a2_fulltext_script"

##get pubmed article(pmcid article)

df_test = litvar_df.head(4).copy()
# df_test = litvar_df

# 🔹 pmcid가 있는 것만 추출
unique_pmcids = df_test["pmcid"].dropna().unique()

# 🔹 pmcid → fulltext 저장할 딕셔너리
pmc_fulltext_dict = {}

def fetch_fulltext_from_pmc(pmcid, max_retry=3):

    if not pmcid:
        return None

    pmc_id = pmcid.replace("PMC", "")

    for attempt in range(max_retry):

        try:
            handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="xml",
                retmode="text"
            )

            xml_content = handle.read()
            handle.close()

            root = ET.fromstring(xml_content)
            body = root.find(".//body")

            if body is not None:
                return "".join(body.itertext()).strip()

            return None

        except Exception as e:
            logging.warning(f"{pmcid} retry {attempt+1}: {e}")
            time.sleep(2)

    return None

# 🔹 fulltext 수집
for pmcid in tqdm(unique_pmcids):
    fulltext = fetch_fulltext_from_pmc(pmcid)
    pmc_fulltext_dict[pmcid] = fulltext
    time.sleep(0.3)  # NCBI rate limit 방지

# 🔹 원래 df에 매핑
df_test["fulltext"] = df_test["pmcid"].map(pmc_fulltext_dict)

# 🔹 저장
df_test.to_csv(r"C:\python\Litvar\output\litvar_df_with_fulltext_added.csv", index=False)



# ##3. get pumed article


# Entrez.email = "your_real_email@example.com"
# # 로깅 설정
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s - %(levelname)s - %(message)s'
# )



# # fulltext 컬럼 초기화
# litvar_df['fulltext'] = None
# litvar_df['fulltext_source'] = None


# class FullTextFetcher:
#     """PMC와 Europe PMC에서 full text를 가져오는 클래스"""

#     def __init__(self):
#         self.session = requests.Session()
#         self.stats = {
#             'pmc_success': 0,
#             'europepmc_success': 0,
#             'abstract_only': 0,
#             'failed': 0
#         }

#     def extract_text_from_xml(self, xml_content):
#         """PMC XML에서 텍스트 추출"""
#         try:
#             root = ET.fromstring(xml_content)

#             # 전체 본문 추출
#             sections = []

#             # Abstract
#             abstract_elem = root.find(".//abstract")
#             if abstract_elem is not None:
#                 abstract_text = "".join(abstract_elem.itertext())
#                 sections.append(f"=== ABSTRACT ===\n{abstract_text}\n")

#             # Body sections
#             body = root.find(".//body")
#             if body is not None:
#                 for sec in body.findall(".//sec"):
#                     title_elem = sec.find("title")
#                     if title_elem is not None and title_elem.text:
#                         section_title = title_elem.text.strip()
#                     else:
#                         section_title = "Section"
#                     section_text = "".join(sec.itertext())
#                     sections.append(f"\n=== {section_title.upper()} ===\n{section_text}\n")

#             return "\n".join(sections) if sections else None

#         except Exception as e:
#             logging.warning(f"XML 파싱 오류: {e}")
#             return None

#     def fetch_from_pmc(self, pmcid):
#         """PMC에서 full text 가져오기"""
#         if not pmcid or pd.isna(pmcid) or pmcid == '':
#             return None

#         try:
#             # PMC prefix 제거 (숫자만 추출)
#             pmc_id = pmcid.replace('PMC', '').strip()

#             # PMC efetch 호출
#             handle = Entrez.efetch(
#                 db="pmc",
#                 id=pmc_id,
#                 rettype="xml",
#                 retmode="text"
#             )
#             xml_content = handle.read()
#             handle.close()

#             # XML에서 텍스트 추출
#             text = self.extract_text_from_xml(xml_content)

#             if text and len(text) > 100:  # 최소 길이 체크
#                 self.stats['pmc_success'] += 1
#                 return text

#             return None

#         except Exception as e:
#             logging.debug(f"PMC {pmcid} 다운로드 실패: {e}")
#             return None

#     def fetch_from_europepmc(self, pmid):
#         """Europe PMC에서 full text 가져오기"""
#         if not pmid or pd.isna(pmid):
#             return None

#         try:
#             # Europe PMC API
#             url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmid}/fullTextXML"
#             response = self.session.get(url, timeout=30)

#             if response.status_code == 200:
#                 text = self.extract_text_from_xml(response.content)
#                 if text and len(text) > 100:
#                     self.stats['europepmc_success'] += 1
#                     return text

#             return None

#         except Exception as e:
#             logging.debug(f"Europe PMC PMID {pmid} 다운로드 실패: {e}")
#             return None

#     def fetch_abstract(self, pmid):
#         """PubMed에서 abstract만 가져오기 (최후의 수단)"""
#         if not pmid or pd.isna(pmid):
#             return None

#         try:
#             handle = Entrez.efetch(
#                 db="pubmed",
#                 id=pmid,
#                 rettype="xml",
#                 retmode="xml"
#             )
#             xml_content = handle.read()
#             handle.close()

#             root = ET.fromstring(xml_content)

#             # Title 추출
#             title_elem = root.find(".//ArticleTitle")
#             title = "".join(title_elem.itertext()).strip() if title_elem is not None else ""

#             # Abstract 추출
#             abstract_elems = root.findall(".//Abstract/AbstractText")
#             abstract = "\n".join(["".join(elem.itertext()) for elem in abstract_elems]).strip()

#             if abstract:
#                 self.stats['abstract_only'] += 1
#                 return f"=== TITLE ===\n{title}\n\n=== ABSTRACT ONLY ===\n{abstract}"

#             return None

#         except Exception as e:
#             logging.debug(f"Abstract 가져오기 실패 PMID {pmid}: {e}")
#             return None

#     def fetch_fulltext_multi_source(self, pmid, pmcid):
#         """여러 소스를 시도하여 full text 가져오기"""

#         # 1순위: PMC (가장 안정적)
#         if pmcid:
#             # 쉼표로 구분된 경우 첫 번째 것만 사용
#             first_pmcid = str(pmcid).split(',')[0].strip()
#             fulltext = self.fetch_from_pmc(first_pmcid)
#             if fulltext:
#                 return fulltext, 'PMC'

#         # 2순위: Europe PMC
#         if pmid:
#             # 쉼표로 구분된 경우 첫 번째 것만 사용
#             first_pmid = str(pmid).split(',')[0].strip()
#             fulltext = self.fetch_from_europepmc(first_pmid)
#             if fulltext:
#                 return fulltext, 'EuropePMC'

#             # 3순위: Abstract만
#             abstract = self.fetch_abstract(first_pmid)
#             if abstract:
#                 return abstract, 'AbstractOnly'

#         self.stats['failed'] += 1
#         return None, 'Failed'


# def main():
#     """메인 실행 함수"""
#     path= r"C:\python\Litvar\output"

#     csv_path = pmcid_output
#     output_path = path+'\litvar_df_with_fulltext.csv'

#     logging.info(f"CSV 파일 로딩: {csv_path}")
#     df = pd.read_csv(csv_path)

#     logging.info(f"총 {len(df)}개의 행 발견")

#     # Full text 저장할 새로운 컬럼 초기화
#     df['fulltext'] = None
#     df['fulltext_source'] = None

#     # Full text fetcher 초기화
#     fetcher = FullTextFetcher()

#     # 각 행에 대해 full text 다운로드
#     logging.info("Full text 다운로드 시작...")

#     for idx, row in tqdm(df.iterrows(), total=len(df), desc="Full text 다운로드"):
#         pmid = row.get('pmid')
#         pmcid = row.get('pmcid')

#         # Full text 다운로드
#         fulltext, source = fetcher.fetch_fulltext_multi_source(pmid, pmcid)

#         if fulltext:
#             df.at[idx, 'fulltext'] = fulltext
#             df.at[idx, 'fulltext_source'] = source
#         else:
#             df.at[idx, 'fulltext'] = 'Not available'
#             df.at[idx, 'fulltext_source'] = 'Failed'

#         # API Rate limit 준수
#         time.sleep(0.4)

#     # 결과 저장
#     logging.info(f"결과 저장 중: {output_path}")
#     df.to_csv(output_path, index=False, encoding='utf-8-sig')

#     # 통계 출력
#     logging.info("\n" + "="*60)
#     logging.info("다운로드 통계:")
#     logging.info(f"  PMC 성공: {fetcher.stats['pmc_success']}개")
#     logging.info(f"  Europe PMC 성공: {fetcher.stats['europepmc_success']}개")
#     logging.info(f"  Abstract만: {fetcher.stats['abstract_only']}개")
#     logging.info(f"  실패: {fetcher.stats['failed']}개")
#     logging.info(f"  총 성공률: {((len(df) - fetcher.stats['failed']) / len(df) * 100):.1f}%")
#     logging.info("="*60)

#     # Full text 있는 행 비율
#     has_fulltext = df[df['fulltext_source'].isin(['PMC', 'EuropePMC'])].shape[0]
#     logging.info(f"\nFull text 확보: {has_fulltext}/{len(df)} ({has_fulltext/len(df)*100:.1f}%)")
#     logging.info(f"\n완료! 결과 파일: {output_path}")


# if __name__ == "__main__":
#     main()

# fetcher = FullTextFetcher()

# for idx, row in tqdm(litvar_df.iterrows(), total=len(litvar_df)):
    
#     pmcid = row.get('pmcid')
    
#     # pmcid만 사용 (요청대로)
#     fulltext = None
#     source = None

#     if pd.notna(pmcid) and str(pmcid).strip() != "":
#         first_pmcid = str(pmcid).split(',')[0].strip()
#         fulltext = fetcher.fetch_from_pmc(first_pmcid)
        
#         if fulltext:
#             source = "PMC"
    
#     if fulltext:
#         litvar_df.at[idx, 'fulltext'] = fulltext
#         litvar_df.at[idx, 'fulltext_source'] = source
#     else:
#         litvar_df.at[idx, 'fulltext'] = None
#         litvar_df.at[idx, 'fulltext_source'] = "Failed"

#     # Rate limit
#     time.sleep(0.4)
# litvar_df['fulltext_source'].value_counts()