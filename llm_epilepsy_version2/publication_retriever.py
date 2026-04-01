import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET
from tqdm import tqdm
import time
import logging

# [제거됨] logging.basicConfig(...) -> main.py에서 총괄 설정함

class PublicationRetriever:
    def __init__(self, email):
        Entrez.email = email
        self.chunk_size = 50

    def _extract_intro_from_xml(self, xml_content):
        """PMC XML에서 Introduction 섹션만 추출"""
        try:
            root = ET.fromstring(xml_content)
            intro_text = []
            
            for sec in root.findall(".//body/sec"):
                title = sec.find("title")
                if title is not None and title.text:
                    t_lower = title.text.lower()
                    if "intro" in t_lower or "background" in t_lower:
                        intro_text.append(f"\n--- {title.text} ---\n")
                        intro_text.append("".join(sec.itertext()))
                        break 
            
            if not intro_text:
                first_sec = root.find(".//body/sec")
                if first_sec is not None:
                    text = "".join(first_sec.itertext())
                    intro_text.append("\n--- Introduction (Assumed) ---\n")
                    intro_text.append(text[:5000]) 

            return "".join(intro_text)
        except Exception:
            return ""

    def fetch_details(self, pmids: list) -> pd.DataFrame:
        logging.info(f"총 {len(pmids)}개 논문의 상세 정보 수집 (Title/Abstract/Intro 태그 처리 강화)...")
        
        all_papers = []
        count_total = len(pmids)
        count_success = 0
        count_skipped_no_pmc = 0
        count_skipped_no_intro = 0
        
        for i in tqdm(range(0, len(pmids), self.chunk_size), desc="Fetching & Filtering"):
            chunk = pmids[i:i+self.chunk_size]
            str_ids = ",".join(map(str, chunk))
            
            try:
                handle = Entrez.efetch(db="pubmed", id=str_ids, rettype="xml", retmode="xml")
                xml_data = handle.read()
                handle.close()
                
                root = ET.fromstring(xml_data)
                
                for article in root.findall(".//PubmedArticle"):
                    pmid_txt = article.findtext(".//PMID")
                    
                    title_node = article.find(".//ArticleTitle")
                    if title_node is not None:
                        title = "".join(title_node.itertext()).strip()
                    else:
                        title = ""

                    journal = article.findtext(".//Journal/Title")
                    pubdate = article.findtext(".//PubDate/Year")
                    
                    abstract_list = article.findall(".//Abstract/AbstractText")
                    abstract_txt = "\n".join(["".join(elem.itertext()) for elem in abstract_list]).strip()
                    
                    pmc_id = None
                    for id_elem in article.findall(".//ArticleIdList/ArticleId"):
                        if id_elem.get("IdType") == "pmc":
                            pmc_id = id_elem.text
                            break
                    
                    if not pmc_id:
                        count_skipped_no_pmc += 1
                        continue
                        
                    intro_content = ""
                    try:
                        time.sleep(0.1)
                        pmc_handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml", retmode="text")
                        pmc_xml = pmc_handle.read()
                        pmc_handle.close()
                        intro_content = self._extract_intro_from_xml(pmc_xml)
                    except Exception:
                        pass 
                    
                    if not intro_content:
                        count_skipped_no_intro += 1
                        continue

                    all_papers.append({
                        "pmid": pmid_txt,
                        "title": title,
                        "journal": journal,
                        "pubdate": pubdate,
                        "abstract": abstract_txt,
                        "introduction": intro_content,
                        "has_fulltext": True
                    })
                    count_success += 1
                    
            except Exception as e:
                logging.error(f"Batch fetch error: {e}")
                continue

        logging.info("="*60)
        logging.info(f"  [Strict Filtering Result]")
        logging.info(f"  - 총 요청 PMID: {count_total}개")
        logging.info(f"  - ✅ 살아남은 논문 (본문 확보): {count_success}개")
        logging.info(f"  - ❌ 삭제됨 (PMC 링크 없음): {count_skipped_no_pmc}개")
        logging.info(f"  - ❌ 삭제됨 (본문 추출 실패): {count_skipped_no_intro}개")
        logging.info("="*60)
                
        return pd.DataFrame(all_papers)