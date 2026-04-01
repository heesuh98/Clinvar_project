import pandas as pd
from Bio import Entrez
from typing import Set, List, Tuple, Dict
import time
import config
from tqdm import tqdm
import requests
import xml.etree.ElementTree as ET 
import logging 
import re 
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# [제거됨] logging.basicConfig(...) -> main.py에서 총괄 설정함

class DataFetcher:
    def __init__(self, gene_name: str):
        Entrez.email = config.ENTREZ_EMAIL
        self.gene_name = gene_name
            
        logging.info(f"DataFetcher가 초기화되었습니다. (Gene: {self.gene_name})")
        retry_strategy = Retry(total=5, status_forcelist=[429, 500, 502, 503, 504], backoff_factor=1)
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session = requests.Session()
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)

    def _fetch_clinvar_data_efetch(self) -> pd.DataFrame:
        logging.info(f"\nClinVar 검색 시작: '{self.gene_name}' [Single Gene]...")
        server = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        # 1. Esearch
        search_term = f"\"{self.gene_name}\"[Gene Name] AND \"single gene\"[properties]"
        clinvar_uids = []
        try:
            logging.info(f"Step 1: 검색 쿼리 -> {search_term}")
            esearch_url = f"{server}esearch.fcgi?db=clinvar&term={search_term}&retmode=json&retmax=10000" 
            r_search = self.session.get(esearch_url)
            r_search.raise_for_status()
            search_data = r_search.json()['esearchresult']
            
            count = int(search_data['count'])
            clinvar_uids = search_data.get('idlist', [])
            
            logging.info(f"-> 검색 결과: 총 {count}개의 변이(UID) 발견.")
            
            if count == 0 or not clinvar_uids:
                logging.warning("-> 검색된 변이가 없습니다.")
                return pd.DataFrame()

        except Exception as e:
            logging.error(f"Esearch Error: {e}")
            return pd.DataFrame()

        # --- 2. Esummary ---
        rcv_ids_to_fetch = set()
        for i in tqdm(range(0, len(clinvar_uids), 200), desc="ClinVar Esummary"):
            id_chunk = clinvar_uids[i:i+200]
            try:
                handle = Entrez.esummary(db="clinvar", id=",".join(id_chunk), retmode="xml")
                xml_text_bytes = handle.read()
                handle.close()
                if not xml_text_bytes: continue
                
                root = ET.fromstring(xml_text_bytes)
                for summary in root.findall(".//DocumentSummary"):
                    desc_elem = summary.find("./germline_classification/description")
                    if desc_elem is None or desc_elem.text is None: continue
                    significance = desc_elem.text.lower()
                    
                    if not any(s in significance for s in ["pathogenic", "likely pathogenic", "uncertain significance"]):
                        continue
                    
                    for rcv_elem in summary.findall(".//supporting_submissions/rcv/string"):
                        if rcv_elem.text: rcv_ids_to_fetch.add(rcv_elem.text.strip())
                time.sleep(0.4)
            except: pass

        if not rcv_ids_to_fetch: return pd.DataFrame()
        rcv_id_list = list(rcv_ids_to_fetch)
        
        # --- 3. Efetch ---
        all_variant_data = [] 
        
        for i in tqdm(range(0, len(rcv_id_list), 100), desc="ClinVar Efetch"):
            id_chunk = rcv_id_list[i:i+100]
            try:
                handle = Entrez.efetch(db="clinvar", id=",".join(id_chunk), rettype="clinvarset", retmode="xml")
                xml_text_bytes = handle.read() 
                handle.close()
                if not xml_text_bytes: continue

                root = ET.fromstring(xml_text_bytes) 
                clinvar_sets = root.findall('.//ClinVarSet') 

                for clinvar_set in clinvar_sets:
                    measure = clinvar_set.find(".//MeasureSet/Measure")
                    if measure is not None:
                        v_type = measure.get("Type", "").lower()
                        if "copy number" in v_type: continue

                    found_gene = False
                    relationships = clinvar_set.findall(".//MeasureSet/Measure/MeasureRelationship")
                    for rel in relationships:
                        symbol_elem = rel.find("Symbol/ElementValue[@Type='Preferred']")
                        if symbol_elem is None or symbol_elem.text != self.gene_name: 
                            continue
                        
                        rel_type = rel.get("Type", "").lower()
                        if "variant in gene" in rel_type or "within single gene" in rel_type:
                            found_gene = True; break
                    
                    if not found_gene: continue

                    set_title = ""
                    title_elem = clinvar_set.find("./Title")
                    if title_elem is not None:
                        set_title = "".join(title_elem.itertext()).strip()
                    
                    if not set_title and measure is not None:
                        name_elem = measure.find("Name/ElementValue[@Type='Preferred']")
                        if name_elem is not None:
                            set_title = "".join(name_elem.itertext()).strip()
                        else:
                            set_title = "No Title Provided"
                    
                    set_title = " ".join(set_title.split())

                    all_significance_text = []
                    ref_sig_elem = clinvar_set.find(".//ReferenceClinVarAssertion/Classifications/GermlineClassification/Description")
                    if ref_sig_elem is not None and ref_sig_elem.text:
                        all_significance_text.append(ref_sig_elem.text.lower())
                    for assert_sig_elem in clinvar_set.findall(".//ClinVarAssertion/Classification/GermlineClassification"):
                        if assert_sig_elem is not None and assert_sig_elem.text:
                            all_significance_text.append(assert_sig_elem.text.lower())
                    
                    if not any(s in sig for sig in all_significance_text for s in ["pathogenic", "likely pathogenic", "uncertain significance"]):
                        continue
                    
                    significance = ref_sig_elem.text if (ref_sig_elem is not None and ref_sig_elem.text) else "Not specified"
                    if (significance == "Not specified" or "benign" in significance.lower()) and all_significance_text:
                        significance = next((s for s in all_significance_text if any(p in s for p in ["pathogenic", "likely pathogenic", "uncertain significance"])), "Not specified").capitalize()
                    
                    hgvsc = "Not specified"
                    hgvsp = "Not specified"
                    rsid = "Not specified"
                    
                    if measure is not None:
                        hgvsc_elem = measure.find("./Name/ElementValue[@Type='HGVS,coding,RefSeq']")
                        if hgvsc_elem is not None: hgvsc = hgvsc_elem.text.strip()
                        hgvsp_elem = measure.find("./Name/ElementValue[@Type='HGVS,protein,RefSeq']")
                        if hgvsp_elem is not None: hgvsp = hgvsp_elem.text.strip()
                        rsid_elem = measure.find("./XRef[@DB='dbSNP']")
                        if rsid_elem is not None: rsid = "rs" + rsid_elem.get('ID').strip()
                    
                    if hgvsp == "Not specified":
                        p_match = re.search(r'(p\.[^)\s]+)', set_title)
                        if p_match: hgvsp = p_match.group(1)
                    if hgvsc == "Not specified":
                        c_match = re.search(r'(c\.[^)\s]+)', set_title)
                        g_match = re.search(r'(g\.[^)\s]+)', set_title)
                        if c_match: hgvsc = c_match.group(1)
                        elif g_match: hgvsc = g_match.group(1)
                    
                    if hgvsc == "Not specified" and hgvsp == "Not specified" and rsid == "Not specified":
                        continue

                    pmids = set()
                    path1 = clinvar_set.findall(".//ReferenceClinVarAssertion/Citation/ID[@Source='PubMed']")
                    path2 = clinvar_set.findall(".//ClinVarAssertion/Citation/ID[@Source='PubMed']")
                    path3 = clinvar_set.findall(".//ClinVarAssertion/Classification/Citation/ID[@Source='PubMed']") 
                    for citation in path1: pmids.add(citation.text)
                    for citation in path2: pmids.add(citation.text)
                    for citation in path3: pmids.add(citation.text)

                    if not pmids: continue 
                    
                    for pmid in pmids:
                        all_variant_data.append({
                            "pmid": int(pmid),
                            "hgvsc": hgvsc, 
                            "hgvsp": hgvsp, 
                            "rsid": rsid,
                            "clinical_significance": significance,
                            "title": set_title
                        })
                time.sleep(0.4) 
            except Exception as e:
                logging.error(f"ClinVar Efetch Error: {e}")

        if not all_variant_data: return pd.DataFrame()
        
        clinvar_df = pd.DataFrame(all_variant_data)
        clinvar_df = clinvar_df.fillna("Not specified")
        
        group_cols = ['pmid', 'hgvsc', 'hgvsp', 'rsid']
        
        def clean_and_merge_sig(series):
            unique_vals = set()
            for val in series:
                parts = str(val).replace(',', '/').split('/')
                for part in parts:
                    clean = part.strip()
                    if "likely" in clean.lower() and "pathogenic" in clean.lower(): clean = "Likely pathogenic"
                    elif "pathogenic" in clean.lower() and "likely" not in clean.lower() and "conflict" not in clean.lower(): clean = "Pathogenic"
                    elif "uncertain" in clean.lower(): clean = "Uncertain significance"
                    
                    if clean and clean.lower() != "not specified":
                        unique_vals.add(clean)
            if not unique_vals: return "Not specified"
            return ', '.join(sorted(list(unique_vals)))

        def pick_best_title(series):
            titles = [str(t).strip() for t in series if t and str(t).lower() != "no title provided"]
            if not titles: return "No Title Provided"
            return max(titles, key=len)

        if not clinvar_df.empty:
            clinvar_df = clinvar_df.groupby(group_cols, as_index=False).agg({
                'title': pick_best_title,
                'clinical_significance': clean_and_merge_sig
            })

        clinvar_df = clinvar_df.drop_duplicates()

        unique_pmid_count = clinvar_df['pmid'].nunique()
        logging.info(f"-> ClinVar 추출 완료: {len(clinvar_df)}행 (PMID: {unique_pmid_count}개)")
        
        return clinvar_df

    def fetch_and_save_data(self) -> Tuple[Set[int], pd.DataFrame, Dict[int, List[str]], Set[int]]:
        logging.info(f"\n--- 1단계: '{self.gene_name}' 변이 정보 및 PMID 수집 (ClinVar Only) ---")
        clinvar_df = self._fetch_clinvar_data_efetch()
        clinvar_pmids = set(clinvar_df['pmid']) if not clinvar_df.empty else set()
        logging.info(f"\n==> 총 {len(clinvar_pmids)}개의 고유 PMID를 수집하여 다음 단계로 전달합니다.")
        return clinvar_pmids, clinvar_df, {}, set(), None, None