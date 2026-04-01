#!/bin/python

import requests
import logging
import time
from tqdm.auto import tqdm
from urllib.parse import quote
import ast
import pandas as pd

print("----run litvar----")

##path##
path= r"C:\python\Litvar\output"
# output_file=path+"\pmid_hgvs.txt"
output_file="/Users/snuh4/python/project/ep_llm/Litvar/output/pmid_hgvs.txt"

class LitVarFetcher:

    def __init__(self, gene_name: str):
        self.gene_name = gene_name
        self.session = requests.Session()
        logging.basicConfig(level=logging.INFO)

    # 1. Gene 기반 variant 목록 가져오기
    def fetch_gene_variants(self):
        base_url = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/search/gene/"
        url = base_url + self.gene_name

        try:
            response = self.session.get(url, timeout=20)
            response.raise_for_status()
        except Exception as e:
            logging.error(f"LitVar gene search API 오류: {e}")
            return []

        variants = []
        lines = response.text.strip().split("\n")

        for line in lines:
            try:
                d = ast.literal_eval(line.strip())
                variants.append(d)
            except Exception:
                continue

        return variants

    # 2. PMID / PMCID 가져오기
    def fetch_pmids_for_variant(self, variant_id: str):
        base_url = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/"
        encoded = quote(variant_id)
        url = f"{base_url}{encoded}/publications"

        try:
            r = self.session.get(url, timeout=10)
            r.raise_for_status()
            j = r.json()

            pmids = j.get("pmids") or []
            pmcids = j.get("pmcids") or []
            return pmids, pmcids

        except Exception as e:
            logging.warning(f"PMID/PMCID 조회 실패 ({encoded}): {e}")
            return [], []

    # 3. HGVS 가져오기
    def fetch_hgvs_for_variant(self, variant_id: str):
        base_url = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/"
        encoded = quote(variant_id)
        url = f"{base_url}{encoded}"

        try:
            r = self.session.get(url, timeout=10)
            r.raise_for_status()
            j = r.json()

            return j.get("hgvs") or []

        except Exception as e:
            logging.warning(f"HGVS 조회 실패 ({encoded}): {e}")
            return []

    # 4. DataFrame 생성
    def build_variant_dataframe(self):

        variants = self.fetch_gene_variants()

        if not variants:
            logging.error("변이 정보를 가져오지 못했습니다.")
            return pd.DataFrame()

        rows = []

        for v in tqdm(variants, desc="LitVar Variant 처리 중"):

            variant_id = v.get("_id")
            if not variant_id:
                continue

            rsid = v.get("rsid")
            gene_list = v.get("gene", [])

            pmids, pmcids = self.fetch_pmids_for_variant(variant_id)
            hgvs = self.fetch_hgvs_for_variant(variant_id)
            pmid_value = ",".join(map(str, pmids)) if pmids else None
            pmcid_value = ",".join(map(str, pmcids)) if pmcids else None
            if isinstance(hgvs, list):
                hgvs_value = str(hgvs) if hgvs else None
            else:
                hgvs_value = str(hgvs) if hgvs else None

            rows.append({
                    "_id": variant_id,
                    "rsid": rsid,
                    "gene": ",".join(gene_list),
                    "pmid": pmid_value,
                    "pmcid": pmcid_value,
                    "hgvs": hgvs if hgvs else None
                })


            time.sleep(0.3)  # API 부하 방지

        df = pd.DataFrame(rows)
        return df


gene = "EEF1A2"
fetcher = LitVarFetcher(gene)
df = fetcher.build_variant_dataframe()
df["gene"] = gene

df.to_csv(output_file, sep="\t", index=False)
