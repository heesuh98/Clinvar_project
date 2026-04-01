#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET
import logging

# ========================
# 설정
# ========================
Entrez.email = "heesuh98@gmail.com"
Entrez.api_key = "17b7a8642acb4874f49ddced130220e1eb09"
Entrez.tool = "eef1a2_fulltext_script"

logging.basicConfig(level=logging.INFO)

# ========================
# 1. PMID 전체 가져오기
# ========================
def search_gene_pubmed_all(gene_name):
    query = f"{gene_name}[All Fields]"

    handle = Entrez.esearch(db="pubmed", term=query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()

    pmids = record["IdList"]
    logging.info(f"총 PMID 수: {len(pmids)}")

    return pmids


# ========================
# 2. PMID → PMCID (개별 + multi 처리)
# ========================
def get_pmcid_from_pmid(pmid):
    try:
        handle = Entrez.elink(
            dbfrom="pubmed",
            db="pmc",
            id=pmid
        )
        record = Entrez.read(handle)
        handle.close()

        pmcids = []

        if record and "LinkSetDb" in record[0]:
            for linksetdb in record[0]["LinkSetDb"]:
                if linksetdb["LinkName"] == "pubmed_pmc":
                    for link in linksetdb["Link"]:
                        pmcids.append("PMC" + link["Id"])

        return pmcids if pmcids else None

    except Exception as e:
        logging.warning(f"PMCID fetch error {pmid}: {e}")
        return None


# ========================
# 3. PMCID → Full text
# ========================
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


# ========================
# 4. RUN
# ========================
def main():

    gene = "EEF1A2"
    output_path = f"{gene}_fulltext_results.csv"

    # 1. PMID 가져오기
    pmids = search_gene_pubmed_all(gene)

    results = []

    # 2. PMID → PMCID → Fulltext (순차 처리)
    for idx, pmid in enumerate(pmids, 1):

        pmcid_list = get_pmcid_from_pmid(pmid)

        # PMCID 없는 경우
        if not pmcid_list:
            results.append({
                "gene": gene,
                "pmid": pmid,
                "pmcid": None,
                "fulltext": None
            })
            logging.info(f"[{idx}/{len(pmids)}] PMCID 없음: {pmid}")
            time.sleep(0.3)
            continue

        # 🔥 여러 PMCID 처리
        for pmcid in pmcid_list:
            fulltext = fetch_fulltext_from_pmc(pmcid)

            results.append({
                "gene": gene,
                "pmid": pmid,
                "pmcid": pmcid,
                "fulltext": fulltext
            })

        logging.info(f"[{idx}/{len(pmids)}] 완료: {pmid}")
        time.sleep(0.5)   # rate limit

    # 3. 저장
    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False, encoding="utf-8-sig")

    logging.info("저장 완료")
    logging.info(f"총 결과 수: {len(df)}")


# ========================
# 실행
# ========================
if __name__ == "__main__":
    main()
