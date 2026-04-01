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
# 1. 연도별 PMID 
#한번에 최대 호출가능한 pmid 갯수가 9999개라서 초과되면 안불러와진다.
# ========================
def fetch_pmids_by_year(query, year, retmax=10000):
    pmids = []
    retstart = 0

    year_query = f'({query}) AND ("{year}"[PDAT])'

    while True:
        handle = Entrez.esearch(
            db="pubmed",
            term=year_query,
            retmax=retmax,
            retstart=retstart
        )
        record = Entrez.read(handle)
        handle.close()

        batch_pmids = record["IdList"]
        pmids.extend(batch_pmids)

        logging.info(f"{year} | {retstart} ~ {retstart + len(batch_pmids)}")

        if len(batch_pmids) < retmax:
            break

        retstart += retmax
        time.sleep(0.3)

    return pmids


def collect_all_pmids(query, start_year=1990, end_year=2026):
    all_pmids = []

    for year in range(start_year, end_year + 1):
        logging.info(f"===== {year} 시작 =====")

        year_pmids = fetch_pmids_by_year(query, year)
        all_pmids.extend(year_pmids)

        logging.info(f"{year} 완료: {len(year_pmids)}")

    # 🔥 중복 제거
    all_pmids = list(dict.fromkeys(all_pmids))

    logging.info(f"최종 PMID 수: {len(all_pmids)}")

    return all_pmids


# ========================
# 2. PMID → PMCID (🔥 1개씩)
# ========================
def get_pmcid_from_pmid(pmid, max_retry=3):
    for attempt in range(max_retry):
        try:
            handle = Entrez.elink(
                dbfrom="pubmed",
                db="pmc",
                id=pmid
            )
            record = Entrez.read(handle)
            handle.close()

            if record and "LinkSetDb" in record[0]:
                for linksetdb in record[0]["LinkSetDb"]:
                    if linksetdb["LinkName"] == "pubmed_pmc":
                        links = linksetdb["Link"]
                        if links:
                            #main article만
                            return "PMC" + links[0]["Id"]

            return None

        except Exception as e:
            logging.warning(f"{pmid} retry {attempt+1}: {e}")
            time.sleep(2)

    return None


# ========================
# 3. PMID → PMCID 전체 변환
# ========================
def map_pmcids(pmids):
    results = []

    for idx, pmid in enumerate(pmids, 1):
        pmcid = get_pmcid_from_pmid(pmid)

        results.append({
            "pmid": pmid,
            "pmcid": pmcid
        })

        logging.info(f"[{idx}/{len(pmids)}] {pmid} → {pmcid}")

        time.sleep(0.3)  #rate limit (중요)

    return pd.DataFrame(results)


# ========================
# 4. 실행
# ========================
if __name__ == "__main__":

    query = '("EEF1A2"[All Fields]) OR ("epilepsy"[All Fields] OR "seizure"[All Fields])'

    # 1. PMID 수집
    pmids = collect_all_pmids(query, start_year=1990, end_year=2026)

    # 2. PMCID 매핑 (🔥 핵심)
    df = map_pmcids(pmids)

    # 3. 저장
    output_file = "EEF1A2_pubmed_add_epilepsy.csv"
    df.to_csv(output_file, index=False, encoding="utf-8-sig")

    print(f"\n저장 완료 → {output_file}")