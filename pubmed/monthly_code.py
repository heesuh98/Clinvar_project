#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import pandas as pd
from Bio import Entrez
import logging

# ========================
# 설정
# ========================
Entrez.email = "heesuh98@gmail.com"
Entrez.api_key = "17b7a8642acb4874f49ddced130220e1eb09"
Entrez.tool = "eef1a2_fulltext_script"

logging.basicConfig(level=logging.INFO)

# ========================
# 공통 fetch 함수
# ========================
def fetch_pmids_by_range(query, retmax=10000):
    pmids = []
    retstart = 0

    while True:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=retmax,
            retstart=retstart
        )
        record = Entrez.read(handle)
        handle.close()

        batch_pmids = record["IdList"]
        pmids.extend(batch_pmids)

        logging.info(f"{retstart} ~ {retstart + len(batch_pmids)}")

        if len(batch_pmids) < retmax:
            break

        retstart += retmax
        time.sleep(0.3)

    return pmids


# ========================
# 연도 → 필요 시 월 분할
# ========================
def fetch_pmids_by_year(query, year):
    pmids = []

    year_query = f'({query}) AND ("{year}"[PDAT])'

    # 🔹 전체 개수 확인
    handle = Entrez.esearch(db="pubmed", term=year_query, retmax=0)
    record = Entrez.read(handle)
    handle.close()

    count = int(record["Count"])
    logging.info(f"{year} 총 논문 수: {count}")

    # ✅ 9999 이하 → 바로 수집
    if count <= 9999:
        return fetch_pmids_by_range(year_query)

    # 🔥 9999 초과 → 월 분할
    logging.info(f"{year} → 월 단위 분할")

    for month in range(1, 13):
        start_date = f"{year}/{month:02d}/01"
        end_date = f"{year}/{month:02d}/31"

        month_query = f'({query}) AND ("{start_date}"[PDAT] : "{end_date}"[PDAT])'

        # 월별 count 확인
        handle = Entrez.esearch(db="pubmed", term=month_query, retmax=0)
        record = Entrez.read(handle)
        handle.close()

        month_count = int(record["Count"])

        # 🔹 월도 9999 초과 → day 분할
        if month_count > 9999:
            logging.info(f"{year}-{month:02d} → day 분할")

            for day in range(1, 32):
                start_day = f"{year}/{month:02d}/{day:02d}"
                end_day = start_day

                day_query = f'({query}) AND ("{start_day}"[PDAT] : "{end_day}"[PDAT])'

                day_pmids = fetch_pmids_by_range(day_query)
                pmids.extend(day_pmids)

                logging.info(f"{year}-{month:02d}-{day:02d}: {len(day_pmids)}")

        else:
            month_pmids = fetch_pmids_by_range(month_query)
            pmids.extend(month_pmids)

            logging.info(f"{year}-{month:02d}: {len(month_pmids)}")

    return pmids


# ========================
# 전체 PMID 수집
# ========================
def collect_all_pmids(query, start_year=1990, end_year=2026):
    all_pmids = []

    for year in range(start_year, end_year + 1):
        logging.info(f"===== {year} 시작 =====")

        year_pmids = fetch_pmids_by_year(query, year)
        all_pmids.extend(year_pmids)

        logging.info(f"{year} 완료: {len(year_pmids)}")

    # 🔥 중복 제거 (순서 유지)
    all_pmids = list(dict.fromkeys(all_pmids))

    logging.info(f"최종 PMID 수 (unique): {len(all_pmids)}")

    return all_pmids


# ========================
# PMID → PMCID (개별 호출)
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
                            return "PMC" + links[0]["Id"]  # main article

            return None

        except Exception as e:
            logging.warning(f"{pmid} retry {attempt+1}: {e}")
            time.sleep(2)

    return None


# ========================
# PMID → PMCID 전체 매핑
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

        time.sleep(0.3)  # rate limit

    return pd.DataFrame(results)


# ========================
# 실행
# ========================
if __name__ == "__main__":

    query = '("EEF1A2"[All Fields]) OR ("epilepsy"[All Fields] OR "seizure"[All Fields])'

    # 1️⃣ PMID 수집
    pmids = collect_all_pmids(query, start_year=1990, end_year=2026)

    # 2️⃣ PMCID 매핑
    df = map_pmcids(pmids)

    # 3️⃣ 저장
    output_file = "EEF1A2_pubmed_add_monthly_epilepsy.csv"
    df.to_csv(output_file, index=False, encoding="utf-8-sig")

    print(f"\n저장 완료 → {output_file}")