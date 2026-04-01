#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
import pandas as pd

# ========================
# 설정
# ========================
Entrez.email = "heesuh98@gmail.com"
Entrez.api_key = "17b7a8642acb4874f49ddced130220e1eb09"
Entrez.tool = "eef1a2_fulltext_script"

# ========================
# 전체 논문 수
# ========================
def get_total_count(query):
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=0
    )
    record = Entrez.read(handle)
    handle.close()

    return int(record["Count"])


# ========================
# 연도별 논문 수
# ========================
def get_yearly_counts(query, start_year=1990, end_year=2026):
    results = []

    for year in range(start_year, end_year + 1):
        year_query = f'({query}) AND ("{year}"[PDAT])'

        handle = Entrez.esearch(
            db="pubmed",
            term=year_query,
            retmax=0
        )
        record = Entrez.read(handle)
        handle.close()

        count = int(record["Count"])

        print(f"{year}: {count}")

        results.append({
            "year": year,
            "count": count
        })

    return pd.DataFrame(results)


# ========================
# 실행
# ========================
if __name__ == "__main__":

    query = '("EEF1A2"[All Fields]) OR ("epilepsy"[All Fields] OR "seizure"[All Fields])'

    # 1️⃣ 전체 count
    total_count = get_total_count(query)
    print(f"\n전체 논문 수: {total_count}")

    # 2️⃣ 연도별 count
    df_year = get_yearly_counts(query, start_year=1990, end_year=2026)

    # 3️⃣ CSV 저장
    output_file = "pubmed_yearly_counts.csv"
    df_year.to_csv(output_file, index=False, encoding="utf-8-sig")

    print(f"\n연도별 결과 저장 완료 → {output_file}")