#!/bin/python

import requests
import pandas as pd
import urllib.parse

path= r"C:\python\Litvar\output"
# output_file=path+"\pmid_hgvs.txt"
output_file="/Users/snuh4/python/project/ep_llm/Litvar/output/pmid_hgvs.txt"

def get_hgvs_from_dbsnp(rsid):
    """
    dbSNP RefSNP API에서 주어진 rsID의 HGVS(cDNA, protein) 리스트 추출
    """
    rs_only = rsid.replace("rs", "")
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rs_only}"

    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {rsid}: {e}")
        return [], []

    data = r.json()

    hgvsc_list = set()
    hgvsp_list = set()

    primary_snapshot_data = data.get("primary_snapshot_data", {})
    placements = primary_snapshot_data.get("placements_with_allele", [])

    for placement in placements:
        seq_id = placement.get("seq_id", "")
        for allele in placement.get("alleles", []):
            hgvs = allele.get("hgvs")
            if not hgvs:
                continue
            # cDNA (coding) : NM_ 또는 NG_ 시작
            if seq_id.startswith("NM_") or seq_id.startswith("NG_"):
                hgvsc_list.add(hgvs)
            # protein : NP_ 시작
            elif seq_id.startswith("NP_"):
                hgvsp_list.add(hgvs)

    return list(hgvsc_list), list(hgvsp_list)


# 사용 예시

file_path = output_file
litvar_df = pd.read_csv(file_path, sep='\t', dtype=str)
litvar_df = litvar_df.dropna(subset=['rsid']).copy()

# dbSNP HGVS 컬럼 초기화
litvar_df['dbsnp_hgvsc'] = None
litvar_df['dbsnp_hgvsp'] = None

for rsid in litvar_df['rsid'].unique():
    # print(f"Fetching {rsid} ...")
    hgvsc_list, hgvsp_list = get_hgvs_from_dbsnp(rsid)

    # print("c.HGVS (hgvsc):", hgvsc_list)
    # print("p.HGVS (hgvsp):", hgvsp_list)
    # print("############")

    # 리스트 → 문자열 변환 (중복 제거 + 안전 처리)
    hgvsc_str = ",".join(sorted(set(hgvsc_list))) if hgvsc_list else None
    hgvsp_str = ",".join(sorted(set(hgvsp_list))) if hgvsp_list else None

    # 해당 rsid 행 전체에 동일 값 넣기
    litvar_df.loc[litvar_df['rsid'] == rsid, 'dbsnp_hgvsc'] = hgvsc_str
    litvar_df.loc[litvar_df['rsid'] == rsid, 'dbsnp_hgvsp'] = hgvsp_str

# ===========================================
# 5) litvar_df에 VariantValidator 정보 채우기
# ===========================================
#dbsnp에서 chromosome position 확인하기
import requests
import urllib.parse
# url2="http://grch37.rest.ensembl.org/variation/human/rs116035550?content-type=application/json"

def get_chr_pos_from_dbsnp(rsid):
    """
    dbSNP RefSNP API에서 주어진 rsID의 chromosome position 리스트 추출
    """
    rs_only = rsid.replace("rs", "")
    url = f"http://grch37.rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {rsid}: {e}")
        return [], []

    data = r.json()
    print(data)
    return data

#get Variant Validator
def get_hgvs_from_variantvalidator(chrom, pos, ref, alt, assembly="hg19"):
    """
    VariantValidator REST API 호출하여 HGVS 정보 반환
    """
    raw = f"{chrom}:{pos}{ref}>{alt}"
    encoded = urllib.parse.quote(raw)
    genome="hg19"
    url = (
        f"https://rest.variantvalidator.org/VariantValidator/variantvalidator/"
        f"{genome}/{encoded}/mane_select?content-type=application/json"
    )

    ##https://rest.variantvalidator.org/VariantValidator/variantvalidator/
    #hg19/chr20%3A62126409C%3ET/mane_select?content-type=application%2Fjson

    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error calling VariantValidator API: {e}")
        return None

    return r.json()
def parse_variantvalidator_result(vv_json):
    # VariantValidator 결과의 최상위 key = "NM_xxxxx:c.xxx"
    top_key = next(k for k in vv_json.keys() if k.startswith("NM_"))
    entry = vv_json[top_key]

    # ─────────────────────────────
    # 기본 정보
    # ─────────────────────────────
    nm_id = top_key.split(":")[0]                         # "NM_001958.5"
    hgvsc = entry.get("hgvs_transcript_variant")          # "NM_001958.5:c.370G>A"

    # HGVSp (단백질) 정보 - slr(short) 우선
    protein_info = entry.get("hgvs_predicted_protein_consequence", {})
    hgvsp = protein_info.get("slr") or protein_info.get("tlr")

    # ─────────────────────────────
    # genomic position 정보: hg19 또는 grch37 기본
    # ─────────────────────────────
    grch37_info = entry["primary_assembly_loci"]["grch37"]  # 또는 entry["primary_assembly_loci"]["hg19"]
    chrom = grch37_info["vcf"]["chr"].replace("chr", "")     # "20"
    pos = int(grch37_info["vcf"]["pos"])                     # 62126409
    ref = grch37_info["vcf"]["ref"]                          # "C"
    alt = grch37_info["vcf"]["alt"]                          # "T"

    # ─────────────────────────────
    # gene info
    # ─────────────────────────────
    gene_symbol = entry["gene_symbol"]                       # "EEF1A2"
    ensembl_gene_id = entry["gene_ids"]["ensembl_gene_id"]   # "ENSG00000101210"

    return {
        "nm_id": nm_id,
        "hgvsc": hgvsc,
        "hgvsp": hgvsp,
        "chromosome": chrom,
        "position": pos,
        "ref": ref,
        "alt": alt,
        "ensembl_gene_id": ensembl_gene_id,
        "gene_symbol": gene_symbol
    }



def rsid_to_hgvs(rsid):
    # 1) dbSNP에서 정보 가져오기
    dbsnp_json = get_chr_pos_from_dbsnp(rsid)
    if dbsnp_json is None:
        return None

    # 2) GRCh38 위치 정보 찾기
    mapping = None
    for m in dbsnp_json.get("mappings", []):
        if m.get("assembly_name") == "GRCh38":
            mapping = m
            break

    if mapping is None:
        print("No GRCh38 mapping found.")
        return None

    chrom = mapping["seq_region_name"]
    pos   = mapping["start"]

    # C/T → REF=C, ALT=T
    allele = mapping.get("allele_string")
    if allele is None:
        print("No allele_string in dbSNP record.")
        return None

    try:
        ref, alt = allele.split("/")
    except ValueError:
        print("Allele parsing error:", allele)
        return None

    # 3) VariantValidator에 질의
    vv = get_hgvs_from_variantvalidator(chrom, pos, ref, alt, assembly="hg19")
    return vv

# 필요한 컬럼 초기화
new_cols = [
    'vv_nm_id', 'vv_hgvsc', 'vv_hgvsp',
    'chromosome', 'position', 'ref', 'alt',
    'ensembl_gene_id', 'gene_symbol'
]

for col in new_cols:
    litvar_df[col] = None

# rsID 목록 순회
for rsid in litvar_df['rsid'].unique():
    print(f"\n[PROCESSING] {rsid}")

    try:
        # 통합 함수 호출: rsid → VariantValidator 결과
        # multi-allele 고려
        dbsnp_json = get_chr_pos_from_dbsnp(rsid)
        if dbsnp_json is None:
            raise ValueError("dbSNP 조회 실패")

        mapping = next((m for m in dbsnp_json.get("mappings", [])
                        if m.get("assembly_name") == "GRCh37"), None)
        if mapping is None:
            raise ValueError("GRCh37 mapping 없음")

        chrom = mapping["seq_region_name"]
        pos = mapping["start"]

        allele_str = mapping.get("allele_string")
        if allele_str is None:
            raise ValueError("allele_string 없음")

        # multi-allele 처리
        alleles = allele_str.split("/")  # ex) G/A/C
        ref = alleles[0]
        alt_list = alleles[1:]           # 나머지는 ALT

        parsed = None
        # 여러 ALT 중 첫 번째 유효한 것만 사용
        for alt in alt_list:
            vv_json = get_hgvs_from_variantvalidator(chrom, pos, ref, alt)
            parsed = parse_variantvalidator_result(vv_json)
            if parsed is not None:
                break

    except Exception as e:
        print(f"[ERROR] {rsid} → {e}")
        parsed = None

    # 동일 rsid 행들 인덱스
    idx_list = litvar_df[litvar_df['rsid'] == rsid].index

    # 파싱 실패한 경우 None 채우기
    if parsed is None:
        for i in idx_list:
            for col in new_cols:
                litvar_df.at[i, col] = None
        continue

    # 파싱 성공 → DF 업데이트
    for i in idx_list:
        litvar_df.at[i, 'vv_nm_id'] = parsed['nm_id']
        litvar_df.at[i, 'vv_hgvsc'] = parsed['hgvsc']
        litvar_df.at[i, 'vv_hgvsp'] = parsed['hgvsp']
        litvar_df.at[i, 'chromosome'] = parsed['chromosome']
        litvar_df.at[i, 'position'] = parsed['position']
        litvar_df.at[i, 'ref'] = parsed['ref']
        litvar_df.at[i, 'alt'] = parsed['alt']
        litvar_df.at[i, 'ensembl_gene_id'] = parsed['ensembl_gene_id']
        litvar_df.at[i, 'gene_symbol'] = parsed['gene_symbol']

    print(f"[SUCCESS] {rsid} 처리 완료")

#normalize column
split_vals = litvar_df["vv_hgvsc"].str.split(":", expand=True)
split_hgvsp = litvar_df["vv_hgvsp"].str.split(":", expand=True)
# print(split_hgvsp)
litvar_df["norm_hgvsc"] = split_vals[1].where(split_vals[1].str.startswith("c."))
litvar_df["norm_hgvsp"] = (
    split_hgvsp[1]
    .where(split_hgvsp[1].str.startswith("p.", na=False))
    .str.replace(r"\(|\)", "", regex=True)
)


path= r"C:\python\Litvar\output"
output_file2=path+"\litvar_vv_gnomad_hgvs_result.csv"
litvar_df.to_csv(output_file2, index=False)
print("output file : ",output_file2)