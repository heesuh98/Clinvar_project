#!/bin/python
import pandas as pd
import requests
from tqdm import tqdm
import math

#path
path= r"C:\python\Litvar\output"
# output_file=path+"\pmid_hgvs.txt"
output_file="/Users/snuh4/python/project/ep_llm/Litvar/output/pmid_hgvs.txt"

GNOMAD_API = "https://gnomad.broadinstitute.org/api"


def clean_value(x):
    """Convert NaN, inf, floats to clean strings."""
    if x is None:
        return ""
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return ""
    return str(x).strip()

def check_dbsnp(rsid):
    if not rsid.startswith("rs"):
        return False
    rs_num = rsid.replace("rs", "")
    if not rs_num.isdigit():
        return False

    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rs_num}"
    try:
        r = requests.get(url, timeout=5)
        return r.status_code == 200
    except:
        return False


def query_gnomad(rsid=None, hgvs=None, ref="GRCh38"):
    """Query gnomAD safely: None filtered, JSON-safe."""
    rsid = rsid if rsid else None
    hgvs = hgvs if hgvs else None

    # If both None, skip
    if rsid is None and hgvs is None:
        return None

    query = """
    query Variant($rsid: String, $hgvs: String, $ref: ReferenceGenomeId!) {
      variant(rsid: $rsid, hgvs: $hgvs, referenceGenome: $ref) {
        variantId
        chrom
        pos
        ref
        alt
        genome { af ac an }
        exome { af ac an }
      }
    }
    """
    variables = {
        "rsid": rsid,
        "hgvs": hgvs,
        "ref": ref
    }

    try:
        r = requests.post(GNOMAD_API, json={"query": query, "variables": variables}, timeout=10)
        return r.json().get("data", {}).get("variant")
    except:
        return None


def vep_translate(hgvs, gene):
    """Convert protein HGVS → canonical HGVS using VEP."""
    hgvs_full = f"{gene}:{hgvs}"  # gene:p.XxxY
    url = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs_full}"
    try:
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
        if r.status_code != 200:
            return None
        return r.json()
    except:
        return None


# -------------------------
# Processing pipeline
# -------------------------

def process_variant(rsid, gene, hgvs):
    rsid = clean_value(rsid)
    gene = clean_value(gene)
    hgvs = clean_value(hgvs)

    dbsnp_exists = False
    gnomad_result = None
    mapped_hgvs = None

    # 1) rsID check
    if rsid.startswith("rs"):
        dbsnp_exists = check_dbsnp(rsid)
        # Try gnomAD with rsID
        gnomad_result = query_gnomad(rsid=rsid)

    # 2) If hgvs exists and no result yet
    if not gnomad_result and hgvs:
        mapped_hgvs = hgvs
        gnomad_result = query_gnomad(hgvs=hgvs)

        # Protein-level HGVS: p.XxxY
        if not gnomad_result and hgvs.startswith("p.") and gene:
            vep_output = vep_translate(hgvs, gene)
            if isinstance(vep_output, list):
                for entry in vep_output:
                    for t in entry.get("transcript_consequences", []):
                        hgvsc = t.get("hgvsc")
                        if hgvsc:
                            mapped_hgvs = hgvsc
                            gnomad_result = query_gnomad(hgvs=hgvsc)
                            if gnomad_result:
                                break

    return {
        "rsid": rsid,
        "gene": gene,
        "input_hgvs": hgvs,
        "mapped_hgvs": mapped_hgvs,
        "dbsnp_exists": dbsnp_exists,
        "gnomad_exists": gnomad_result is not None,
        "gnomad_variantId": gnomad_result["variantId"] if gnomad_result else None,
        "gnomad_af_genome": gnomad_result["genome"]["af"] if (gnomad_result and gnomad_result["genome"]) else None,
        "gnomad_af_exome": gnomad_result["exome"]["af"] if (gnomad_result and gnomad_result["exome"]) else None
    }


# -------------------------
# Run
# -------------------------

INPUT_FILE = output_file
OUTPUT_FILE = path+"\litvar_gnomad_variants.tsv"

df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str)
df = df.fillna("") 

results = []

for _, row in tqdm(df.iterrows(), total=len(df)):
    result = process_variant(
        rsid=row.get("rsid", ""),
        gene=row.get("gene", ""),
        hgvs=row.get("hgvs", "")
    )
    results.append(result)

outdf = pd.DataFrame(results)
outdf.to_csv(OUTPUT_FILE, sep="\t", index=False)

print("\n=== 검증 완료 ===")
print(f"\n결과 저장됨 → {OUTPUT_FILE}")