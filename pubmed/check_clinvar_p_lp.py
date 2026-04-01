#!/bin/python

import pandas as pd
import os

target_gene = "EEF1A2"
os.system("mkdir clinvar_results")
# os.system("python /content/drive/MyDrive/BMI_LAB/연구/Clinvar 뇌전증 연구(임병찬 교수님)/script/llm_epilepsy_version_2/main.py EEF1A2")
#clinvar_result.txt는 web에서 다운로드 받은결과.
df=pd.read_csv('/Users/snuh4/python/project/ep_llm/clinvar/clinvar_result.txt',sep='\t')
df

# transcript: NM_ 으로 시작해서 ( 전까지
df["transcript"] = df["Name"].str.extract(r'(NM_[^\(]+)')
# hgvsc: : 뒤에 c. 로 시작하는 부분
df["hgvsc"] = df["Name"].str.extract(r':\s*(c\.[^\s]+)')
# hgvsp: (p. ... ) 괄호 안
df["hgvsp"] = df["Name"].str.extract(r'\((p\.[^)]+)\)')
df["Gene"] = df["Name"].str.extract(r"\(([^)]+)\):")


clinvar_df=df[["Gene","transcript", "hgvsc", "hgvsp", "Germline classification"]]
# print(clinvar_df.head)

########
##match clinvar with pubmed gpt result
########

#get pubmed gpt result
p_df=pd.read_csv("/Users/snuh4/python/project/ep_llm/Litvar/output/pubmed_hgvsp_exist_EEF1A2_phenotypes.csv")
print(p_df.columns)
# print(p_df.shape)


#p_df에서 hgvsp와 clinvar_df의 hgvsp가 일치하는 행 필터링
merged_df = pd.merge(clinvar_df, p_df, left_on="hgvsp", right_on="variant_p", how="inner")
print(merged_df.shape)
print(merged_df[["hgvsc", "hgvsp", "variant_p", "Germline classification"]])
merged_df.to_csv("./merged_df.csv", index=False, encoding="utf-8-sig")