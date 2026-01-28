# This code was used to download the GWAS summary of imaging-derived phenotypes in East Asian participants
# Writted by Yang Xiao, PKU, 2025
# @: xiaoyang9604@gmail.com

import pandas as pd
import os
import urllib.request
from urllib.error import URLError

excel_path = 'IDPGWAS_used.xlsx' 
output_dir = '' 

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

df = pd.read_excel(excel_path, header=None)

raw_urls = df.iloc[:, 5].dropna().tolist()

total = len(raw_urls)
for i, url in enumerate(raw_urls):
    base = str(url).strip().rstrip('/')
    gcst_id = base.split('/')[-1]
    full_url = f"{base}/harmonised/{gcst_id}.h.tsv.gz"
    file_name = f"{gcst_id}.h.tsv.gz"
    save_path = os.path.join(output_dir, file_name)
    if os.path.exists(save_path):
        print(f"[{i+1}/{total}] skipping: {file_name}")
        continue
    print(f"[{i+1}/{total}] downloading: {file_name} ...")
    try:
        urllib.request.urlretrieve(full_url, save_path)
    except URLError as e:
        print(f"fail: {full_url}, error: {e.reason}")
    except Exception as e:
        print(f"fail: {full_url}, error: {e}")

print(f"\ndone: saved in {os.path.abspath(output_dir)}")