import pandas as pd
import requests
from io import BytesIO
import gzip
import numpy as np

#   Extract CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz
#
rna_seq_url = "https://storage.googleapis.com/depmap-external-downloads/ccle/ccle_2019/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz?GoogleAccessId=depmap-external-downloads%40broad-achilles.iam.gserviceaccount.com&Expires=1700798657&Signature=KBe3RDYk76KfAHgr5NU9DDnBx%252BPd7od%252BkzAw%252F4B%252BDvmTAvlVO%252BcfWJSygu7gRFOMAEGBWyWLFk0E2kUG94JqOflFV%252B191qOoWE9CA0XGsS3hEfZAzuLZdKuGSVfz5hwDonh%252B3uaO0vBGaawHW1DY1Ar9Ve14eLxln4AJnOFGw87vX7%252FqSca%252Bq7X%252FyRmYqwvVuyMTJ2CvRvhRny%252F3SdJnuOWuVABp9IasTic5iarbu1s3vrB7uMMsaqSF8AgB3z22eKLkBgQGT69NSIpply4c3g5ifRG0cYCu%252FkP9mTuHjKkh%252F4EHLDfm2KlcYa3vAcsbTIwkio1e9BeqOJC%252F0C7qYw%3D%3D&userProject=broad-achilles"
rna_seq_response = requests.get(rna_seq_url)
if rna_seq_response.status_code != 200:
    print(f"HTTP error: {rna_seq_response.status_code}")
rna_seq_data = pd.read_csv(BytesIO(gzip.decompress(rna_seq_response.content)), delimiter='\t')

# Extract Cell_lines_annotations_20181226.txt
metadata_url = "https://storage.googleapis.com/depmap-external-downloads/ccle/ccle_2019/Cell_lines_annotations_20181226.txt?GoogleAccessId=depmap-external-downloads%40broad-achilles.iam.gserviceaccount.com&Expires=1700798813&Signature=Zez9H9VLTbc%252BC1IIxp3mEzQOjDRZGgX3R4LrS7Lhe6DeSRY3l%252FfiUgsZNzpYNOe%252FgV7dOto19RJBKaDoXhjzsyL%252FLOE0fYZTIIPdkXcTxMWicnZT9GpjUGp8nYj5nyHVQPL42SS65jPPQ0PkTMKleDI1XRTqdEB0Fr4EC%252BXp316yuXwzoS7mnJ4bSyawPK3mltv5lYmyq0yreH66wzstcMc481yxNrZ34KLvAwJykA4kwqnLREr6jZzocTmOZtseZzqt5ziPWu9pe16YJI49eGqMT6%252F3P8rZ5o693H8IfumSZ8rYH3nDJM7SIshgUl2DbcT89x6Td69hHqfC3Ti4RA%3D%3D&userProject=broad-achilles"
metadata_response = requests.get(metadata_url)
metadata_data = pd.read_csv(BytesIO(metadata_response.content), delimiter='\t')

# Transform : Load with variable names

rnaseq_tpm = rna_seq_data.copy()
rnaseq_metadata = metadata_data.copy()

# Check for column-wise missing values in rnaseq_metadata and drop columns
missing_values = rnaseq_metadata.isnull().sum()
columns_to_drop = missing_values[missing_values > 700].index
rnaseq_metadata = rnaseq_metadata.drop(columns=columns_to_drop)

# Load rnaseq_tpm
rnaseq_tpm = rnaseq_tpm.drop(columns='transcript_ids')
rnaseq_tpm = rnaseq_tpm.apply(lambda x: np.log2(x + 0.001) if np.issubdtype(x.dtype, np.number) else x)

# Subset rnaseq_metadata based on common cell line names
common_cell_lines = set(rnaseq_metadata['CCLE_ID']).intersection(set(rnaseq_tpm.columns))
rnaseq_metadata = rnaseq_metadata[rnaseq_metadata['CCLE_ID'].isin(common_cell_lines)]

# Check and reorder rnaseq_metadata dataframe
rnaseq_metadata_order = list(rnaseq_tpm.columns)
#rnaseq_metadata = rnaseq_metadata[['CCLE_ID'] + rnaseq_metadata_order]

# Load :Save transformed data into CSV files
rnaseq_tpm.to_csv('rnaseq_tpm_transformed.csv', index=False)
rnaseq_metadata.to_csv('rnaseq_metadata_transformed.csv', index=False)
