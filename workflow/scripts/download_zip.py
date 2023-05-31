import pandas as pd
from pathlib import Path
import xnat
import os

include=snakemake.config['suffix_lut'][snakemake.wildcards.suffix]['include']
exclude=snakemake.config['suffix_lut'][snakemake.wildcards.suffix]['exclude']


include_pattern = '|'.join(include)
exclude_pattern = '|'.join(exclude)

df = pd.read_table(snakemake.input.tsv,dtype={'subject':str})
filtered_df = df[df['subject'].eq(snakemake.wildcards.subject) 
                    & df['session'].eq(snakemake.wildcards.session)
                    & df['site'].eq(snakemake.wildcards.site)
                    & df['scan_name'].str.contains(include_pattern) 
                    & ~df['scan_name'].str.contains(exclude_pattern)]

selected_uris = filtered_df['scan_uri'].to_list()

if len(selected_uris) > 1:
    print('multiple files matched for {snakemake.wildcards}')
    print(filtered_df)
    print('picking first entry')
   
uri = selected_uris[0]

with xnat.connect(snakemake.config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS']) as xnat_connection:
    scanobj = xnat_connection.create_object(uri)
    scanobj.download(snakemake.output.zip_file)

