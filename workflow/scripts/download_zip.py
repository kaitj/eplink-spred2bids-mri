import pandas as pd
import xnat
import os
df = pd.read_table(snakemake.input.tsv,dtype={'subject':str}).query(snakemake.params.query)
selected_uris = df['scan_uri'].to_list()
uri=selected_uris[0]

with xnat.connect(snakemake.config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS']) as xnat_connection:
    scanobj = xnat_connection.create_object(uri)
    scanobj.download(snakemake.output.zip_file)

