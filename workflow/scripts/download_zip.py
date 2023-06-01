import xnat
import os
"""
from helpers import get_uri

uri = get_uri(snakemake.input.tsv,
            snakemake.wildcards.subject,
            snakemake.wildcards.session,
            snakemake.wildcards.site,
            snakemake.wildcards.suffix,
            snakemake.config)[0]
"""
uri = snakemake.params.uri[0]

with xnat.connect(snakemake.config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS']) as xnat_connection:
    scanobj = xnat_connection.create_object(uri)
    scanobj.download(snakemake.output.zip_file)

