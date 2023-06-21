import pandas as pd

df = pd.read_table(snakemake.input.tsv,dtype={'subject':str})
df = df.query(snakemake.params.query)
df.sort_values(by=['site','subject','session']).to_csv(snakemake.output.tsv,sep='\t',index=False)


