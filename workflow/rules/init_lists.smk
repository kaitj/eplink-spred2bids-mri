import pandas as pd


rule get_subject_by_suffix_table:
    input:
        tsv='resources/subjects.tsv'
    params:
        query = lambda wildcards: config['suffix_lut'][wildcards.suffix]['query']
    output:
        tsv='resources/subjects_{suffix}.tsv'
    script: 
        '../scripts/get_subject_by_suffix_table.py'



rule get_subjects_table:
    output:
        tsv='resources/subjects.tsv'
    script: 
        '../scripts/get_subjects_table.py'

rule get_bad_subjects_from_log:
    input:
        tsv='resources/subjects_{suffix}.tsv',
        log_dir='logs/convert_dcm_to_bids' #make sure all zero-sized files gone..
    output:
        tsv='resources/failed_subjects_{suffix}.tsv'
    run:
        from glob import glob
        from pathlib import Path
        bad_logs = [Path(bad_log).name for bad_log in sorted(glob(f'{input.log_dir}/sub-*_ses-*_{wildcards.suffix}.txt'))]
        sessions = [bad_log.split('_')[1].split('-')[1] for bad_log in bad_logs]
        sites = [bad_log.split('_')[0].split('-')[1][:3] for bad_log in bad_logs]
        subjects = [bad_log.split('_')[0].split('-')[1][3:] for bad_log in bad_logs]
        pd.DataFrame({'site':sites,'subject':subjects,'session':sessions}).to_csv(output.tsv,sep='\t',index=False)
            

rule get_nofailed_subjects_list:
    input:
        failed_tsv='resources/failed_subjects_{suffix}.tsv',
        all_tsv='resources/subjects_{suffix}.tsv',
    output:
        tsv='resources/nofailed_subjects_{suffix}.tsv'
    run:
        df_all=pd.read_table(input.all_tsv,dtype={'subject':str})
        df_failed=pd.read_table(input.failed_tsv,dtype={'subject':str})
        
        # Merge on the three columns with an indicator
        merged_df = pd.merge(df_all, df_failed, on=['subject', 'site', 'session'], how='outer', indicator=True)

        # Filter rows where the merge is only present in the left dataframe (df_all)
        df1_filtered = merged_df[merged_df['_merge'] == 'left_only']

        # Remove the indicator column
        df1_filtered = df1_filtered.drop(columns=['_merge'])
        df1_filtered.to_csv(output.tsv,sep='\t',index=False)

