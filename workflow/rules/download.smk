import pandas as pd


envvars:
    'SPRED_USER',
    'SPRED_PASS'

wildcard_constraints:
    site='[a-zA-Z]{3}'

localrules:  download_zip



#this requires that the subject suffix zip list tables are created (all_zip_lists)
zip_lists={}
for suffix in config['suffix_lut'].keys():
    zip_lists[suffix] = pd.read_table(f'resources/nofailed_subjects_{suffix}.tsv',dtype={'subject':str})[['subject','session','site']].to_dict(orient='list')


rule download_zip:
    input:
        tsv='resources/subjects_{suffix}.tsv'
    params:
        query="subject=='{subject}' and session=='{session}' and site=='{site}'"
    output:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_{suffix}.zip'
    threads: 32
    script:
        '../scripts/download_zip.py'


for suffix in config['suffix_lut']:
    rule:
        name: f'convert_dcm_to_bids_{suffix}'
        input:
            zip_file=f'dicom_zips/site-{{site}}_subject-{{subject}}_ses-{{session}}_{suffix}.zip'
        output:
            **config['suffix_lut'][suffix]['outputs']
        log: f'logs/convert_dcm_to_bids/sub-{{site}}{{subject}}_ses-{{session}}_{suffix}.txt'
        group: 'convert'
        script: 
            '../scripts/convert_dcm_to_bids.py'

#this hack is needed because the same filenames are required for two diff phase maps (HSC uses twophase, others use phasediff)
ruleorder:  convert_dcm_to_bids_fmapMagImages > convert_dcm_to_bids_fmapTwoPhase 


rule create_bval_bvec_pepolar:
    input:
        nii='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-{acq}_dir-PA_dwi.nii.gz',
    output:
        bvec='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-{acq}_dir-PA_dwi.bvec',
        bval='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-{acq}_dir-PA_dwi.bval',
    group: 'convert'
    script: '../scripts/create_bval_bvec_pepolar.py'

  
rule create_dataset_json:
    input:
        **{ suffix:expand(config['suffix_lut'][suffix]['outputs']['nii'],zip,**zip_lists[suffix]) for suffix in config['suffix_lut'].keys()},
        dd_json = 'resources/dataset_description_template.json',
        bidsignore = 'resources/bids_root_files/bidsignore',
        rest_json = 'resources/bids_root_files/task-rest_bold.json',
        movie_json = 'resources/bids_root_files/task-movie_bold.json',
    output:
        dd_jsons = expand('bids_{session}/dataset_description.json',session=config['session_lut'].keys()),
        bidsignores= expand('bids_{session}/.bidsignore',session=config['session_lut'].keys()),
        rest_jsons= expand('bids_{session}/task-rest_bold.json',session=config['session_lut'].keys()),
        movie_jsons= expand('bids_{session}/task-movie_bold.json',session=config['session_lut'].keys())

    run: 
        for out_dd_json in output.dd_jsons:
            shell('cp {input.dd_json} {out_dd_json}')
        for out_rest_json in output.rest_jsons:
            shell('cp {input.rest_json} {out_rest_json}')
        for out_movie_json in output.movie_jsons:
            shell('cp {input.movie_json} {out_movie_json}')
        for out_bidsignore in output.bidsignores:
            shell('cp {input.bidsignore} {out_bidsignore}')
 



