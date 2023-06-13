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

rule convert_t1:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_T1w.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/anat/sub-{site}{subject}_T1w.nii.gz',
        json='bids_{session}/sub-{site}{subject}/anat/sub-{site}{subject}_T1w.json'
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_T1w.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_dwi_multishell_AP:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwimultishell{dir}.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-{dir,AP}_dwi.nii.gz',
        json='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-{dir,AP}_dwi.json',
        bval='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-{dir,AP}_dwi.bval',
        bvec='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-{dir,AP}_dwi.bvec'
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_dwimultishell{dir,AP}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'

rule convert_dwi_multishell_PA:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwimultishell{dir}.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-{dir,PA}_dwi.nii.gz',
        json='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-{dir,PA}_dwi.json',
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_dwimultishell{dir,PA}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'



rule convert_dwi_singleshell_AP:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwisingleshell{dir}.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-{dir,AP}_dwi.nii.gz',
        json='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-{dir,AP}_dwi.json',
        bval='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-{dir,AP}_dwi.bval',
        bvec='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-{dir,AP}_dwi.bvec'
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_dwisingleshell{dir,AP}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'

rule convert_dwi_singleshell_PA:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwisingleshell{dir}.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-{dir,PA}_dwi.nii.gz',
        json='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-{dir,PA}_dwi.json',
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_dwisingleshell{dir,PA}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'



rule convert_dwi_pepolar:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwiPepolar.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-pepolar_dir-PA_dwi.nii.gz',
        json='bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-pepolar_dir-PA_dwi.json',
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_dwiPepolar.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_rest:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_rest.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/func/sub-{site}{subject}_task-rest_bold.nii.gz',
        json='bids_{session}/sub-{site}{subject}/func/sub-{site}{subject}_task-rest_bold.json'
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_rest.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'

rule convert_movie:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_movie.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/func/sub-{site}{subject}_task-movie_bold.nii.gz',
        json='bids_{session}/sub-{site}{subject}/func/sub-{site}{subject}_task-movie_bold.json'
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_movie.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_fmapMagImages:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_fmapMagImages.zip'
    output:
        nii=['bids_{session}/sub-{site,LHS|TWH}{subject}/fmap/sub-{site}{subject}_magnitude1.nii.gz',
                'bids_{session}/sub-{site,LHS|TWH}{subject}/fmap/sub-{site}{subject}_magnitude2.nii.gz'],
        json=['bids_{session}/sub-{site,LHS|TWH}{subject}/fmap/sub-{site}{subject}_magnitude1.json',
                'bids_{session}/sub-{site,LHS|TWH}{subject}/fmap/sub-{site}{subject}_magnitude2.json']
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site,LHS|TWH}{subject}_fmapMagImages.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_fmapPhaseDiff:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_fmapPhaseDiff.zip'
    output:
        nii='bids_{session}/sub-{site}{subject}/fmap/sub-{site}{subject}_phasediff.nii.gz',
        json='bids_{session}/sub-{site}{subject}/fmap/sub-{site}{subject}_phasediff.json',
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site}{subject}_fmapPhaseDiff.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_fmapTwoPhase:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_fmapTwoPhase.zip'
    output:
        nii=['bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_magnitude1.nii.gz',
                'bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_phase1.nii.gz',
                'bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_magnitude2.nii.gz',
                'bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_phase2.nii.gz'],
        json=['bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_magnitude1.json',
                'bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_phase1.json',
                'bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_magnitude2.json',
                'bids_{session}/sub-{site,HSC}{subject}/fmap/sub-{site}{subject}_phase2.json']
    log: 'logs_{session}/convert_dcm_to_bids/sub-{site,HSC}{subject}_fmapTwoPhase.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'



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
        t1=expand('bids_{session}/sub-{site}{subject}/anat/sub-{site}{subject}_T1w.nii.gz',zip,**zip_lists['T1w']),
        rest=expand('bids_{session}/sub-{site}{subject}/func/sub-{site}{subject}_task-rest_bold.nii.gz',zip,**zip_lists['rest']),
        movie=expand('bids_{session}/sub-{site}{subject}/func/sub-{site}{subject}_task-movie_bold.nii.gz',zip,**zip_lists['movie']),
        dwi_multishell_ap=expand('bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-AP_dwi.bvec',zip,**zip_lists['dwimultishellAP']),
        dwi_multishell_pa=expand('bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-multishell_dir-PA_dwi.bvec',zip,**zip_lists['dwimultishellPA']),
        dwi_singleshell_ap=expand('bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-AP_dwi.bvec',zip,**zip_lists['dwisingleshellAP']),
        dwi_singleshell_pa=expand('bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-singleshell_dir-PA_dwi.bvec',zip,**zip_lists['dwisingleshellPA']),
        dwi_pepolar=expand('bids_{session}/sub-{site}{subject}/dwi/sub-{site}{subject}_acq-pepolar_dir-PA_dwi.bvec',zip,**zip_lists['dwiPepolar']),
        fmap_magimages=expand('bids_{session}/sub-{site}{subject}/fmap/sub-{site}{subject}_magnitude1.nii.gz',zip,**zip_lists['fmapMagImages']),
        fmap_phdiff=expand('bids_{session}/sub-{site}{subject}/fmap/sub-{site}{subject}_phasediff.nii.gz',zip,**zip_lists['fmapPhaseDiff']),
        fmap_twophase=expand('bids_{session}/sub-{site}{subject}/fmap/sub-{site}{subject}_magnitude1.nii.gz',zip,**zip_lists['fmapTwoPhase']),
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
 



