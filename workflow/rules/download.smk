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
        nii='bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_ses-{session}_T1w.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_ses-{session}_T1w.json'
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_T1w.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_dwi_multishell_AP:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwimultishell{dir}.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-{dir,AP}_dwi.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-{dir,AP}_dwi.json',
        bval='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-{dir,AP}_dwi.bval',
        bvec='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-{dir,AP}_dwi.bvec'
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_dwimultishell{dir}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'

rule convert_dwi_multishell_PA:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwimultishell{dir}.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-{dir,PA}_dwi.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-{dir,PA}_dwi.json',
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_dwimultishell{dir}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'



rule convert_dwi_singleshell_AP:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwisingleshell{dir}.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-{dir,AP}_dwi.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-{dir,AP}_dwi.json',
        bval='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-{dir,AP}_dwi.bval',
        bvec='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-{dir,AP}_dwi.bvec'
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_dwisingleshell{dir}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'

rule convert_dwi_singleshell_PA:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwisingleshell{dir}.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-{dir,PA}_dwi.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-{dir,PA}_dwi.json',
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_dwisingleshell{dir}.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'



rule convert_dwi_pepolar:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwiPepolar.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-pepolar_dir-PA_dwi.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-pepolar_dir-PA_dwi.json',
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_dwiPepolar.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_rest:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_rest.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_ses-{session}_task-rest_bold.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_ses-{session}_task-rest_bold.json'
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_rest.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'

rule convert_movie:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_movie.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_ses-{session}_task-movie_bold.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_ses-{session}_task-movie_bold.json'
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_movie.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_fmapMagImages:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_fmapMagImages.zip'
    output:
        nii=['bids/sub-{site,LHS|TWH}{subject}/ses-{session}/fmap/sub-{site,LHS|TWH}{subject}_ses-{session}_magnitude1.nii.gz',
                'bids/sub-{site,LHS|TWH}{subject}/ses-{session}/fmap/sub-{site,LHS|TWH}{subject}_ses-{session}_magnitude2.nii.gz'],
        json=['bids/sub-{site,LHS|TWH}{subject}/ses-{session}/fmap/sub-{site,LHS|TWH}{subject}_ses-{session}_magnitude1.json',
                'bids/sub-{site,LHS|TWH}{subject}/ses-{session}/fmap/sub-{site,LHS|TWH}{subject}_ses-{session}_magnitude2.json']
    log: 'logs/convert_dcm_to_bids/sub-{site,LHS|TWH}{subject}_ses-{session}_fmapMagImages.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_fmapPhaseDiff:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_fmapPhaseDiff.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/fmap/sub-{site}{subject}_ses-{session}_phasediff.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/fmap/sub-{site}{subject}_ses-{session}_phasediff.json',
    log: 'logs/convert_dcm_to_bids/sub-{site}{subject}_ses-{session}_fmapPhaseDiff.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_fmapTwoPhase:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_fmapTwoPhase.zip'
    output:
        nii=['bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_magnitude1.nii.gz',
                'bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_phase1.nii.gz',
                'bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_magnitude2.nii.gz',
                'bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_phase2.nii.gz'],
        json=['bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_magnitude1.json',
                'bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_phase1.json',
                'bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_magnitude2.json',
                'bids/sub-{site,HSC}{subject}/ses-{session}/fmap/sub-{site,HSC}{subject}_ses-{session}_phase2.json']
    log: 'logs/convert_dcm_to_bids/sub-{site,HSC}{subject}_ses-{session}_fmapTwoPhase.txt'
    group: 'convert'
    script: 
        '../scripts/convert_dcm_to_bids.py'



rule create_bval_bvec_pepolar:
    input:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-{acq}_dir-PA_dwi.nii.gz',
    output:
        bvec='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-{acq}_dir-PA_dwi.bvec',
        bval='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-{acq}_dir-PA_dwi.bval',
    group: 'convert'
    script: '../scripts/create_bval_bvec_pepolar.py'
  

rule create_dataset_json:
    input:
        t1=expand('bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_ses-{session}_T1w.nii.gz',zip,**zip_lists['T1w']),
        rest=expand('bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_ses-{session}_task-rest_bold.nii.gz',zip,**zip_lists['rest']),
        movie=expand('bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_ses-{session}_task-movie_bold.nii.gz',zip,**zip_lists['movie']),
        dwi_multishell_ap=expand('bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-AP_dwi.bvec',zip,**zip_lists['dwimultishellAP']),
        dwi_multishell_pa=expand('bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-multishell_dir-PA_dwi.bvec',zip,**zip_lists['dwimultishellPA']),
        dwi_singleshell_ap=expand('bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-AP_dwi.bvec',zip,**zip_lists['dwisingleshellAP']),
        dwi_singleshell_pa=expand('bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-singleshell_dir-PA_dwi.bvec',zip,**zip_lists['dwisingleshellPA']),
        dwi_pepolar=expand('bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_ses-{session}_acq-pepolar_dir-PA_dwi.bvec',zip,**zip_lists['dwiPepolar']),
        fmap_magimages=expand('bids/sub-{site}{subject}/ses-{session}/fmap/sub-{site}{subject}_ses-{session}_magnitude1.nii.gz',zip,**zip_lists['fmapMagImages']),
        fmap_phdiff=expand('bids/sub-{site}{subject}/ses-{session}/fmap/sub-{site}{subject}_ses-{session}_phasediff.nii.gz',zip,**zip_lists['fmapPhaseDiff']),
        fmap_twophase=expand('bids/sub-{site}{subject}/ses-{session}/fmap/sub-{site}{subject}_ses-{session}_magnitude1.nii.gz',zip,**zip_lists['fmapTwoPhase']),
        json = 'resources/dataset_description_template.json'
    output:
        json = 'bids/dataset_description.json'
    shell: 'cp {input.json} {output.json}'


