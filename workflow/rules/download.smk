import xnat
import pandas as pd


envvars:
    'SPRED_USER',
    'SPRED_PASS'

wildcard_constraints:
    site='[a-zA-Z]{3}'

localrules: get_subjects_table, download_zip

      



rule get_subjects_table:
    output:
        tsv='resources/subjects.tsv'
    script: 
        '../scripts/get_subjects_table.py'

rule download_zip:
    input:
        tsv='resources/subjects.tsv'
    params:
        uri=lambda wildcards, input: get_uri(input.tsv,wildcards.subject,wildcards.session,wildcards.site,wildcards.suffix,config)
    output:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_{suffix}.zip'
    threads: 32
    script:
        '../scripts/download_zip.py'

rule convert_t1:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_T1w.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_T1w.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_T1w.json'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_dwi:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_dwi{dwiacq}.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.json',
        bval='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.bval',
        bvec='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.bvec'
    script: 
        '../scripts/convert_dcm_to_bids.py'


rule convert_rest:
    input:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_rest.zip'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_task-rest_bold.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_task-rest_bold.json'
    script: 
        '../scripts/convert_dcm_to_bids.py'

rule create_dataset_json:
    input:
        t1=expand('bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_T1w.nii.gz',zip,**get_zip_lists_bysuffix('resources/subjects.tsv',suffix='T1w')),
        rest=expand('bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_task-rest_bold.nii.gz',zip,**get_zip_lists_bysuffix('resources/subjects.tsv',suffix='rest')),
        dwi_ap=expand('bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-AP_dwi.nii.gz',zip,**get_zip_lists_bysuffix('resources/subjects.tsv',suffix='dwiAP')),
        dwi_pa=expand('bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-PA_dwi.nii.gz',zip,**get_zip_lists_bysuffix('resources/subjects.tsv',suffix='dwiPA')),

        json = 'resources/dataset_description_template.json'
    output:
        json = 'bids/dataset_description.json'
    shell: 'cp {input.json} {output.json}'


