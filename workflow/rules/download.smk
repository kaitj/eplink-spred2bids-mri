import xnat
import pandas as pd

envvars:
    'SPRED_USER',
    'SPRED_PASS'

wildcard_constraints:
    site='[a-zA-Z]{3}'

def get_subjects(site):
    if os.path.exists(f'resources/subjects_{site}.txt'):
        with open(f'resources/subjects_{site}.txt','r') as f:
            return [s.replace('\n','') for s in f.readlines()]

localrules: get_subjects_table, download_mri_zip


        

# create table with scans for each subject

# then have rules for each image type


rule get_subjects_table:
    output:
        tsv='resources/subjects.tsv'
    script: 
        '../scripts/get_subjects_table.py'

rule download_zip:
    input:
        tsv='resources/subjects.tsv'
    output:
        zip_file='dicom_zips/site-{site}_subject-{subject}_ses-{session}_{suffix}.zip'
    script:
        '../scripts/download_zip.py'

def get_zips():
    df = pd.read_table('resources/subjects.tsv',dtype={'subject':str})
    zips=[]
    for site in config['sites']:
        for session in config['session_lut'].keys():
            for suffix in config['suffix_lut'].keys():
                #get subjects:
                subjects = df.query("site==@site and session==@session").subject.to_list()
                zips.extend(expand('dicom_zips/site-{site}_subject-{subject}_ses-{session}_{suffix}.zip',site=site,session=session,subject=subjects, suffix=suffix))
    return zips


rule all_zips:
    input:
        zips_dir=get_zips()


rule convert_t1:
    input:
        zip_file=directory('dicom_zips/site-{site}_subject-{subject}_ses-{session}_T1w.zip')
    params:
        tmpdir=lambda wildcards, resources: os.path.join(resources.tmpdir,'_'.join(wildcards),'T1w'),
        outdir='bids/sub-{site}{subject}/ses-{session}/anat/'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_T1w.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/anat/sub-{site}{subject}_T1w.json'
    shadow: 'minimal'
    shell:
        'mkdir -p {params.tmpdir} && '
        'unzip -d {params.tmpdir} {input.zip_file} && '
        'dcm2niix -d 9 -z y -f output {params.tmpdir} && '
        'cp {params.tmpdir}/output.nii.gz {output.nii} && '
        'cp {params.tmpdir}/output.json {output.json} && '
        'rm -rf {params.tmpdir}'

rule convert_dwi:
    input:
        zip_file=directory('dicom_zips/site-{site}_subject-{subject}_ses-{session}_{dwiacq}.zip')
    params:
        tmpdir=lambda wildcards, resources: os.path.join(resources.tmpdir,'_'.join(wildcards),'dwi'),
        outdir='bids/sub-{site}{subject}/ses-{session}/dwi/'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.json',
        bval='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.bval',
        bvec='bids/sub-{site}{subject}/ses-{session}/dwi/sub-{site}{subject}_acq-{dwiacq}_dwi.bvec'
    shadow: 'minimal'
    shell:
        'mkdir -p {params.tmpdir} && '
        'unzip -d {params.tmpdir} {input.zip_file} && '
        'dcm2niix -d 9 -z y -f output {params.tmpdir} && '
        'cp {params.tmpdir}/output.nii.gz {output.nii} && '
        'cp {params.tmpdir}/output.json {output.json} && '
        'cp {params.tmpdir}/output.bvec {output.bvec} && '
        'cp {params.tmpdir}/output.bval {output.bval} && '
        'rm -rf {params.tmpdir}'



rule convert_rest:
    input:
        zip_file=directory('dicom_zips/site-{site}_subject-{subject}_ses-{session}_rest.zip')
    params:
        tmpdir=lambda wildcards, resources: os.path.join(resources.tmpdir,'_'.join(wildcards),'rest'),
        outdir='bids/sub-{site}{subject}/ses-{session}/func/'
    output:
        nii='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_task-rest_bold.nii.gz',
        json='bids/sub-{site}{subject}/ses-{session}/func/sub-{site}{subject}_task-rest_bold.json'
    shadow: 'minimal'
    shell:
        'mkdir -p {params.tmpdir} && '
        'unzip -d {params.tmpdir} {input.zip_file} && '
        'dcm2niix -d 9 -z y -f output {params.tmpdir} && '
        'cp {params.tmpdir}/output.nii.gz {output.nii} && '
        'cp {params.tmpdir}/output.json {output.json} && '
        'rm -rf {params.tmpdir}'


rule make_dicom_tar:
    input:
        zipfile = 'raw/site-{site}/sub-{subject}/mri.zip'
    params:
        file_match = '*/scans/*/resources/DICOM/files/*',
        temp_dir = os.path.join(config['tmp_download'],'raw/site-{site}/sub-{subject}/mri_unzip')
    output:
        tar = 'raw/site-{site}/sub-{subject}/mri/sub-{subject}.tar'
    group: 'dl'
    shell:
        'mkdir -p {params.temp_dir} && '
        'unzip -d {params.temp_dir} {input.zipfile} {params.file_match} && '
        'tar -cvf {output.tar} {params.temp_dir} && '
        'rm -rf {params.temp_dir}'

       
rule tar_to_bids:
    input:
        tar = 'raw/site-{site}/sub-{subject}/mri/sub-{subject}.tar',
        heuristic = lambda wildcards: config['tar2bids'][wildcards.site],
        container = 'resources/singularity/tar2bids.sif'
    params:
        temp_bids_dir = 'raw/site-{site}/sub-{subject}/mri/temp_bids',
        heudiconv_tmpdir = os.path.join(config['tmp_download'],'{site}','{subject}')
    output:
        dir = directory('bids/site-{site}/sub-{subject}')
    group: 'dl'
    shell: 
        "mkdir -p {params.heudiconv_tmpdir} && "
        "singularity run -e {input.container} -o {params.temp_bids_dir} "
        " -w {params.heudiconv_tmpdir} -h {input.heuristic} -N 1 -T 'sub-{{subject}}' {input.tar} || true && " #force clean exit for tar2bids
        "mkdir -p {output.dir} && rsync -av {params.temp_bids_dir}/sub-{wildcards.subject}/ {output.dir} && "
        "rm -rf {params.heudiconv_tmpdir}"


rule create_dataset_json:
    input:
        dir = lambda wildcards: expand('bids/site-{site}/sub-{subject}',site=wildcards.site,subject=get_subjects(wildcards.site)),
        json = 'resources/dataset_description_template.json'
    output:
        json = 'bids/site-{site}/dataset_description.json'
    shell: 'cp {input.json} {output.json}'



