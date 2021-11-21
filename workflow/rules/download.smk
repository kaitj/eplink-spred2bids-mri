import xnat

envvars:
    'SPRED_USER',
    'SPRED_PASS'


def get_subjects(site):
    if os.path.exists(f'resources/subjects_{site}.txt'):
        with open(f'resources/subjects_{site}.txt','r') as f:
            return [s.replace('\n','') for s in f.readlines()]

localrules: get_subject_list, download_mri_zip

rule get_subject_list:
    params:
        site_id = 'EPL31_{site}'
    output:
        subj_list = 'resources/subjects_{site}.txt'
    run:
        session = xnat.connect(config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'])
        subjects = [row[0] for row in session.projects[params.site_id].subjects.tabulate(columns=['label'])]
        subjects_with_mr = list()
        for subject in subjects:
            try:
                exp = session.create_object(f'/data/projects/{params.site_id}/experiments/{subject}_01_SE01_MR')
                subjects_with_mr.append(subject.split('_')[2]) #strip off all but numeric part of ID
            except:
                print(f'{subject} does not have mri')

        with open(output.subj_list, "w") as out:
            for s in subjects_with_mr:
                out.write(s+'\n') 
        session.disconnect()


rule download_mri_zip:
    params:
        remote_path = config['remote_path_mri']
    output:
        zipfile = 'raw/site-{site}/sub-{subject}/mri.zip'
    run:
        session = xnat.connect(config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'])
        experiment = session.create_object(params.remote_path)
        experiment.download(output.zipfile)
        session.disconnect()

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
        container = '../resources/singularity/tar2bids.sif'

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


rule make_root_bids_files:
    input:
        dir = lambda wildcards: expand('bids/site-{site}/sub-{subject}',site=wildcards.site,subject=get_subjects(wildcards.site)),
        json = 'resources/dataset_description_template.json'
    output:
        json = 'bids/site-{site}/dataset_description.json'
    shell: 'cp {input.json} {output.json}'
