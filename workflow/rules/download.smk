import xnat

envvars:
    'SPRED_USER',
    'SPRED_PASS'


rule get_subject_list:
    params:
        site_id = 'EPL31_{site}'
    output:
        subj_list = 'resources/subjects_{site}.txt'
    run:
        session = xnat.connect(config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'])
        subjects = [row[0] for row in session.projects[params.site_id].subjects.tabulate(columns=['label'])]
        with open(output.subj_list, "w") as out:
            for s in subjects:
                out.write(s.split('_')[2]+'\n') #pull out EPL31_LHS_XXXX
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
        temp_dir = 'raw/site-{site}/sub-{subject}/mri_unzip'
    output:
        tar = 'raw/site-{site}/sub-{subject}/mri/sub-{subject}.tar'
    shell:
        'unzip -d {params.temp_dir} {input.zipfile} {params.file_match} && '
        'tar -cvf {output.tar} {params.temp_dir}  && rm -rf {params.temp_dir}'

       
rule tar_to_bids:
    input:
        tar = 'raw/site-{site}/sub-{subject}/mri/sub-{subject}.tar',
        heuristic = lambda wildcards: config['tar2bids'][wildcards.site],
        container = '../resources/singularity/tar2bids.sif'
    params:
        temp_bids_dir = 'raw/site-{site}/sub-{subject}/mri/temp_bids'
    output:
        dir = directory('bids/site-{site}/sub-{subject}')
    shell: 
        "singularity run -e {input.container} -o {params.temp_bids_dir} "
        "  -h {input.heuristic} -N 1 -T 'sub-{{subject}}' {input.tar} || true && " #force clean exit for tar2bids
        "mkdir -p {output.dir} && rsync -av {params.temp_bids_dir}/sub-{wildcards.subject}/ {output.dir}"


