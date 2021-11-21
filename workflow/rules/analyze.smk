import snakebids
import json
import pandas as pd
from snakebids import bids
from os.path import join




sites = config['bids_dir'].keys()


#bids parsing for fmriprep 
if os.path.exists('snakebids.json'):

    #read it in:
    with open('snakebids.json') as f:
        pybids = json.load(f)

else:
    pybids = { site: dict() for site in sites }

    for site in sites:
        pybids[site].update(
            snakebids.generate_inputs(
                bids_dir=config["bids_dir"][site],
                pybids_inputs=config["pybids_inputs"],
                derivatives=config["fmriprep_dir"][site],
            )
        )
        #this adds constraints to the bids naming
        wildcard_constraints:  **snakebids.get_wildcard_constraints(\
            config["pybids_inputs"]\
        )


        # write pybids parsing output to json for easier retrieval later..
        with open('snakebids.json', "w") as f:
            # write either as JSON or YAML
            import json
            json.dump(pybids, f, indent=4)

          
wildcard_constraints: 
    site = '[a-zA-Z0-9]+',
    denoise = '[a-zA-Z0-9]+',
    fwhm = '[a-zA-Z0-9]+'


 
def get_targets():
    targets = list()
    for site in sites:
        targets.extend(\
            expand(bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',atlas='{atlas}',suffix='conn.txt'),\
                subject=pybids[site]['input_lists']['preproc_bold']['subject'],
                site=site,
                task='rest',
                denoise=config['denoise'].keys(),
                fwhm='5',
                space='MNI152NLin6Asym',
                atlas='schaefer')
        )
    return targets
       
        
rule all_analyze:
    input: get_targets()
                


rule gen_aparc_stats_table:
    input: 
        fs_dir = lambda wildcards: config['fs_dir'][wildcards.site],
        subjlist = lambda wildcards: config['subject_list'][wildcards.site]
    output: 
        tsv = 'results/site-{site}_hemi-{hemi}_desc-{meas}_stats.tsv'
    shell: 
        'SUBJECTS_DIR={input.fs_dir} aparcstats2table -m {wildcards.meas}  --subjects `cat {input.subjlist}` --hemi {wildcards.hemi} -t {output}'

rule gen_aseg_stats_table:
    input: 
        fs_dir = lambda wildcards: config['fs_dir'][wildcards.site],
        subjlist = lambda wildcards: config['subject_list'][wildcards.site]
    output: 
        tsv = 'results/site-{site}_desc-vol_asegstats.tsv'
    shell: 
        'SUBJECTS_DIR={input.fs_dir} asegstats2table  --subjects `cat {input.subjlist}`  -t {output}'

rule add_spike_regressors:
    """ this isn't totally necessary, as fmriprep gives motion_outliers_XX with fd_th=0.5, dvars_th=1.5"""
    input:
        confounds_tsv = lambda wildcards: pybids[wildcards.site]['input_path']['confounds']
    params:
        spike_fd_th = 0.5,
        spike_dvars_th = 3 
    output:
        confounds_tsv = bids(root='results',subject='{subject}',site='{site}',task='{task}',desc='confoundswithspikes',suffix='timeseries.tsv')
    script: 'scripts/add_spike_regressors.py'




rule denoise:
    input: 
        nii = lambda wildcards: pybids[wildcards.site]['input_path']['preproc_bold'],
        json = lambda wildcards: re.sub('.nii.gz','.json',pybids[wildcards.site]['input_path']['preproc_bold']),
        confounds_tsv = lambda wildcards: pybids[wildcards.site]['input_path']['confounds'],
        mask_nii = lambda wildcards: pybids[wildcards.site]['input_path']['preproc_mask']
    params:
        denoise_params = lambda wildcards: config['denoise'][wildcards.denoise],
    output: 
        nii = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',space='{space}',suffix='bold.nii.gz'),
        json = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',space='{space}',suffix='bold.json')
    group: 'subj'
    script: 'scripts/denoise.py'

  
rule smooth:
    input:
        nii = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',suffix='bold.nii.gz'),
        json = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',suffix='bold.json')
    params:
        fwhm = lambda wildcards: float(wildcards.fwhm)
    output:
        nii = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',fwhm='{fwhm}',suffix='bold.nii.gz'),
        json = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',fwhm='{fwhm}',suffix='bold.json')
    group: 'subj'
    script: 'scripts/smooth.py'


rule schaefer_connectivity:
    input:
        nii = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',fwhm='{fwhm}',suffix='bold.nii.gz'),
    params:
        n_rois = 300,
        yeo_networks = 7,
        data_dir = 'resources',
        conn_measure = 'correlation'
    output:
        txt = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',fwhm='{fwhm}',atlas='schaefer',suffix='conn.txt'),
        png = bids(root='results',subject='{subject}',site='{site}',task='{task}',denoise='{denoise}',fwhm='{fwhm}',atlas='schaefer',suffix='conn.png'),
    script: 'scripts/connectivity_matrix.py'
     


