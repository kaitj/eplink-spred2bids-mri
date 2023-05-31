import xnat
import pandas as pd
import os

xnat_connection = xnat.connect(snakemake.config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'])
project_id = snakemake.config['project_id']

df = pd.DataFrame(columns=['site','subject','session','scan_name','scan_uri'])

for session in snakemake.config['session_lut'].keys():

    mr_session = snakemake.config['session_lut'][session]
    for site in snakemake.config['sites']:

        site_id = f'{project_id}_{site}'
        subjects = [row[0] for row in xnat_connection.projects[site_id].subjects.tabulate(columns=['label'])]

        for subject in subjects:

            subject_id = subject.split('_')[2] #strip off all but numeric part of ID

            try:
                exp_uri=f'/data/projects/{site_id}/experiments/{subject}_{mr_session}_SE01_MR'
                print(exp_uri)
                exp = xnat_connection.create_object(f'/data/projects/{site_id}/experiments/{subject}_{mr_session}_SE01_MR')
            except:
                print(f'exception: {subject} does not have mri')
                continue

            #now get the scans
            scans = exp.scans
            for scan in scans.values():
                scan_uri = scan.uri
                scan_name = scan.type
                df = pd.concat([df,pd.DataFrame({'subject':[subject_id],'site':[site],'session':[session],'scan_name':[scan_name],'scan_uri':[scan_uri]})],ignore_index=True)

df.to_csv(snakemake.output.tsv,sep='\t',index=False)

xnat_connection.disconnect()


