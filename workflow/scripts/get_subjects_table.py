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
                #print(f'exception: {subject} does not have mri')
                continue

            #now get the scans
            scans = exp.scans
            for scan in scans.values():
                scan_uri = scan.uri
#                scan_name = scan.type
                scan_metadata = scan.data
                series_description = scan_metadata.get('series_description','')

                if 'parameters/imageType' in scan_metadata:
                    image_type = scan_metadata['parameters/imageType'].split('\\\\')
                    image_types = {f'image_type_{i}':[image_type[i]] for i in range(len(image_type))}
                else:
                    image_types = {}

#   {'series_description': 'gre_field_mapping', 'scanner/manufacturer': 'SIEMENS', 'image_session_ID': 'spred_E33978', 'type': 'gre_field_mapping', 'xnat_imageScanData_id': 103831, 'parameters/voxelRes/z': 3, 'xnat_imagescandata_id': 103831, 'parameters/voxelRes/x': 3, 'parameters/voxelRes/y': 3, 'scanner': 'MRC35368', 'startTime': '15:36:25', 'parameters/imageType': 'ORIGINAL\\\\PRIMARY\\\\M\\\\ND', 'ID': '17', 'parameters/flip': 60, 'parameters/seqVariant': 'SP', 'parameters/pixelBandwidth': 290, 'frames': 92, 'parameters/fov/y': 80, 'scanner/model': 'Prisma_fit', 'parameters/fov/x': 80, 'parameters/tr': 500, 'quality': 'unknown', 'UID': '1.3.12.2.1107.5.2.43.67007.2019103015345663164066300.0.0.0', 'parameters/scanSequence': 'GR', 'parameters/sequence': '*fm2d2r', 'parameters/orientation': 'Tra', 'fieldStrength': '3.0', 'parameters/acqType': '2D'}
 

                df = pd.concat([df,pd.DataFrame({'subject':[subject_id],'site':[site],'session':[session],**image_types,'series_description':[series_description],'scan_uri':[scan_uri]})],ignore_index=True)

df.to_csv(snakemake.output.tsv,sep='\t',index=False)

xnat_connection.disconnect()


