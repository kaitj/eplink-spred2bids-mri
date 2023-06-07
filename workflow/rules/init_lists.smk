import pandas as pd

rule get_subject_by_suffix_table:
    input:
        tsv='resources/subjects.tsv'
    output:
        tsv='resources/subjects_{suffix}.tsv'
    run:
        in_zip_lists = get_zip_lists_bysuffix(input.tsv,wildcards.suffix)
        pd.DataFrame(in_zip_lists).sort_values(by=['site','subject','session']).to_csv(output.tsv,sep='\t',index=False)





def get_zip_lists(tsv):
    df = pd.read_table(tsv,dtype={'subject':str})
    in_zip_lists={'site':[],'subject':[],'session':[],'suffix':[]}

    for site in config['sites']:
        for session in config['session_lut'].keys():
            for suffix in config['suffix_lut'].keys():
                #get subjects:
                for subject in df.query("site==@site and session==@session").subject.to_list():
                    #check if the subject actually has this scan:
                    uri = get_uri(tsv,subject,session, site, suffix, config)
                    if len(uri) > 0:
                        in_zip_lists['subject'].append(subject)
                        in_zip_lists['session'].append(session)
                        in_zip_lists['site'].append(site)
                        in_zip_lists['suffix'].append(suffix)
    return in_zip_lists

def get_zip_lists_bysuffix(tsv,suffix):
    df = pd.read_table(tsv,dtype={'subject':str})
    in_zip_lists={'site':[],'subject':[],'session':[]}

    for site in config['sites']:
        for session in config['session_lut'].keys():
            #get subjects:
            for subject in list(set(df.query("site==@site and session==@session").subject.to_list())):
                #check if the subject actually has this scan:
                uri = get_uri(tsv,subject,session, site, suffix, config)
                if len(uri) > 0:
                    in_zip_lists['subject'].append(subject)
                    in_zip_lists['session'].append(session)
                    in_zip_lists['site'].append(site)
    return in_zip_lists

  

rule get_subjects_table:
    output:
        tsv='resources/subjects.tsv'
    script: 
        '../scripts/get_subjects_table.py'


