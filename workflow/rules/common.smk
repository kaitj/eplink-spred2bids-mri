def get_uri(tsv,subject,session, site, suffix, config):

    include_str=config['suffix_lut'][suffix]['include']
    exclude_str=config['suffix_lut'][suffix]['exclude']


    include_pattern = '|'.join(include_str)
    exclude_pattern = '|'.join(exclude_str)

    df = pd.read_table(tsv,dtype={'subject':str})
    filtered_df = df[df['subject'].eq(subject) 
                        & df['session'].eq(session)
                        & df['site'].eq(site)
                        & df['scan_name'].str.contains(include_pattern) 
                        & ~df['scan_name'].str.contains(exclude_pattern)]

    selected_uris = filtered_df['scan_uri'].to_list()

       
    return selected_uris

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
            for subject in df.query("site==@site and session==@session").subject.to_list():
                #check if the subject actually has this scan:
                uri = get_uri(tsv,subject,session, site, suffix, config)
                if len(uri) > 0:
                    in_zip_lists['subject'].append(subject)
                    in_zip_lists['session'].append(session)
                    in_zip_lists['site'].append(site)
    return in_zip_lists

  

