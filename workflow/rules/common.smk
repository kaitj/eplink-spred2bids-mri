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


