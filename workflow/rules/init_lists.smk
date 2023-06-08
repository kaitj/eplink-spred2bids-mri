import pandas as pd

rule get_subject_by_suffix_table:
    input:
        tsv='resources/subjects.tsv'
    params:
        query = lambda wildcards: config['suffix_lut'][wildcards.suffix]['query']
    output:
        tsv='resources/subjects_{suffix}.tsv'
    script: 
        '../scripts/get_subject_by_suffix_table.py'



rule get_subjects_table:
    output:
        tsv='resources/subjects.tsv'
    script: 
        '../scripts/get_subjects_table.py'


