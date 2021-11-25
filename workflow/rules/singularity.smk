rule get_container:
    """ generic rule for pulling container into resources, so it can be used""" 
    """with singularity run instead of singularity exec"""
    params:
        uri = lambda wildcards: config['containers'][wildcards.container]
    output: 'resources/singularity/{container}.sif'
    shell: 'singularity pull {output} {params.uri}'
 
