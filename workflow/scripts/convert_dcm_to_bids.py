from snakemake.shell import shell
from glob import glob
import tempfile
import shutil
from pathlib import Path

#get folder name
out_dir=Path(snakemake.output.nii).parent

#get filename without folder or .nii.gz
out_name=str(Path(snakemake.output.nii).name).split('.')[0] 

#get list of extensions expected
out_exts=[''.join(Path(f).suffixes) for f in snakemake.output]

with tempfile.TemporaryDirectory() as tmp_dir:
 
    #unzip to tmp_dir
    shell('unzip -d {tmp_dir} {snakemake.input.zip_file}')

    #run dcm2niix, outputting to the tmp_dir root:
    shell('./deps/dcm2niix -d 9 -z y -f output {tmp_dir} 2> {snakemake.log} ')

    #this nominally creates an output.nii.gz file (along with .json), 
    # but dcm2niix may add an unknown suffix to the output file too, so we need to 
    # glob, grab first nii.gz, then get the basename

    name_to_copy = str(Path(glob(f'{tmp_dir}/output*.nii.gz')[0]).name).split('.')[0] 

    #now, copy each {tmp_dir}/{out_name}{ext} to {out_dir}/{out_name}{out_ext}

    out_dir.mkdir(exist_ok=True)
    for ext in out_exts:
        shutil.copy(Path(tmp_dir) / (name_to_copy + ext), out_dir / (out_name + ext))

