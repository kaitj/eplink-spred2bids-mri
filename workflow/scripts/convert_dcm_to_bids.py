from snakemake.shell import shell
from glob import glob
import tempfile
import shutil
from pathlib import Path

if isinstance(snakemake.output.nii,str):
    nii_list=[snakemake.output.nii]
else:
    nii_list=list(snakemake.output.nii)

#get list of extensions expected
out_exts=list(set([''.join(Path(f).suffixes) for f in snakemake.output]))
print(out_exts)
print(nii_list)


with tempfile.TemporaryDirectory() as tmp_dir:
 
    #unzip to tmp_dir
    shell('unzip -d {tmp_dir} {snakemake.input.zip_file}')

    #run dcm2niix, outputting to the tmp_dir root:
    shell('./deps/dcm2niix -d 9 -z y -f output {tmp_dir} 2> {snakemake.log} ')

    #this nominally creates an output.nii.gz file (along with .json), 
    # but dcm2niix may add an unknown suffix to the output file too, so we need to 
    # glob, grab list of nii.gz files, then get the basename
    base_names = [ str(Path(nii).name).split('.')[0] for nii in sorted(glob(f'{tmp_dir}/output*.nii.gz'))]

    print(base_names)
    print(nii_list)
    #now, copy each {tmp_dir}/{out_name}{ext} to {out_dir}/{out_name}{out_ext}

    for base_name, nii in zip(base_names,nii_list):

        out_dir=Path(nii).parent

        #get filename without folder or .nii.gz
        out_name=str(Path(nii).name).split('.')[0] 


        out_dir.mkdir(exist_ok=True)
        for ext in out_exts:
            src=Path(tmp_dir) / (base_name + ext)
            dest=out_dir / (out_name + ext)
            print(f'copying {src} to {dest}')
            shutil.copy(src,dest)

