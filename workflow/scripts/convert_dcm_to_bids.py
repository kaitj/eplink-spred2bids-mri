from snakemake.shell import shell
from glob import glob
import tempfile
import shutil
from pathlib import Path

out_dir=Path(snakemake.output.nii).parent
out_name=str(Path(snakemake.output.nii).name.rstrip('.nii.gz'))
out_exts=[''.join(Path(f).suffixes) for f in snakemake.output]

with tempfile.TemporaryDirectory() as tmpdirname:
 
    #unzip to tmpdir
    shell('unzip -d {tmpdirname} {snakemake.input.zip_file}')

    #run dcm2niix, outputting to the tmpdir root:
    shell('dcm2niix -d 9 -z y -f output {tmpdirname}')

    #this nominally creates an output.nii.gz file (along with .json), 
    # but dcm2niix may add an unknown suffix to the output file too, so we need to glob

    nii_to_copy = glob(f'{tmpdirname}/output*.nii.gz')[0]

    basename_to_copy=nii_to_copy.rstrip('.nii.gz')

    #now, copy each {tmpdirname}/output*.{ext} to {out_dir}/{out_name}{out_ext}
    out_dir.mkdir(exist_ok=True)
    for ext in out_exts:
        shutil.copy(basename_to_copy+ext, out_dir / (out_name + ext))

