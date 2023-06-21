import nibabel as nib

def create_bval_bvec_files(input_file, output_bval, output_bvec):
    img = nib.load(input_file)
    if len(img.shape) == 4:
        N=img.shape[-1]
    else:
        N=1
    print(N)
    with open(output_bval, 'w') as f:
        for _ in range(N):
            f.write('0 ')
        f.write('\n')

    with open(output_bvec, 'w') as f:
        for _ in range(3):
            for _ in range(N):
                f.write('0 ')
            f.write('\n')

# Usage:
create_bval_bvec_files(snakemake.input.nii, snakemake.output.bval, snakemake.output.bvec)

