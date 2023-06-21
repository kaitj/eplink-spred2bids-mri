# eplink-spred2bids-mri
Snakemake workflow for downloading and converting EpLink MRI data to BIDS. See http://github.com/khanlab/eplink-spred2bids-eeg for EEG workflow.

## Prerequisites:

You need to have BrainCODE SPRED (XNAT) access to the Eplink EPL31 project. 
Put your username and password in environment variables named `SPRED_USER` and `SPRED_PASS` respectively.

## Instructions:

1. Clone this repository and install with `pip install <path_to_cloned_repo>`

2. Update the config file to change the `tmp_download` folder to a local disk with large enough space. 

3. Run snakemake with a dry-run first:
```
snakemake -np
```

4. If everything looks fine, run with the specified number of parallel cores, e.g. 4 as below:
```
snakemake --cores 4
```

 


