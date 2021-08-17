# scriptsForAmandioEtAl2021

All scripts necessary to build figures from raw data in Amandio et al. 2021

The analyses rely on a conda environments. To create it use:

```bash
conda env create -f conda_env.yml
```

Each of these directories contain sbatch scripts or command lines for the corresponding analyses: `cHi-C`, `ChIP`, `ChIPmentation`.

`RNAseq` contains both sbatch scripts and Rscripts and results from detailed analysis.

`scripts` contains python or R scripts.

`plots` contains 2 bash scripts with all commands to generate annotations and configuration files and output files of all figures generated by [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks).

Finally, `figures` contains symlink with the figures as used in the paper before processing with illustator.