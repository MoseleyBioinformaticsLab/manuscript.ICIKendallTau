
# ICI-Kendall-Tau Manuscript

[![DOI](https://zenodo.org/badge/429873882.svg)](https://zenodo.org/badge/latestdoi/429873882)

All the code for generating the manuscript:

Information-Content-Informed Kendall-tau Correlation: Utilizing Missing Values, Robert M Flight, Praneeth S Bhatt, Hunter NB Moseley [biorxiv](https://doi.org/10.1101/2022.02.24.481854) is contained within this repository.

The current versions can be obtained directly from this repository:

  * [Manuscript html](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTau/ici_kt_manuscript.html)
  * [Supplemental html](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTau/supplemental_materials.html)
  * [Manuscript word](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTau/ici_kt_manuscript.docx)
  * [Supplemental word](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTau/supplemental_materials.docx)

## License

The contents of this work are licensed under a CC-BY license.
If you use any content, you must give attribution to this original work.

## R Packages Needed

The repository used R 4.1.0, and {renv} 1.0.0
For all of the other packages needed, see the file `renv.lock`.

To setup to be able to rerun everything here, you can clone the repo from [github](https://github.com/MoseleyBioinformaticsLab/manuscript.ICIKendallTau) or download it from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.6309187), and then from within that folder:

```
# clone from github
git clone https://github.com/MoseleyBioinformaticsLab/manuscript.ICIKendallTau.git
```

```
# download from zenodo
wget https://zenodo.org/records/10573598/files/MoseleyBioinformaticsLab/manuscript.ICIKendallTau-draft_v4.1_2024-01-26.zip?download=1 --output-document=manuscript.ICIKendallTau.zip
unzip manuscript.ICIKendallTau.zip
```

```r
# make sure renv is installed
install.packages("renv")
# restore the packages
renv::restore()
```

The {ICIKendallTau} package on GitHub is now different than the one used for this manuscript, you should install the one archived on Zenodo (v 0.3.20).

```
wget https://zenodo.org/records/10580528/files/ICIKendallTau_0.3.20.tar.gz?download=1 --output-document=ICIKendallTau_0.3.20.tar.gz
tar -xzvf ICIKendallTau_0.3.20.tar.gz
cd ICIKendallTau
```

Then start an R session from within the directory, and install it.

```r
remotes::install_local()
```

## Rerun All Underlying Analyses

If you want to rerun everything from the beginning, you can do `tar_make()` on the project, and it should just run it.
However, some of the calculations require a **lot** of compute resources, or will just take a long time, even with multiple cores.
I would definitely recommend at the very least 64 GB of RAM, and you may need more depending on how many cores you are using.
The compute node we ran on had 80 cores, and 1 TB of RAM, and we regularly used 500 GB of RAM for some of the calculations.

```r
targets::tar_make()
```


## Obtaining targets Cache

To make it easier to at least generate the manuscript, there is a copy of the {targets} cache on Zenodo, in two parts [here](https://zenodo.org/doi/10.5281/zenodo.10570285) and [here](https://zenodo.org/doi/10.5281/zenodo.10570255).
Make sure you have lots of room, there are 68 GB worth of files (and those are already compressed R data files).
This is in three parts, and can be downloaded using `wget`.

Having this cache, you have the state of the computations when we submitted the manuscript, and can examine almost any of the outputs you want by loading them using `tar_load(object)`.

```
# make sure you are wherever you want the targets cache, like wherever
# you cloned the github repo to, at the top level of the directory.
wget https://zenodo.org/records/10570286/files/manuscript.ICIKendallTau.targets_cache.parts0.tar?download=1 --output-document=manuscript.ICIKendallTau.targets_cache.parts0.tar
wget https://zenodo.org/records/10570286/files/manuscript.ICIKendallTau.targets_cache.parts1.tar?download=1 --output-document=manuscript.ICIKendallTau.targets_cache.parts1.tar
wget https://zenodo.org/records/10570256/files/manuscript.ICIKendallTau.targets_cache.parts2.tar?download=1 --output-document=manuscript.ICIKendallTau.targets_cache.parts2.tar

# now untar them. Note the "k", this keeps the directory from getting overwritten by subsequent un-tar
tar -xkvf manuscript.ICIKendallTau.targets_cache.parts0.tar
tar -xkvf manuscript.ICIKendallTau.targets_cache.parts1.tar
tar -xkvf manuscript.ICIKendallTau.targets_cache.parts2.tar
```
