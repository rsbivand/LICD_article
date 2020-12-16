# LICD_article

[![CODECHECK](https://codecheck.org.uk/img/codeworks-badge.svg)](https://doi.org/10.5281/zenodo.4279275)

Reproduction materials for "The application of Local Indicators for Categorical Data (LICD) to explore spatial dependence in archaeological spaces" by Francesco Carrer, Tomasz M. Kossowski, Justyna Wilk, Michał B. Pietrzak and Roger Bivand

Script `LICD-Archaeo-HLC.R` runs the Torridge example, downloading data from given URLs, so requires internet access; it sources local_JC0.R. It generates `Torridge_*.jpeg` and `jc_out.csv`.

Script `LICD-Archaeo-Grid.R` runs the Barmose example, using data from the **archdata** package. It generates `Barmose_*.jpeg` and `barmose_jc_out.csv`.

Both scripts are run using `source(..., echo=TRUE)`, with output to `*_output.Rout` respectively. This also contains the output  of `sessionInfo()` for the output files included in this repository.

File `HLC_map.zip` is a compressed stand-alone interactive webmap produced in running `LICD-Archaeo-HLC.R`; to use, uncompress and open `index.html`. It was made using the development version of **mapview** (see https://github.com/r-spatial/mapview/issues/336) as noted on line 151 of the script.
