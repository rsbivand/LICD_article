# LICD_article

Reproduction materials for "The application of Local Indicators for Categorical Data (LICD) to explore spatial dependence in archaeological spaces" by Francesco Carrer, Tomasz M. Kossowski, Justyna Wilk, Michał B. Pietrzak and Roger Bivand

Script `LICD-Archaeo-HLC.R` runs the Torridge example, downloading data from given URLs, so requires internet access; it sources local_JC0.R. It generates `Torridge_*.jpeg` and `jc_out.csv`.

Script `LICD-Archaeo-Grid.R` runs the Barmose example, using data from the **archdata** package. It generates `Barmose_*.jpeg` and `barmose_jc_out.csv`.

Both scripts are run using `source(..., echo=TRUE)`, with output to `*_output.Rout` respectively. This also contains the output  of `sessionInfo()` for the output files included in this repository.

