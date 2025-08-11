# SNP-Slice: Bayesian Nonparametric Resolution of Multi-Strain Infections
This library provides an installable implementation of the original `snp-slice` [model](https://github.com/nianqiaoju/snp-slice). We have no affiliation with this author and provide this implementation as is.

## Overview

SNP-Slice is a Bayesian nonparametric method for resolving multi-strain infections using slice sampling with stick-breaking construction. The algorithm simultaneously unveils strain haplotypes and links them to hosts from sequencing data.

## Installation

```r
# Install from GitHub (when available)
devtools::install_github("m-murphy/snp.slice")
```

## Features

- **Multiple Observation Models**: Supports categorical, Poisson, binomial, and negative binomial models
- **Flexible Data Input**: Accepts read count data, categorical data, or file paths
- **Convergence Diagnostics**: Built-in monitoring and early stopping
- **Comprehensive Results**: Returns allocation matrix, dictionary matrix, and diagnostics

## Documentation

- `?snp_slice` - Main function documentation
- `vignette("introduction")` - Getting started guide
- `vignette("model_comparison")` - Comparing different models

## Citation

If you use SNP-Slice in your research, please cite:

> SNP-Slice Resolves Mixed Infections: Simultaneously Unveiling Strain Haplotypes and Linking Them to Hosts
> [Bioinformatics Article](https://academic.oup.com/bioinformatics/article/40/6/btae344/7695237)

## License

This package is licensed under the GPL-3 License.
