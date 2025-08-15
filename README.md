# SNP-Slice: Bayesian Nonparametric Resolution of Multi-Strain Infections
This library provides an installable implementation of the original `snp-slice` [model](https://github.com/nianqiaoju/snp-slice). We have no affiliation with this author and provide this implementation as is.

## Overview

SNP-Slice is a Bayesian nonparametric method for resolving multi-strain infections using slice sampling. The algorithm simultaneously estimates strain haplotypes and links them to hosts from SNP data.

## Installation

```r
# Install from our r-universe
install.packages(
    'snp.slicer', 
    repos = c('https://plasmogenepi.r-universe.dev', 'https://cloud.r-project.org')
)
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

```
@article{Ju_Liu_He_2024,
    title={SNP-slice resolves mixed infections: simultaneously unveiling strain haplotypes and linking them to hosts},
    volume={40},
    ISSN={1367-4811},
    DOI={10.1093/bioinformatics/btae344},
    number={6},
    journal={Bioinformatics},
    author={Ju, Nianqiao and Liu, Jiawei and He, Qixin},
    year={2024},
    month=june,
    pages={btae344}
}
```
> [Bioinformatics Article](https://academic.oup.com/bioinformatics/article/40/6/btae344/7695237)

## License

This package is licensed under the GPL-3 License.
