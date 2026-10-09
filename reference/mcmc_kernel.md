# MCMC kernel seam

The kernel module owns the hot slice-sampler loops. Callers cross this
seam via `model_obj$kernel`; tests can swap the R reference adapter for
the compiled adapter without changing
[`slice_iter()`](https://plasmogenepi.github.io/snp.slicer/reference/slice_iter.md).

## Details

For count/categorical models with the compiled backend, log-likelihood
layout caches (`kernel_loglik_const`, `kernel_obs_code`) are built
automatically on first use (`resolve_mcmc_kernel`, `slice_iter`).
