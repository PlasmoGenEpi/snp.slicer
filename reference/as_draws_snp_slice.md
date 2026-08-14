# Assemble a posterior draws array from a SNP-Slice run

Converts the retained MCMC samples of every chain into a
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
(dimensions iteration \\\times\\ chain \\\times\\ variable) of
permutation-invariant scalar summaries. This is the reusable primitive
behind
[`convergence_diagnostics`](https://plasmogenepi.github.io/snp.slicer/reference/convergence_diagnostics.md):
because it returns a standard `draws_array`, any function in the
posterior package can be applied to it directly.

Diagnostics are computed on invariant summaries rather than the raw
`A`/`D` matrices because strain labels are not identifiable across
iterations or chains, which makes element-wise R-hat / ESS on `A` or `D`
uninterpretable.

## Usage

``` r
as_draws_snp_slice(
  results,
  pars = c("logpost", "n_strains", "kstar", "ktrunc"),
  additional_burnin = 0
)
```

## Arguments

- results:

  A `snp_slice_results` object (with `store_mcmc = TRUE`).

- pars:

  Character vector of parameters to include. One of `"logpost"`,
  `"n_strains"`, `"kstar"`, `"ktrunc"`, which are all scalar per sample,
  or `"coi"` (per-host complexity of infection), which is a vector that
  expands to `coi[1]..coi[N]`).

- additional_burnin:

  Number of additional stored samples to discard from the start of each
  retained chain before assembling the array.

## Value

A
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html).

## Details

All chains are pooled into one array so that R-hat and ESS are computed
for the whole run. Without early stopping every chain retains the same
number of samples; when `gap` early-stopping fired, chains can differ in
length, in which case every chain is truncated to the shortest and a
warning is issued.

## Examples

``` r
result <- load_example_results()
draws <- as_draws_snp_slice(result)
posterior::summarise_draws(draws)
#> # A tibble: 4 × 10
#>   variable    mean  median     sd    mad      q5     q95  rhat ess_bulk ess_tail
#>   <chr>      <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <dbl>    <dbl>    <dbl>
#> 1 logpost  -7.51e4 -75010. 123.   109.   -75254. -74903.  2.19     4.01    41.7 
#> 2 n_strai…  5.11e1     52    2.16   1.48     48      53   7.59     3.32     3.52
#> 3 kstar     1.13e2    113    0      0       113     113  NA       NA       NA   
#> 4 ktrunc    1.34e2    134    0      0       134     134  NA       NA       NA   
```
