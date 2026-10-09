# Ensure obs cache vectors exist on `model_obj` (idempotent, in-place safe)

Called automatically before compiled kernel updates; callers do not need
to invoke this explicitly.

## Usage

``` r
ensure_kernel_obs_cache(model_obj)

attach_kernel_obs_cache(model_obj)
```
