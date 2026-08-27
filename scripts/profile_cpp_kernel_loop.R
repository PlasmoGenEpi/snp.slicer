#!/usr/bin/env Rscript

# Long fused-kernel loop for macOS `sample`. Prints PID and loops until killed.
# Usage: Rscript scripts/profile_cpp_kernel_loop.R

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(file_arg) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."))
} else {
  getwd()
}

if (requireNamespace("pkgbuild", quietly = TRUE)) {
  pkgbuild::compile_dll(root, debug = FALSE)
}
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, compile = FALSE, quiet = TRUE)
} else {
  devtools::load_all(root, quiet = TRUE)
}

data("example_snp_data", package = "snp.slicer")
processed <- preprocess_data(example_snp_data, "negative_binomial")
model_obj <- create_model("negative_binomial", processed)
model_obj <- setup_mcmc_kernel(model_obj, backend = "cpp")
set.seed(123)
state <- model_obj$initialize_state(model_obj, threshold = 0.001)
state <- slice_init(state, model_obj)
for (i in 1:20) {
  state <- slice_iter(state, model_obj)
}

cat("pid=", Sys.getpid(), " N=", model_obj$N, " P=", model_obj$P,
    " kplus=", state$kplus, " mixed=", length(state$mixed), "\n", sep = "")
flush.console()

repeat {
  state <- slice_iter(state, model_obj)
}
