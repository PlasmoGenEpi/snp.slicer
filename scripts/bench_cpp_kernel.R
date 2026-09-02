#!/usr/bin/env Rscript

# Repeatable C++ kernel bench. Writes TSV + JSON under profile/.
# Usage: Rscript scripts/bench_cpp_kernel.R [warmup] [iters] [label]

args <- commandArgs(trailingOnly = TRUE)
warmup_n <- if (length(args) >= 1) as.integer(args[[1]]) else 20L
iters_n <- if (length(args) >= 2) as.integer(args[[2]]) else 40L
label <- if (length(args) >= 3) args[[3]] else "baseline"

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(file_arg) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."))
} else {
  getwd()
}
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  root <- getwd()
}

suppressPackageStartupMessages({
  if (!requireNamespace("pkgload", quietly = TRUE) &&
      !requireNamespace("devtools", quietly = TRUE)) {
    stop("Need pkgload or devtools to load snp.slicer from source")
  }
})

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

elapsed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}

for (i in seq_len(warmup_n)) {
  state <- slice_iter(state, model_obj)
}

obs <- kernel_obs_args(model_obj)
llik_tab <- kernel_llik_tab(model_obj)
model_type <- model_type_id(model_obj$name)
kernel <- model_obj$kernel

fused_slice <- numeric(iters_n)
fused_cpp <- numeric(iters_n)
loglik_t <- numeric(iters_n)
logprior_t <- numeric(iters_n)
update_s_t <- numeric(iters_n)
update_a_t <- numeric(iters_n)
update_d_t <- numeric(iters_n)
update_mu_t <- numeric(iters_n)

state_fused <- state
for (i in seq_len(iters_n)) {
  fused_slice[[i]] <- elapsed(state_fused <- slice_iter(state_fused, model_obj))
}

state_cpp <- state
for (i in seq_len(iters_n)) {
  fused_cpp[[i]] <- elapsed({
    result <- cpp_slice_iter(
      A = state_cpp$A,
      D = state_cpp$D,
      mu = state_cpp$mu,
      mixed = as.integer(state_cpp$mixed),
      kplus = as.integer(state_cpp$kplus),
      kstar = as.integer(state_cpp$kstar),
      kmin = as.integer(state_cpp$kmin),
      ktrunc = as.integer(state_cpp$ktrunc),
      y = model_obj$y,
      r = model_obj$r,
      rho = model_obj$rho,
      alpha = model_obj$alpha,
      N = as.integer(model_obj$N),
      P = as.integer(model_obj$P),
      model_type = model_type,
      llik_tab = llik_tab,
      loglik_const = obs$loglik_const,
      obs_code = obs$obs_code
    )
    state_cpp$A <- result$A
    state_cpp$D <- result$D
    state_cpp$mu <- result$mu
    state_cpp$kplus <- result$kplus
    state_cpp$kstar <- result$kstar
    state_cpp$ktrunc <- result$ktrunc
  })
}

state_parts <- state
for (i in seq_len(iters_n)) {
  update_s_t[[i]] <- elapsed(state_parts <- kernel$update_s(state_parts, model_obj))
  update_a_t[[i]] <- elapsed(state_parts <- kernel$update_a(state_parts, model_obj))
  update_d_t[[i]] <- elapsed(state_parts <- kernel$update_d(state_parts, model_obj))
  update_mu_t[[i]] <- elapsed(state_parts <- kernel$update_mu(state_parts, model_obj))
  loglik_t[[i]] <- elapsed({
    state_parts$loglik <- as.numeric(cpp_loglik_total(
      A = state_parts$A,
      D = state_parts$D,
      y = model_obj$y,
      r = model_obj$r,
      model_type = model_type,
      llik_tab = llik_tab,
      loglik_const = obs$loglik_const,
      obs_code = obs$obs_code
    ))
  })
  logprior_t[[i]] <- elapsed({
    state_parts$logpost <- as.numeric(state_parts$loglik) +
      logprior_a(state_parts$A, state_parts$mu, model_obj$alpha) +
      logprior_mu(state_parts$mu, model_obj$alpha, model_obj$N) +
      logprior_d(state_parts$D, model_obj$rho)
  })
}

median_n <- function(x) stats::median(x)
sum_n <- function(x) sum(x)

rows <- list(
  list(op = "slice_iter", median_s = median_n(fused_slice), total_s = sum_n(fused_slice), n = iters_n),
  list(op = "cpp_slice_iter", median_s = median_n(fused_cpp), total_s = sum_n(fused_cpp), n = iters_n),
  list(op = "update_s", median_s = median_n(update_s_t), total_s = sum_n(update_s_t), n = iters_n),
  list(op = "update_a", median_s = median_n(update_a_t), total_s = sum_n(update_a_t), n = iters_n),
  list(op = "update_d", median_s = median_n(update_d_t), total_s = sum_n(update_d_t), n = iters_n),
  list(op = "update_mu", median_s = median_n(update_mu_t), total_s = sum_n(update_mu_t), n = iters_n),
  list(op = "cpp_loglik_total", median_s = median_n(loglik_t), total_s = sum_n(loglik_t), n = iters_n),
  list(op = "logprior", median_s = median_n(logprior_t), total_s = sum_n(logprior_t), n = iters_n)
)

out_dir <- file.path(root, "profile")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
tsv_path <- file.path(out_dir, sprintf("%s-%s.tsv", label, stamp))
json_path <- file.path(out_dir, sprintf("%s-%s.json", label, stamp))
meta_path <- file.path(out_dir, sprintf("%s-%s.meta.tsv", label, stamp))

tsv <- do.call(rbind, lapply(rows, function(r) {
  data.frame(op = r$op, median_s = r$median_s, total_s = r$total_s, n = r$n,
             stringsAsFactors = FALSE)
}))
write.table(tsv, tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)

meta <- data.frame(
  key = c("label", "N", "P", "n_mixed", "kplus", "kstar", "ktrunc", "warmup", "iters", "model"),
  value = c(
    label,
    as.character(model_obj$N),
    as.character(model_obj$P),
    as.character(length(state$mixed)),
    as.character(state$kplus),
    as.character(state$kstar),
    as.character(state$ktrunc),
    as.character(warmup_n),
    as.character(iters_n),
    model_obj$name
  ),
  stringsAsFactors = FALSE
)
write.table(meta, meta_path, sep = "\t", row.names = FALSE, quote = FALSE)

json_rows <- paste(vapply(rows, function(r) {
  sprintf(
    '{"op":"%s","median_s":%.8f,"total_s":%.8f,"n":%d}',
    r$op, r$median_s, r$total_s, r$n
  )
}, character(1)), collapse = ",")
writeLines(
  sprintf(
    '{"label":"%s","N":%d,"P":%d,"n_mixed":%d,"kplus":%d,"kstar":%d,"ktrunc":%d,"warmup":%d,"iters":%d,"model":"%s","ops":[%s]}',
    label, model_obj$N, model_obj$P, length(state$mixed),
    as.integer(state$kplus), as.integer(state$kstar), as.integer(state$ktrunc),
    warmup_n, iters_n, model_obj$name, json_rows
  ),
  json_path
)

cat("Wrote ", tsv_path, "\n", sep = "")
print(tsv)
cat("\nmeta:\n")
print(meta)
cat("JSON: ", json_path, "\n", sep = "")
