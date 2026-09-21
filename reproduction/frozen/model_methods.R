# Self-contained experiment-only interfaces for the frozen lean-seq kernels.

load_lean_kernels <- function(experiment_dir = Sys.getenv(
  "SSQR_EXPERIMENT_DIR", "experiments/hpc_u11_rep500")) {
  build_dir <- file.path(experiment_dir, "build")
  list(
    llqr = dyn.load(file.path(build_dir, "llqr_seq_lean_sortskip.so")),
    tvcqr = dyn.load(file.path(build_dir, "tvcqr_seq_lean_nohistory.so"))
  )
}

run_llqr_lean <- function(dll, x, y, z, tau, h, case) {
  n <- length(y)
  ne <- length(z)
  .Fortran(
    "llqr_seq_fortran", PACKAGE = dll[["name"]],
    x = as.double(x), y = as.double(y), z = as.double(z),
    m = as.integer(n), nvar = 1L, rounds = as.integer(ne),
    tau = as.double(tau), h = as.double(h), tol = 1e-14,
    maxit = 1000000L, case_int = as.integer(case), bland_int = 0L,
    ll_est = double(ne), d_ll_est = double(ne), it_num = integer(ne),
    residual_est = double(n), H_mat = integer(2L * ne)
  )
}

run_tvcqr_lean <- function(dll, x, y, tau) {
  n <- nrow(x)
  nvar <- ncol(x)
  q <- 2L * (nvar + 1L)
  .Fortran(
    "tvcqr_seq_fortran", PACKAGE = dll[["name"]],
    x = as.double(x), y = as.double(y), m = as.integer(n),
    nvar = as.integer(nvar), tau = as.double(tau), h = 0.0,
    tol = 1e-14, maxit = 1000000L, bland_int = 0L,
    theta_ll_est = double(n * (nvar + 1L)), it_num = integer(n),
    residual_est = double(1L), H_mat = integer(n * q)
  )
}
