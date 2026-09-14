# New Case 2 only; sample size and evaluation-grid size are separate objects.
llqr_interior_grid <- function(x) {
  selected <- sort(x[x >= 0.1 & x <= 0.9])
  list(z = if (length(selected)) selected else 0.5,
       n_interior = length(selected), placeholder = !length(selected))
}

generate_logistic_case2 <- function(n, seed) {
  set.seed(seed)
  x <- runif(n)
  error <- rlogis(n, location = 0, scale = sqrt(3) / pi)
  grid <- llqr_interior_grid(x)
  c(list(x = x, y = 1 + 2 * x^2 + error), grid)
}
