#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

project <- normalizePath(Sys.getenv("SSQR_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(project)
run_tag <- Sys.getenv("SSQR_RUN_TAG", "")
output_root <- Sys.getenv("SSQR_OUTPUT_ROOT", "")
stopifnot(nzchar(run_tag), nzchar(output_root))
run_dir <- file.path(output_root, run_tag)
tables <- file.path(run_dir, "tables")
gate <- readLines(file.path(tables, "run_summary.txt"), warn = FALSE)
stopifnot(any(gate == "paper_staging_ready: TRUE"))

summary <- read.csv(file.path(tables, "config_method_summary.csv"), check.names = FALSE)
ratios <- read.csv(file.path(tables, "time_ratios.csv"), check.names = FALSE)
discrepancy <- read.csv(file.path(tables, "numerical_discrepancy.csv"), check.names = FALSE)
screening <- read.csv(file.path(tables, "paper_screening_tau05_staging.csv"), check.names = FALSE)
retained <- read.csv(file.path(tables, "retained_size_replication_metrics.csv"), check.names = FALSE)
stage <- file.path(run_dir, "paper_staging")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)

paper <- summary[summary$method %in% c("direct_baseline", "unified_u11"), ]
paper$raw_method <- ifelse(paper$method == "direct_baseline",
                           ifelse(paper$model == "llqr", "llqr", "tvc_rq"),
                           "unified_u11")
paper$paper_method <- ifelse(paper$method == "direct_baseline",
                             "direct local fit", "screen-seq")
paper$Mm.factor <- ifelse(paper$method == "unified_u11",
                          ifelse(paper$model == "llqr", 0.1, 1e-5), NA_real_)
paper$time_range_sec <- paper$time_max_sec - paper$time_min_sec
paper$max_average_relative_bias <- paper$max_average_relative_discrepancy
paper$paper_status <- ifelse(
  paper$solver_stable & paper$path_accepted & paper$fallback_count == 0L,
  ifelse(paper$method == "unified_u11", "all_true_ppro", "complete"),
  "invalid"
)
write.csv(paper, file.path(stage, "rep500_seed2025_method_summary_candidate.csv"), row.names = FALSE)
write.csv(screening, file.path(stage, "rep500_seed2025_screening_tau05_candidate.csv"), row.names = FALSE)
write.csv(ratios, file.path(stage, "runtime_ratios_complete.csv"), row.names = FALSE)
write.csv(discrepancy, file.path(stage, "numerical_discrepancy_complete.csv"), row.names = FALSE)
write.csv(retained, file.path(stage, "retained_size_replication_metrics.csv"), row.names = FALSE)
write.csv(summary[summary$method == "lean_seq", ],
          file.path(stage, "lean_seq_ablation_candidate.csv"), row.names = FALSE)

tau05 <- paper[paper$tau == 0.5, ]
tau05$case_label <- paste0("Case ", tau05$paper_case)
write.csv(tau05, file.path(stage, "runtime_figure_data_tau05.csv"), row.names = FALSE)

draw_runtime <- function(device, path) {
  device(path)
  op <- par(mfrow = c(2, 2), mar = c(3.4, 3.7, 2.2, 0.8), mgp = c(2.2, 0.7, 0),
            family = "serif")
  on.exit({ par(op); dev.off() }, add = FALSE)
  for (case_id in 1:4) {
    d <- tau05[tau05$paper_case == case_id, ]
    d <- d[order(d$n, d$method), ]
    ns <- sort(unique(d$n))
    direct <- d[d$method == "direct_baseline", ]; direct <- direct[match(ns, direct$n), ]
    u11 <- d[d$method == "unified_u11", ]; u11 <- u11[match(ns, u11$n), ]
    z <- rbind(direct$time_mean_sec, u11$time_mean_sec)
    ymax <- max(direct$time_max_sec, u11$time_max_sec) * 1.08
    mids <- barplot(z, beside = TRUE, col = c("grey45", "#EE7733"), border = "black",
                    names.arg = ns, ylim = c(0, ymax), ylab = "Time (s)", xlab = "n",
                    main = paste0("Case ", case_id))
    segments(mids[1, ], direct$time_min_sec, mids[1, ], direct$time_max_sec)
    segments(mids[2, ], u11$time_min_sec, mids[2, ], u11$time_max_sec)
    segments(mids[1, ] - .08, direct$time_min_sec, mids[1, ] + .08, direct$time_min_sec)
    segments(mids[1, ] - .08, direct$time_max_sec, mids[1, ] + .08, direct$time_max_sec)
    segments(mids[2, ] - .08, u11$time_min_sec, mids[2, ] + .08, u11$time_min_sec)
    segments(mids[2, ] - .08, u11$time_max_sec, mids[2, ] + .08, u11$time_max_sec)
    if (case_id == 1L) legend("topleft", c("direct local fit", "screen-seq U11"),
                              fill = c("grey45", "#EE7733"), bty = "n", cex = .78)
  }
}
draw_runtime(function(path) png(path, width = 2200, height = 1600, res = 250),
             file.path(stage, "runtime_tau05_preview.png"))
draw_runtime(function(path) pdf(path, width = 7.4, height = 5.2, family = "serif"),
             file.path(stage, "runtime_tau05_preview.pdf"))

fmt <- function(x) ifelse(abs(x) >= 10, sprintf("%.2f", x), sprintf("%.4f", x))
runtime_lines <- c(
  "\\begin{tabular}{lllrrrr}", "\\toprule",
  "Case & $\\tau$ & Method & $n=1000$ & $n=2000$ & $n=5000$ & $n=10000$ \\\\",
  "\\midrule"
)
for (pc in 1:4) for (tau in c(.2, .5, .8)) for (method in c("direct_baseline", "unified_u11")) {
  d <- paper[paper$paper_case == pc & paper$tau == tau & paper$method == method, ]
  d <- d[order(d$n), ]
  stopifnot(nrow(d) == 4L)
  label <- if (method == "direct_baseline") "direct local fit" else "\\texttt{screen-seq}"
  runtime_lines <- c(runtime_lines, sprintf(
    "%s & %.1f & %s & %s \\\\", if (method == "direct_baseline") paste0("Case ", pc) else "",
    tau, label, paste(fmt(d$time_mean_sec), collapse = " & ")
  ))
}
runtime_lines <- c(runtime_lines, "\\bottomrule", "\\end{tabular}")
writeLines(runtime_lines, file.path(stage, "runtime_complete_candidate.tex"))

for (tau in c(.2, .5, .8)) {
  d <- paper[paper$tau == tau & paper$method == "unified_u11", ]
  d <- d[order(d$paper_case, d$n), ]
  lines <- c("\\begin{tabular}{lrrrr}", "\\toprule",
             "Case & $n=1000$ & $n=2000$ & $n=5000$ & $n=10000$ \\\\", "\\midrule")
  for (pc in 1:4) {
    z <- d[d$paper_case == pc, ]; z <- z[order(z$n), ]
    lines <- c(lines, sprintf("Case %d & %s \\\\", pc,
      paste(sprintf("%.2e", z$max_average_relative_discrepancy), collapse = " & ")))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(stage, sprintf("discrepancy_tau%02d_candidate.tex", round(100 * tau))))
}

writeLines(c(
  "Paper staging generated only after the 48-configuration stability gate passed.",
  "Canonical paper/data, manuscript sources, and existing figures were not modified.",
  "Cases 1--2 are LLQR; Cases 3--4 are TVCQR model cases 1--2.",
  "Numerical discrepancy is direct-relative; lean_seq is staged separately as an ablation."
), file.path(stage, "README.txt"))
cat(sprintf("paper staging written to %s\n", stage))
