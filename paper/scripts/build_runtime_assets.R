#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0L) {
  stop("This script should be run via Rscript.")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
paper_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
data_path <- file.path(paper_dir, "data", "runtime_summary.csv")
art_dir <- file.path(paper_dir, "art")
generated_dir <- file.path(paper_dir, "generated")

dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

runtime <- read.csv(
  data_path,
  stringsAsFactors = FALSE,
  colClasses = rep("character", 7L)
)
runtime$case_id <- as.integer(runtime$case_id)
runtime$tau <- as.numeric(runtime$tau)
runtime$n <- as.integer(runtime$n)
runtime$mean_str <- runtime$mean
runtime$range_str <- runtime$range
runtime$mean_value <- as.numeric(runtime$mean_str)
runtime$range_value <- as.numeric(runtime$range_str)

case_map <- c(
  "1" = "Case 1: LLQR (Gaussian design)",
  "2" = "Case 2: LLQR (Uniform design)",
  "3" = "Case 3: TVCQR (i.i.d. design)",
  "4" = "Case 4: TVCQR (dependent design)"
)

panel_map <- c(
  "1" = "(a) Case 1\nLLQR (Gaussian design)",
  "2" = "(b) Case 2\nLLQR (Uniform design)",
  "3" = "(c) Case 3\nTVCQR (i.i.d. design)",
  "4" = "(d) Case 4\nTVCQR (dependent design)"
)

baseline_map <- c(
  "1" = "llqr",
  "2" = "llqr",
  "3" = "rq",
  "4" = "rq"
)

all_methods <- c(
  "llqr", "rq", "seq", "seq(ft)",
  "screen-seq(1)", "screen-seq(2)",
  "screen-seq(1)(ft)", "screen-seq(2)(ft)"
)

expected_rows <- 4L * 3L * 7L * 5L
if (nrow(runtime) != expected_rows) {
  stop(sprintf("Expected %d runtime rows, found %d.", expected_rows, nrow(runtime)))
}

build_main_figure <- function(dat) {
  selected_n <- c(500L, 2000L, 5000L)
  figure_rows <- dat[dat$tau == 0.5 & dat$n %in% selected_n, ]
  keep <- mapply(
    function(case_id, method) {
      method %in% c(baseline_map[[as.character(case_id)]], "seq(ft)", "screen-seq(2)(ft)")
    },
    figure_rows$case_id,
    figure_rows$method
  )
  figure_rows <- figure_rows[keep, ]

  figure_rows$method_display <- ifelse(
    figure_rows$method == baseline_map[as.character(figure_rows$case_id)],
    "baseline",
    figure_rows$method
  )
  figure_rows$method_display <- factor(
    figure_rows$method_display,
    levels = c("baseline", "seq(ft)", "screen-seq(2)(ft)")
  )
  figure_rows$panel <- factor(
    panel_map[as.character(figure_rows$case_id)],
    levels = unname(panel_map)
  )
  figure_rows$n_label <- factor(
    paste0("n=", figure_rows$n),
    levels = paste0("n=", selected_n)
  )

  palette <- c(
    "baseline" = "#4A4A4A",
    "seq(ft)" = "#9B9B9B",
    "screen-seq(2)(ft)" = "#B7C3D0"
  )

  p <- ggplot(
    figure_rows,
    aes(x = n_label, y = mean_value, fill = method_display)
  ) +
    geom_col(
      position = position_dodge(width = 0.72),
      width = 0.64,
      colour = "black",
      linewidth = 0.25
    ) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Mean computation time (s)") +
    theme_bw(base_size = 9, base_family = "serif") +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 8.2),
      legend.key.size = grid::unit(10, "pt"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
      strip.text = element_text(size = 8.4, face = "plain", margin = margin(2, 2, 2, 2)),
      axis.title.y = element_text(size = 8.8),
      axis.text = element_text(size = 8),
      axis.ticks = element_line(linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.border = element_rect(linewidth = 0.35),
      plot.margin = margin(4, 8, 2, 2)
    )

  ggsave(
    filename = file.path(art_dir, "runtime_tau05_main.pdf"),
    plot = p,
    width = 7.2,
    height = 5.2,
    units = "in",
    device = grDevices::pdf
  )
  ggsave(
    filename = file.path(art_dir, "runtime_tau05_main.png"),
    plot = p,
    width = 7.2,
    height = 5.2,
    units = "in",
    dpi = 600
  )
}

emit_longtable <- function(dat) {
  output_path <- file.path(generated_dir, "runtime_summary_longtable.tex")
  con <- file(output_path, open = "w")
  on.exit(close(con), add = TRUE)

  writeLines(c(
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2.8pt}",
    "\\renewcommand{\\arraystretch}{0.95}",
    "\\begin{longtable}{lllrrrrrrrrrr}",
    "\\caption{Complete runtime summaries for the four simulation settings. Case~1 = LLQR (Gaussian design), Case~2 = LLQR (Uniform design), Case~3 = TVCQR (i.i.d. design), and Case~4 = TVCQR (dependent design). Entries report the Monte Carlo mean computation time and range (max--min), in seconds, over $1000$ replications.}\\label{tab:runtime_full}\\\\",
    "\\toprule",
    "Case & $\\tau$ & Method & \\multicolumn{2}{c}{$n=200$} & \\multicolumn{2}{c}{$n=500$} & \\multicolumn{2}{c}{$n=1000$} & \\multicolumn{2}{c}{$n=2000$} & \\multicolumn{2}{c}{$n=5000$}\\\\",
    "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\\cmidrule(lr){12-13}",
    " &  &  & Mean & Range & Mean & Range & Mean & Range & Mean & Range & Mean & Range\\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\multicolumn{13}{l}{\\tablename\\ \\thetable\\ (continued)}\\\\",
    "\\toprule",
    "Case & $\\tau$ & Method & \\multicolumn{2}{c}{$n=200$} & \\multicolumn{2}{c}{$n=500$} & \\multicolumn{2}{c}{$n=1000$} & \\multicolumn{2}{c}{$n=2000$} & \\multicolumn{2}{c}{$n=5000$}\\\\",
    "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\\cmidrule(lr){12-13}",
    " &  &  & Mean & Range & Mean & Range & Mean & Range & Mean & Range & Mean & Range\\\\",
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{13}{r}{Continued on next page}\\\\",
    "\\midrule",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot"
  ), con)

  tau_order <- c(0.2, 0.5, 0.8)
  for (case_id in 1:4) {
    case_block <- dat[dat$case_id == case_id, ]
    baseline <- baseline_map[[as.character(case_id)]]
    method_order <- if (case_id <= 2) {
      c("llqr", "seq", "seq(ft)", "screen-seq(1)", "screen-seq(2)", "screen-seq(1)(ft)", "screen-seq(2)(ft)")
    } else {
      c("rq", "seq", "seq(ft)", "screen-seq(1)", "screen-seq(2)", "screen-seq(1)(ft)", "screen-seq(2)(ft)")
    }

    first_case_row <- TRUE
    for (tau in tau_order) {
      tau_block <- case_block[case_block$tau == tau, ]
      first_tau_row <- TRUE
      for (method in method_order) {
        method_block <- tau_block[tau_block$method == method, ]
        method_block <- method_block[order(method_block$n), ]
        if (nrow(method_block) != 5L) {
          stop(sprintf("Incomplete block for case %d, tau %.1f, method %s.", case_id, tau, method))
        }

        case_cell <- if (first_case_row) sprintf("Case %d", case_id) else ""
        tau_cell <- if (first_tau_row) sprintf("%.1f", tau) else ""
        numeric_cells <- unlist(lapply(seq_len(nrow(method_block)), function(i) {
          c(method_block$mean_str[i], method_block$range_str[i])
        }))
        line <- paste(
          c(case_cell, tau_cell, method, numeric_cells),
          collapse = " & "
        )
        writeLines(paste0(line, "\\\\"), con)

        first_case_row <- FALSE
        first_tau_row <- FALSE
      }
      writeLines("\\midrule", con)
    }
  }

  writeLines(c("\\end{longtable}", "\\endgroup"), con)
}

build_main_figure(runtime)
emit_longtable(runtime)
