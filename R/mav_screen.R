# MAV Screen Function
# Post: model mean-variance relationship and identify overdispersed features.
# Parameter: path: single input .feather file containing quantile-normalized count matrix (first column must be "pos").
#            fitting_model: Model to fit mean-variance trend: "gam" (default) or "loess".
#            spline_basis: For "gam": number of basis functions (k = 40); for "loess": smoothing span (0.5).
#            seed: Random seed for reproducibility (default = 42).
#            font_size: Font size for plot labels and axes (default = 10).
#            nrow_sample_per: Percentage for number of rows to randomly sample for the density plot (default = 0.2).
#            out_dir: Output Directory (Optional)
# Output: Annotated dataframe saved as a .feather file
#         Two diagnostic plots saved as .PNG in "mav_screen/" folder.
mav_screen <- function(path, fitting_model = "gam", spline_basis = NULL, seed = 42, font_size = 10, nrow_sample_per = 0.2, out_dir = NULL) {
  # load library
  suppressPackageStartupMessages({
    library(arrow)
    library(matrixStats)
    library(mgcv)
    library(ggplot2)
    library(cowplot)
    library(ggExtra)
    library(ggpointdensity)
    library(tools)
  })

  # Default Setting
  if (fitting_model == "gam" && (is.null(spline_basis) || length(spline_basis) == 0)) {
    spline_basis <- 40
  } else if (fitting_model == "loess" && (is.null(spline_basis) || length(spline_basis) == 0)) {
    spline_basis <- 0.5
  }
  # return a list of 40


  import_name <- basename(tools::file_path_sans_ext(path))
  output_name <- paste0(import_name, "_mav_screen") # new file name

  if (is.null(out_dir)) {
    out_dir <- "mav_screen"
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Set up mav screen
  wgc_qnorm = read_feather(path)
  pos_list <- wgc_qnorm[, "pos"]
  wgc_qnorm = wgc_qnorm[, -1]

  nrow_sample <- nrow(wgc_qnorm) * nrow_sample_per
  data <- as.matrix(wgc_qnorm)

  clip_max <- sqrt(ncol(data))
  print(paste0("using max clipping value: ",clip_max))

  gene_mean_all <- rowMeans(data)
  gene_var_all <- rowVars(data)

  data_filter <- data

  gene_mean <- rowMeans(data_filter)
  gene_var <- rowVars(data_filter)
  data_fit <- data.frame(X=gene_mean, Y=gene_var)

  # regression model
  k <- spline_basis
  span <- spline_basis
  set.seed(21) # no need
  if (fitting_model == "gam") {
    fit_model <- gam(formula = log2(x=Y) ~ s(log2(x=X), k=k), data = data_fit)
  } else if (fitting_model == "loess") {
    fit_model <- loess(formula = log2(x=Y) ~ log2(x=X), data=data_fit, span=span)
  }

  # manipulation
  gene_var_expect <- 2^(fit_model$fitted)
  gene_sd_expect <- sqrt(gene_var_expect)

  gene_var_norm <- (data_filter - gene_mean)/gene_sd_expect
  gene_var_norm_mean = rowMeans(gene_var_norm)
  gene_var_norm_var = rowVars(gene_var_norm)
  gene_hyper_var <- rowSums(gene_var_norm^2)/(ncol(data_filter) - 1)
  gene_var_norm_clip = gene_var_norm

  # record the clipped rows
  clip_or_not = matrix(0, nrow = dim(gene_var_norm)[1], ncol = dim(gene_var_norm)[2])
  clip_or_not[which(gene_var_norm > clip_max)] <- 1
  clip_or_not_vec = rowSums(clip_or_not)
  gene_var_norm_clip[which(gene_var_norm > clip_max)] <- clip_max

  gene_var_norm_mean_clip = rowMeans(gene_var_norm_clip)
  gene_var_norm_var_clip = rowVars(gene_var_norm_clip)
  gene_hyper_var_clip <- rowSums(gene_var_norm_clip^2)/(ncol(data_filter) - 1)

  result <- data.frame(pos=pos_list, clip=clip_or_not_vec, mean_orig=gene_mean_all, mean=gene_mean, var=gene_var, norm_mean=gene_var_norm_mean, norm_var=gene_var_norm_var,
                       var_expect=gene_var_expect, hypervar=gene_hyper_var, norm_mean_clip=gene_var_norm_mean_clip, norm_var_clip=gene_var_norm_var_clip, hypervar_clip=gene_hyper_var_clip,
                       log2Mean = log2(gene_mean), log2Var = log2(gene_var))

  path_feather <- file.path(out_dir, paste0(output_name, ".feather"))
  write_feather(result, path_feather)

  # Plot
  p1 <- ggplot(result, aes(log2(mean), log2(var))) + geom_point(alpha = 1/20) +
    geom_point(data=result,aes(log2(mean),log2(var_expect)), color="red", size=0.1) +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))

  p2 <- ggplot(result, aes(log2(mean), hypervar)) + geom_point(alpha = 1/20) +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))

  combined_plot <- plot_grid(p1, p2, labels = c('A', 'B'))
  print(combined_plot)

  path_main_plot <- file.path(out_dir, paste0(output_name, ".png"))
  ggsave(path_main_plot, plot = combined_plot)


  set.seed(seed)
  nrow_sample = min(nrow(result), nrow_sample)
  result_sample = result[sample(nrow(result), nrow_sample), ]

  p3 <- ggplot(result_sample, aes(log2(mean), log2(var))) +
    geom_pointdensity() + scale_color_viridis_c() +
    geom_point(data=result,aes(log2(mean),log2(var_expect)), color="red", size=0.1) +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  print(p3)


  density_fig_filename = paste0(output_name, "_smooth-", format(nrow_sample, scientific=FALSE), "_seed-", format(seed, scientific=FALSE), ".png")
  density_fig_dir_filename = file.path(out_dir, density_fig_filename)
  ggsave(density_fig_dir_filename, plot = p3)

}
