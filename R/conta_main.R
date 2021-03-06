# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Run conta.
#'
#' Reads a given counts file, processes it and calculates
#' contamination likelihood for a set of ranges.
#'
#' @param tsv_file input tsv file
#' @param metrics_file input tsv metrics file in long format
#' @param filename_prefix experiment name (basename for outputs)
#' @param save_dir output folder
#' @param lr_th min avg. likelihood ratio per SNP to make a call, this
#'     number is highly dependent on the data type used and should be optimized
#'     by the user based on a sensitivity - specificity trade-off.
#' @param sim_level if non-zero, a contaminant at this level will
#'     be added to the current sample. The contaminant will be randomly
#'     generated from minor allele frequencies of the SNPs in the population.
#' @param baseline baseline file for blacklist or noise model
#' @param min_depth minimum depth for a SNP to be considered
#' @param max_depth maximum depth for a SNP to be considered
#' @param loh_lr_cutoff minimum likelihood ratio to call a region as LOH
#' @param loh_delta_cutoff minimum delta (het deviation) to call a region as LOH
#' @param loh_min_snps minimum number of SNPs in a region to consider LOH
#' @param loh_max_snps maxmimum number of SNPs in a region to use for LOH, if
#'     there are more SNPs, they are subsampled to this number
#' @param cf_correction cf correction calculated from empirical data
#' @param min_maf minimum minor allele frequency to include a SNP
#' @param min_cf minimum contamination fraction to call
#' @param blackswan blackswan term for maximum likelihood estimation
#' @param outlier_frac fraction of outlier SNPs (based on depth) to remove
#' @param tsv_rev_file input tsv file for reverse strand reads, if this option
#'     is provided then first tsv_file is considered as positive strand reads
#' @param cores number of cores to be used for parallelization
#' @param subsample Either NA (use all SNPs) or number of SNPs to subsample to
#' @param context_mode whether to run with errors calculated in 3-base context
#' @param default_het_mean default for heterozygote mean allele frequency
#' @param default_het_sum default heterozygote allele frequency alpha + beta
#'     from beta binomial alternative formulation
#' @param error_quantile_filter remove SNPs with error rates higher than this
#'     quantile. Default 1.0 does not remove any SNPs. The idea is to remove
#'     SNPs with high mean error rates, since they might also contain high
#'     variance and thus introduce spurious contamination likelihoods.
#' @param min_lr_to_cf_ratio do not make a conta call if likelihood ratio to
#'     contamination fraction is equal to or below this ratio. In general
#'     conta likelihood ratio is greater than and is correlated with
#'     the contamination fraction. In cases where this ratio is greatly violated
#'     such as 10 times or less likelihood ratio compared to contamination
#'     fraction, it is very likely that is is a spurious call generated due to
#'     LOH or other artifacts.
#' @param seed random seed
#'
#' @return none
#'
#' @importFrom utils packageVersion
#' @importFrom utils write.table
#' @export
conta_main <- function(tsv_file, sample_id, save_dir, filename_prefix = sample_id,
                       metrics_file = "",
                       lr_th = 0.005, sim_level = 0, baseline_file = NA,
                       min_depth = 15, max_depth = 10000, loh_lr_cutoff = 0.001,
                       loh_delta_cutoff = 0.2,
                       loh_min_snps = 20, loh_max_snps = 1000,
                       min_maf = 0.001, subsample = NA,
                       cf_correction = 0, min_cf = 0.002, blackswan = 0.2,
                       outlier_frac = 0, tsv_rev_file = NA,
                       cores = 2, context_mode = FALSE,
                       chr_y_male_threshold = 0.0005,
                       default_het_mean = 0.5,
                       default_het_sum = 200,
                       error_quantile_filter = 0.75,
                       min_lr_to_cf_ratio = 0.2,
                       seed = 1359) {

  options("digits" = 8)
  options("mc.cores" = cores)
  set.seed(seed)

  # Create output folder if it doesn't exist
  dir.create(save_dir, showWarnings = FALSE)

  # Write out a default TSV result file in case we exit early from fail_test
  empty_result <- empty_result(sample_id)
  out_file <- file.path(save_dir, paste(filename_prefix, "conta.tsv", sep = "."))
  write.table(empty_result, file = out_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  # Read baseline if available
  baseline <- NA
  if (!is.na(baseline_file)) {
    baseline <- read_data_table(baseline_file)
  }

  # Prep snp counts
  dat <- read_and_prep(tsv_file, tsv_rev_file, baseline, default_het_mean,
                       default_het_sum)

  # Set dat as a subset of dat if it exceeds a pre-determined size
  if (!is.na(subsample) & nrow(dat) > subsample) {
    dat <- dat[sort(sample(1:nrow(dat), subsample)), ]
  }

  # Simulate contamination if sim_level is non-0
  dat <- sim_conta(dat, sim_level)

  # Original depth by chr
  plot_depth_by_chr(dat, save_dir, filename_prefix, min_depth = min_depth,
                    ext_plot = "depth.png", min_maf = min_maf)

  # Add in more useful fields and filter based on depth and outlier fractions
  dat <- annotate_and_filter(dat, min_depth = min_depth, max_depth = max_depth,
                             out_frac = outlier_frac)

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat)

  # Depth by chr after filters
  plot_depth_by_chr(dat, save_dir, filename_prefix, min_depth = min_depth,
                    ext_plot = "filtered.depth.png", min_maf = min_maf)

  # Calculate substitution rates per base and add them to SNP data table
  EE <- calculate_error_model(dat, save_dir, filename_prefix,
                              context_mode = context_mode)
  dat <- add_error_rates(dat, EE, context_mode)

  # Remove low maf positions (they were kept for error model calculations)
  dat <- dat %>%
    dplyr::filter(maf > min_maf, maf < (1 - min_maf))

  # Optionally remove high error positions
  dat <- dat %>%
    dplyr::filter(er <= quantile(dat$er, error_quantile_filter))

  # Add experimental bayesian genotype estimation.
  dat <- bayesian_genotype(dat)

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat)

  # Obtain results and plot lr without LOH filter applied yet
  result <- optimize_likelihood(dat, lr_th, save_dir, filename_prefix,
                                loh = FALSE, blackswan, min_cf, cf_correction,
                                outlier_frac = outlier_frac)

  # Remove loss of heterozygosity regions
  bin_stats <- get_per_bin_loh(dat, save_dir, filename_prefix, min_lr = loh_lr_cutoff,
                               blackswan = blackswan, conta_cf = result$cf,
                               min_loh = loh_delta_cutoff,
                               min_snps = loh_min_snps,
                               max_snps = loh_max_snps)
  dat_loh <- exclude_high_loh_regions(dat, bin_stats)

  # Plot minor allele ratio plot (.vr) with LOH
  plot_minor_ratio(dat, dat_loh, save_dir, filename_prefix, max_snps = 100000)

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat_loh)

  # Calculate likelihood, max likelihood, and whether to call it
  result <- optimize_likelihood(dat_loh, lr_th, save_dir,
                                filename_prefix, loh = TRUE, blackswan, min_cf,
                                cf_correction, outlier_frac = outlier_frac,
                                min_lr_to_cf_ratio = min_lr_to_cf_ratio)

  # Calculate chr X and Y counts from metrics file if provided
  chr_y_stats <- chr_stats(metrics_file, "chrY")
  chr_x_stats <- chr_stats(metrics_file, "chrX")
  x_het_snps <- nrow(dat[gt == "0/1" & (chrom == "X" | chrom == "chrX"), ])

  # Make a sex/gender call for the host
  sex <- get_sex_call(chr_y_stats, chr_y_male_threshold, result)

  # Make a pregnancy call for the contaminant
  pregnancy <- get_pregnancy_call(result, sex,
                                  chr_y_stats, chr_y_male_threshold)

  # Final results table
  max_result <- data.table(conta_version = packageVersion("conta"),
                           sample = sample_id,
                           sex = sex,
                           format(result, digits = 5, trim = TRUE),
                           y_count = round(chr_y_stats$count, 4),
                           y_norm_count = round(chr_y_stats$normalized_count, 6),
                           y_fraction_covered = round(chr_y_stats$fraction_covered, 4),
                           x_count = round(chr_x_stats$count, 4),
                           x_norm_count = round(chr_x_stats$normalized_count, 6),
                           x_fraction_covered = round(chr_x_stats$fraction_covered, 4),
                           x_het_snps = x_het_snps,
                           pregnancy = pregnancy,
                           excluded_regions = bin_stats[loh == TRUE, .N],
                           error_rate = round(mean(EE$er), 7),
                           round(t(data.frame(EE$er,
                                 row.names = rownames(EE))), digits = 7))
  # Remove sample ID column if not provided or empty.
  if (is.na(sample_id) || is.null(sample_id) || sample_id == "") {
    max_result$sample <- NULL
  }

  # Write conta results
  write.table(max_result, file = out_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  # Write genotype files for all SNPs
  write_gt_file(dat, max_result, blackswan, outlier_frac,
                file.path(save_dir, paste(filename_prefix, "gt.tsv",
                                          sep = ".")))

  # Write genotype files for LOH-excluded SNPs
  write_gt_file(dat_loh, max_result, blackswan, outlier_frac,
                file.path(save_dir, paste(filename_prefix, "gt.loh.tsv",
                                          sep = ".")))
}
