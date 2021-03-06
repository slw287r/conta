# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' @useDynLib conta
#' @importFrom Rcpp sourceCpp
#' @import data.table parallel ggplot2
NULL

#' Return nucleotide bases
#' @param context_mode whether to get bases from 3-base context
#'
#' @export
get_bases <- function(context_mode = FALSE) {
  if (!context_mode)
    return(c("A", "T", "G", "C"))
  else {
    bases <- expand.grid(get_bases(), get_bases(), get_bases())
    return(paste(bases$Var1, bases$Var2, bases$Var3, sep = ""))
  }
}

#' Test whether the two contexts are compatible with a one base variance in the
#' middle base
#' @param context1 context of the first case
#' @param context2 context of the second case
#'
#' @export
is_compatible_context <- function(context1, context2) {
  if (substring(context1, 1, 1) == substring(context2, 1, 1) &&
      substring(context1, 3, 3) == substring(context2, 3, 3))
    return(TRUE)
  else
    return(FALSE)
}

#' Return initial test range
#'
#' @export
get_initial_range <- function() {
  return(c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3,
           5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1))
}

#' Return an empty results table to output when quitting early
#' @param sample_id sample ID to use in empty result.
#' @importFrom utils packageVersion
#' @export
empty_result <- function(sample_id = "") {
  vals <- list(conta_version = as.character(packageVersion("conta")),
               sample = sample_id)
  na_names <- c("sex", "conta_call", "cf", "sum_log_lr", "avg_log_lr", "hom_snps",
                "depth", "pos_lr_all", "pos_lr_x", "pos_lr_chr_cv", "y_count",
                "y_norm_count", "y_fraction_covered", "x_count", "x_norm_count",
                "x_fraction_covered", "x_het_snps", "pregnancy",
                "excluded_regions", "error_rate", "T>A", "G>A", "C>A", "A>T",
                "G>T", "C>T", "A>G", "T>G", "C>G", "A>C", "T>C", "G>C")
  na_vals <- rep(NA_character_, length(na_names))
  names(na_vals) <- na_names
  empty_result <- do.call(data.table, c(vals, na_vals))
  if (sample_id == "") {
    empty_result[, sample := NULL]
  }
  return(empty_result)
}

#' Fail if certain conditions are not met
#'
#' @param dat data.table with SNP counts and genotypes
#'
#' @export
fail_test <- function(dat) {
  if ( is.null(dat) || nrow(dat) == 0 || is.null(dat$gt) ||
       (dat %>% dplyr::filter(gt == "0/0") %>% nrow()) == 0 ||
       (dat %>% dplyr::filter(gt == "1/1") %>% nrow()) == 0 ||
       (dat %>% dplyr::filter(gt == "0/1") %>% nrow()) == 0) {
    msg <- "Either no SNPs passed filters or not all genotypes were called."
    if (interactive()) {
      stop(msg)
    } else {
      message(msg)
      quit(save = "no", status = 2, runLast = FALSE)
    }
  }
}

#' Return expected contamination fraction
#'
#' Average ratio of heterozygotes to homzoygotes for contaminant SNPs between
#' 100 individual samples was determined to be 80%. For these SNPs, only one
#' allele (half the fraction) contributes to the contamination. For the
#' remaining 20% of SNPs which are homozygotes, both alleles contribute to the
#' contamination. This 80% number is an approximation and may be changed based
#' on prior information.
#'
#' @param cf numeric contamination fraction
#' @param het_rate ratio of heterozygote to homozygote contaminated alleles
#' @return expected contamination fraction
#'
#' @export
get_exp_cf <- function(cf, het_rate = 0.8) {
  return(het_rate * (cf / 2) + (1 - het_rate) * cf)
}

#' Simulate contamination
#'
#' It works by adding in each allele proportional to specified
#' contamination fraction and maf. So it actually simulates a
#' sample randomly ignoring linkage disequilibrium.
#'
#' @param dat data.frame containing counts and metrics per SNP
#' @param cf numeric contamination fraction
#'
#' @return data.frame containing counts and metrics per SNP
#'
#' @importFrom stats rpois runif
#' @export
sim_conta <- function(dat, cf) {

  if (cf == 0) return(dat)

  # Calculate allele depths for each chromosome
  allele1 <- ifelse(runif(nrow(dat)) > dat$maf, dat$major, dat$minor)
  allele2 <- ifelse(runif(nrow(dat)) > dat$maf, dat$major, dat$minor)

  # Add in contamination
  bases <- get_bases()
  for (b in bases) {
    locs1 <- which(allele1 == b)
    allele1_probs <- c(dat[locs1, "depth"] * cf / 2)
    allele1_draws <- rpois(length(allele1_probs), allele1_probs)
    dat[locs1, b] <- dat[locs1, b] + allele1_draws

    locs2 <- which(allele2 == b)
    allele2_probs <- c(dat[locs2, "depth"] * cf / 2)
    allele2_draws <- rpois(length(allele2_probs), allele2_probs)
    dat[locs2, b] <- dat[locs2, b] + allele2_draws
  }

  # Recalculate depth, ratio and counts
  dat <- ratio_and_counts(dat)

  return(dat)
}

#' Simulate LOH and/or contamination
#'
#' This simulation returns a simple data structure with metrics required to
#' test LOH functions. It works by simulating 4 chromosomes (2 for host, and
#' 2 for contaminant). It also simulates the ratio of the host chromosomes which
#' is analogous to loss of heterozygosity.
#'
#' @param n number of SNPs to simulate
#' @param min_maf minimum maf to simulate (and maximum is 1 - maf)
#' @param dp_min min epth to simulate
#' @param dp_max max depth to simulate
#' @param er_min minimum error rate to simulate
#' @param er_max maximum error rate to simulate
#' @param delta LoH deviation from heterozygosity for host
#' @param max_delta maximum allowed delta, if larger value set to this
#' @param het_mean_binomial_size higher value has het_means tighter across SNPs
#' @param delta_cont LoH for contaminant
#' @param default_het_mean default expected mean for beta binomial
#' @param default_het_size default expected size for beta binomial
#' @param alpha contamination level
#' @param seed numeric set seed for simulation
#'
#' @return data.frame containing counts and metrics per SNP
#'
#' @export
simulate_loh_conta <- function(n, min_maf, dp_min, dp_max, er_min, er_max,
                               delta, alpha, delta_cont = 0,
                               default_het_mean = 0.43,
                               max_delta = 0.49,
                               het_mean_binomial_size = 50,
                               default_het_size = 200,
                               seed = 1359) {

  # Delta should not exceed 0.49 for this simulation
  if (delta > max_delta) {
    delta <- max_delta
  }

  set.seed(seed)

  # Simulate maf for each SNP (at range as provied)
  maf <- runif(n, min = min_maf, max = 1 - min_maf)

  # Simulate error rate for each SNP randomly at specified range
  er <- runif(n, min = er_min, max = er_max)

  # Simulate depth for each SNP randomly around specified target
  dp <- rpois(n, runif(n, dp_min, dp_max))

  # Simulate mean expected allele frequency and its variance based on alternative
  # beta binomial formulation, size sets the shape of binomial
  het_mean <- rbinom(n, size = het_mean_binomial_size,
                    prob = default_het_mean) / het_mean_binomial_size

  # Simulate alleles for host, 1 = reference, 0 = minor allele
  a1 <- runif(n) > maf
  a2 <- runif(n) > maf

  # Simulate alleles for contaminant, 1 = reference, 0 = minor allele
  a3 <- runif(n) > maf
  a4 <- runif(n) > maf

  # Simulate depth for host
  dp_host <- rpois(n, dp * (1 - alpha))
  dp_cont <- rpois(n, dp * alpha)

  # Draw minor allele and reference depths for host
  ad1 <- rpois(n, dp_host * ifelse(a1, er, 1 - er) * nloh(het_mean, -delta))
  ref1 <- rpois(n, dp_host * ifelse(a1, 1 - er, er) * nloh(het_mean, -delta))
  ad2 <- rpois(n, dp_host * ifelse(a2, er, 1 - er) * nloh(het_mean, delta))
  ref2 <- rpois(n, dp_host * ifelse(a2, 1 - er, er) * nloh(het_mean, delta))

  # Draw minor allele and reference depths for contaminant
  ad3 <- rpois(n, dp_cont * ifelse(a3, er, 1 - er) * nloh(het_mean, -delta_cont))
  ref3 <- rpois(n, dp_cont * ifelse(a3, 1 - er, er) * nloh(het_mean, -delta_cont))
  ad4 <- rpois(n, dp_cont * ifelse(a4, er, 1 - er) * nloh(het_mean, delta_cont))
  ref4 <- rpois(n, dp_cont * ifelse(a4, 1 - er, er) * nloh(het_mean, delta_cont))

  # Make a mix of the two simulated samples
  ad <- ad1 + ad2 + ad3 + ad4
  dp <- ad + ref1 + ref2 + ref3 + ref4

  return(data.frame(depth = dp, minor_count = ad, maf = maf, er = er,
                    het_mean = het_mean, het_sum = default_het_size))
}

#' Load conta files
#'
#' For a given record, mt and snps, read vfn or gt file, and return it.
#'
#' @param conta_loc location of conta vfn or gt file
#' @param snps targeted dbsnp data.table
#' @return data.table with conta genotypes
#'
#' @export
load_conta_file <- function(conta_loc, snps = NULL) {

  if ((!is.null(snps)) & !is.na(conta_loc) &
      (file.exists(conta_loc) ||
       (requireNamespace("grails3r") && grails3r::s3_file_exists(conta_loc)) ||
       as.logical(do.call(aws.s3::head_object,
                          c(file, aws.signature::locate_credentials()))))) {

    # Read targeted genotypes and rsid from vfn file format
    conta_loc_dt <- read_data_table(conta_loc, sep = ",", showProgress = FALSE)

    # Add rsid if missing.
    if (!is.null(conta_loc_dt))
      conta_loc_dt[, rsid := snps$ID]

  } else if (!is.na(conta_loc) &
             (file.exists(conta_loc) ||
              (requireNamespace("grails3r") &&
               grails3r::s3_file_exists(conta_loc)) ||
              as.logical(do.call(aws.s3::head_object,
                                 c(file, locate_credentials()))))) {

    # Read conta genotypes and rsid from gt file format
    conta_loc_dt <- read_data_table(conta_loc, showProgress = FALSE)

    # Vf is the variant allele frequency, needs to be set of wgs files
    conta_loc_dt[, vf := vr / dp]

  }
  else {

    conta_loc_dt <- NULL

  }

  return(conta_loc_dt)
}

#' Read conta genotype file
#'
#' @param gt_file input gt file
#' @return data frame containing gt file contents transformed in a way
#'     useful to downstream methods
#' @export
read_gt_file <- function(gt_file) {

  dat <- read.table(gt_file, header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(dplyr::everything(), chrom = chr, depth = dp) %>%
    factorize_chroms() %>%
    arrange(chrom, pos) %>%
    add_chunks()

  return(dat)
}

#' Calculate concordance between two samples' genotypes. Annotate the
#' type of genotype swaps per SNP.
#'
#'
#' @param dt1 data.table of genotype set 1
#' @param dt2 data.table of genotype set 2
#'
#' @importFrom dplyr left_join group_by count bind_rows mutate
#' @importFrom tidyr spread
#' @export
genotype_concordance <- function(dt1, dt2) {
  # Merge the two gt data tables to put each SNP in a single row.
  # Annotate the type of genotype swap per SNP.
  merged_gt_df <- dplyr::inner_join(dt1, dt2, by = "rsid") %>%
    dplyr::mutate(swap_type = ifelse(gt.x == gt.y,
                                     "No_Change",
                                     ifelse(gt.x == "0/0" & gt.y == "1/1",
                                            "Hom_Ref_to_Hom_Alt",
                                            ifelse(gt.x == "1/1" & gt.y == "0/0",
                                                   "Hom_Alt_to_Hom_Ref",
                                                   "Het_Change"))))

  # Count the number of SNPs per swap type
  swap_type_counts <- merged_gt_df %>%
    dplyr::group_by(swap_type) %>%
    dplyr::count() %>%
    tidyr::spread(swap_type, n)

  # Calculate total number of intersecting snps between samples
  total_snps <- rowSums(swap_type_counts)

  # Set blank dataframe, bind 'swap_type_counts' to blank template, set NAs to 0
  swap_template <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                            c("Hom_Ref_to_Hom_Alt", "Hom_Alt_to_Hom_Ref",
                              "Het_Change", "No_Change"))
  swap_metrics <- dplyr::bind_rows(swap_template, swap_type_counts)
  swap_metrics[is.na(swap_metrics)] <- 0

  # Calculate gt concordance which is the fraction of matching genotypes
  swap_metrics <- swap_metrics %>%
    dplyr::mutate(Concordance =
                    1-(Hom_Ref_to_Hom_Alt + Hom_Alt_to_Hom_Ref + Het_Change)/total_snps)
  return(swap_metrics)
}

#' Get folders under a specified s3 directory
#'
#' Returns the set of folders under a given s3 path. get_bucket_df normally
#' returns all folders and files recursively, this method parses get_bucket_df
#' output to return only folders directly under the specified dir.
#'
#' @param s3_path to display the contents
#' @export
get_s3_folders <- function(s3_path) {

  # Stop if this is not an s3 path
  if (!startsWith(s3_path, "s3://"))
    stop("get_s3_ls() requires a valid s3 path")

  # Retireve bucket text
  bucket_text <- strsplit(s3_path, "/")[[1]][3]

  # Retrieve prefix text
  pre_text <- paste(strsplit(s3_path, "/")[[1]][c(-1, -2, -3)],
                    sep = "/", collapse = "/")

  # Get all files (and folders) under specified bucket/prefix
  bc_all_files <- aws.s3::get_bucket_df(bucket = bucket_text,
                                        prefix = pre_text, max = Inf)

  # Number of folders under prefix
  pl <- length(strsplit(pre_text, "/")[[1]])

  # Get all paths directly under prefix (omit files and folder in deeper levels)
  paths <- unique(sapply(strsplit(bc_all_files$Key, "/"),
                  function(x) {
                    if (length(x) > (pl + 1)) x[[pl + 1]]
                    }))

  # Add a slash in the end before returning to denote a folder
  return(paste(unlist(paths), "/", sep = ""))
}

#' Set numeric equivalent of chromosomes for sorting purposes
#'
#' @section TODO: Handle chromosomes outside X, Y and M/MT
#'
#' @param dat data.table to add numeric chromosome
#' @export
set_numeric_chrs <- function(dat) {

  suppressMessages(require(dplyr))

  datt <- dat %>%
    mutate(
      chrom_int = substring(chrom, 4)) %>%
    mutate(
      chrom_int = dplyr::case_when(
        chrom_int == "X" ~ 23,
        chrom_int == "Y" ~ 24,
        chrom_int == "M" ~ 25,
        chrom_int == "MT" ~ 25,
        TRUE ~ as.numeric(chrom_int))
    )
  return(data.table(datt))
}

#' Set factor equivalent of chromosomes for sorting purposes
#' This function also removes Chr or chr prefixes.
#'
#' @param dat data.table to add factorized chromosome
#' @export
factorize_chroms <- function(dat) {

  dat$chrom <- gsub("chr", "", dat$chrom, ignore.case = TRUE)
  dat$chrom <- factor(dat$chrom, levels = as.character(c(1:22, "X", "Y")))
  return(dat)
}

#' Make a sex call for the host
#'
#' @section TODO: Add a joint caller for chrX and chrY
#' @section TODO: Consider adding chrX heterozygous SNP counts to the model
#' @section TODO: Report sex chromosome aneuploidies
#'
#' @param chr_y_stats data frame containing read counts for chrY
#' @param chr_y_male_threshold threshold to call male on chrY data
#' @param result data frame with conta results
#' @param chr_x_stats data frame containing read counts for chrX
#' @param chr_x_female_threshold threshold to call female on chrX reads
#' @param ftm female threshold multiplier minimum fraction of chrY reads to call
#'    This multiplier defines a gray zone for the sex call.
#' @param cf_threshold contamination fraction to skip sex call in case there
#'     is a host switch (used only when there are enough chrY reads). i.e if
#'     both host and contaminant are females, then we can still make a sex call.
#'
#' @return sex estimated sex call, i.e. FEMALE, MALE or NO CALL
#' @export
get_sex_call <- function(chr_y_stats, chr_y_male_threshold, result,
                         chr_x_stats = NA, chr_x_female_threshold = NA,
                         ftm = 0.75, cf_threshold = 0.2) {

  # Base logic:
  # NO CALL if conta call and contamination fraction is high
  # MALE if normalized chrY count is more than the threshold
  # FEMALE if chrY count is less than specified fraction of the threshold
  # NO CALL otherwise (which is the grey zone, could be lower percentages of
  # contamination or mosaic aneuploidies)
  sex <- dplyr::case_when(
    result$conta_call & (result$cf >= cf_threshold) &
      chr_y_stats$normalized_count > (ftm * chr_y_male_threshold) ~ "No_call",
    chr_y_stats$normalized_count >= chr_y_male_threshold ~ "Male",
    chr_y_stats$normalized_count <= (ftm * chr_y_male_threshold) ~ "Female",
    TRUE ~ "No_call"
  )

  return(sex)
}

#' Make a pregnancy call if there is contamination and host is female
#'
#' Pregnancy may be defined as a female host with a contamination call that has
#' less signal on the X chromosome. In general if the child is a male, there
#' should be no contamination signal on the X chromosome, presence of Y
#' fragments, and an overall contamination signal. For a female child, there
#' should be half the expected contamination signal on the X chromosome, and
#' no presence of Y fragments.
#'
#' @param result data frame with conta results
#' @param sex estimated sex call for the host
#' @param chr_y_stats data frame containing read counts for chrY
#' @param chr_y_male_threshold threshold to call male on chrY data
#' @param chr_y_female_baseline expected chrY normalized reads for a female
#' @param chr_y_male_baseline expected chrY normalized reads for a male
#' @param conta_variance expected variance in case of male contamination
#' @param min_pregnancy_level minimum contamination fraction for pregnancy
#' @param mpm male pregnancy multiplier minimum fraction conta likelihood
#'     expected in the case of a male pregnancy
#' @param fpm female_pregnancy_multiplier minimum fraction conta likelihood
#'     expected in the case of a female pregnancy
#' @return pregnancy estimated pregnancy call, i.e. FEMALE, MALE or NA
#'     (NA is reported in the case of no contamination or host sex not female)
#' @export
get_pregnancy_call <- function(result, sex, chr_y_stats, chr_y_male_threshold,
                               chr_y_female_baseline = 6 * 10^-6,
                               chr_y_male_baseline = 0.001,
                               conta_variance = 0.75,
                               min_pregnancy_level = 0.2,
                               mpm = 0.2, fpm = 0.6) {

  no_X_conta <- result$pos_lr_x <= (mpm * result$pos_lr_all)
  one_X_conta <- result$pos_lr_x <= (fpm * result$pos_lr_all)
  observed_Y <- chr_y_stats$normalized_count >= chr_y_female_baseline +
    conta_variance * result$cf * chr_y_male_baseline
  pregnancy <- dplyr::case_when(
    result$conta_call & result$cf >= min_pregnancy_level & (sex != "MALE") &
      no_X_conta & observed_Y ~ "Male",
    result$conta_call & result$cf >= min_pregnancy_level & (sex != "MALE") &
      one_X_conta & !observed_Y ~ "Female",
    TRUE ~ "NA"
  )

  return(pregnancy)
}

#' Return snps with the following filters:
#'    - loh call absent or false
#'    - valid rsid
#'    - maf >= maf_threshold
#'    - dp >= depth_threshold
#'
#' @param dt data frame with conta snps
#' @param maf_threshold
#' @param depth_threshold
#' @return filtered data table
#  @export
filter_gt_file <- function(dt, maf_threshold, depth_threshold) {
  # TODO(jrosenfield): remove this fix once rsids are formatted properly upstream
  result <- dplyr::mutate(dt, rsid = as.character(rsid)) %>%
    dplyr::mutate(rsid = ifelse(startsWith(rsid, "rs"), rsid, paste("rs", rsid, sep="")))
  result <- dplyr::filter(result, rsid!="rs0", maf>=maf_threshold, dp>=depth_threshold)
  if ("loh" %in% colnames(result)) {
    result <- dplyr::filter(result, !loh)
  }
  return(result)
}
