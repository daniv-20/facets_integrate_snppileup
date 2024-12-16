#' SNP Pileup Analysis
#'
#' Runs the SNP pileup analysis with the given arguments.
#'
#' @param vcffile Input file path (e.g., VCF file).
#' @param output Output file path.
#' @param bamfiles Vector of bamfile paths.
#' @param  count_orphans Do not discard anomalous read paris.
#' @param max_depth Sets the maximum depth. Default is 4000.
#' @param gzipped Compresses the output file with BGZF. If false, output file is CSV.
#' @param progress If true, shows a progress bar. WARNING: requires additional time to calculate number of SNPs, and will take longer than normal.
#' @param psuedo_snps Every n positions, if there is no SNP, insert a blank record with the total count at the position.
#' @param min_map_quality Sets the minimum threshold for mapping quality. The default is 0.
#' @param min_base_quality Sets the minimum threshold for base quality. Default is 0.
#' @param min_read_counts Comma separated list of minimum read counts for position to be output. Default is 0.
#' @param verbose If true, shows detailed messages.
#' @param ignore_overlaps If true, disables read-pair overlap detection.
#' @param debug_mode If true, shows many detailed messages.
#' @return None. Results are written to the output file.
#' @export
run_snp_pileup <- function(
  vcffile, ## vcf file path
  output, ## output file path
  bamfiles, ## vector of file paths
  count_orphans = FALSE,
  gzipped = FALSE,
  ignore_overlaps = FALSE,
  min_base_quality = 0,
  min_map_quality = 0,
  min_read_counts = 0,
  max_depth = 4000,
  progress = FALSE,
  psuedo_snps = 0,
  verbose = FALSE,
  debug_mode = FALSE) {
if (!file.exists(vcffile)) stop("Input file does not exist: ", vcffile)
if (file.exists(output)) stop("Output file already exists: ", output)
if (!gzipped) {
  if (!grepl("\\.csv$", output)) {
    output <- paste0(output, ".csv")
    warning("Output file does not end in .csv. Adding .csv. \n New output file: ", output)
  }
}
for (bamfile in bamfiles) {
  if (!file.exists(bamfile)) stop("Bam file does not exist: ", bamfile)
}
if (class(min_read_counts) != "numeric") {
  stop(paste0("Error: Min_read_counts is ", class(min_read_counts), ". Min_read_counts must be numeric."))
}
user_values <- list(
  count_orphans = count_orphans,
  gzipped = gzipped,
  ignore_overlaps = ignore_overlaps,
  min_base_quality = min_base_quality,
  min_map_quality = min_map_quality,
  min_read_counts = min_read_counts,
  max_depth = max_depth,
  progress = progress,
  psuedo_snps = psuedo_snps,
  verbose = verbose,
  debug_mode = debug_mode
)

# Categorize variables by expected type
numeric_vars <- c("min_base_quality", "min_map_quality", "min_read_counts", "max_depth", "psuedo_snps")
boolean_vars <- c("count_orphans", "gzipped", "ignore_overlaps", "progress", "verbose", "debug_mode")

# Loop to validate user values
for (var_name in names(user_values)) {
  value <- user_values[[var_name]]

  if (var_name %in% numeric_vars && !is.numeric(value)) {
    stop(sprintf("Error: %s is not numeric (type: %s).\n", var_name, class(value)))
  }

  if (var_name %in% boolean_vars && !is.logical(value)) {
    stop(sprintf("Error: %s is not boolean (type: %s).\n", var_name, class(value)))
  }
}


qual_args <- c(
  "-Q", min_base_quality,
  "-q", min_map_quality,
  "-r", min_read_counts,
  "-d", max_depth,
  "-P", psuedo_snps
)

if (count_orphans) {
  qual_args <- c(qual_args, "-A")
}
if (gzipped) {
  qual_args <- c(qual_args, "-g")
}
if (ignore_overlaps) {
  qual_args <- c(qual_args, "-x")
}
if (progress) {
  qual_args <- c(qual_args, "-p")
}
if (verbose) {
  qual_args <- c(qual_args, "-v")
}


args <- c(qual_args, vcffile, output, bamfiles)
.Call("_facets_run_snp_pileup_logic", args)
}
