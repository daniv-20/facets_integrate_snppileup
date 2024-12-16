#include <Rcpp.h>
#include "../inst/include/snp-pileup-rev.h"
#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib> // For exit()
#include <cstdio>  // For printf, FILE
#include <sstream>
#include <ctime>

using namespace Rcpp;
// #include <htslib/synced_bcf_reader.h>
// #include <htslib/sam.h>
// #include <htslib/bgzf.h>
// #include <htslib/vcf.h>
// #include <htslib/hts.h>

// Define the arguments structure
#ifndef SNP_PILEUP_REV_H_DEFINED
#define SNP_PILEUP_REV_H_DEFINED
struct arguments
{
  std::vector<std::string> args;
  bool count_orphans = false;
  bool gzipped = false;
  bool ignore_overlaps = false;
  int min_base_quality = 0;
  int min_map_quality = 0;
  std::vector<int> min_read_counts;
  int max_depth = 4000;
  int rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  void (*outFunc)(arguments, std::string, FILE *) = nullptr; // set default to print output
  bool progress = false;
  int pseudo_snps = 0;
  bool verbose = false;
  BGZF *gzippedPointer = nullptr; // Added
  bool debug_mode = false;        // added debug mode for my own sanity
};
#endif

// Function to read and parse the header
bool read_vcf_header(BGZF *bgzf, std::vector<std::string> &columns)
{
  kstring_t line = {0, 0, NULL};
  bool header_found = false;
  // size_t len = 0;

  Rcpp::Rcout << "Debug: Entering read_vcf_header function" << std::endl;

  while (bgzf_getline(bgzf, '\n', &line) >= 0)
  {
    // Check if it's a header line
    Rcpp::Rcout << "Debug: Read line: " << line.s << std::endl;
    if (line.s[0] == '#')
    {
      // Check for the #CHROM line (last header line)
      Rcpp::Rcout << "Debug: Line is a header line" << std::endl;
      if (strncmp(line.s, "#CHROM", 6) == 0)
      {
        // Parse the #CHROM line into columns
        std::string header_line(line.s);
        Rcpp::Rcout << "Debug: Found #CHROM line" << std::endl;
        // Debug: Print the header line
        Rcpp::Rcout << "Debug: Header line content: " << header_line << std::endl;

        std::istringstream iss(header_line);
        std::string column;

        while (std::getline(iss, column, '\t'))
        {
          columns.push_back(column);
          Rcpp::Rcout << "Debug: Parsed column: " << column << std::endl;
        }

        if (columns.size() != 4)
        {
          Rcpp::Rcerr << "Error: Expected 4 columns, but parsed " << columns.size() << " columns." << std::endl;
          return false;
        }
        header_found = true;
        break;

        // Debug: Check the number of columns parsed

        free(line.s); // Clean up
        return true;  // Header successfully read
      }
    }
    else
    {
      break; // Stop at the first non-header line
    }
  }

  if (!header_found)
  {
    // File does not have a #CHROM line, check the first non-header line for format
    std::string first_line(line.s);
    std::istringstream iss(first_line);
    std::string column;
    int column_count = 0;

    while (std::getline(iss, column, '\t'))
    {
      columns.push_back(column);
      column_count++;
    }

    if (column_count != 4)
    {
      Rcpp::Rcerr << "Error: File does not have a header and the first line is not in 4-column format." << std::endl;
      free(line.s);
      return false;
    }

    // Warn and continue
    Rcpp::Rcout << "Warning: File does not have a header but appears to be in correct 4-column format. Continuing." << std::endl;
    free(line.s); // Clean up
    return true;  // Header successfully read
  }

  free(line.s);
  return false; // No valid #CHROM line found
}

// Parse arguments from R inputs
// Argument parser
void parse_arguments(const std::vector<std::string> &input_args, arguments &args)
{
  // Rcpp::Rcout << "Debug: Entered parse args" << std::endl;
  if (input_args.size() < 2)
  {
    Rcpp::stop("Usage: <vcf file> <output file> <sequence files...>");
  }

  for (size_t i = 0; i < input_args.size(); ++i)
  {
    const std::string &arg = input_args[i];
    if (arg == "--count-orphans" || arg == "-A")
    {
      args.count_orphans = true;
    }
    else if (arg == "--gzip" || arg == "-g")
    {
      args.gzipped = true;
    }
    else if (arg == "--ignore-overlaps" || arg == "-x")
    {
      args.ignore_overlaps = true;
    }
    else if (arg == "--min-base-quality" || arg == "-Q")
    {
      if (++i < input_args.size())
      {
        // Rcpp::Rcout << "Debug: mbc" << input_args[i] << std::endl;
        args.min_base_quality = std::stoi(input_args[i]);
      }
      else
      {
        Rcpp::stop("Missing value for --min-base-quality");
      }
    }
    else if (arg == "--min-map-quality" || arg == "-q")
    {
      if (++i < input_args.size())
      {
        // Rcpp::Rcout << "Debug: mmc" << input_args[i] << std::endl;
        args.min_map_quality = std::stoi(input_args[i]);
      }
      else
      {
        Rcpp::stop("Missing value for --min-map-quality");
      }
    }
    else if (arg == "--min-read-counts" || arg == "-r")
    {
      if (++i < input_args.size())
      {
        std::stringstream ss(input_args[i]);
        std::string token;
        while (std::getline(ss, token, ','))
        {
          args.min_read_counts.push_back(std::stoi(token));
        }
      }
      else
      {
        Rcpp::stop("Missing value for --min-read-counts");
      }
    }
    else if (arg == "--max-depth" || arg == "-d")
    {
      if (++i < input_args.size())
      {
        // Rcpp::Rcout << "Debug: md " << input_args[i] << std::endl;
        args.max_depth = std::stoi(input_args[i]);
      }
      else
      {
        Rcpp::stop("Missing value for --max-depth");
      }
    }
    else if (arg == "--progress" || arg == "-p")
    {
      args.progress = true;
    }
    else if (arg == "--pseudo-snps" || arg == "-P")
    {
      if (++i < input_args.size())
      {
        // Rcpp::Rcout << "Debug: ps" << input_args[i] << std::endl;
        args.pseudo_snps = std::stoi(input_args[i]);
      }
      else
      {
        Rcpp::stop("Missing value for --pseudo-snps");
      }
    }
    else if (arg == "--verbose" || arg == "-v")
    {
      args.verbose = true;
    }
    else if (arg == "--debug")
    {
      args.debug_mode = true;
      Rcpp::Rcout << "Debug: will print debug statements" << std::endl;
    }
    else
    {
      args.args.push_back(arg);
    }
  }

  if (args.args.size() < 2)
  {
    Rcpp::stop("Error: Must provide at least a VCF file and an output file.");
  }

  // Ensure min_read_counts is initialized properly
  size_t num_sequence_files = args.args.size() - 2; // Number of sequence files
  if (args.min_read_counts.empty())
  {
    Rcpp::Rcout << "Warning: No min-read-counts provided. Defaulting to 0 for all files." << std::endl;
    args.min_read_counts.assign(num_sequence_files, 0); // Default to 0 for all files
  }
  else if (args.min_read_counts.size() == 1)
  {
    Rcpp::Rcout << "Single min-read-count provided. Using it for all sequence files." << std::endl;
    args.min_read_counts.assign(num_sequence_files, args.min_read_counts[0]); // Use the single value for all files
  }
  else if (args.min_read_counts.size() < num_sequence_files)
  {
    Rcpp::Rcout << "Warning: Fewer min-read-counts provided than sequence files. Filling with defaults." << std::endl;
    size_t missing_count = num_sequence_files - args.min_read_counts.size();
    args.min_read_counts.insert(args.min_read_counts.end(), missing_count, 0); // Fill with 0
  }
  else if (args.min_read_counts.size() > num_sequence_files)
  {
    Rcpp::Rcout << "Warning: More min-read-counts provided than sequence files. Truncating excess values." << std::endl;
    args.min_read_counts.resize(num_sequence_files); // Truncate to match number of files
  }
}

// A helper function for conditional debugging
void debugPrint(const std::string &message, bool debug)
{
  if (debug)
  {
    Rcpp::Rcout << message << std::endl;
  }
}

// Prints plain output to a file
void print_output(arguments arguments, std::string str, FILE *fp)
{
  fputs(str.c_str(), fp);
}

// Writes gzipped output to a file
void gzip_output(arguments arguments, std::string str, FILE *fp)
{
  if (bgzf_write(arguments.gzippedPointer, str.c_str(), str.length()) < 0)
  {
    Rcpp::stop("Failed to write to file, terminating.");
  }
}

// Checks if a string ends with a certain substring
inline bool ends_with(std::string const &value, std::string const &ending)
{
  if (ending.size() > value.size())
    return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

static int mplp_func(void *data, bam1_t *b)
{
  // it seems to me that this function is run once every read, and lets you control if the read should be skipped or not
  // here is where you would put quality checks and things like that
  //	printf("hello\n");
  // char *ref;
  mplp_aux_t *ma = (mplp_aux_t *)data;
  int ret, skip = 0; // ref_len;
  do
  {
    // int has_ref;
    ret = sam_read1(ma->fp, ma->h, b);
    if (ret < 0)
      break;

    /*if (b->core.flag & READ_FAILED_QUALITY_CHECKS || b->core.flag & READ_SECONDARY_ALIGNMENT || b->core.flag & READ_PCR_OPTICAL_DUPLICATE || b->core.flag & READ_UNMAPPED) {
        skip = 1;
        continue;
    }*/

    if (ma->conf->rflag_filter && ma->conf->rflag_filter & b->core.flag)
    {
      skip = 1;
      continue;
    }

    if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP))
    { // exclude unmapped reads
      skip = 1;
      continue;
    }

    /*if (ma->conf->bed && ma->conf->all == 0) { // test overlap
        skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
        if (skip) continue;
    }*/

    skip = 0;
    /*if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
    if (has_ref && ma->conf->capQ_thres > 10) {
        int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
        if (q < 0) skip = 1;
        else if (b->core.qual > q) b->core.qual = q;
        } */
    if (b->core.qual < ma->conf->min_map_quality)
    {
      skip = 1;
    }
    else if (b->core.qual == 0)
    {
      skip = 1;
    }
    else if (!ma->conf->count_orphans && (b->core.flag & BAM_FPAIRED) && !(b->core.flag & BAM_FPROPER_PAIR))
    {
      skip = 1;
    }
  } while (skip);
  return ret;
}

static int vcf_chr_to_bam(char *chromosome, char **bam_chrs, int32_t n_targets)
{
  // try to remove chr prefix if it exists
  // Rcpp::Rcout << "Debug: remove chr prefix if exists" << std::endl;
  if (!strncmp(chromosome, "chr", 3))
  {
    chromosome += 3;
  }
  bool bam_chr_prefix = (!strncmp(bam_chrs[0], "chr", 3));
  for (int i = 0; i < n_targets; i++)
  {
    if (!bam_chr_prefix)
    {
      if (!strcmp(bam_chrs[i], chromosome))
      {
        return i;
      }
    }
    else
    {
      if (!strcmp(bam_chrs[i] + 3, chromosome))
      {
        return i;
      }
    }
  }
  return -1;
}

uint64_t get_snp_count(char *file)
{
  // Rcpp::Rcout << "Debug: entered get snp count" << std::endl;
  uint64_t count = 0;
  // Open BGZF compressed vcf file in read mode
  BGZF *bgzf = bgzf_open(file, "r");
  if (!bgzf)
  {
    std::cerr << "Failed to open BGZF file: " << file << std::endl;
    return -1; // return -1 if there is an error
  }

  // Buffer for reading lines
  kstring_t line = {0, 0, NULL}; // Initialize kstring_t
  // size_t len = 0;

  // read lines from BGZF file
  while (bgzf_getline(bgzf, '\n', &line) >= 0)
  { // reaD A SINGLE LINE FROM BGZF FILE UNTIL newline character, store line in provided buffer 'line'
    // ignore header lines (start with #)
    if (line.s[0] != '#')
    {
      count++;
    }
  }
  // clean up
  free(line.s);     // free the allocated line buffer
  bgzf_close(bgzf); // close the bgzf handle

  return count;
}

int program_main(arguments arguments)
{
  debugPrint("Debug: enter program main", arguments.debug_mode);

  clock_t start = clock();

  // initialize engine
  int i = 0;
  int n = arguments.args.size() - 2; // n is the number of files. currently hardcoded at 2
  bam_mplp_t iter;
  const bam_pileup1_t **plp = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t *));
  int *n_plp = (int *)calloc(n, sizeof(int));
  mplp_ref_t mp_ref = MPLP_REF_INIT;
  struct arguments *conf = &arguments;
  hts_verbose = 1;

  // Load VCF file
  debugPrint("Debug: About to load vcf file", arguments.debug_mode);

  // open file
  BGZF *bgzf = bgzf_open(arguments.args[0].c_str(), "r");
  if (!bgzf)
  {
    std::cerr << "Failed to open BGZF file: " << arguments.args[0] << std::endl;
    return 1;
  }

  // Read and validate the VCF header
  std::vector<std::string> header_columns;
  if (!read_vcf_header(bgzf, header_columns))
  {

    Rcpp::Rcout << "Attempted to read header and failed: " << read_vcf_header(bgzf, header_columns) << std::endl;

    std::cerr << "Invalid or missing VCF header in file: " << arguments.args[0] << std::endl;
    bgzf_close(bgzf);
    return 1;
  }

  debugPrint("Debug: Loaded vcf file", arguments.debug_mode);

  // calculate snp count for progress tracking
  uint64_t count = 0;
  if (arguments.progress)
  {
    std::cout << "Calculating SNP count..." << std::flush;
    count = get_snp_count(const_cast<char *>(arguments.args[0].c_str()));
    std::cout << "done." << std::endl;
  }

  debugPrint("Debug: Construct BAM Data to pass to pileup engine", arguments.debug_mode);
  // Initialize BAM file data
  mplp_aux_t **data = (mplp_aux_t **)calloc(n, sizeof(mplp_aux_t *));
  for (i = 0; i < n; ++i)
  {
    data[i] = (mplp_aux_t *)calloc(1, sizeof(mplp_aux_t));

    // open file
    hFILE *hfp = hopen(arguments.args[i + 2].c_str(), "r");
    if (!hfp)
    {
      std::cerr << "Failed to open BAM file: " << arguments.args[i + 2] << std::endl;
      return 1;
    }

    // ommitting detect file format because ... ???

    data[i]->fp = hts_hopen(hfp, arguments.args[i + 2].c_str(), "rb");
    if (!data[i]->fp)
    {
      // Error handling: Couldn't open the file
      std::cerr << "Couldn't open sequence file: " << arguments.args[i + 2]
                << " (" << strerror(errno) << ")" << std::endl;
      // hclose(hfp); // Clean up the file handle
      if (hclose(hfp) != 0)
      {
        // Handle the error, e.g., log an error message or throw an exception.
        std::cerr << "Error: Failed to close file handle." << std::endl;
      }
      return 1;
    }

    hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0);

    data[i]->conf = conf;
    data[i]->h = sam_hdr_read(data[i]->fp);
    data[i]->ref = &mp_ref;
  }

  // Read BAM header
  bam_hdr_t *hdr = data[0]->h;

  // Start pileup engine
  debugPrint("Debug: Start pileup engine", arguments.debug_mode);
  iter = bam_mplp_init(n, mplp_func, (void **)data);
  if (!arguments.ignore_overlaps)
  {
    bam_mplp_init_overlaps(iter);
  }
  bam_mplp_set_maxcnt(iter, arguments.max_depth);

  if (arguments.verbose)
  {
    printf("Max per-file depth set to %d.\n", arguments.max_depth);
  }

  // setup output file

  string fname = string(arguments.args[1]);
  if (arguments.gzipped)
  {
    if (!ends_with(fname, ".gz"))
    {
      fname += ".gz";
    }
    arguments.outFunc = gzip_output;
  }
  else
  {
    arguments.outFunc = print_output; // Default to plain text output
  }
  // check if output exists
  {
    std::ostringstream oss;
    oss << "Debug: Checking file existence for: " << fname;
    debugPrint(oss.str(), arguments.debug_mode);
  }

  // Rcpp::Rcout << "Debug: check if output exists." << std::endl;
  FILE *test_output = fopen(fname.c_str(), "r");
  if (test_output)
  {
    std::cerr << "Output file already exists: "
              << arguments.args[1]
              << std::endl;
    fclose(test_output);
    return 1;
  }
  debugPrint("Debug: Ready to open output", arguments.debug_mode);
  // Rcpp::Rcout << "Debug: checked output." << std::endl;
  //  DON'T CLOSE test_output HERE because if you are here, test_output is null
  FILE *output_file = NULL;
  if (arguments.gzipped)
  {
    arguments.gzippedPointer = bgzf_open(fname.c_str(), "w+");
    if (!arguments.gzippedPointer)
    {
      Rcpp::Rcerr << "Error: Failed to open gzipped output file." << std::endl;
      return 1;
    }
    debugPrint("Debug: opened gzipped file", arguments.debug_mode);
  }
  else
  {

    // open output file
    output_file = fopen(fname.c_str(), "w+");
    debugPrint("Debug: opened output", arguments.debug_mode);
    // Rcpp::Rcout << "Debug: opened output." << std::endl;

    if (!output_file)
    {
      printf("Failed to open output file for writing: %s\n", strerror(errno));
      Rcpp::Rcerr << "Error: Failed to open output file for writing: " << strerror(errno) << std::endl;
      debugPrint("Debug: output file could not be opened", arguments.debug_mode);
      return 1;
    }
    else
    {
      debugPrint("Debug: File opened successfully", arguments.debug_mode);
    }
  }

  // Rcpp::Rcout << "Debug: File opened successfully." << std::endl;
  ostringstream output;

  // output header to file
  output << "Chromosome,Position,Ref,Alt";
  for (int i = 0; i < n; ++i)
  {
    output << ",File" << (i + 1) << "R";
    output << ",File" << (i + 1) << "A";
    output << ",File" << (i + 1) << "E";
    output << ",File" << (i + 1) << "D";
  }
  output << "\n";
  // Debug header content
  // Rcpp::Rcout << "Debug: Header generated: " << output.str() << std::endl;
  {
    std::ostringstream oss;
    oss << "Debug: Header generated: " << output.str();
    debugPrint(oss.str(), arguments.debug_mode);
  }

  // Write header to file -> does this need to be different for gzip???
  std::string header = output.str();

  // check on outfunc
  if (arguments.outFunc != nullptr)
  {
    (arguments.outFunc)(arguments, output.str(), output_file);
  }
  else
  {
    Rcpp::Rcerr << "Error: outFunc is not initialized. No output can be written" << std::endl;
  }

  // Process VCF records
  debugPrint("Debug: go thru it", arguments.debug_mode);
  int ret;
  int tid, pos;
  vector<file_info> f_info = vector<file_info>();
  bool first = true;
  kstring_t line = {0, 0, NULL};
  uint64_t current_count = 0;
  float last_progress = 100.0;
  bool have_snpped = false;
  debugPrint("Debug: ready to loop", arguments.debug_mode);
  while (bgzf_getline(bgzf, '\n', &line) >= 0)
  {
    // Skip header lines
    // debugPrint("Debug: in the loop, no parse yet", arguments.debug_mode);
    if (line.s[0] == '#')
      continue;

    // Parse the record
    std::istringstream iss(line.s);
    std::vector<std::string> fields;
    std::string field;

    // debugPrint("Debug: did parse", arguments.debug_mode);

    while (std::getline(iss, field, '\t'))
    {
      fields.push_back(field);
      // debugPrint("Debug: in field while", arguments.debug_mode);
      // debugPrint(field, arguments.debug_mode);
    }

    // Validate record
    if (fields.size() < 4)
    {
      std::cerr << "Malformed VCF record: " << line.s << std::endl;
      continue;
    }

    // Extract fields
    // debugPrint("Debug: extracting fields", arguments.debug_mode);

    //         {
    //     std::ostringstream oss;
    //     oss << "Debug: Fields: ";

    //     // Check if fields is a container
    //     for (const auto &field : fields) {
    //         oss << field << " ";
    //     }

    //     debugPrint(oss.str(), arguments.debug_mode);
    // }

    std::string chrom = fields[0];
    // debugPrint(chrom, arguments.debug_mode);
    // Rcpp::Rcout << "Debug int vcf_pos: " << fields[1] << std::endl;
    int vcf_pos = std::stoi(fields[1]);
    std::string ref = fields[2];
    // debugPrint(ref, arguments.debug_mode);
    std::string alt = fields[3];
    // debugPrint(alt, arguments.debug_mode);

    // debugPrint("Debug: extracted fields", arguments.debug_mode);

    // Filter for SNPs (single nucleotide variants only)
    if (ref.size() != 1 || alt.size() != 1)
    {
      debugPrint("Skipping non-SNP record: " + std::string(line.s), arguments.debug_mode);
      continue;
    }
    // Map VCF chromosome and position to BAM target indices
    char chrom_buffer[chrom.size() + 1];
    std::strcpy(chrom_buffer, chrom.c_str()); // copy string into mutable buffer

    int vcf_tid = vcf_chr_to_bam(chrom_buffer, hdr->target_name, hdr->n_targets);
    current_count++;

    // track progress
    float progress = ((float)current_count / count);
    if (arguments.progress && (int)(last_progress * 100) != (int)(progress * 100))
    {
      int barWidth = 70;
      std::cout << "[";
      int pos = barWidth * progress;
      for (int i = 0; i < barWidth; ++i)
      {
        if (i < pos)
          std::cout << "=";
        else if (i == pos)
          std::cout << ">";
        else
          std::cout << " ";
      }
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();
      last_progress = progress;
    }

    // plp processing loop
    while ((ret = bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0)
    {
      if (!first && tid > vcf_tid)
        continue; // Skip ahead if pileup is beyond VCF tid
      if (!first && tid == vcf_tid && pos < vcf_pos)
        continue; // Skip if position is earlier

      // moving forward....
      if (tid == vcf_tid && pos == vcf_pos)
      {
        debugPrint("Debug: Found matching SNP position!", arguments.debug_mode);

        //  {
        //         std::ostringstream oss;
        //         oss << "Debug: Match found - Chrom=" << chrom
        //             << ", Pos=" << vcf_pos
        //             << ", Ref=" << ref
        //             << ", Alt=" << alt;
        //         debugPrint(oss.str(), arguments.debug_mode);
        //     }

        // Process the pileup data for each file
        // output data to output file
        bool is_not_zero = false;
        bool fails_min = false;
        for (i = 0; i < n; ++i)
        {

          {
            std::ostringstream oss;
            oss << "Debug: Processing file" << (i + 1) << " of " << n;
            debugPrint(oss.str(), arguments.debug_mode);
          }
          file_info this_file = {0, 0, 0, 0}; // Initialize counts
          if (n_plp[i] >= arguments.min_read_counts[i])
          {
            //         {
            //   std::ostringstream oss;
            //   oss << "Debug: Sufficient reads for file " << (i + 1);
            //   debugPrint(oss.str(), arguments.debug_mode);
            // }

            for (int j = 0; j < n_plp[i]; ++j)
            {
              const bam_pileup1_t *p = plp[i] + j;
              char base = bam_seqi(bam_get_seq(p->b), p->qpos);
              if (p->is_del)
              {
                this_file.deletions++;
              }
              else if (seq_nt16_str[(unsigned int)base] == ref[0])
              {
                this_file.refs++;
              }
              else if (seq_nt16_str[(unsigned int)base] == alt[0])
              {
                this_file.alts++;
              }
              else
              {
                this_file.errors++;
              }
            }
          }
          else
          {
            {
              std::ostringstream oss;
              oss << "Debug: Insufficient reads for file " << i + 1;
              debugPrint(oss.str(), arguments.debug_mode);
            }
            fails_min = true;
          }
          if (this_file.refs != 0 || this_file.alts != 0 || this_file.errors != 0)
          {
            is_not_zero = true;
          }
          f_info.push_back(this_file);
        }

        debugPrint("Debug: end for loop, start print output", arguments.debug_mode);

        // Write the SNP information to the output
        if (is_not_zero && !fails_min)
        {
          output.str("");
          output.clear();
          output << chrom << "," << vcf_pos << "," << ref << "," << alt;
          for (i = 0; i < n; ++i)
          {
            output << "," << f_info[i].refs
                   << "," << f_info[i].alts
                   << "," << f_info[i].errors
                   << "," << f_info[i].deletions;
          }
          output << "\n";
          if (arguments.gzipped)
          {
            if (bgzf_write(arguments.gzippedPointer, output.str().c_str(), output.str().size()) < 0)
            {
              std::cerr << "Error writing to BGZF output file" << std::endl;
              return 1;
            }
          }
          else
          {
            fwrite(output.str().c_str(), sizeof(char), output.str().size(), output_file);
          }
        }
        f_info.clear();
        break;
      }
      if (arguments.pseudo_snps && !have_snpped && (((pos + 1) % arguments.pseudo_snps) == 0))
      {
        bool is_not_zero = false;
        bool fails_min = false;
        for (i = 0; i < n; i++)
        {
          file_info this_file = file_info();
          this_file.refs = 0;
          this_file.alts = 0;
          this_file.errors = 0;
          this_file.deletions = 0;
          if (n_plp[i] >= arguments.min_read_counts[i])
          {
            for (int j = 0; j < n_plp[i]; ++j)
            {
              const bam_pileup1_t *p = plp[i] + j;
              int c = p->qpos < p->b->core.l_qseq
                          ? bam_get_qual(p->b)[p->qpos]
                          : 0;
              if (c == 0)
              {
                continue; // no
              }
              if (c < arguments.min_base_quality)
              {
                continue; // skip anything with quality below threshold
              }
              this_file.refs++;
            }
          }
          else
          {
            fails_min = true;
          }
          if (this_file.refs != 0)
          {
            is_not_zero = true;
          }
          f_info.push_back(this_file);
        }
        // add a pseudo snp!
        if (is_not_zero && !fails_min)
        {
          output.str("");
          output.clear();
          output << hdr->target_name[tid];
          output << ",";
          output << (pos + 1);
          output << ",.,.";
          for (i = 0; i < n; i++)
          {
            output << ",";
            output << f_info[i].refs;
            output << ",0,0,0";
          }
          output << "\n";
          (arguments.outFunc)(arguments, output.str(), output_file);
        }
        f_info.clear();
      }
      // Break the loop if the BAM position surpasses the VCF position
      if (tid > vcf_tid || (tid == vcf_tid && pos > vcf_pos))
      {
        break;
      }
    }
  }

  // Cleanup
  free(line.s);
  bgzf_close(bgzf);
  bam_mplp_destroy(iter);

  for (i = 0; i < n; ++i)
  {
    if (data[i]->fp)
      hts_close(data[i]->fp);
    free(data[i]);
  }
  free(data);

  debugPrint("Debug: Completed program_main", arguments.debug_mode);
  printf("Finished in %f seconds.\n", (clock() - start) / (double)CLOCKS_PER_SEC);
  if (arguments.gzipped)
  {
    bgzf_close(arguments.gzippedPointer);
  }
  else
  {
    fclose(output_file);
  }
  return 0;
}

// Expose to Rcpp

//' Run Snp Pileup
//'
//' This function runs snp_pileup
//'
//' @param debug A logical Value. If TRUE will print 'debug' statements.
//' @return None.
//' @export
// [[Rcpp::export]]
void run_snp_pileup_logic(const std::vector<std::string> &input_args)
{
  arguments args;

  // Parse the arguments
  parse_arguments(input_args, args);
  debugPrint("Debug: args parsed", args.debug_mode);

  // Run the main program
  // Rcpp::Rcout << "Debug: Running main program" << std::endl;
  debugPrint("Debug: running main program", args.debug_mode);

  int status = program_main(args);
  // Rcpp::Rcout << "Debug: Ran main program" << std::endl;
  debugPrint("Debug: Ran main program", args.debug_mode);
  if (status != 0)
  {
    Rcpp::stop("Program terminated with errors.");
  }
}

//' Check Htslib Version
//'
//' This function prints the version of HTSLIB being used for snp-pileup
//'
//' @return None.
//' @export
// [[Rcpp::export]]
void htslib_version()
{
  // Rcpp::Rcout << "Htslib linked successfully - snp-pileup!" << std::endl;
  Rcpp::Rcout << "HTSlib version:" << hts_version() << std::endl;
}
