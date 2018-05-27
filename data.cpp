/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 */


#include "data.h"

Data::Data()
{
   N = 0;
   p = 0;
   K = 0;
   nsnps = 0;
   visited = NULL;
   genovec = NULL;
   // dosage_present = NULL;
   // dosage_main = NULL;
   avg = NULL;
   verbose = false;
   use_preloaded_maf = false;
   pgfi_alloc = nullptr;
   pgr_alloc = nullptr;
   plink2::PreinitPgfi(&pgfi);
   plink2::PreinitPgr(&pgr);
}

Data::~Data()
{
   if(visited)
      delete[] visited;
   plink2::aligned_free_cond(genovec);
   if(avg)
      delete[] avg;
   plink2::CleanupPgr(&pgr);
   plink2::CleanupPgfi(&pgfi);
   plink2::aligned_free_cond(pgfi_alloc);
   plink2::aligned_free_cond(pgr_alloc);
   in.close();
}

/*
 *                   plink BED           sparsnp
 * minor homozyous:  00 => numeric 0     10 => numeric 2
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * major homozygous: 11 => numeric 3     00 => numeric 0
 * missing:          01 => numeric 1     11 => numeric 3
 *
 *
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says,
 * The bytes in plink are read backwards HGFEDCBA, not GHEFCDAB, but we read
 * them forwards as a character (a proper byte)
 *
 * By default, plink usage dosage of the *major* allele, since allele A1 is
 * usually the minor allele and the code "1" refers to the second allele A2,
 * so that "11" is A2/A2 or major/major.
 *
 * We always use minor allele dosage, to be consistent with the output from
 * plink --recodeA which used minor allele dosage by default.
 *
 * out: array of genotypes
 * in: array of packed genotypes (bytes)
 * n: number of bytes in input
 *
 */

#define PLINK2_NA 65535
#define PLINK2_SCALE 16384
#define PLINK2_INVSCALE (1.0 / PLINK2_SCALE)

void decode_plink2_hc(const uintptr_t* genovec, const uint32_t sample_ct, const uint32_t bidx, double mean, double sd, MatrixXd* outp) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / plink2::kBitsPerWordD2;
  double lookup[4];
  if (sd > VAR_TOL) {
    const double inv_sd = 1.0 / sd;
    lookup[0] = (0 - mean) * inv_sd;
    lookup[1] = (1 - mean) * inv_sd;
    lookup[2] = (2 - mean) * inv_sd;
    lookup[3] = 0;  // impute to average
    // std::cout << lookup[0] << "," << lookup[1] << "," << lookup[2] << std::endl;
  } else {
    lookup[0] = 0;
    lookup[1] = 0;
    lookup[2] = 0;
    lookup[3] = 0;
  }
  uint32_t loop_len = plink2::kBitsPerWordD2;
  uint32_t widx = 0;
  while (1) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        return;
      }
      loop_len = 1 + (sample_ct - 1) % plink2::kBitsPerWordD2;
    }
    uintptr_t cur_word = genovec[widx];
    const uint32_t offset = widx * plink2::kBitsPerWordD2;
    for (uint32_t uii = 0; uii < loop_len; ++uii) {
      const uintptr_t cur_plink1_geno = cur_word & 3;
      (*outp)(uii + offset, bidx) = lookup[cur_plink1_geno];
      cur_word >>= 2;
    }
    ++widx;
  }
}

void Data::get_size()
{
   verbose && STDOUT << timestamp() << "Analyzing BED/PGEN file '"
      << geno_filename << "'";
   plink2::PgenHeaderCtrl header_ctrl;
   char errstr_buf[plink2::kPglErrstrBufBlen];
   uintptr_t cur_alloc_cacheline_ct;
   plink2::PglErr reterr = PgfiInitPhase1(geno_filename, 0xffffffffU, N, 0, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, errstr_buf);
   if (reterr) {
      throw std::runtime_error(errstr_buf);
   }
   if (plink2::cachealigned_malloc(cur_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
      throw std::runtime_error("Out of memory.");
   }
   nsnps = pgfi.raw_variant_ct;
   uint32_t max_vrec_width;
   reterr = PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, nsnps, &max_vrec_width, &pgfi, pgfi_alloc, &cur_alloc_cacheline_ct, errstr_buf);
   if (reterr) {
      throw std::runtime_error("Out of memory.");
   }
   if (plink2::cachealigned_malloc(cur_alloc_cacheline_ct * plink2::kCacheline, &pgr_alloc)) {
      throw std::runtime_error("Out of memory.");
   }
   reterr = PgrInit(geno_filename, max_vrec_width, &pgfi, &pgr, pgr_alloc);
   if (reterr) {
      throw std::runtime_error("Out of memory.");
   }
}

// Prepare input stream etc before reading in SNP blocks
void Data::prepare()
{
  if (plink2::cachealigned_malloc(plink2::RoundUpPow2((N + 3) / 4, plink2::kCacheline), &genovec)) {
      throw std::runtime_error("Out of memory.");
   }
   // todo: dosage_present, dosage_main

   avg = new double[nsnps]();
   visited = new bool[nsnps]();
   X_meansd = MatrixXd::Zero(nsnps, 2); // TODO: duplication here with avg

   // scaled_geno_lookup = ArrayXXd::Zero(4, nsnps);

   verbose && STDOUT << timestamp() << "Detected BED file: "
      << geno_filename << " with " << (len + 3)
      << " bytes, " << N << " samples, " << nsnps
      << " SNPs." << std::endl;
}

// Reads a _contiguous_ block of SNPs [start, stop] at a time.
// The block will contain standardised genotypes already, no need to
// standardise them again.
//
// If resize is false, then the calling code is responsible for ensuring that
// X is handled accordingly with the block size (X may be bigger than the
// block).
void Data::read_snp_block(unsigned int start_idx, unsigned int stop_idx,
                          bool transpose, bool resize)
{
   in.seekg(3 + np * start_idx);

   unsigned int actual_block_size = stop_idx - start_idx + 1;

   // Resize matrix, e.g., with final block that may be smaller than
   // $blocksize$
   if(transpose)
   {
      if(X.rows() == 0 || (resize && X.rows() != actual_block_size))
      {
         verbose && STDOUT << timestamp()
	    << "Reallocating memory: " << X.rows() << " -> " <<
	    actual_block_size << std::endl;
         if(X.rows() > actual_block_size)
	 {
            X = MatrixXd(actual_block_size, N);
	 }
      }
   }
   else if(X.cols() == 0 || (resize && X.cols() != actual_block_size))
   {
      verbose && STDOUT << timestamp()
	 << "Reallocating memory: " << X.cols() << " -> " <<
	 actual_block_size << std::endl;
      X = MatrixXd(N, actual_block_size);
   }

   for(unsigned int j = 0; j < actual_block_size; j++)
   {
      unsigned int k = start_idx + j;

      // read raw genotypes
      // replace this with PgrGetD(nullptr, nullptr, N, k, &pgr, genovec, dosage_present, dosage_main, &dosage_ct)
      plink2::PglErr reterr = PgrGet(nullptr, nullptr, N, k, &pgr, genovec);
      if (reterr) {
         throw std::runtime_error("File read/decode failure.");
      }

      // Compute average per SNP, excluding missing values
      double snp_avg = 0;
      unsigned int ngood = 0;

      // We've seen this SNP, don't need to compute its average again
      if(!visited[k])
      {
	 double P, sd;

	 if(!use_preloaded_maf)
	 {
            std::array<uint32_t, 4> genocounts;
            plink2::ZeroTrailingQuaters(N, genovec);
            plink2::GenovecCountFreqsUnsafe(genovec, N, genocounts);
            ngood = genocounts[0] + genocounts[1] + genocounts[2];
            snp_avg = genocounts[1] + 2 * genocounts[2];
      	    snp_avg /= ngood;

	    // Store the 4 possible standardised genotypes for each SNP
	    P = snp_avg / 2.0;
	    if(stand_method_x == STANDARDISE_BINOM)
	       sd = sqrt(P * (1 - P));
	    else if(stand_method_x == STANDARDISE_BINOM2)
	       sd = sqrt(2.0 * P * (1 - P));
	    else
	    {
	       std::string err = std::string("unknown standardisation method: ")
		  + std::to_string(stand_method_x);
	       throw std::runtime_error(err);
	    }

	    X_meansd(k, 0) = snp_avg;
	    X_meansd(k, 1) = sd;
	 }
	 else
	 {
	    snp_avg = X_meansd(k, 0);
	    sd = X_meansd(k, 1);
	 }

	 // scaled genotyped initialised to zero
	 visited[k] = true;
      }

      // Unpack the genotypes, but don't convert to 0/1/2/NA, keep in
      // original format (see comments for decode_plink).
      // There is a bit of waste here in the first time the SNP is visited, as
      // we unpack the data twice, once with decoding and once without.
      decode_plink2_hc(genovec, N, j, X_meansd(k, 0), X_meansd(k, 1), &X);
   }
}

// Reads an entire bed file into memory
// Expects PLINK bed in SNP-major format
void Data::read_bed() {
   X = MatrixXd(N, nsnps);

   unsigned int md = nsnps / 50;

   // iterate over all SNPs
   for(unsigned int j = 0 ; j < nsnps; j++)
   {
      // read raw genotypes
      // replace this with PgrGetD(nullptr, nullptr, N, k, pgrp, genovec, dosage_present, dosage_main, &dosage_ct)
      plink2::PglErr reterr = PgrGet(nullptr, nullptr, N, j, &pgr, genovec);
      if (reterr) {
         throw std::runtime_error("File read/decode failure.");
      }

      // Compute average per SNP, excluding missing values
      std::array<uint32_t, 4> genocounts;
      plink2::ZeroTrailingQuaters(N, genovec);
      plink2::GenovecCountFreqsUnsafe(genovec, N, genocounts);
      uint32_t ngood = genocounts[0] + genocounts[1] + genocounts[2];
      double snp_avg = genocounts[1] + 2 * genocounts[2];
      snp_avg /= ngood;
      avg[j] = snp_avg;

      // Impute using average per SNP
      decode_plink2_hc(genovec, N, j, snp_avg, 1.0, &X);

      if(verbose && j % md == md - 1)
	 STDOUT << timestamp() << "Reading genotypes, "
	    << roundl(((double)j / nsnps) * 100) << "% done"
	    << std::endl;
   }

   p = X.cols();

   verbose && STDOUT << timestamp() << "Loaded genotypes: "
      << N << " samples, " << p << " SNPs" << std::endl;
}

void Data::read_pheno(const char *filename, unsigned int firstcol)
{
    NamedMatrixWrapper M = read_text(filename, firstcol);
    Y = M.X;
    N = M.X.rows();
}

// Reads PLINK phenotype files:
// FID IID pheno1 pheno2 ...
// Need to be able to read continuous phenotypes
//
// firstcol is _one-based_, 3 for pheno file, 6 for FAM file (ignoring sex),
// 5 for FAM file (with gender)
NamedMatrixWrapper read_text(const char *filename,
   unsigned int firstcol, unsigned int nrows, unsigned int skip,
   bool verbose)
{
   NamedMatrixWrapper M;

   unsigned int line_num = 0;

   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::string err = std::string("Error reading file '")
	 + filename + "': " + strerror(errno);
      throw std::runtime_error(err);
   }
   std::vector<std::string> lines;

   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof() && (nrows == -1 || line_num < nrows))
      {
	 if(line_num >= skip)
	    lines.push_back(line);
	 line_num++;
      }
   }

   verbose && STDOUT << timestamp() << "Detected text file " <<
      filename << ", " << lines.size() << " rows" << std::endl;

   in.close();

   unsigned int numtok = 0, numfields, numfields_1st = 0;

   M.X = MatrixXd(0, 0);

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);

      numtok = tokens.size();
      numfields = numtok - firstcol + 1;
      if(i == 0)
      {
	 M.X.resize(lines.size(), numfields);
	 numfields_1st = numfields;
      }
      else if(numfields_1st != numfields)
      {
	 std::string err = std::string("Error reading file '")
	    + filename + "': inconsistent number of columns";
	 throw std::runtime_error(err);
      }

      VectorXd y(numfields);
      char* err;
      errno = 0;
      for(unsigned int j = 0 ; j < numfields ; j++)
      {
	 //y(j) = std::atof(tokens[j + firstcol - 1].c_str());
	 double m = std::strtod(tokens[j + firstcol - 1].c_str(), &err);
	 if(*err != '\0' || errno != 0)
	 {
	    std::string err = std::string("Error reading file '")
	       + filename + "', line " + std::to_string(i + 1)
	       + ": '" + tokens[j + firstcol - 1] + "'"
	       + " cannot be parsed as a number";
	    throw std::runtime_error(err);
	 }
	 y(j) = m;
      }
      M.X.row(i) = y;
   }

   return M;
}

void Data::read_plink_bim(const char *filename)
{
   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::string err = std::string("Error reading file ")
	 + filename;
      throw std::runtime_error(err);
   }
   std::vector<std::string> lines;

   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof())
	 lines.push_back(line);
   }

   verbose && STDOUT << timestamp() << "Detected bim file " <<
      filename << ", " << lines.size() << " SNPs" << std::endl;
   in.close();

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);
      snp_ids.push_back(tokens[1]);
      ref_alleles.push_back(tokens[4]);
      alt_alleles.push_back(tokens[5]);

      char* err;
      errno = 0;
      unsigned long long m = std::strtol(tokens[3].c_str(), &err, 10);
      if(*err != '\0' || errno != 0)
      {
	 std::string err = std::string("Error reading file '")
	    + filename + "', line " + std::to_string(i + 1)
	    + ": '" + tokens[3] + "' cannot be parsed as a number";
	 throw std::runtime_error(err);
      }
      bp.push_back(m);
   }
}

void Data::read_plink_fam(const char *filename)
{
   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::string err = std::string(
	 "[Data::read_plink_fam] Error reading file ") + filename;
      throw std::runtime_error(err);
   }
   std::vector<std::string> lines;

   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof())
	 lines.push_back(line);
   }

   in.close();

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);
      fam_ids.push_back(tokens[0]);
      indiv_ids.push_back(tokens[1]);
   }
}

std::string Data::tolower(const std::string& v)
{
   std::string r = v;
   std::transform(r.begin(), r.end(), r.begin(), ::tolower);
   return r;
}
