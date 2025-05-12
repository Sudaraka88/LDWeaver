#include <Rcpp.h>
//#include "kseq.h"
#include "kseq2.h" //with gz reading
#include <iostream>
#include <fcntl.h>
#include <stdio.h>
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

  using namespace Rcpp;

// [[Rcpp::export(name = '.extractAlnParam')]]
List extractAlnParam(std::string file, int filter, double gap_thresh, double maf_thresh) {
  // filter = 0, default in spydrpick
  // filter = 1, relaxed
  // default gap_thresh = 0.15, maf_thresh = 0.01
  Rcout << "Step 1: Checking alignment: ";
  int n = 0;
  int l = 0;
  // kseq seq;

  gzFile fp_, fp;
  kseq_t *seq;

  const char * f = file.c_str();
  // int fp_ = open(f, O_RDONLY);
  // FunctorRead r;
  // kstream<int, FunctorRead> ks_(fp_, r);
  // l = ks_.read(seq);
  // int seq_length = seq.seq.length(); // This is the first sequence, while loop reads from the second, re-run below?
  // close(fp_);

  fp_ = gzopen(f, "r");
  seq = kseq_init(fp_);
  l = kseq_read(seq);
  int seq_length = strlen(seq->seq.s); // This is the first sequence, while loop reads from the second, re-run below?
  Rcout << seq_length << " bp per seq \n";
  kseq_destroy(seq);
  gzclose(fp_);

  Rcout << "Step 2: Extracting seq info: ";
  // int fp = open(f, O_RDONLY);
  // kstream<int, FunctorRead> ks(fp, r);
  fp = gzopen(f, "r");
  seq = kseq_init(fp);

  NumericMatrix allele_counts(5, seq_length);
  Rcpp::StringVector seq_names;
  // while((l = ks.read(seq)) >= 0) {
  while((l = kseq_read(seq)) >= 0) {
    int s_len_t = strlen(seq->seq.s);
    seq_names.push_back(seq->name.s);
    // if(seq.seq.length()!=seq_length) {
    if(s_len_t!=seq_length) {
      return List::create(Named("seq.length") = -1); // failed due to mismatching seq. lengths
      break;
    }
    for(int j=0; j<seq_length; j++){
      if((seq->seq.s[j]=='a') || (seq->seq.s[j]=='A')){
        allele_counts(0, j) += 1;
      } else if((seq->seq.s[j]=='c') || (seq->seq.s[j]=='C')){
        allele_counts(1, j) += 1;
      } else if((seq->seq.s[j]=='g') || (seq->seq.s[j]=='G')){
        allele_counts(2, j) += 1;
      } else if((seq->seq.s[j]=='t') || (seq->seq.s[j]=='T')){
        allele_counts(3, j) += 1;
      } else {
        allele_counts(4, j) += 1;
      }
    }
    // for(int j=0; j<seq_length; j++){
    //   if((seq.seq[j]=='a') || (seq.seq[j]=='A')){
    //     allele_counts(0, j) += 1;
    //   } else if((seq.seq[j]=='c') || (seq.seq[j]=='C')){
    //     allele_counts(1, j) += 1;
    //   } else if((seq.seq[j]=='g') || (seq.seq[j]=='G')){
    //     allele_counts(2, j) += 1;
    //   } else if((seq.seq[j]=='t') || (seq.seq[j]=='T')){
    //     allele_counts(3, j) += 1;
    //   } else { // 'n', 'N', '-', etc.
    //     allele_counts(4, j) += 1;
    //   }
    // }
    n++;
  }
  Rcout << n << " seqs found \n";
  // close(fp);
  kseq_destroy(seq);
  gzclose(fp);

  Rcout << "Step 3: Filtering SNPs: ";
  // Rcout << n  << '\n';
  // allele_counts contain info required to filter...
  IntegerVector idx_ngss = IntegerVector::create(0, 1, 2, 3);
  // LogicalVector filt(seq_length);
  int chk_flag;
  std::vector<int> POS; // holder for indices
  //std::vector<int> consensus; // here we don't really care about the consensus
  POS.reserve(2*10000);
  //consensus.reserve(2*10000);

  int n_snp = 0;
  // to avoid repeated checking of filter technique, repeat whole chunk
  if(filter==0){ // spydrpick default
    int min_maf = n*maf_thresh;
    // Rcout << min_maf << '\n';
    for(int j = 0; j<seq_length; j++){
      if ((POS.capacity() - 2*n_snp)<100){
        POS.reserve( 2*POS.capacity() );
        //    consensus.reserve(2*consensus.capacity());
      }

      chk_flag = 0;
      for(int k = 0; k<4; k++){
        if(allele_counts(k, j) > 0){ // we need at least one non-gap allele
          chk_flag += 1;
          if(chk_flag > 1){ // seems polymorphic
            if(allele_counts(4, j)/n < gap_thresh){ // now check gap content
              NumericVector snp = allele_counts(_, j);
              snp = snp[idx_ngss];
              std::sort(snp.begin(), snp.end());
              if(snp[2] > min_maf){ // second largest non-gap element
                POS.push_back(j+1);
                //          consensus.push_back(which_max(allele_counts(_,j)));
                // filt[j] = true;
                n_snp += 1; // keep count of snps
              }
            }
            chk_flag = 0;
            break;
          }
        }
      }
    }
  } else{ // relaxed
    int min_maf = n*(1-maf_thresh);
    for(int j = 0; j<seq_length; j++){
      chk_flag = 0;
      for(int k = 0; k<4; k++){
        if ((POS.capacity() - 2*n_snp)<100){
          POS.reserve( 2*POS.capacity() );
          //  consensus.reserve(2*consensus.capacity());
        }

        if(allele_counts(k, j) > 0){ // we need at least one non-gap allele
          chk_flag += 1;
          if(chk_flag > 1){ // seems polymorphic
            if(allele_counts(4, j)/n < gap_thresh){ // now check gap content
              // NumericVector snp = allele_counts(_, j);
              // snp = snp[idx_ngss];
              // std::sort(snp.begin(), snp.end());
              //Rcout << j << ":" <<  (1-(max(allele_counts(_,j))/n)) << ' ';
              if( max(allele_counts(_,j)) <= min_maf){
                POS.push_back(j+1);
                //      consensus.push_back(which_max(allele_counts(_,j)));
                // filt[j] = true;
                n_snp += 1;
              }
            }
            chk_flag = 0;
            break;
          }
        }
      }
    }
  }

  Rcout << n_snp << " sites retained \n";
  List params = List::create(Named("num.seqs") = n,
                             Named("num.snps") = n_snp,
                             Named("seq.length") = seq_length,
                             Named("seq.names") = seq_names,
                             Named("pos") = POS);
  return(params);

}

// [[Rcpp::export(name = '.extractSNPs')]]
List extractSNPs(std::string file, int n_seq, int n_snp, std::vector<int> POS) {
  gzFile fp2;
  kseq_t *seq;
  const char * f = file.c_str();
  // Prepare the data for sparse matrix
  Rcout << "Step 4: Building matrices:";

  int n = 0;
  int k; // indexing the positions
  std::vector<int> m_i_A; m_i_A.reserve(2*n_snp);
  std::vector<int> m_j_A; m_j_A.reserve(2*n_snp);
  std::vector<int> m_x_A; m_x_A.reserve(2*n_snp);

  std::vector<int> m_i_C; m_i_C.reserve(2*n_snp);
  std::vector<int> m_j_C; m_j_C.reserve(2*n_snp);
  std::vector<int> m_x_C; m_x_C.reserve(2*n_snp);

  std::vector<int> m_i_G; m_i_G.reserve(2*n_snp);
  std::vector<int> m_j_G; m_j_G.reserve(2*n_snp);
  std::vector<int> m_x_G; m_x_G.reserve(2*n_snp);

  std::vector<int> m_i_T; m_i_T.reserve(2*n_snp);
  std::vector<int> m_j_T; m_j_T.reserve(2*n_snp);
  std::vector<int> m_x_T; m_x_T.reserve(2*n_snp);

  std::vector<int> m_i_N; m_i_N.reserve(2*n_snp);
  std::vector<int> m_j_N; m_j_N.reserve(2*n_snp);
  std::vector<int> m_x_N; m_x_N.reserve(2*n_snp);

  NumericMatrix ACGTN_table(5,n_snp);

  // int fp2 = open(f, O_RDONLY);
  // kstream<int, FunctorRead> ks2(fp2, r);
  fp2 = gzopen(f, "r");
  seq = kseq_init(fp2);
  char temp_char;

  // std::string nt ("AaCcGgTt");
  // we need to add the functionality to find *R* and *UQE* required for MI computation

  // while((l = ks2.read(seq)) >= 0) {
  Rcpp::StringVector seq_names(n_seq);
  int l = 0;
  while((l = kseq_read(seq)) >= 0) {
    // Record sequence names
    // seq_names[n] = seq.name;
    seq_names[n] = seq->name.s;
    k = 1;
    // Rcout << seq_names[n]   << '\n';

    for(int j : POS){
      // temp_char = seq.seq[j-1];
      temp_char = seq->seq.s[j-1];
      // Rcout << j << ':' << temp_char << ' ';
      if((temp_char=='A') || (temp_char=='a')) { // && (consensus[k-1]!=0))){ // A = 0
        m_i_A.push_back(n+1); // n+1 for base 1 indexing in R
        // m_j.push_back(j); // j's are from POS - based 1 R idx already applied
        m_j_A.push_back(k);
        m_x_A.push_back(1); // 'A' but different from reference (is this what we need?)
        ++ACGTN_table(0, (k-1));
      } else if((temp_char=='C') || (temp_char=='c')) {// && (consensus[k-1]!=1))){ // C = 1
        m_i_C.push_back(n+1);
        // m_j.push_back(j);
        m_j_C.push_back(k);
        m_x_C.push_back(2);
        ++ACGTN_table(1, (k-1));
      } else if((temp_char=='G') || (temp_char=='g')) { //}&& (consensus[k-1]!=2))){ //G = 2
        m_i_G.push_back(n+1);
        // m_j.push_back(j);
        m_j_G.push_back(k);
        m_x_G.push_back(3);
        ++ACGTN_table(2, (k-1));
      } else if((temp_char=='T') || (temp_char=='t')) {//}&& (consensus[k-1]!=3))){ // T = 3
        m_i_T.push_back(n+1);
        // m_j.push_back(j);
        m_j_T.push_back(k);
        m_x_T.push_back(4);
        ++ACGTN_table(3, (k-1));
      } else { //}if((temp_char=='N') || (temp_char=='n')){ // not possible
        m_i_N.push_back(n+1);
        // m_j.push_back(j);
        m_j_N.push_back(k);
        m_x_N.push_back(5);
        ++ACGTN_table(4, (k-1));
      }
      k += 1; // SNP indexes
    }
    n += 1;
  }
  // close(fp2);
  kseq_destroy(seq);
  gzclose(fp2);
  // Rcout << n  << '\n';
  Rcout << " Done! \n";

  return List::create(Named("seq.names") = seq_names,
                      Named("ACGTN_table") = wrap(ACGTN_table),
                      Named("i_A") = wrap(m_i_A),
                      Named("j_A") = wrap(m_j_A),
                      Named("x_A") = wrap(m_x_A),
                      Named("i_C") = wrap(m_i_C),
                      Named("j_C") = wrap(m_j_C),
                      Named("x_C") = wrap(m_x_C),
                      Named("i_G") = wrap(m_i_G),
                      Named("j_G") = wrap(m_j_G),
                      Named("x_G") = wrap(m_x_G),
                      Named("i_T") = wrap(m_i_T),
                      Named("j_T") = wrap(m_j_T),
                      Named("x_T") = wrap(m_x_T),
                      Named("i_N") = wrap(m_i_N),
                      Named("j_N") = wrap(m_j_N),
                      Named("x_N") = wrap(m_x_N));
}


// [[Rcpp::export(name = '.extractRef')]]
List extractRef(std::string file) {
  Rcout << "Checking reference... ";
  gzFile fp_;
  kseq_t *seq;

  int l = 0;

  const char * f = file.c_str();
  Rcpp::StringVector seq_names;

  fp_ = gzopen(f, "r");
  seq = kseq_init(fp_);
  l = kseq_read(seq);
  int seq_length = strlen(seq->seq.s); // This is the first sequence, while loop reads from the second, re-run below?
  Rcout << seq_length << " bp seq found!\n";

  Rcpp::String ref = seq->seq.s;


  seq_names.push_back(seq->name.s);
  // seq->seq.s

  kseq_destroy(seq);
  gzclose(fp_);

  return List::create(Named("seq.name") = seq_names,
                      Named("seq.length") = seq_length,
                      Named("seq") = ref);
}
