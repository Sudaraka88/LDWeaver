#include <Rcpp.h>
#include "kseq2.h" //with gz reading
#include <iostream>
#include <fcntl.h>
#include <stdio.h>
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

  using namespace Rcpp;
// [[Rcpp::export(name = '.readFasta')]]
List readFasta(std::string file, int pos_len) {
  //Rcout << "Step 1: Checking alignment: ";
  int n = 0;
  int l = 0;

  gzFile fp_, fp;
  kseq_t *seq;

  const char * f = file.c_str();

  fp_ = gzopen(f, "r");
  seq = kseq_init(fp_);
  l = kseq_read(seq);
  int seq_length = strlen(seq->seq.s); // This is the first sequence, while loop reads from the second, re-run below?

  if(seq_length != pos_len) { // pos_len must match the seq_length
    return List::create(Named("seq.length") = -1);
  }

//  Rcout << seq_length << " bp per seq \n";
  kseq_destroy(seq);
  gzclose(fp_);


  std::vector<std::string> char_vec; char_vec.reserve(100*seq_length);

  //Rcout << "Step 2: Extracting seq info: ";
  // int fp = open(f, O_RDONLY);
  // kstream<int, FunctorRead> ks(fp, r);
  fp = gzopen(f, "r");
  seq = kseq_init(fp);



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
    char_vec.emplace_back(seq->seq.s);
    n++;
  }
  //Rcout << n << " seqs found \n";
  // close(fp);
  kseq_destroy(seq);
  gzclose(fp);

  List params = List::create(Named("num.seqs") = n,
                             Named("seq.length") = seq_length,
                             Named("seq.names") = seq_names,
                             Named("seq.s") = char_vec);
  return(params);

}
