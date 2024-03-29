/*
    PredictHaplo-Paired: a program for estimating haplotypes from "next
   generation sequencing" reads Copyright (C) 2014 Volker Roth

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// TODO: where does one actually get this definition from?
#define R_finite(x) std::isfinite(x)

#include "phaplo/Exception.hpp"
#include "phaplo/ParseStatistics.hpp"
#include "phaplo/binary.hpp"
#include "phaplo/is_alignment_match.hpp"
#include "phaplo/is_cigar_code.hpp"
#include "phaplo/median.hpp"
#include "phaplo/statistics_to_messages.hpp"
#include "phaplo/to_vector_of_ints.hpp"
#include "scythestat/distributions.h"
#include "scythestat/ide.h"
#include "scythestat/la.h"
#include "scythestat/lapack.h"
#include "scythestat/matrix.h"
#include "scythestat/optimize.h"
#include "scythestat/rng.h"
#include "scythestat/rng/mersenne.h"
#include "scythestat/smath.h"
#include "scythestat/stat.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <unordered_map>

using namespace scythe;
using namespace std;

template <class T> struct index_cmp {
  index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const {
    return arr[a] > arr[b];
  }
  const T arr;
};

bool myfunction(int i, int j) { return (i > j); }

int visualization_level = 1;
int local_window_size;

vector<string> tokenize(const string &str, const string &delimiters) {
  vector<string> tokens;

  string::size_type lastPos = 0, pos = 0;
  int count = 0;

  if (str.length() < 1)
    return tokens;

  // skip delimiters at beginning.
  lastPos = str.find_first_not_of(delimiters, 0);

  if ((str.substr(0, lastPos - pos).length()) > 0) {
    count = str.substr(0, lastPos - pos).length();

    for (int i = 0; i < count; i++)
      tokens.push_back("");

    if (string::npos == lastPos)
      tokens.push_back("");
  }

  // find first \"non-delimiter\".
  pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // skip delimiters.  Note the \"not_of\"
    lastPos = str.find_first_not_of(delimiters, pos);

    if ((string::npos != pos) &&
        (str.substr(pos, lastPos - pos).length() > 1)) {
      count = str.substr(pos, lastPos - pos).length();

      for (int i = 0; i < count - 1; i++)
        tokens.push_back("");
    }

    pos = str.find_first_of(delimiters, lastPos);
  }

  return tokens;
}

int qual_subtract;
bool include_deletions = true;
int levels;

int parse_sam_line(const vector<string> &tokens, string &used_qual,
                   vector<int> &seq_b, double &a_score, const char gap_quality,
                   int &indels, bool &have_quality_scores,
                   int &al_start) noexcept {

  used_qual.clear();
  seq_b.clear();
  a_score = std::numeric_limits<double>::quiet_NaN();
  indels = 0;
  have_quality_scores = true;
  al_start = atoi(tokens[3].c_str());

  // cout <<"tokens.size(): " <<  tokens.size() << endl;
  for (int j = 11; j < tokens.size(); j++) {
    vector<string> tokens_as = tokenize(tokens[j].c_str(), ":");
    if (tokens_as[0] == "AS") {
      a_score = atof(tokens_as[2].c_str());
      break;
    }
  }

  vector<int> isCode(tokens[5].size(), 1);
  for (int i = 0; i < tokens[5].size(); i++) {
    if (!phaplo::is_cigar_code(tokens[5][i]))
      isCode[i] = 0;
  }

  vector<int> sub_length_vec;
  vector<char> symbols;
  int sub_length = 0;
  for (int i = 0; i < tokens[5].size(); i++) {
    if (isCode[i] == 1) {
      sub_length_vec.push_back(
          atoi(tokens[5].substr(i - sub_length, sub_length).c_str()));
      symbols.push_back(tokens[5][i]);
      sub_length = 0;
    }
    if (isCode[i] == 0)
      sub_length++;
  }

  string qual_s = tokens[10];
  if (qual_s.size() == 1 && qual_s[0] == '*')
    have_quality_scores = false;
  int c = 0;

  for (int i = 0; i < sub_length_vec.size(); i++) {

    if (symbols[i] == 'S')
      for (int j = 0; j < sub_length_vec[i]; j++) {
        c++;
      }

    if (phaplo::is_alignment_match(symbols[i]))
      for (int j = 0; j < sub_length_vec[i]; j++) {
        int k = 0;
        if (tokens[9][c] == 'A' || tokens[9][c] == 'a')
          k = 1;
        else if (tokens[9][c] == 'C' || tokens[9][c] == 'c')
          k = 2;
        else if (tokens[9][c] == 'G' || tokens[9][c] == 'g')
          k = 3;
        else if (tokens[9][c] == 'T' || tokens[9][c] == 't')
          k = 4;
        seq_b.push_back(k);
        if (have_quality_scores) {
          // cout << int(qual_s[c] - qual_subtract) <<':'<<pow(10,-
          // int(qual_s[c] - qual_subtract)/10.0)<<' '<<flush; //newCout

          int q =
              int(41 * (1 - pow(10, -int(qual_s[c] - qual_subtract) / 10.0)) +
                  qual_subtract);
          qual_s[c] = char(q);
          used_qual.append(1, qual_s[c]);
        } else
          used_qual.append(1, 'I');
        // cout<< tokens[9][c];
        c++;
      }
    else if (symbols[i] == 'I')
      for (int j = 0; j < sub_length_vec[i]; j++) {

        // cout<< tokens[9][c];
        c++;
        // indels++;
      }
    else if (symbols[i] == 'D')
      for (int j = 0; j < sub_length_vec[i]; j++) {

        used_qual.append(1, gap_quality);

        seq_b.push_back(0);
        indels++;
        // cout<< '*';
      }
    else if (symbols[i] == 'P') {
      // for(int j =0; j< sub_length_vec[i];j++){
      // cout<< '*';
      // }
    }
  }

  return 0;
}

phaplo::ParseStatistics
parseSAMpaired(string al, double max_gap_fraction,
               double min_align_score_fraction, double min_qual, int min_length,
               char gap_quality, double &mean_length, int &min_seq_start,
               int &max_seq_start, int &max_sequence_stop,
               int &min_sequence_stop, bool &have_quality_scores,
               vector<string> &quality_scores, vector<int> &strand,
               vector<vector<int>> &Reads, vector<int> &Positions_Start,
               vector<string> &IDs) {

  string line;

  have_quality_scores = true;

  ifstream inf6(al.c_str(), ios::in);

  vector<int> seq_b, seq_b_pairs;
  string used_qual, used_qual_pairs;
  double a_score, a_score_pairs;

  int al_start, al_start_pairs;

  int indels, indels_pairs;

  auto statistics = phaplo::ParseStatistics();

  std::unordered_map<std::string, std::vector<std::string>> candidates;
  while (getline(inf6, line, '\n')) {
    if (line[0] == '@') {
      continue;
    }

    statistics.reads.total_count++;

    const auto tokens = tokenize(line, "\t");
    const auto RC = static_cast<unsigned int>(atoi(tokens[1].c_str()));
    const auto id = tokens[0];

    const auto sRC = phaplo::binary(RC);

    if (phaplo::used_bits(sRC) > 8) {
      statistics.reads.unsupported_attribute++;
      continue;
    }

    const auto is_unmapped = bool{sRC[2] || sRC[3]};
    const auto is_paired_read = sRC[0];

    if (is_unmapped) {
      statistics.reads.unmapped++;
      continue;
    }

    if (!is_paired_read) {
      statistics.reads.unpaired++;
      continue;
    }

    const auto found = candidates.find(id);

    if (found != candidates.end()) {
      const auto is_read1 = bool{sRC[6]};
      const auto tokens_1 =
          is_read1 ? std::move(tokens) : std::move(found->second);
      const auto tokens_2 =
          !is_read1 ? std::move(tokens) : std::move(found->second);
      candidates.erase(found);

      statistics.pairs.total_count++;

      parse_sam_line(tokens_1, used_qual, seq_b, a_score, gap_quality, indels,
                     have_quality_scores, al_start);
      parse_sam_line(tokens_2, used_qual_pairs, seq_b_pairs, a_score_pairs,
                     gap_quality, indels_pairs, have_quality_scores,
                     al_start_pairs);

      if (!(seq_b.size() > min_length && seq_b_pairs.size() > min_length)) {
        statistics.pairs.sequence_too_short++;
        continue;
      }

      if (!(double(indels) / seq_b.size() < max_gap_fraction &&
            double(indels_pairs) / seq_b_pairs.size() < max_gap_fraction)) {
        statistics.pairs.gap_fraction_too_high++;
        continue;
      }

      if (!((std::isnan(a_score) ||
             a_score / seq_b.size() > min_align_score_fraction) &&
            (std::isnan(a_score_pairs) ||
             a_score_pairs / seq_b_pairs.size() > min_align_score_fraction))) {
        statistics.pairs.align_score_fraction_too_low++;
        continue;
      }

      int StartPos = al_start;
      vector<int> SEQ_combined = seq_b;
      string Qscores = used_qual;
      bool is_gap = false;
      int Nlength = 0;

      if (al_start_pairs == al_start) {
        // of_gaps << 0<<',';
      }

      if (al_start_pairs < al_start) {
        StartPos = al_start_pairs;

        is_gap = false;
        if (al_start - (al_start_pairs + (int)seq_b_pairs.size()) > 0) {
          is_gap = true;
        }

        if (is_gap) {

          Nlength = al_start - (al_start_pairs + (int)seq_b_pairs.size());

          vector<int> Ns(Nlength, 7);

          seq_b_pairs.insert(seq_b_pairs.end(), Ns.begin(), Ns.end());
          seq_b_pairs.insert(seq_b_pairs.end(), seq_b.begin(), seq_b.end());
          SEQ_combined = seq_b_pairs;

          string Nstr(Nlength, '_');

          Qscores = used_qual_pairs + Nstr + used_qual;

        } else {

          // of_gaps << 0<<',';
          vector<int> first_part(seq_b_pairs.begin(),
                                 seq_b_pairs.begin() +
                                     (al_start - al_start_pairs));
          first_part.insert(first_part.end(), seq_b.begin(), seq_b.end());
          SEQ_combined = first_part;

          Qscores =
              used_qual_pairs.substr(0, al_start - al_start_pairs) + used_qual;
        }
      } else {
        if (al_start_pairs > al_start) {

          is_gap = false;
          if (al_start_pairs - (al_start + (int)seq_b.size()) > 0) {
            is_gap = true;
          }

          if (is_gap) {
            Nlength = al_start_pairs - (al_start + (int)seq_b.size());
            // of_gaps << Nlength<<',';

            vector<int> Ns(Nlength, 7);

            seq_b.insert(seq_b.end(), Ns.begin(), Ns.end());
            seq_b.insert(seq_b.end(), seq_b_pairs.begin(), seq_b_pairs.end());
            SEQ_combined = seq_b;

            string Nstr(Nlength, 'Z');
            Qscores = used_qual + Nstr + used_qual_pairs;
          } else {
            // of_gaps << 0<<',';
            vector<int> first_part(seq_b.begin(),
                                   seq_b.begin() + (al_start_pairs - al_start));

            first_part.insert(first_part.end(), seq_b_pairs.begin(),
                              seq_b_pairs.end());
            SEQ_combined = first_part;

            Qscores = used_qual.substr(0, al_start_pairs - al_start) +
                      used_qual_pairs;
          }
        }
      }

      if (!is_gap || Nlength < 200) {
        Positions_Start.push_back(StartPos);
        Reads.push_back(SEQ_combined);
        IDs.push_back(id);
        mean_length += SEQ_combined.size();

        string qNeqStr(SEQ_combined.size(), 'Z');

        quality_scores.push_back(Qscores);
        strand.push_back(RC);

        if (al_start + seq_b.size() > max_sequence_stop)
          max_sequence_stop = al_start + seq_b.size();

        if (al_start + seq_b.size() < min_sequence_stop)
          min_sequence_stop = al_start + seq_b.size();

        if (al_start < min_seq_start)
          min_seq_start = al_start;

        if (al_start > max_seq_start)
          max_seq_start = al_start;
      }
    } else {
      candidates[id] = tokens;
    }
  }
  return statistics;
}

int MultiNomialDPMReadsSemiEntropy(
    const int &nSample,

    const vector<vector<int>> &Reads,

    const vector<string> &quality_scores, const vector<int> &Positions_Start,

    const vector<int> &reads_in_window,

    vector<vector<double>> Assignments, // supervision information
    const int &WindowStart, const int &WindowStop,

    const int &K,

    const double &alphaMN, const double &alpha,

    Matrix<double> &MNprob, //(GD,K*levels) log-scale (+1e-6)
    Matrix<double> &MNprobMEAN,
    Matrix<int> &maxTab, // (GD,K);
    Matrix<double> &pi,  //(K,1)
    Matrix<double> &piMean, vector<int> &C,
    Matrix<int> &CMAX, // (n,K)

    const int &lseed, const double &burnin_fraction,

    const vector<int> &refSeq,

    const vector<int> &entropy_select,

    bool local) {

  int n = reads_in_window.size();
  int GD = WindowStop - WindowStart + 1;

  int burnin = int(nSample * burnin_fraction);

  int idx, idx_global;

  vector<int> C_old(n, 0);

  mersenne myrng;
  myrng.initialize(lseed);

  piMean = 0;
  MNprobMEAN = 0;
  CMAX = 0;

  Matrix<double> alphaProb(K, 1);
  vector<double> bins(K);
  Matrix<double> cprobL(K, 1);
  Matrix<double> cprob(K, n);

  // vector<double> logPi(K,log(1.0/K));

  Matrix<double> logPi(n, K);

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < K; k++) {
      logPi(j, k) = log(1.0 / K);
    }
  }

  vector<int> nK_vec(K, 0);

  vector<long> lev(levels, 0);
  vector<vector<long>> GD_lev(GD, lev);
  vector<vector<vector<long>>> counts(K, GD_lev); // dims = (K,GD,levels);

  Matrix<double> probs(GD, levels);
  Matrix<double> tot_counts(GD, levels);
  Matrix<double> diriMN(levels, 1);

  for (int j = 0; j < n; j++) {
    nK_vec[C[j]] += 1;
  }

  int nuc;
  // add one "reference" read...

  for (int k = 0; k < K; k++) {
    for (int l = 0; l < GD; l++) {
      idx_global = l + WindowStart - 1;
      if (entropy_select[idx_global] == 1) {

        nuc = refSeq[l];
        if (nuc < 5)
          counts[k][l][nuc] += 1;
      }
    }
    nK_vec[k] += 1;
  }

  vector<int> is_clamped(n, 0);

  int qual;
  for (int o = 0; o < n; o++) {
    for (int r = 0; r < Reads[reads_in_window[o]].size(); r++) {

      idx_global = r + Positions_Start[reads_in_window[o]] - 1;
      idx = r + Positions_Start[reads_in_window[o]] - WindowStart;

      if (Positions_Start[reads_in_window[o]] - WindowStart >
              local_window_size &&
          Positions_Start[reads_in_window[o]] +
                  Reads[reads_in_window[o]].size() - WindowStart <
              GD - local_window_size)
        is_clamped[o] = 1;

      if (idx >= 0 && idx < GD && entropy_select[idx_global] == 1) {

        nuc = Reads[reads_in_window[o]][r];
        qual = int(quality_scores[reads_in_window[o]][r]) - qual_subtract;
        // cout << qual <<' '<< flush;
        if (!include_deletions && nuc == 0)
          continue;

        if (nuc < 5) {

          counts[C[o]][idx][nuc] += qual;
        }
      }
    }
  }

  cout << "Now running MCMC iterations... " << endl;

  for (int i = 0; i < nSample; i++) {

    C_old = C;

    for (int k = 0; k < K; k++) {

      for (int l = 0; l < GD; l++) {

        idx_global = l + WindowStart - 1;
        if (entropy_select[idx_global] == 1) {
          for (int m = 0; m < levels; m++)
            probs(l, m) = counts[k][l][m] + 1e-3;
          probs(l, _) = probs(l, _) * nK_vec[k] / sum(probs(l, _));
          // if( nK_vec[k] == 1)
          // cout << k <<' '<< probs(l,_)<<' '<< nK_vec[k] << endl;
        }
      }

      for (int l = 0; l < GD; l++) {
        idx_global = l + WindowStart - 1;
        if (entropy_select[idx_global] == 1) {
          diriMN = myrng.rdirich(t(probs(l, _)) + alphaMN);

          for (int p = 0; p < levels; p++) {
            if (i > burnin)
              MNprobMEAN(l, k * levels + p) =
                  MNprobMEAN(l, k * levels + p) + diriMN[p];
            MNprob(l, k * levels + p) = log(diriMN[p]);
          }
        }
      }
    }

    for (int o = 0; o < n; o++) {

      if (is_clamped[o])
        continue;

      vector<double> lMN(K, 0);

      for (int r = 0; r < Reads[reads_in_window[o]].size(); r++) {

        idx_global = r + Positions_Start[reads_in_window[o]] - 1;
        idx = r + Positions_Start[reads_in_window[o]] - WindowStart;

        if (idx >= 0 && idx < GD && entropy_select[idx_global] == 1) {

          nuc = Reads[reads_in_window[o]][r];

          if (!include_deletions && nuc == 0)
            continue;

          if (nuc < 5)
            for (int k = 0; k < K; k++)
              lMN[k] += MNprob(idx, k * levels + nuc);
        }
      }
      for (int k = 0; k < K; k++)
        cprob(k, o) = lMN[k] + logPi(o, k);
    }

    for (int l = 0; l < n; l++) {

      if (is_clamped[l]) {
        if (i > burnin) {
          CMAX(l, C[l]) += 1;
        }
        continue;
      }
      double lMax = max(cprob(_, l));
      cprobL = exp(cprob(_, l) - lMax);

      cprobL = cprobL / sum(cprobL);

      bins[0] = cprobL[0];

      for (int k = 1; k < K; k++) {
        bins[k] = bins[k - 1] + cprobL[k];
      }

      double unif = myrng.runif() * bins[K - 1];
      bins[K - 1] += 1e-6;
      for (int k = 0; k < K; k++) {
        if (unif < bins[k]) {
          C[l] = k;
          if (i > burnin) {
            CMAX(l, k) += 1;
          }

          if (C[l] != C_old[l]) {

            for (int r = 0; r < Reads[reads_in_window[l]].size(); r++) {

              idx_global = r + Positions_Start[reads_in_window[l]] - 1;
              idx = r + Positions_Start[reads_in_window[l]] - WindowStart;

              if (idx >= 0 && idx < GD && entropy_select[idx_global] == 1) {

                nuc = Reads[reads_in_window[l]][r];
                qual =
                    int(quality_scores[reads_in_window[l]][r]) - qual_subtract;
                if (!include_deletions && nuc == 0)
                  continue;

                if (nuc < 5) {
                  // downdate couts in  C_old[l]
                  counts[C_old[l]][idx][nuc] -= qual;
                  // update couts in  C[l]
                  counts[C[l]][idx][nuc] += qual;
                }
              }
            }

            nK_vec[C_old[l]] -= 1;

            nK_vec[C[l]] += 1;
          }
          break;
        }
      }
    }

    if (local) {

      for (int k = 0; k < K; k++) {

        alphaProb[k] = alpha / K + nK_vec[k];
      }

      pi = myrng.rdirich(alphaProb);

      for (int l = 0; l < n; l++) {
        if (is_clamped[l]) {
          continue;
        }
        for (int k = 0; k < K; k++) {
          logPi(l, k) = log(pi[k]);
        }
      }

    } else {

      for (int l = 0; l < n; l++) {

        if (is_clamped[l]) {
          continue;
        }

        for (int k = 0; k < K; k++) {

          alphaProb[k] = alpha * Assignments[l][k] + nK_vec[k];
        }

        pi = myrng.rdirich(alphaProb);

        for (int k = 0; k < K; k++) {
          logPi(l, k) = log(pi[k]);
        }
      }
    }

    if (i > burnin) {
      piMean = piMean + pi;
    }
    if (i % 50 == 0) {
      cout << i << ' ' << flush;
    }
  } // end sampling with entropy thresholding

  cout << endl;

  // recompute MNprobMEAN for positions not selected

  for (int k = 0; k < K; k++)
    counts[k] = GD_lev;

  for (int o = 0; o < n; o++) {

    for (int r = 0; r < Reads[reads_in_window[o]].size(); r++) {

      idx_global = r + Positions_Start[reads_in_window[o]] - 1;
      idx = r + Positions_Start[reads_in_window[o]] - WindowStart;

      if (idx >= 0 && idx < GD && entropy_select[idx_global] == 0) {

        nuc = Reads[reads_in_window[o]][r];
        qual = int(quality_scores[reads_in_window[o]][r]) - qual_subtract;
        if (!include_deletions && nuc == 0)
          continue;

        if (nuc < 5) {

          counts[C[o]][idx][nuc] += qual;
        }
      }
    }
  }

  // add one "reference" read...

  for (int k = 0; k < K; k++) {
    for (int l = 0; l < GD; l++) {
      idx_global = l + WindowStart - 1;
      if (entropy_select[idx_global] == 0) {

        nuc = refSeq[l];
        if (nuc < 5)
          counts[k][l][nuc] += 1;
      }
    }
  }

  for (int k = 0; k < K; k++) {
    for (int l = 0; l < GD; l++) {
      idx_global = l + WindowStart - 1;
      if (entropy_select[idx_global] == 0) {
        for (int m = 0; m < levels; m++)
          probs(l, m) = counts[k][l][m] + 1e-3;
        probs(l, _) = probs(l, _) / sum(probs(l, _));

        for (int p = 0; p < levels; p++) {
          MNprobMEAN(l, k * levels + p) = probs(l, p);
        }
      }
    }
  }

  for (int k = 0; k < K; k++) {

    piMean[k] = piMean[k] / (nSample - burnin + 1);
  }

  for (int l = 0; l < GD; l++) {
    idx_global = l + WindowStart - 1;
    if (entropy_select[idx_global] == 1) {
      MNprobMEAN(l, _) = MNprobMEAN(l, _) / (nSample - burnin + 1);
    }

    for (int k = 0; k < K; k++) {

      double mt = -1e100;

      for (int p = 0; p < levels; p++) {
        if (MNprobMEAN(l, k * levels + p) > mt) {
          mt = MNprobMEAN(l, k * levels + p);
          maxTab(l, k) = p;
        }
      }
    }
  }

  return 0;
}

/**
 * @return The path to the FASTA file that the reconstructed haplotypes were
 * written to.
 */
std::string visualizeAlignments(
    int K, int WindowStart, int GD, const vector<int> &Positions_Start,
    const vector<vector<int>> &Reads, const vector<int> &reads_in_window,
    const vector<string> &IDs_in_window, const vector<string> &IDs,
    const vector<int> &strand, string prefix, string SQ_string,
    const vector<int> &entropy_selection, const vector<int> &C_index,
    const vector<int> &L, const Matrix<int> &maxTab, bool have_true_haplotypes,
    const vector<vector<int>> &Haplos, vector<int> Haplo_Positions_Start,
    const vector<string> &trueHaplo_IDs,
    const vector<vector<int>> &reconstructedHaplos,
    const vector<vector<double>> &reconstructedConf, vector<int> &match,
    vector<int> &match_costs, bool make_pgm,
    vector<double> reconstructedFrequency,
    const vector<vector<int>> &reconstructed_overlaps, bool plot_reads,
    const vector<vector<int>> &error_positions) {

  bool SNV_only = true;
  if (visualization_level > 1)
    SNV_only = false;

  if (visualization_level > 0) {
    vector<bool> is_SNV(GD);

    if (have_true_haplotypes) {

      for (int i = 0; i < GD; i++) {

        vector<int> hk(Haplos.size(), 0);
        for (int k = 0; k < Haplos.size(); k++) {
          int iH = i + WindowStart - Haplo_Positions_Start[k];
          if (iH >= 0 && iH < Haplos[k].size())
            hk[k] = Haplos[k][iH];
        }

        sort(hk.begin(), hk.end());

        vector<int>::iterator new_end = unique(hk.begin(), hk.end());

        hk.erase(new_end, hk.end());
        is_SNV[i] = false;
        if (hk.size() > 1) {
          is_SNV[i] = true;
        }
      }
    }

    ostringstream sVis1;
    sVis1 << prefix << "visuAlign_" << WindowStart << "_"
          << WindowStart + GD - 1 << ".html";
    string sVis = sVis1.str();

    ofstream ofAW(sVis.c_str(), ios::out);
    ofAW << fixed << setprecision(4);
    ofAW << "<html lang=\"eng\" xml:lang=\"eng\" "
            "xmlns=\"http://www.w3.org/1999/xhtml\">"
         << endl;
    ofAW << "<head>" << endl;
    ofAW << "<link href=\"nuc-colors.css\" rel=\"stylesheet\" "
            "type=\"text/css\" />"
         << endl;
    ofAW << "</head>" << endl;
    ofAW << "<body id=\"top\">" << endl;

    ofAW << "<pre id=\"alignmentContent\" xml:space=\"preserve\"  "
            "style=\"cursor:default\">"
         << endl;

    // consensus: entropy_selection
    ofAW << '\t';
    for (int i = 0; i < GD; i++) {
      if (SNV_only) {
        if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
          continue;
      }

      char nuc = SQ_string[i + WindowStart - 1];

      if (entropy_selection[i + WindowStart - 1]) {
        ofAW << "<span class=\"entropySel\">" << nuc << "</span>";
      } else
        ofAW << nuc;
    }

    ofAW << '\t' << "reference: entropy_selection" << '\t';
    ofAW << endl;

    if (have_true_haplotypes) {
      // consensus: true mutations

      ofAW << '\t';
      for (int i = 0; i < GD; i++) {

        if (SNV_only) {
          if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
            continue;
        }

        vector<int> hk(Haplos.size(), 0);
        for (int k = 0; k < Haplos.size(); k++) {
          int iH = i + WindowStart - Haplo_Positions_Start[k];
          if (iH >= 0 && iH < Haplos[k].size())
            hk[k] = Haplos[k][iH];
        }

        sort(hk.begin(), hk.end());

        vector<int>::iterator new_end = unique(hk.begin(), hk.end());

        hk.erase(new_end, hk.end());
        bool is_true = false;
        if (hk.size() > 1) {
          is_true = true;
        }
        char nuc = SQ_string[i + WindowStart - 1];

        if (is_true) {
          ofAW << "<span class=\"trueMut\">" << nuc << "</span>";
        } else
          ofAW << nuc;
      }

      ofAW << '\t' << "reference: true mutations" << '\t';
      ofAW << endl;

      ofAW << '\t';
      for (int i = 0; i < GD; i++) {
        if (SNV_only) {
          if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
            continue;
        }
        ofAW << "-";
      }

      ofAW << '\t';
      ofAW << endl;
    }
    // the consensus:
    ofAW << '\t';
    for (int i = 0; i < GD; i++) {

      if (SNV_only) {
        if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
          continue;
      }

      char nuc = SQ_string[i + WindowStart - 1];

      if (nuc == '-') {
        ofAW << "<a title=\"" << i + WindowStart << "\">-</a>";
      } else {
        if (nuc == 'A') {
          ofAW << "<span class=\"A\"><a title=\"" << i + WindowStart
               << "\">A</a></span>";
        } else {
          if (nuc == 'C') {
            ofAW << "<span class=\"C\"><a title=\"" << i + WindowStart
                 << "\">C</a></span>";
          } else {
            if (nuc == 'G') {
              ofAW << "<span class=\"G\"><a title=\"" << i + WindowStart
                   << "\">G</a></span>";
            } else {
              if (nuc == 'T') {
                ofAW << "<span class=\"T\"><a title=\"" << i + WindowStart
                     << "\">T</a></span>";
              } else {
                if (nuc == 'N') {
                  ofAW << "<span class=\"N\"><a title=\"" << i + WindowStart
                       << "\">N</a></span>";
                } else {
                  ofAW << "<a title=\"" << i + WindowStart << "\">.</a>";
                }
              }
            }
          }
        }
      }
    }

    ofAW << '\t' << "reference" << '\t' << endl;

    ofAW << '\t';
    for (int i = 0; i < GD; i++) {

      if (SNV_only) {
        if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
          continue;
      }

      ofAW << "-";
    }
    ofAW << '\t' << endl;

    ofAW << '\t';
    for (int i = 0; i < GD; i++) {
      if (SNV_only) {
        if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
          continue;
      }
      ofAW << "-";
    }
    ofAW << '\t' << endl;

    if (have_true_haplotypes) {
      // true haplotypes

      for (int k = 0; k < Haplos.size(); k++) {

        ofAW << '\t';

        for (int i = 0; i < GD; i++) {

          if (SNV_only) {
            if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
              continue;
          }

          int iH = i + WindowStart - Haplo_Positions_Start[k];
          if (iH < 0 || iH >= Haplos[k].size())
            ofAW << ".";
          else {
            int nuc = Haplos[k][iH];
            if (nuc == 0) {
              ofAW << "-";
            } else {
              if (nuc == 1) {
                ofAW << "<span class=\"A\">A</span>";
              } else {
                if (nuc == 2) {
                  ofAW << "<span class=\"C\">C</span>";
                } else {
                  if (nuc == 3) {
                    ofAW << "<span class=\"G\">G</span>";
                  } else {
                    if (nuc == 4) {
                      ofAW << "<span class=\"T\">T</span>";
                    } else {
                      if (nuc == 5) {
                        ofAW << "<span class=\"N\">N</span>";
                      } else {
                        ofAW << ".";
                      }
                    }
                  }
                }
              }
            }
          }
        }
        ofAW << '\t' << "true_" << k << ": " << trueHaplo_IDs[k] << '\t'
             << endl;

        ofAW << '\t';
        for (int i = 0; i < GD; i++) {
          if (SNV_only) {
            if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
              continue;
          }
          ofAW << "-";
        }

        ofAW << '\t' << endl;
      }

      ofAW << '\t';
      for (int i = 0; i < GD; i++) {
        if (SNV_only) {
          if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
            continue;
        }
        ofAW << "-";
      }

      ofAW << '\t' << endl;
    }

    for (int k = 0; k < reconstructedHaplos.size(); k++) {

      ofAW << '\t';
      for (int i = 0; i < GD; i++) {
        if (SNV_only) {
          if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
            continue;
        }
        int nuc = reconstructedHaplos[k][i];
        if (nuc == 0) {
          ofAW << "-";
        } else {
          if (nuc == 1) {
            ofAW << "<span class=\"A\">A</span>";
          } else {
            if (nuc == 2) {
              ofAW << "<span class=\"C\">C</span>";
            } else {
              if (nuc == 3) {
                ofAW << "<span class=\"G\">G</span>";
              } else {
                if (nuc == 4) {
                  ofAW << "<span class=\"T\">T</span>";
                } else {
                  if (nuc == 5) {
                    ofAW << "<span class=\"N\">N</span>";
                  } else {
                    ofAW << ".";
                  }
                }
              }
            }
          }
        }
      }

      ofAW << '\t' << "reconstructed_" << k
           << " Freq: " << reconstructedFrequency[k] << '\t';

      ofAW << "  Overlaps:";

      for (int i = 0; i < reconstructed_overlaps[k].size(); i++) {
        ofAW << i + 1 << ":";
        if (i == reconstructed_overlaps[k].size() - 1 &&
            reconstructed_overlaps[k][i] >= 10) {
          ofAW.width(3);
          ofAW << left << "<span class=\"G\">" << reconstructed_overlaps[k][i]
               << "</span>"
               << "|";
        } else {
          ofAW.width(3);
          ofAW << left << reconstructed_overlaps[k][i] << "|";
        }
      }
      if (have_true_haplotypes) {
        ofAW << '\t' << "Best match: true_" << match[k] << '\t'
             << "Costs: " << match_costs[k] << '\t';
        ofAW << "Error positions: ";
        if (error_positions[k].size() < 15) {
          for (int i = 0; i < error_positions[k].size(); i++)
            ofAW << error_positions[k][i] << "|";
          ofAW << '\t';
        } else
          ofAW << "> 15 errors" << '\t';
      }

      ofAW << '\t' << endl;

      ///////

      ofAW << '\t';
      for (int i = 0; i < GD; i++) {
        if (SNV_only) {
          if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
            continue;
        }
        double conf = reconstructedConf[k][i];
        if (conf > 0.99) {
          ofAW << "<span class=\"c99\">&hearts;</span>";
        } else {
          if (conf > 0.95) {
            ofAW << "<span class=\"c95\">&hearts;</span>";
          } else {
            if (conf > 0.9) {
              ofAW << "<span class=\"c90\">&diams;</span>";
            } else {
              if (conf > 0.85) {
                ofAW << "<span class=\"c85\">&diams;</span>";
              } else {
                if (conf > 0.8) {
                  ofAW << "<span class=\"c80\">*</span>";
                } else {
                  ofAW << "<span class=\"less80\">*</span>";
                }
              }
            }
          }
        }
      }
      ofAW << '\t' << endl;
      ///////////////////

      ofAW << '\t';
      for (int i = 0; i < GD; i++) {
        if (SNV_only) {
          if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
            continue;
        }
        ofAW << "-";
      }

      ofAW << '\t' << endl;
    }

    if (plot_reads) {
      int cc = 0;
      for (int k = 0; k < K; k++) {
        if (C_index[k] == 1) {
          ofAW << '\t';
          for (int i = 0; i < GD; i++) {
            if (SNV_only) {
              if (!entropy_selection[i + WindowStart - 1] && !is_SNV[i])
                continue;
            }
            ofAW << "-";
          }

          ofAW << '\t';
          ofAW << endl;
          for (int o = 0; o < reads_in_window.size(); o++) {
            string re(GD, '.');
            if (L[o] == k) {
              ofAW << '\t';

              for (int r = 0; r < Reads[reads_in_window[o]].size(); r++) {
                int nuc = Reads[reads_in_window[o]][r];
                int idx = r + Positions_Start[reads_in_window[o]] - WindowStart;

                if (idx >= 0 && idx < GD) {
                  stringstream ss(stringstream::out);
                  ss << nuc;

                  re[idx] = ss.str()[0];
                }
              }

              for (int l = 0; l < re.size(); l++) {

                if (SNV_only) {
                  if (!entropy_selection[l + WindowStart - 1] && !is_SNV[l])
                    continue;
                }

                if (re[l] == '0') {
                  ofAW << "-";
                } else {
                  if (re[l] == '1') {
                    ofAW << "<span class=\"A\">A</span>";
                  } else {
                    if (re[l] == '2') {
                      ofAW << "<span class=\"C\">C</span>";
                    } else {
                      if (re[l] == '3') {
                        ofAW << "<span class=\"G\">G</span>";
                      } else {
                        if (re[l] == '4') {
                          ofAW << "<span class=\"T\">T</span>";
                        } else {
                          if (re[l] == '5') {
                            ofAW << "<span class=\"N\">N</span>";
                          } else {
                            ofAW << ".";
                          }
                        }
                      }
                    }
                  }
                }
              }
              ofAW << '\t';
              ofAW << IDs_in_window[o] << '\t' << ", reconstructed_" << cc
                   << '\t';
              ofAW << endl;
            }
          }
          cc++;
        }
      }
    }

    ofAW << "</pre>" << endl;

    ofAW << " </body>" << endl;
    ofAW << "</html>" << endl;
    ofAW.close();
  }

  if (make_pgm) {
    ostringstream sVis2;
    sVis2 << prefix << "visuAlign_" << WindowStart << "_"
          << WindowStart + GD - 1 << ".pgm";
    string s = sVis2.str();
    int n = reads_in_window.size();
    int m = GD;
    ofstream of2(s.c_str(), ios::out);
    of2 << "P2" << endl;
    of2 << "# dmat.pgm " << endl;
    of2 << m << ' ' << n << endl;
    of2 << 255 << endl;
    int counter = 0;
    for (int k = 0; k < K; k++) {
      if (C_index[k] == 1) {
        counter++;

        for (int o = 0; o < reads_in_window.size(); o++) {
          string re(GD, '.');
          if (L[o] == k) {

            for (int r = 0; r < Reads[reads_in_window[o]].size(); r++) {
              int nuc = Reads[reads_in_window[o]][r];
              int idx = r + Positions_Start[reads_in_window[o]] - WindowStart;

              if (idx >= 0 && idx < GD) {
                stringstream ss(stringstream::out);
                ss << nuc;

                re[idx] = ss.str()[0];
              }
            }

            for (int l = 0; l < re.size(); l++) {
              if (re[l] != '.') {
                of2 << 0 << ' ';
              } else {
                if (counter % 2 == 0)
                  of2 << 255 << ' ';
                else
                  of2 << 200 << ' ';
              }
            }
            of2 << endl;
          }
        }
      }
    }

    of2.close();
  }

  // export reconstructed haplotypes in fasta format

  const auto fasta_output_path = prefix + std::to_string(WindowStart) + "_" +
                                 std::to_string(WindowStart + GD - 1) + ".fas";

  ofstream ofFASTA(fasta_output_path.c_str(), ios::out);

  for (int k = 0; k < reconstructedHaplos.size(); k++) {

    ofFASTA << ">"
            << "reconstructed_" << k << endl;
    ofFASTA << ";Freq:" << reconstructedFrequency[k] << endl;
    ofFASTA << ";Overlap quality scores (5,10):";
    ofFASTA << reconstructed_overlaps[k][4] << ","
            << reconstructed_overlaps[k][9] << endl;
    if (have_true_haplotypes) {
      ofFASTA << ";Best match:true_" << match[k] << endl;
      ofFASTA << ";Costs:" << match_costs[k] << endl;
    }
    ofFASTA << ";Confidence scores" << endl;
    ofFASTA << ";";
    for (int i = 0; i < GD; i++) {
      if (i > 0 && i % 69 == 0)
        ofFASTA << endl << ";";
      double rconf = reconstructedConf[k][i];

      int conf = floor(1.0 / 9.0 * (pow(10, rconf) - 1) * (126 - 33)) + 33;

      char myChar;
      myChar = (char)conf;
      ofFASTA << myChar;
    }
    ofFASTA << endl;
    ofFASTA << ";EndOfComments" << endl;

    for (int i = 0; i < GD; i++) {
      if (i > 0 && i % 70 == 0)
        ofFASTA << endl;
      int nuc = reconstructedHaplos[k][i];
      if (nuc == 0) {
        ofFASTA << "-";
      } else {
        if (nuc == 1) {
          ofFASTA << "A";
        } else {
          if (nuc == 2) {
            ofFASTA << "C";
          } else {
            if (nuc == 3) {
              ofFASTA << "G";
            } else {
              if (nuc == 4) {
                ofFASTA << "T";
              } else {
                if (nuc == 5) {
                  ofFASTA << "N";
                } else {
                  ofFASTA << ".";
                }
              }
            }
          }
        }
      }
    }
    ofFASTA << endl;
  }
  ofFASTA.close();

  return fasta_output_path;
}

int local_Analysis(
    int min_overlap, int K, int nSample, int nT, int max_reads_in_window,
    int WindowStart, int WindowIncrement, vector<vector<int>> &WindowStartStop,
    int max_seq_length, int min_seq_start, const vector<int> &Positions_Start,
    const vector<vector<int>> &Reads, bool have_quality_scores,
    const vector<string> &quality_scores, const vector<string> &IDs,
    const vector<int> &strand, string prefix, string SQ_string,
    bool have_true_haplotypes, const vector<vector<int>> &trueHaplos,
    const vector<int> &trueHaplo_Positions_Start,
    const vector<string> &trueHaplo_IDs, vector<int> &foundClusters,
    double entropy_threshold, double entropy_fraction,
    const vector<int> &entropy_select, double alphaMN) {

  mersenne myrng;
  myrng.initialize(time(NULL));

  vector<int> reads_in_window;
  vector<string> IDs_in_window;

  int WindowStop = WindowStart + WindowIncrement;

  bool last_run = false;
  int iter = 0;

  if (WindowStop > max_seq_length) {
    WindowStop = max_seq_length;
    last_run = true;
  }

  while (WindowStop <= max_seq_length) {

    vector<int> startStop(2);
    startStop[0] = WindowStart;
    startStop[1] = WindowStop;
    WindowStartStop.push_back(startStop);

    int GD = WindowStop - WindowStart + 1;

    cout << "Local analysis in window from " << WindowStart << " to "
         << WindowStop << endl;

    reads_in_window.clear();
    IDs_in_window.clear();

    vector<int> coverage(WindowStop - WindowStart + 1, 0);
    int i_count = 0;
    for (int i = 0; i < nT; i++) {

      if (Positions_Start[i] < WindowStop - min_overlap &&
          Positions_Start[i] + Reads[i].size() > WindowStart + min_overlap) {
        reads_in_window.push_back(i);

        IDs_in_window.push_back(IDs[i]);

        for (int j = WindowStart; j < Positions_Start[i] + Reads[i].size(); j++)

          if (j - WindowStart < coverage.size())
            coverage[j - WindowStart] += 1;
      }
    }

    // cout << "coverage in current window: "<<endl;
    int min_cov = coverage[0];
    for (int i = 0; i < coverage.size(); i++) {
      // cout << coverage[i]<<"|";
      if (coverage[i] < min_cov)
        min_cov = coverage[i];
    }

    cout << endl << "Minimum coverage in current window: " << min_cov << endl;

    vector<int> v(reads_in_window.size());

    for (int k = 0; k < reads_in_window.size(); k++) {

      v[k] = reads_in_window[k];
    }

    std::shuffle(v.begin(), v.end(), std::mt19937{std::random_device{}()});

    int number_of_reads_in_window = reads_in_window.size();

    cout << "Number of reads in window: " << number_of_reads_in_window
         << " Selected upper bound: " << max_reads_in_window << endl;
    if (number_of_reads_in_window > max_reads_in_window) {
      number_of_reads_in_window = max_reads_in_window;

      reads_in_window.clear();
      IDs_in_window.clear();

      for (int j = 0; j < number_of_reads_in_window; j++) {
        reads_in_window.push_back(v[j]);
        IDs_in_window.push_back(IDs[v[j]]);
      }
    }

    int n = reads_in_window.size(); // number of reads in window

    double alpha = K / log(double(n));

    Matrix<double> MNprob(GD, K * levels); // log-scale (+1e-6)
    Matrix<double> MNprobMEAN(GD, K * levels);

    Matrix<int> maxTab(GD, K);

    Matrix<double> pi(K, 1);
    Matrix<double> piMean(K, 1);
    vector<int> C(n, 0);

    for (int i = 0; i < n; i++)
      C[i] = rand() % K;

    Matrix<int> CMAX(n, K);
    vector<double> uniform(K, 1.0 / K);
    vector<vector<double>> Assignments(n, uniform);

    int lseed = rand() % 123456;
    double burnin_fraction = 0.4;

    vector<int> refSeq(GD);

    for (int i = 0; i < GD; i++) {
      char nuc = SQ_string[i + WindowStart - 1];
      // cout << nuc;
      if (nuc == '-')
        refSeq[i] = 0;
      else if (nuc == 'A')
        refSeq[i] = 1;
      else if (nuc == 'C')
        refSeq[i] = 2;
      else if (nuc == 'G')
        refSeq[i] = 3;
      else if (nuc == 'T')
        refSeq[i] = 4;
      else
        refSeq[i] = 5;
    }
    // cout << endl;
    MNprobMEAN = 0;
    MultiNomialDPMReadsSemiEntropy(nSample,

                                   Reads,

                                   quality_scores, Positions_Start,

                                   reads_in_window,

                                   Assignments, WindowStart, WindowStop,

                                   K,

                                   alphaMN, alpha,

                                   MNprob, //(GD,K*levels) log-scale (+1e-6)
                                   MNprobMEAN,
                                   maxTab, //(GD,2);
                                   pi,     //(K,1)
                                   piMean, C,
                                   CMAX, // (n,K)

                                   lseed, burnin_fraction, refSeq,

                                   entropy_select, true);

    Matrix<double> maxConf(GD, K);

    double lconf_scale = 0.75;
    vector<vector<int>> reconstructedHaplos;
    vector<vector<double>> reconstructedConf;
    vector<double> reconstructedFrequency;

    vector<int> C_index(K, 0);
    vector<int> L(n);
    vector<double> L_conf(n);
    for (int i = 0; i < n; i++) {
      int cm = -1;
      double c_sum = 0;
      for (int k = 0; k < K; k++) {
        if (CMAX(i, k) > cm) {
          cm = CMAX(i, k);
          L[i] = k;
        }

        c_sum += CMAX(i, k);
      }
      C_index[L[i]] = 1;

      L_conf[i] = lconf_scale * cm / c_sum;
    }

    //////////////////////////////////////

    Matrix<double> freq_M(K, levels);
    freq_M = 0;
    vector<Matrix<double>> freq(GD, freq_M);
    vector<int> nK(K, 0);
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < n; i++) {
        if (L[i] == k) { // the reads belonging to cluster k
          nK[k] += 1;
          for (int p = 0; p < GD; p++) {

            int idx = p - Positions_Start[reads_in_window[i]] + WindowStart;

            if (idx >= 0 && idx < Reads[reads_in_window[i]].size()) {
              int l = Reads[reads_in_window[i]][idx];
              if (include_deletions) {
                if (l < 5)
                  freq[p](k, l) += 1;
              } else {
                if (l > 0 && l < 5)
                  freq[p](k, l) += 1;
              }
            }
          }
        }
      }
      if (nK[k] == 0)
        nK[k] = 1;
      for (int p = 0; p < GD; p++) {
        freq[p](k, _) /= sum(freq[p](k, _)); // nK[k];
      }
    }

    for (int k = 0; k < K; k++) {
      for (int l = 0; l < GD; l++) {
        maxConf(l, k) = freq[l](k, maxTab(l, k));
      }
    }

    // quality scoring
    vector<vector<int>> reconstructed_overlaps;
    bool count_only_selected_positions = true;

    for (int k = 0; k < K; k++) {
      if (C_index[k] == 1) {

        int overlap = 20;
        if (GD - 2 * overlap < 1)
          overlap = 0.5 * (GD - 1);

        Matrix<> over_counts(10, GD - 2 * overlap);
        vector<int> overlap_rank(over_counts.rows(), 0);

        over_counts = 0;

        for (int p = overlap; p < GD - overlap; p++) {

          vector<int> over_vec;
          for (int o = 0; o < n; o++) {
            if (L[o] == k) {

              int idx_start = Positions_Start[reads_in_window[o]] - WindowStart;
              int idx_stop = Reads[reads_in_window[o]].size() - 1 +
                             Positions_Start[reads_in_window[o]] - WindowStart;

              if (idx_start <= p && idx_stop >= p) {

                int over = p - idx_start;
                if (idx_stop - p < over)
                  over = idx_stop - p;

                /////////////////
                if (count_only_selected_positions) {
                  int over_1 = 0;
                  for (int q = idx_start; q <= p; q++) {
                    int idx_global = q + WindowStart - 1;
                    if (idx_global >= 0) {
                      if (entropy_select[idx_global] == 1) {
                        over_1++;
                      }
                    }
                  }
                  int over_2 = 0;
                  for (int q = p; q <= idx_stop; q++) {
                    int idx_global = q + WindowStart - 1;
                    if (idx_global < max_seq_length) {
                      if (entropy_select[idx_global] == 1) {
                        over_2++;
                      }
                    }
                  }
                  over = over_1;
                  if (over_2 < over_1)
                    over = over_2;
                }

                //////////////////

                over_vec.push_back(over);
              }
            }
          }
          sort(over_vec.begin(), over_vec.end(), myfunction);
          for (int i = 0; i < over_counts.rows(); i++)
            if (i < over_vec.size())
              over_counts(i, p - overlap) = over_vec[i];
        }
        for (int i = 0; i < over_counts.rows(); i++) {
          overlap_rank[i] = min(over_counts(i, _));
          // cout << i<<':'<<overlap_rank[i]<<'|';
        }
        // cout << endl;

        reconstructed_overlaps.push_back(overlap_rank); // min_count_fraction);
      }
    }
    ///////////////////////////////

    double sum_freq = 0.0;
    for (int k = 0; k < K; k++) {
      if (C_index[k] == 1) {

        vector<double> conf(GD);
        vector<int> rec(GD);
        for (int i = 0; i < GD; i++) {
          rec[i] = maxTab(i, k);
          conf[i] = maxConf(i, k);
        }
        reconstructedHaplos.push_back(rec);
        reconstructedConf.push_back(conf);
        sum_freq += piMean(k, 0);
        reconstructedFrequency.push_back(piMean(k, 0));
      }
    }
    for (int k = 0; k < reconstructedFrequency.size(); k++) {
      reconstructedFrequency[k] = reconstructedFrequency[k] / sum_freq;
    }

    ostringstream s1;
    s1 << prefix << "hiv_labels_" << WindowStart << ".lab";
    string s2 = s1.str();

    ofstream ofLab(s2.c_str(), ios::out);
    for (int j = 0; j < n; j++) {
      ofLab << L[j] << "|" << L_conf[j] << endl;
    }
    ofLab.close();

    ostringstream s3;
    s3 << prefix << "hiv_reads_" << WindowStart << ".reads";
    string s4 = s3.str();
    ofstream ofReads(s4.c_str(), ios::out);
    for (int j = 0; j < reads_in_window.size(); j++) {
      ofReads << reads_in_window[j] << endl;
    }
    ofReads.close();

    vector<double> strand_K(reconstructedHaplos.size(), 0);

    int good_clusters = 0;
    int cc = 0;
    for (int k = 0; k < K; k++) {
      int kcount = 0;
      if (C_index[k] == 1) {

        for (int o = 0; o < reads_in_window.size(); o++) {

          if (L[o] == k) {

            strand_K[cc] += strand[reads_in_window[o]];
            kcount++;
          }
        }

        strand_K[cc] /= (16.0);

        strand_K[cc] = strand_K[cc] / (1e-5 + kcount - strand_K[cc]);
        if (reconstructed_overlaps[cc][reconstructed_overlaps[cc].size() - 1] >
            20) {
          good_clusters++;
        } else {
          good_clusters--;
        }

        if (reconstructed_overlaps[cc][reconstructed_overlaps[cc].size() - 1] <
            5) {
          good_clusters -= 2;
        }

        cc++;
      }
    }

    foundClusters.push_back(good_clusters);

    vector<int> match(reconstructedHaplos.size(), -1);
    vector<int> match_costs(reconstructedHaplos.size(), 30000);

    vector<vector<int>> error_positions;
    vector<vector<int>> diff_vec;
    if (have_true_haplotypes) {

      for (int k = 0; k < reconstructedHaplos.size(); k++) {

        double minD = 1e100;
        vector<int> dv(trueHaplos.size());
        // cout << k<<": " ;
        for (int l = 0; l < trueHaplos.size(); l++) {

          int diff = 0;

          for (int i = 0; i < GD; i++) {
            int iH = WindowStart - trueHaplo_Positions_Start[l] + i;
            if (iH >= 0 && iH < trueHaplos[l].size() && trueHaplos[l][iH] < 5 &&
                reconstructedHaplos[k][i] < 5 &&
                reconstructedHaplos[k][i] != trueHaplos[l][iH])
              diff++;
          }

          dv[l] = diff;
          // cout << diff <<' ';
        }
        // cout << endl;
        diff_vec.push_back(dv);
      }

      for (int l = 0; l < reconstructedHaplos.size(); l++) {

        int min_match = diff_vec[l][0];
        int min_i = 0;
        for (int k = 1; k < trueHaplos.size(); k++) {

          if (diff_vec[l][k] < min_match) {
            min_match = diff_vec[l][k];
            min_i = k;
          }
        }
        match[l] = min_i;
        match_costs[l] = min_match;

        cout << "reconstructed_" << l << '\t'
             << "Freq: " << reconstructedFrequency[l] << '\t' << "Overlaps:";

        for (int i = 0; i < reconstructed_overlaps[l].size(); i++) {
          cout << i + 1 << ":";
          cout.width(3);
          cout << left << reconstructed_overlaps[l][i] << "|";
        }

        cout << '\t' << "Best match:" << min_i << "  Costs:" << min_match
             << '\t';

        vector<int> e_positions;
        for (int i = 0; i < GD; i++) {
          int iH = WindowStart - trueHaplo_Positions_Start[min_i] + i;
          if (iH >= 0 && iH < trueHaplos[min_i].size() &&
              reconstructedHaplos[l][i] < 5 && trueHaplos[min_i][iH] < 5 &&
              reconstructedHaplos[l][i] != trueHaplos[min_i][iH]) {
            cout << i + WindowStart << ": " << reconstructedHaplos[l][i]
                 << trueHaplos[min_i][iH] << "|";
            e_positions.push_back(i + WindowStart);
          }
        }
        error_positions.push_back(e_positions);
        cout << endl;
      }

    } else {

      for (int l = 0; l < reconstructedHaplos.size(); l++) {

        cout << "reconstructed_" << l << '\t'
             << "Freq: " << reconstructedFrequency[l] << '\t' << "Overlaps:";

        for (int i = 0; i < reconstructed_overlaps[l].size(); i++) {
          cout << i + 1 << ":";
          cout.width(3);
          cout << left << reconstructed_overlaps[l][i] << "|";
        }
        cout << endl;
      }
    }

    string this_prefix = prefix + "local_";

    visualizeAlignments(
        K, WindowStart, GD, Positions_Start, Reads, reads_in_window,
        IDs_in_window, IDs, strand, this_prefix, SQ_string, entropy_select,
        C_index, L, maxTab, have_true_haplotypes, trueHaplos,
        trueHaplo_Positions_Start, trueHaplo_IDs, reconstructedHaplos,
        reconstructedConf, match, match_costs, false, reconstructedFrequency,
        reconstructed_overlaps, true, error_positions);

    cout << "=========================================================="
         << endl;

    //  double increment =  myrng.runif()*30 + 0.175*WindowIncrement;

    double increment = myrng.runif() * 30 + 0.4 * WindowIncrement;

    WindowStart = WindowStart + int(increment);
    WindowStop = WindowStop + int(increment);

    if (last_run)
      break;

    if (WindowStop >= max_seq_length) {
      WindowStop = max_seq_length;
      last_run = true;
    }

    if (WindowStart > max_seq_length - 200) {
      WindowStart = max_seq_length - 200;
    }

    if (WindowStart <= min_seq_start)
      WindowStart = min_seq_start;

    iter++;
  }

  return 0;
}

int reconstruct_global(
    int min_overlap, int K, int nSample, vector<vector<int>> WindowStartStop,
    int WindowIncrement, int select_start, int nT, int max_seq_length,
    int min_seq_start, const vector<int> &Positions_Start,
    const vector<vector<int>> &Reads, bool have_quality_scores,
    const vector<string> &quality_scores, vector<string> &IDs,
    const vector<int> &strand, string prefix, string SQ_string,
    bool have_true_haplotypes, const vector<vector<int>> &trueHaplos,
    const vector<int> &trueHaplo_Positions_Start,
    const vector<string> &trueHaplo_IDs, double entropy_threshold,
    double entropy_fraction, const vector<int> &entropy_select,
    const std::optional<std::string> &reconstructed_haplos_path) {

  string line;
  vector<string> tokens;

  // read in all used reads
  vector<int> readsIndicator(nT, 1);

  string isInWindow(nT, 'n');

  int WindowStart = WindowStartStop[select_start][0];
  int WindowStop = WindowStartStop[select_start][1];

  ostringstream s1;
  s1 << prefix << "hiv_labels_" << WindowStart << ".lab";
  string s2 = s1.str();

  ifstream ifLab(s2.c_str(), ios::in);

  vector<int> L;
  vector<double> L_conf;
  while (getline(ifLab, line, '\n')) {
    tokens = tokenize(line, "|");
    L.push_back(atoi(tokens[0].c_str()));
    L_conf.push_back(atoi(tokens[1].c_str()));
  }
  ifLab.close();

  ostringstream s3;
  s3 << prefix << "hiv_reads_" << WindowStart << ".reads";
  string s4 = s3.str();
  ifstream ifReads(s4.c_str(), ios::in);

  vector<int> reads_in_window;
  vector<string> IDs_in_window;
  while (getline(ifReads, line, '\n')) {
    int i = atoi(line.c_str());
    reads_in_window.push_back(i);
    IDs_in_window.push_back(IDs[i]);
    isInWindow[i] = 'y';
  }
  ifReads.close();

  int n = reads_in_window.size();
  int n_init = n;

  vector<int> nK_vec(K, 0);
  for (int j = 0; j < n; j++) {
    nK_vec[L[j]] = nK_vec[L[j]] + 1;
  }

  vector<int> K_match(K);

  int K_count = 0;
  for (int k = 0; k < K; k++) {

    if (nK_vec[k] == 0) {
      K_match[k] = -1;
      continue;
    }
    cout << k << ' ' << nK_vec[k] << ' ';
    cout << K_count << endl;
    K_match[k] = K_count;
    K_count++;
  }

  vector<int> L_new(n);
  for (int j = 0; j < n; j++) {
    L_new[j] = K_match[L[j]];
  }
  L = L_new;
  K = K_count;

  vector<double> almostZero(K, 1e-1 / K);
  vector<vector<double>> Assignments(n, almostZero);
  for (int j = 0; j < n; j++) {
    Assignments[j][L[j]] = 1 - (K - 1) * almostZero[0];
  }

  int GD = WindowStop - WindowStart + 1;

  double alphaMN = 1;
  double alpha = n / K;

  Matrix<double> pi(K, 1);
  Matrix<double> piMean(K, 1);
  vector<int> C = L;

  Matrix<int> CMAX(n, K);

  int lseed = rand() % 123456;
  double burnin_fraction = 0.4;

  //////////////// increase window ///////////////

  bool make_pgm = false;
  mersenne myrng;
  myrng.initialize(lseed);

  int first = 1;

  const auto continue_increasing_window = [&]() {
    return WindowStart > min_seq_start || WindowStop < max_seq_length - 1;
  };

  while (continue_increasing_window()) {
    int n_old = n;

    double increment = myrng.runif() * 50 + 0.4 * WindowIncrement;

    if (!first) {
      WindowStart = WindowStart - int(increment);
      WindowStop = WindowStop + int(increment);
    }
    first = 0;
    if (WindowStop >= max_seq_length) {
      WindowStop = max_seq_length - 1;
    }

    if (WindowStart <= min_seq_start)
      WindowStart = min_seq_start;

    if (WindowStop == max_seq_length - 1 && WindowStart == min_seq_start) {
      // make_pgm = true;
    }

    cout << "================================================================"
         << endl;
    cout << "Global analysis in window from " << WindowStart << " to "
         << WindowStop << endl;

    vector<double> uniform(K, 1.0 / K);

    vector<int>::iterator it;

    for (int i = 0; i < nT; i++) {
      if (readsIndicator[i] == 1) {
        if (Positions_Start[i] < WindowStop - min_overlap &&
            Positions_Start[i] + Reads[i].size() > WindowStart + min_overlap) {

          if (isInWindow[i] == 'n') {
            // it = find (reads_in_window.begin(), reads_in_window.end(), i);
            // if(it ==  reads_in_window.end()){
            // is new...

            reads_in_window.push_back(i);

            isInWindow[i] = 'y';

            IDs_in_window.push_back(IDs[i]);

            L.push_back(-1);
            L_conf.push_back(0);
            Assignments.push_back(uniform);
          }
        }
      }
    }

    n = reads_in_window.size(); // number of reads in window

    C = L;
    for (int i = n_old; i < n; i++)
      C[i] = rand() % K;
    ;

    GD = WindowStop - WindowStart + 1;

    Matrix<double> MNprob(GD, K * levels); // log-scale (+1e-6)
    Matrix<double> MNprobMEAN(GD, K * levels);

    Matrix<int> maxTab(GD, K);

    Matrix<int> CMAX(n, K);

    int lseed = rand() % 123456;

    vector<int> refSeq(GD);

    for (int i = 0; i < GD; i++) {
      char nuc = SQ_string[i + WindowStart - 1];
      if (nuc == '-')
        refSeq[i] = 0;
      else if (nuc == 'A')
        refSeq[i] = 1;
      else if (nuc == 'C')
        refSeq[i] = 2;
      else if (nuc == 'G')
        refSeq[i] = 3;
      else if (nuc == 'T')
        refSeq[i] = 4;
      else
        refSeq[i] = 5;
    }

    MultiNomialDPMReadsSemiEntropy(nSample,

                                   Reads,

                                   quality_scores, Positions_Start,

                                   reads_in_window,

                                   Assignments, WindowStart, WindowStop,

                                   K,

                                   alphaMN, alpha,

                                   MNprob, //(GD,K*levels) log-scale (+1e-6)
                                   MNprobMEAN,
                                   maxTab, //(GD,2);
                                   pi,     //(K,1)
                                   piMean, C,
                                   CMAX, // (n,K)

                                   lseed, burnin_fraction,

                                   refSeq,

                                   entropy_select, false);

    ////////////////
    Matrix<double> maxConf(GD, K);

    Matrix<double> freq_M(K, levels);
    freq_M = 0;
    vector<Matrix<double>> freq(GD, freq_M);
    vector<int> nK(K, 0);
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < n; i++) {
        if (L[i] == k) { // the reads belonging to cluster k
          nK[k] += 1;
          for (int p = 0; p < GD; p++) {

            int idx = p - Positions_Start[reads_in_window[i]] + WindowStart;

            if (idx >= 0 && idx < Reads[reads_in_window[i]].size()) {
              int l = Reads[reads_in_window[i]][idx];
              if (include_deletions) {
                if (l < 5)
                  freq[p](k, l) += 1;
              } else {
                if (l > 0 && l < 5)
                  freq[p](k, l) += 1;
              }
            }
          }
        }
      }
      if (nK[k] == 0)
        nK[k] = 1;
      for (int p = 0; p < GD; p++) {
        freq[p](k, _) /= sum(freq[p](k, _)); // nK[k];
      }
    }

    for (int k = 0; k < K; k++) {
      for (int l = 0; l < GD; l++) {
        maxConf(l, k) = freq[l](k, maxTab(l, k));
      }
    }

    double assign_regularizer = 1e-5; // (1- burnin_fraction)*nSample*1e-1;

    for (int i = 0; i < n; i++) { // new
      int cm = -1;
      double c_sum = 0;
      for (int k = 0; k < K; k++) {
        if (CMAX(i, k) > cm) {
          cm = CMAX(i, k);
          L[i] = k;
        }
        c_sum += CMAX(i, k);
      }
      L_conf[i] = 0.6 * cm / c_sum;

      for (int k = 0; k < K; k++) {
        Assignments[i][k] = (CMAX(i, k) + assign_regularizer) /
                            (K * assign_regularizer + c_sum);
      }
    }

    n_init = n;

    vector<int> C_index(K, 0);

    for (int i = 0; i < n; i++) {
      C_index[L[i]] = 1;
    }

    // quality scoring
    vector<vector<int>> reconstructed_overlaps;

    for (int k = 0; k < K; k++) {
      if (C_index[k] == 1) {

        int overlap = 0.5 * WindowIncrement;

        if (GD - 2 * overlap < 1)
          overlap = 0.5 * (GD - 1);

        Matrix<> over_counts(10, GD - 2 * overlap);
        vector<int> overlap_rank(over_counts.rows(), 0);

        bool count_only_selected_positions = true;
        over_counts = 0;

        for (int p = overlap; p < GD - overlap; p++) {

          vector<int> over_vec;
          for (int o = 0; o < n; o++) {
            if (L[o] == k) {

              int idx_start = Positions_Start[reads_in_window[o]] - WindowStart;
              int idx_stop = Reads[reads_in_window[o]].size() - 1 +
                             Positions_Start[reads_in_window[o]] - WindowStart;

              if (idx_start <= p && idx_stop >= p) {

                int over = p - idx_start;
                if (idx_stop - p < over)
                  over = idx_stop - p;

                /////////////////
                if (count_only_selected_positions) {
                  int over_1 = 0;
                  for (int q = idx_start; q <= p; q++) {
                    int idx_global = q + WindowStart - 1;
                    if (idx_global >= 0) {
                      if (entropy_select[idx_global] == 1) {
                        over_1++;
                      }
                    }
                  }
                  int over_2 = 0;
                  for (int q = p; q <= idx_stop; q++) {
                    int idx_global = q + WindowStart - 1;
                    if (idx_global < max_seq_length) {
                      if (entropy_select[idx_global] == 1) {
                        over_2++;
                      }
                    }
                  }
                  over = over_1;
                  if (over_2 < over_1)
                    over = over_2;
                }

                //////////////////

                over_vec.push_back(over);
              }
            }
          }
          sort(over_vec.begin(), over_vec.end(), myfunction);
          for (int i = 0; i < over_counts.rows(); i++)
            if (i < over_vec.size())
              over_counts(i, p - overlap) = over_vec[i];
        }
        for (int i = 0; i < over_counts.rows(); i++) {
          overlap_rank[i] = min(over_counts(i, _));
        }

        reconstructed_overlaps.push_back(overlap_rank); // min_count_fraction);
      }
    }

    vector<double> reconstructedFrequency;
    vector<vector<int>> reconstructedHaplos;
    vector<vector<double>> reconstructedConf;
    double sum_freq = 0.0;
    for (int k = 0; k < K; k++) {
      if (C_index[k] == 1) {

        vector<int> rec(GD);
        vector<double> conf(GD);
        for (int i = 0; i < GD; i++) {
          rec[i] = maxTab(i, k);
          conf[i] = maxConf(i, k);
        }

        reconstructedConf.push_back(conf);

        reconstructedHaplos.push_back(rec);
        reconstructedFrequency.push_back(piMean(k, 0));
        sum_freq += piMean(k, 0);
      }
    }
    for (int k = 0; k < reconstructedFrequency.size(); k++) {
      reconstructedFrequency[k] = reconstructedFrequency[k] / sum_freq;
    }

    vector<int> match(reconstructedHaplos.size(), -1);
    vector<int> match_costs(reconstructedHaplos.size(), 30000);

    vector<vector<int>> error_positions;
    vector<vector<int>> diff_vec;
    if (have_true_haplotypes) {

      for (int k = 0; k < reconstructedHaplos.size(); k++) {

        double minD = 1e100;
        vector<int> dv(trueHaplos.size());
        for (int l = 0; l < trueHaplos.size(); l++) {

          int diff = 0;

          for (int i = 0; i < GD; i++) {
            int iH = WindowStart - trueHaplo_Positions_Start[l] + i;
            if (iH >= 0 && iH < trueHaplos[l].size() && trueHaplos[l][iH] < 5 &&
                reconstructedHaplos[k][i] < 5 &&
                reconstructedHaplos[k][i] != trueHaplos[l][iH])
              diff++;
          }

          dv[l] = diff;
        }
        diff_vec.push_back(dv);
      }

      for (int l = 0; l < reconstructedHaplos.size(); l++) {

        int min_match = diff_vec[l][0];
        int min_i = 0;
        for (int k = 1; k < trueHaplos.size(); k++) {

          if (diff_vec[l][k] < min_match) {
            min_match = diff_vec[l][k];
            min_i = k;
          }
        }
        match[l] = min_i;
        match_costs[l] = min_match;
        cout << "reconstructed_" << l << "  Freq: " << reconstructedFrequency[l]
             << "  Overlaps:";

        for (int i = 0; i < reconstructed_overlaps[l].size(); i++) {
          cout << i + 1 << ":";
          cout.width(3);
          cout << left << reconstructed_overlaps[l][i] << "|";
        }

        cout << '\t' << "Best match:" << min_i << "  Costs:" << min_match
             << '\t';

        vector<int> e_positions;
        for (int i = 0; i < GD; i++) {
          int iH = WindowStart - trueHaplo_Positions_Start[min_i] + i;
          if (iH >= 0 && iH < trueHaplos[min_i].size() &&
              reconstructedHaplos[l][i] < 5 && trueHaplos[min_i][iH] < 5 &&
              reconstructedHaplos[l][i] != trueHaplos[min_i][iH]) {
            cout << i + WindowStart << ": " << reconstructedHaplos[l][i]
                 << trueHaplos[min_i][iH] << "|";
            e_positions.push_back(i + WindowStart);
          }
        }
        error_positions.push_back(e_positions);
        cout << endl;
      }

    } else {

      for (int l = 0; l < reconstructedHaplos.size(); l++) {

        cout << "reconstructed_" << l << "  Freq: " << reconstructedFrequency[l]
             << "  Overlaps:";

        for (int i = 0; i < reconstructed_overlaps[l].size(); i++) {
          cout << i + 1 << ":";
          cout.width(3);
          cout << left << reconstructed_overlaps[l][i] << "|";
        }
        cout << endl;
      }
    }

    const auto this_prefix = prefix + "global_";

    const auto fasta_output_path = visualizeAlignments(
        K, WindowStart, GD, Positions_Start, Reads, reads_in_window,
        IDs_in_window, IDs, strand, this_prefix, SQ_string, entropy_select,
        C_index, L, maxTab, have_true_haplotypes, trueHaplos,
        trueHaplo_Positions_Start, trueHaplo_IDs, reconstructedHaplos,
        reconstructedConf, match, match_costs, make_pgm, reconstructedFrequency,
        reconstructed_overlaps, false, error_positions);

    if (reconstructed_haplos_path && !continue_increasing_window()) {
      try {
        std::filesystem::copy_file(fasta_output_path,
                                   *reconstructed_haplos_path);
      } catch (const std::filesystem::filesystem_error &error) {
        std::cerr << "Warning: Copying the reconstructed haplotypes to '"
                  << *reconstructed_haplos_path << "' failed ('" << error.what()
                  << "')." << std::endl;
      }
    }
  }
  return 0;
}

int main(int argc, char *argv[]) {
  try {
    qual_subtract = 33;

    srand(time(NULL));

    string prefix = "predicthaplo_output/ph_";
    string cons;

    std::optional<std::string> FASTAreads;
    bool have_true_haplotypes = true;
    string FASTAhaplos;
    std::optional<std::string> reconstructed_haplos_path;

    bool do_local_Analysis = true;
    double entropy_threshold = 4e-2;
    double entropy_fraction = 0.25;

    string line, line_stats, line_ID;
    string tok = ":";
    vector<string> tokens, tokens2, tokens_as;
    int max_reads_in_window = 10000;

    string prefix_extension = "";

    int reconstruction_start, reconstruction_stop;
    double alpha_MN_local =
        25; // prior parameter for multinomial tables over the nucleotides
    double max_gap_fraction = 0.05;         // relative to alignment length
    double min_align_score_fraction = 0.35; // relative to read length

    double min_qual = 30;
    int min_length = 220;

    double min_overlap_factor =
        0.85; // reads must have an overlap with the local reconstruction window
              // of at least this factor times the window size
    double local_window_size_factor =
        0.7; // size of  local reconstruction window relative to the median of
             // the read lengths

    double mismatch = 3.0, gap_open = 7.0, gap_extension = 3.0;
    int K = 25; // max number of clusters (in the truncated Dirichlet process)
    int nSample = 501; // MCMC iterations

    vector<vector<string>> arg_buffer;
    int count = 0;
    vector<string> arg_buffer_sub;

    // commandline interface
    static struct option longopts[] = {
        {"prefix", required_argument, NULL, 0},
        {"reference", required_argument, NULL, 0},
        {"visualization_level", required_argument, NULL, 0},
        {"sam", required_argument, NULL, 0},
        {"have_true_haplotypes", required_argument, NULL, 0},
        {"true_haplotypes", required_argument, NULL, 0},
        {"do_local_Analysis", required_argument, NULL, 0},
        {"max_reads_in_window", required_argument, NULL, 0},
        {"entropy_threshold", required_argument, NULL, 0},
        {"reconstruction_start", required_argument, NULL, 0},
        {"reconstruction_stop", required_argument, NULL, 0},
        {"min_qual", required_argument, NULL, 0},
        {"min_length", required_argument, NULL, 0},
        {"max_gap_fraction", required_argument, NULL, 0},
        {"min_align_score_fraction", required_argument, NULL, 0},
        {"alpha_MN_local", required_argument, NULL, 0},
        {"min_overlap_factor", required_argument, NULL, 0},
        {"local_window_size_factor", required_argument, NULL, 0},
        {"cluster_number", required_argument, NULL, 0},
        {"nSample", required_argument, NULL, 0},
        {"include_deletions", required_argument, NULL, 0},
        {"reconstructed_haplotypes", required_argument, NULL, 0},
        {"help", no_argument, NULL, 0},
        {NULL, 0, NULL, 0}};

    int ch;
    int longindex = -1;
    while ((ch = getopt_long(argc, argv, "", longopts, &longindex)) != -1) {
      if (ch == '?') {
        throw phaplo::Error(phaplo::ErrorCode::unsupported_flag);
      }

      string choice = longopts[longindex].name;
      if (choice == "prefix") {
        prefix = optarg;
      } else if (choice == "reference") {
        cons = optarg;
      } else if (choice == "visualization_level") {
        visualization_level = atoi(optarg);
      } else if (choice == "sam") {
        if (FASTAreads.has_value()) {
          throw phaplo::Error(phaplo::ErrorCode::multiple_sam_files);
        }
        FASTAreads = optarg;
      } else if (choice == "have_true_haplotypes") {
        have_true_haplotypes = atoi(optarg);
      } else if (choice == "true_haplotypes") {
        FASTAhaplos = optarg;
      } else if (choice == "do_local_Analysis") {
        do_local_Analysis = atoi(optarg);
      } else if (choice == "max_reads_in_window") {
        max_reads_in_window = atoi(optarg);
      } else if (choice == "entropy_threshold") {
        entropy_threshold = atof(optarg);
      } else if (choice == "reconstruction_start") {
        reconstruction_start = atoi(optarg);
      } else if (choice == "reconstruction_stop") {
        reconstruction_stop = atoi(optarg);
      } else if (choice == "min_qual") {
        min_qual = atoi(optarg);
      } else if (choice == "min_length") {
        min_length = atoi(optarg);
      } else if (choice == "max_gap_fraction") {
        max_gap_fraction = atof(optarg);
      } else if (choice == "min_align_score_fraction") {
        min_align_score_fraction = atof(optarg);
      } else if (choice == "alpha_MN_local") {
        alpha_MN_local = atof(optarg);
      } else if (choice == "min_overlap_factor") {
        min_overlap_factor = atof(optarg);
      } else if (choice == "local_window_size_factor") {
        local_window_size_factor = atof(optarg);
      } else if (choice == "cluster_number") {
        K = atoi(optarg);
      } else if (choice == "nSample") {
        nSample = atoi(optarg);
      } else if (choice == "include_deletions") {
        include_deletions = atoi(optarg);
      } else if (choice == "reconstructed_haplotypes") {
        reconstructed_haplos_path = optarg;
      } else if (choice == "help") {
        cout << "Usage: " << argv[0] << " [OPTIONS]\n"
             << "\n"
             << "  This software aims at reconstructing haplotypes from "
                "next-generation sequencing data.\n"
             << "\n"
             << "Options:\n"
             << "  --sam FILE                        Filename of the aligned "
                "reads "
                "(sam format).\n"
             << "  --reference FILE                  Filename of reference "
                "sequence "
                "(FASTA).\n"
             << "  --prefix STR                      Prefix of output files.\n"
             << "  --visualization_level INT         do_visualize (1 = true, 0 "
                "= "
                "false).\n"
             << "  --have_true_haplotypes INT        have_true_haplotypes (1 = "
                "true, 0 "
                "= false).\n"
             << "  --true_haplotypes FILE            Filename of the true "
                "haplotypes "
                "(MSA in "
                "FASTA format) (fill in any dummy filename if there is no "
                "\"true\" haplotypes).\n"
             << "  --do_local_Analysis INT           do_local_analysis (1 = "
                "true, 0 = "
                "false) "
                "(must be 1 in the first run).\n"
             << "  --max_reads_in_window INT         ...\n"
             << "  --entropy_threshold FLOAT         ...\n"
             << "  --reconstruction_start INT        ...\n"
             << "  --reconstruction_stop INT         ...\n"
             << "  --min_qual INT                    ...\n"
             << "  --min_length INT                  ...\n"
             << "  --max_gap_fraction FLOAT          Relative to alignment "
                "length.\n"
             << "  --min_align_score_fraction FLOAT  Relative to read length.\n"
             << "  --alpha_MN_local FLOAT            Prior parameter for "
                "multinomial tables over the nucleotides.\n"
             << "  --min_overlap_factor FLOAT        Reads must have an "
                "overlap with the local reconstruction window of at least this "
                "factor times the window size.\n"
             << "  --local_window_size_factor FLOAT  Size of  local "
                "reconstruction window relative to the median of the read "
                "lengths.\n"
             << "  --cluster_number INT              Max number of clusters "
                "(in the truncated Dirichlet process).\n"
             << "  --nSample INT                     MCMC iterations.\n"
             << "  --include_deletions INT           Include deletions (0 = "
                "no, 1 = yes).\n"
             << "  --reconstructed_haplotypes FILE   Results of global "
                "haplotype reconstruction are copied to FILE (FASTA format)."
             << "  --help                            Show this message and "
                "exit.\n"
             << endl;
        return 0;
      } else {
        throw phaplo::Error(phaplo::ErrorCode::unsupported_flag);
      }
    }

    if (!FASTAreads.has_value()) {
      throw phaplo::Error(phaplo::ErrorCode::no_sam_file);
    }
    if (cons == "") {
      throw phaplo::Error(phaplo::ErrorCode::no_reference_file);
    }

    cout << "Configuration:\n"
         << "  prefix = " << prefix << "\n"
         << "  cons = " << cons << "\n"
         << "  visualization_level = " << visualization_level << "\n"
         << "  FASTAreads = " << FASTAreads.value() << "\n"
         << "  have_true_haplotypes = " << have_true_haplotypes << "\n"
         << "  FASTAhaplos = " << FASTAhaplos << "\n"
         << "  do_local_Analysis = " << do_local_Analysis << "" << endl;

    // let's go!
    include_deletions = false;

    prefix = prefix + prefix_extension;

    // TODO: this can only go wrong at some point
    int mkdir_ret = system(("mkdir -p $(dirname " + prefix + ")").c_str());

    string fas_rm = "rm -f " + prefix + "*.fas";
    string fas_global_rm = "rm -f " + prefix + "global*.fas";
    string lab_rm = "rm -f " + prefix + "*.lab";
    string reads_rm = "rm -f " + prefix + "*.reads";
    string html_rm = "rm -f " + prefix + "*.html";
    string pgm_rm = "rm -f " + prefix + "*.pgm";
    string html_global_rm = "rm -f " + prefix + "global*.html";
    bool do_rm = true; // false;

    if (do_rm) {
      if (do_local_Analysis) {
        int i = system(fas_rm.c_str());
        i = system(lab_rm.c_str());
        i = system(reads_rm.c_str());
        i = system(html_rm.c_str());
        i = system(pgm_rm.c_str());
      } else {
        int i = system(fas_global_rm.c_str());
        i = system(html_global_rm.c_str());
        i = system(pgm_rm.c_str());
      }
    }

    string haplos_fstr;

    cout << fixed << setprecision(3);

    string SQ_string;

    ifstream infRef(cons.c_str(), ios::in);

    getline(infRef, line, '\n');
    while (getline(infRef, line, '\n')) {

      SQ_string = SQ_string + line;
    }
    infRef.close();

    string searchString("+");
    string searchString_2("-");
    string replaceString("N");
    string::size_type pos = 0;
    while ((pos = SQ_string.find(searchString, pos)) != string::npos) {
      SQ_string.replace(pos, searchString.size(), replaceString);
      pos++;
    }

    pos = 0;
    while ((pos = SQ_string.find(searchString_2, pos)) != string::npos) {
      SQ_string.replace(pos, searchString_2.size(), replaceString);
      pos++;
    }

    vector<vector<int>> Reads;
    vector<int> Positions_Start;
    vector<string> IDs;

    bool have_quality_scores = true;

    vector<string> quality_scores;
    vector<int> strand;

    char gap_quality = '*'; //'I';

    double mean_length = 0.;

    int max_sequence_stop = 0;
    int min_seq_start = 100000;

    int min_sequence_stop = 100000;
    int max_seq_start = 0;

    for (auto &&message : phaplo::statistics_to_messages(parseSAMpaired(
             FASTAreads.value(), max_gap_fraction, min_align_score_fraction,
             min_qual, min_length, gap_quality, mean_length, min_seq_start,
             max_seq_start, max_sequence_stop, min_sequence_stop,
             have_quality_scores, quality_scores, strand, Reads,
             Positions_Start, IDs))) {
      std::cerr << "Warning: " << message << std::endl;
    }

    if (Reads.empty()) {
      throw phaplo::Error(phaplo::ErrorCode::no_valid_reads);
    }

    cout << "After parsing the reads in file " << FASTAreads.value()
         << ": average read length= " << mean_length / Reads.size() << ' '
         << Reads.size() << endl;

    if (count > 10) {
      if (reconstruction_start > min_seq_start &&
          reconstruction_start < max_sequence_stop)
        min_seq_start = reconstruction_start;
      if (reconstruction_stop < max_sequence_stop &&
          reconstruction_stop > min_seq_start)
        max_sequence_stop = reconstruction_stop;
    }

    cout << "First read considered in the analysis starts at position "
         << min_seq_start << ". Last read ends at position "
         << max_sequence_stop << endl;
    cout << "There are " << Reads.size() << " reads" << endl;

    vector<int> si(Reads.size());
    for (int i = 0; i < Reads.size(); i++) {
      si[i] = Reads[i].size();
    }

    mean_length = phaplo::median(std::move(si)); // mean_length/Reads.size();
    si.clear();

    local_window_size = int(local_window_size_factor * mean_length);
    int min_overlap = int(min_overlap_factor * local_window_size);

    cout << "Median of read lengths: " << mean_length << endl;
    cout << "Local window size: " << local_window_size << endl;
    cout << "Minimum overlap of reads to local analysis windows: "
         << min_overlap << endl;

    vector<vector<int>> trueHaplos;
    vector<int> trueHaplo_Positions_Start;
    vector<string> trueHaplo_IDs;

    if (have_true_haplotypes) {

      ofstream HAP("myREFERENCE_MSA.txt", ios::out);

      ifstream inf6(FASTAhaplos.c_str(), ios::in);

      string SQ_string = "";
      getline(inf6, line, '\n');
      HAP << line << endl;
      HAP << "dummy line..." << endl;
      while (getline(inf6, line, '\n')) {

        if (line[0] == '>') {

          int count = 0;
          for (int i = 0; i < SQ_string.size(); i++) {

            if (SQ_string[i] != '-') {
              break;
            }
            count++;
          }
          for (int i = 0; i < count; i++) {
            SQ_string[i] = 'N';
          }
          count = SQ_string.size() - 1;
          for (int i = SQ_string.size() - 1; i >= 0; i--) {

            if (SQ_string[i] != '-') {
              break;
            }
            count--;
          }
          for (int i = SQ_string.size() - 1; i >= count; i--) {
            SQ_string[i] = 'N';
          }

          HAP << "1:";
          for (int i = 0; i < SQ_string.size(); i++) {
            if (SQ_string[i] == '-') {
              HAP << "0";
            } else {
              if (SQ_string[i] == 'A' || SQ_string[i] == 'a') {
                HAP << "1";
              } else {
                if (SQ_string[i] == 'C' || SQ_string[i] == 'c') {
                  HAP << "2";
                } else {
                  if (SQ_string[i] == 'G' || SQ_string[i] == 'g') {
                    HAP << "3";
                  } else {
                    if (SQ_string[i] == 'T' || SQ_string[i] == 't') {
                      HAP << "4";
                    } else {
                      HAP << "5";
                    }
                  }
                }
              }
            }
          }
          HAP << endl;
          HAP << line << endl;
          HAP << "dummy line..." << endl;
          SQ_string = "";
        } else {
          SQ_string = SQ_string + line;
        }
      }

      int count = 0;
      for (int i = 0; i < SQ_string.size(); i++) {

        if (SQ_string[i] != '-') {
          break;
        }
        count++;
      }
      for (int i = 0; i < count; i++) {
        SQ_string[i] = 'N';
      }
      count = SQ_string.size() - 1;
      for (int i = SQ_string.size() - 1; i >= 0; i--) {

        if (SQ_string[i] != '-') {
          break;
        }
        count--;
      }
      for (int i = SQ_string.size() - 1; i >= count; i--) {
        SQ_string[i] = 'N';
      }

      HAP << 1 << ":";

      for (int i = 0; i < SQ_string.size(); i++) {
        if (SQ_string[i] == '-') {
          HAP << "0";
        } else {
          if (SQ_string[i] == 'A' || SQ_string[i] == 'a') {
            HAP << "1";
          } else {
            if (SQ_string[i] == 'C' || SQ_string[i] == 'c') {
              HAP << "2";
            } else {
              if (SQ_string[i] == 'G' || SQ_string[i] == 'g') {
                HAP << "3";
              } else {
                if (SQ_string[i] == 'T' || SQ_string[i] == 't') {
                  HAP << "4";
                } else {
                  HAP << "5";
                }
              }
            }
          }
        }
      }
      HAP << endl;
      HAP.close();
      inf6.close();
      haplos_fstr = "myREFERENCE_MSA.txt";

      /// haplos

      string line, line_ID, line_stats;
      vector<string> tokens;
      string tok = ":";

      ifstream inf4(haplos_fstr.c_str(), ios::in);

      int haplo_max_start = -1;
      int haplo_min_stop = 100000;
      while (getline(inf4, line_ID, '\n')) {

        getline(inf4, line_stats, '\n');

        getline(inf4, line, '\n');

        trueHaplo_IDs.push_back(line_ID);

        tokens = tokenize(line, tok);
        trueHaplo_Positions_Start.push_back(atoi(tokens[0].c_str()));
        if (atoi(tokens[0].c_str()) > haplo_max_start) {
          haplo_max_start = atoi(tokens[0].c_str());
        }

        trueHaplos.push_back(phaplo::to_vector_of_ints(tokens[1]));

        if (tokens[1].size() < haplo_min_stop)
          haplo_min_stop = atoi(tokens[0].c_str()) + tokens[1].size();
      }

      inf4.close();

      if (haplo_min_stop < max_sequence_stop)
        max_sequence_stop = haplo_min_stop;

      if (haplo_max_start > min_seq_start) {
        min_seq_start = haplo_max_start;
      }
    }

    cout << "Reconstruction starts at position " << min_seq_start
         << " and stops at position " << max_sequence_stop << endl;

    int nT = Reads.size();

    int WindowStart = min_seq_start;

    vector<vector<int>> WindowStartStop;

    vector<int> foundClusters;

    // entopy selection

    levels = 5;
    int GD = max_sequence_stop;
    Matrix<double> tot_counts(GD, levels);

    tot_counts = 1;

    for (int o = 0; o < nT; o++) {
      for (int r = 0; r < Reads[o].size(); r++) {
        int nuc = Reads[o][r];

        int qual = int(quality_scores[o][r]) - qual_subtract;

        // cout << qual << endl;

        int idx = r + Positions_Start[o] - 1;

        if (idx < GD && nuc < 5) {

          tot_counts(idx, nuc) += qual;
        }
      }
    }

    vector<int> is_gap(GD, 0);
    for (int i = WindowStart - 1; i < GD; i++) {

      if (tot_counts(i, 0) > 0.05 * sum(tot_counts(i, _)))
        is_gap[i] = 1;

      if (!include_deletions)
        tot_counts(i, 0) = 0;

      tot_counts(i, _) = tot_counts(i, _) + 30;
      tot_counts(i, _) = tot_counts(i, _) / sum(tot_counts(i, _));
    }

    vector<double> entropy(GD, 0);
    for (int i = WindowStart - 1; i < GD; i++) {
      double sum = 0.;
      for (int j = 0; j < levels; j++) {
        sum -= tot_counts(i, j) * log(tot_counts(i, j));
      }
      entropy[i] = sum;
      if (!include_deletions && is_gap[i]) {
        entropy[i] = 0;
      }
    }

    // cout << "yep" << endl;

    vector<int> e_index(entropy.size());
    for (int i = 0; i < entropy.size(); i++)
      e_index[i] = i;

    sort(e_index.begin(), e_index.end(), index_cmp<vector<double> &>(entropy));

    double eThresh = -(1 - entropy_threshold) * log(1 - entropy_threshold) -
                     entropy_threshold * log(entropy_threshold);

    vector<int> entropy_select(GD, 0);
    for (int i = 0; i < int(GD * entropy_fraction); i++) {
      if (entropy[e_index[i]] > eThresh) {

        entropy_select[e_index[i]] = 1;
      }
    }

    entropy.clear();
    e_index.clear();
    is_gap.clear();
    tot_counts.resize(0, 0);

    ////////////////

    ////// quality scoring
    string qstr = prefix + "overlap.txt";
    ofstream ofQ(qstr.c_str());
    int overlap = 201;
    if (GD - 2 * overlap - WindowStart < 10) {
      overlap = 0.5 * (GD - WindowStart - 10);
    }

    Matrix<> over_counts(10, GD - 2 * overlap - WindowStart);

    over_counts = 0;

    for (int p = WindowStart + overlap; p < GD - overlap; p++) {

      if (entropy_select[p] == 0)
        continue;

      vector<int> over_vec;
      for (int o = 0; o < nT; o++) {

        int idx_start = Positions_Start[o];
        int idx_stop = Reads[o].size() - 1 + Positions_Start[o];

        if (idx_start <= p && idx_stop >= p) {

          int over = p - idx_start;
          if (idx_stop - p < over)
            over = idx_stop - p;

          over_vec.push_back(over);
        }
      }
      sort(over_vec.begin(), over_vec.end(), myfunction);
      for (int i = 0; i < over_counts.rows(); i++)
        if (i < over_vec.size())
          over_counts(i, p - overlap - WindowStart) = over_vec[i];

      ofQ << p;
      // cout << p ;
      int i;
      for (i = 0; i < over_counts.rows(); i++) {
        if (over_counts(i, p - overlap - WindowStart) < overlap - 1) {
          ofQ << '\t' << i;
          // cout <<'\t'<< i;
          break;
        }
      }
      if (i == over_counts.rows()) {
        ofQ << '\t' << over_counts.rows();
        // cout <<'\t'<<">= "<<  over_counts.rows();
      }
      for (i = 0; i < over_counts.rows(); i++) {
        if (over_counts(i, p - overlap - WindowStart) < 0.75 * (overlap - 1)) {
          ofQ << '\t' << i;
          // cout <<'\t'<< i;
          break;
        }
      }
      if (i == over_counts.rows()) {
        ofQ << '\t' << over_counts.rows();
        //  cout <<'\t'<<">= "<<  over_counts.rows();
      }

      ofQ << endl;
      // cout << endl;
    }
    ofQ.close();
    // cout << endl;

    ///////////////////////////////

    if (do_local_Analysis) {

      local_Analysis(
          min_overlap, K, nSample, nT, max_reads_in_window,
          WindowStart + floor(0.5 * local_window_size), local_window_size,
          WindowStartStop, max_sequence_stop, min_seq_start, Positions_Start,
          Reads, have_quality_scores, quality_scores, IDs, strand, prefix,
          SQ_string, have_true_haplotypes, trueHaplos,
          trueHaplo_Positions_Start, trueHaplo_IDs, foundClusters,
          entropy_threshold, entropy_fraction, entropy_select, alpha_MN_local);

      string startStopStr = prefix + "StartStop.txt";

      ofstream outS(startStopStr.c_str(), ios::out);
      for (int i = 0; i < WindowStartStop.size(); i++) {
        outS << i << ' ' << WindowStartStop[i][0] << ' '
             << WindowStartStop[i][1] << ' ' << foundClusters[i] << endl;
      }
      outS.close();
    }

    string startStopStr = prefix + "StartStop.txt";
    int max_clusters = -1;
    int max_index;
    WindowStartStop.clear();
    ifstream inS(startStopStr.c_str(), ios::in);
    vector<int> clusterNums;
    while (getline(inS, line, '\n')) {
      tokens = tokenize(line, " ");
      vector<int> wss(2);
      wss[0] = atoi(tokens[1].c_str());
      wss[1] = atoi(tokens[2].c_str());
      WindowStartStop.push_back(wss);

      int clnum = atoi(tokens[3].c_str());
      clusterNums.push_back(clnum);

      if (clnum > max_clusters) {
        max_clusters = clnum;
        max_index = atoi(tokens[0].c_str());
      }
    }
    inS.close();

    if (max_clusters < 0) {
      throw phaplo::Error(phaplo::ErrorCode::no_maximum_cluster_number_window);
    }

    double min_d = 1e100;
    int min_i;
    for (int i = 0; i < clusterNums.size(); i++) {
      if (clusterNums[i] == max_clusters) {
        double dist = fabs(i - 0.5 * clusterNums.size());
        if (dist < min_d) {
          min_d = dist;
          min_i = i;
        }
      }
    }

    max_index = min_i;

    int select_start = max_index;
    reconstruct_global(
        min_overlap, K, nSample, WindowStartStop, local_window_size,
        select_start, nT, max_sequence_stop + 1, min_seq_start, Positions_Start,
        Reads, have_quality_scores, quality_scores, IDs, strand, prefix,
        SQ_string, have_true_haplotypes, trueHaplos, trueHaplo_Positions_Start,
        trueHaplo_IDs, entropy_threshold, entropy_fraction, entropy_select,
        reconstructed_haplos_path);
  } catch (const phaplo::Error &error) {
    std::cerr << "Error: " << error.what() << std::endl;
    return error.id();
  }
}
