// Copyright 2014 Hamed S. Najafabadi

/********************************************************************

This file is part of RCOpt.

RCOpt is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RCOpt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RCOpt.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/

#ifndef _H_DECLARATIONS // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_DECLARATIONS // this macro indicates that this file is now included


#include "structures.h" // this header file is required for definition of structures


double normal_cum_p( double z );
double cumhyper(int i, int s1, int s2, int N);
bool read_arguments( int argc, char* argv[] );
void print_arguments();
void print_commandline_format();

int find_motif( const char *name, s_motif *motifs[], int num_motifs );
char* extract_phrase( const char *line, int phrase_index, char *phrase, const char *delimiters );

bool open_output( ofstream &ofs, const char *path, const char *extension );
bool open_rndforest( const char *FASTA_file, s_motif *motifs[], int *num_motifs, char const *filter );
bool open_weights( const char *weight_file, s_motif *motifs[], int num_motifs );
bool open_FASTA( const char *file, s_seq *seqs[], int *num_seqs );

bool read_rndforest( ifstream &ifs, s_motif *motifs[], int *num_motifs, char const *filter, char this_EOL );
bool read_weights( ifstream &ifs, s_motif *motifs[], int num_motifs, char this_EOL );
bool read_FASTA( ifstream &ifs, s_seq *seqs[], int *num_seqs, char this_EOL );

bool generate_PWMs( s_motif *motifs[], int num_motifs, double enrichment_fold );
void generate_metaPFM( s_motif *motifs[], int num_motifs, s_motif *metaPFM, int metaPFM_width );

void scan_sequences( s_seq *seqs[], int num_seqs, double *PFM[], int PFM_width, int direction );
int calculate_initial_enrichments( s_seq *seqs[], int num_seqs, s_motif *motifs[], int num_motifs );
int calculate_optimized_enrichments( s_seq *seqs[], int num_seqs, s_motif *motifs[], int num_motifs );

double get_correlation( double *PWM1[], double *PWM2[], int PWM_width, double *mean, double *stdev );

int optimize_PFMs( s_seq *seqs[], int num_seqs, s_motif *motifs[], int num_motifs, int num_iterations, int reoptimization );

bool write_scores( ofstream &ofs, s_seq *seqs[], int num_seqs, s_motif *motifs[], int num_motifs, const char *experiment_name );
bool write_metaPFM_scores( ofstream &ofs, s_seq *seqs[], int num_seqs, int metaPFM_width );
bool write_PFMs( ofstream &ofs, s_motif *motifs[], int num_motifs, const char *experiment_name, bool opt_only );
void write_opt_MEME( ofstream &ofs, s_motif *motif, const char *experiment_name );
void write_report( ofstream &ofs, s_motif *motifs[], int num_motifs, const char *experiment_name );
void write_graphics( ofstream &ofs, s_motif *motifs[], int num_motifs );


//*******************************************************************
#endif // this is to make sure that this file will not be included twice
