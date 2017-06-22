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


// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

extern const char *__rf_file;
extern const char *__filter;
extern const char *__weight_file;
extern const char *__fasta_file;
extern const double *__enrichment_fold;
extern const char *__experiment;
extern const char *__output_file;


void welcome_message()
{
	cout << endl
		<< "****************************** RCOpt version 1.01 *****************************" << endl
		<< "*.............................................................................*" << endl
		<< "* Copyright 2014 Hamed S. Najafabadi .........................................*" << endl
		<< "*.............................................................................*" << endl
		<< "*.............................................................................*" << endl
		<< endl;
}


int main( int argc, char* argv[] )
{	
	welcome_message();
	
	srand ( time(NULL) );

	if( argc <= 1 )
	// no argument is profided
	{
		print_commandline_format();
		return 0;
	}
	
	if( !read_arguments( argc, argv ) )
		return 1;
		
	print_arguments();

	//******************* open output

	ofstream ofs_PFM;
	if( !open_output( ofs_PFM, __output_file, ".PFM.txt" ) )
		return 1;

	ofstream ofs_opt_PFM;
	if( !open_output( ofs_opt_PFM, __output_file, ".opt.PFM.txt" ) )
		return 1;

	ofstream ofs_ps;
	if( !open_output( ofs_ps, __output_file, ".ps" ) )
		return 1;

	ofstream ofs_report;
	if( !open_output( ofs_report, __output_file, ".report.txt" ) )
		return 1;

	ofstream ofs_scores;
	if( !open_output( ofs_scores, __output_file, ".PFM.scores.txt" ) )
		return 1;

	ofstream ofs_metaPFM_scores;
	if( !open_output( ofs_metaPFM_scores, __output_file, ".metaPFM.profiles.txt" ) )
		return 1;

	//******************* open the result file

	s_motif *motifs[ MAX_MOTIFS ];
	int num_motifs = 0;

	cout << "Opening the randomForest output file..." << endl;
	if( !open_rndforest( __rf_file, motifs, &num_motifs, __filter ) )
		return 1;
		
	//******************* open the weight file

	if( strcmp( __weight_file, "None" ) != 0 ) // a weight file has been provided
	{
		cout << "Opening the weight file..." << endl;
		if( !open_weights( __weight_file, motifs, num_motifs ) )
			return 1;
	}
	
	//******************* open the FASTA file

	s_seq *seqs[ MAX_SEQS ];
	int num_seqs = 0;
	cout << "Opening the FASTA file..." << endl;
	if( !open_FASTA( __fasta_file, seqs, &num_seqs ) )
		return 1;
		
	//******************* calculations
	cout << "Generating PWMs..." << endl;
	if( !generate_PWMs( motifs, num_motifs, *__enrichment_fold ) )
		return 1;
		
	cout << "Finding enriched PFMs..." << endl;
	int num_enriched = calculate_initial_enrichments( seqs, num_seqs, motifs, num_motifs );
	if( !num_enriched )
	{
		cout << "ERROR: No motifs are enriched." << endl;
		return 1;
	}
	
	int metaPFM_width;
	int i;
	for( i = 0; ; i ++ )
	{
		cout << "(Re-)optimizing " << num_enriched << " PFMs with significant enrichment at FDR<0.01 ..." << endl;
		metaPFM_width = optimize_PFMs( seqs, num_seqs, motifs, num_enriched, 30, i );

		// calculate optimized enrichments
		int new_num_enriched = calculate_optimized_enrichments( seqs, num_seqs, motifs, num_enriched );
		
		if( new_num_enriched == num_enriched )
			break;
			
		if( !new_num_enriched )
		{
			cout << "ERROR: No motifs are enriched." << endl;
			return 1;
		}
		
		num_enriched = new_num_enriched;
	}
		
	//******************* write the output files
	cout << "Writing to output..." << endl;
	
	// write scores
	write_scores( ofs_scores, seqs, num_seqs, motifs, num_enriched, __experiment );
	write_metaPFM_scores( ofs_metaPFM_scores, seqs, num_seqs, metaPFM_width );
	
	motifs[ num_enriched ] = new s_motif;
	generate_metaPFM( motifs, num_enriched, motifs[ num_enriched ], metaPFM_width );

	// write all PFMs
	if( !write_PFMs( ofs_PFM, motifs, num_enriched+1, __experiment, false ) )
		return 1;

	// write the optimized PFMs
	if( !write_PFMs( ofs_opt_PFM, motifs, num_enriched+1, __experiment, true ) )
		return 1;

	// write the postscript
	write_graphics( ofs_ps, motifs, num_enriched+1 );

	// write report
	write_report( ofs_report, motifs, num_enriched+1, __experiment );

	cout << endl << "Job finished successfully." << endl;

	return 0;
}
