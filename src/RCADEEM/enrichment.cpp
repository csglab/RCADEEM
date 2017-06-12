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


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

////////////////////////////////////////////////////////////////////////////////
int _compare_p_hyper( const void *a, const void *b )
{
	s_motif *motif1 = *( (s_motif **) a );
	s_motif *motif2 = *( (s_motif **) b );

	if( motif1 ->p_hyper > motif2 ->p_hyper )
		return 1;
	else if( motif1 ->p_hyper < motif2 ->p_hyper )
		return -1;
	
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
int _compare_p_correl( const void *a, const void *b )
{
	s_motif *motif1 = *( (s_motif **) a );
	s_motif *motif2 = *( (s_motif **) b );

	if( motif1 ->p_correl > motif2 ->p_correl )
		return 1;
	else if( motif1 ->p_correl < motif2 ->p_correl )
		return -1;
	
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
int calculate_initial_enrichments(
	s_seq *seqs[],
	int num_seqs,
	s_motif *motifs[],
	int num_motifs )
{
	// count the total number of positive sequences
	int K = 0;
	int i, j;
	for( i = 0; i < num_seqs; i ++ )
		if( seqs[ i ] ->positive )
			K ++;

	cout << "Motif" << char(9) << "k" << char(9) << "n" << char(9) << "K" << char(9) << "N" << char(9) << "P-value" << endl;
	for( i = 0; i < num_motifs; i ++ )
	{
		scan_sequences( seqs, num_seqs, motifs[ i ] ->PFM, motifs[ i ] ->PWM_width, BOTH_STRANDS );

		// count the number of positive sequences among the top 100 that have the highest score for this PFM
		int k = 0;
		for( j = 0; j < 100; j ++ )
			if( seqs[ j ] ->positive )
				k ++;
		
		motifs[ i ] ->p_hyper = cumhyper( k, 100, K, num_seqs );
		cout << motifs[ i ] ->name << char(9) << k << char(9) << 100 << char(9) << K
			<< char(9) << num_seqs << char(9) << motifs[ i ] ->p_hyper << endl;
	}

	// sort the motifs by ascending order of hypergeometric p-value
	qsort( motifs, num_motifs, sizeof(s_motif*), _compare_p_hyper );

	// mark the motifs that pass FDR <= 0.01
	
	for( i = num_motifs - 1; i >= 0; i -- )
		if( motifs[ i ] ->p_hyper * num_motifs / (i+1.0) <= 0.01 ) // this is the motif that passes the cutoff
		{
			for( j = 0; j <= i; j ++ )
				motifs[ j ] ->PFM_optimized = 1;
				
			break;
		}
	
	// return the number of motifs that were marked for optimization	
	return i+1;
}


///////////////////////////////////////////////////////////////////////////////////////////
int calculate_optimized_enrichments(
	s_seq *seqs[],
	int num_seqs,
	s_motif *motifs[],
	int num_motifs )
{
	int initial_num_motifs = num_motifs;
	s_motif *motif_order[ MAX_MOTIFS ];
	memcpy( motif_order, motifs, sizeof(s_motif*) * num_motifs );
	
	// count the total number of positive sequences
	int K = 0;
	int i, j;
	for( i = 0; i < num_seqs; i ++ )
		if( seqs[ i ] ->positive )
			K ++;

	///////////////////////////////////// find motifs that are still similar to their seed		
	static double *PWM1[ 4 ];
	static double *PWM2[ 4 ];
	static int initialized = 0;
	
	if( !initialized )
	{
		for( i = 0; i < 4; i ++ )
		{
			PWM1[ i ] = new double[ MAX_ZF_PER_GENE*3 + 1 ];
			PWM2[ i ] = new double[ MAX_ZF_PER_GENE*3 + 1 ];
		}
		
		initialized = 1;
	}
			

	for( i = 0; i < num_motifs; i ++ )
		if( motifs[ i ] ->PFM_optimized )
		{
			int n, x;
			for( n = 0; n < NUM_N_LETTERS; n ++ )
				for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
				{
					PWM1[ n ][ x ] = log10( motifs[ i ] ->PFM[ n ][ x ] * 4 );
					PWM2[ n ][ x ] = log10( motifs[ i ] ->opt_PFM[ n ][ x ] * 4 );
				}
			
			double mean_correl = 0, stdev_correl = 1;
			double correl = get_correlation( PWM1, PWM2, motifs[ i ] ->PWM_width, &mean_correl, &stdev_correl );
			double z_correl = ( correl - mean_correl ) / stdev_correl;
			double p_correl = normal_cum_p( z_correl );
			
			motifs[ i ] ->correl = correl;
			motifs[ i ] ->p_correl = p_correl;
		}


	// sort the motifs by ascending order of correlation p-value
	qsort( motifs, num_motifs, sizeof(s_motif*), _compare_p_correl );

	// mark the motifs that pass FDR <= 0.01
	
	for( i = num_motifs - 1; i >= 0; i -- )
		if( motifs[ i ] ->p_correl * num_motifs / (i+1.0) <= 0.01 ) // this is the motif that passes the cutoff
			break;

	num_motifs = i+1; // i+1 motifs pass this filter

	///////////////////////////////////// find motifs that are enriched among positive sequences
	cout << "Motif" << char(9) << "k" << char(9) << "n" << char(9) << "K" << char(9) << "N" << char(9) << "P-value" << endl;
	for( i = 0; i < num_motifs; i ++ )
	{
		motifs[ j ] ->PFM_optimized = 0; // reset the PFM_optimized status to zero
		scan_sequences( seqs, num_seqs, motifs[ i ] ->opt_PFM, motifs[ i ] ->PWM_width, BOTH_STRANDS );

		// count the number of positive sequences among the top 100 that have the highest score for this PFM
		int k = 0;
		for( j = 0; j < 100; j ++ )
			if( seqs[ j ] ->positive )
				k ++;
		
		motifs[ i ] ->p_hyper = cumhyper( k, 100, K, num_seqs );
		cout << motifs[ i ] ->name << char(9) << k << char(9) << 100 << char(9) << K
			<< char(9) << num_seqs << char(9) << motifs[ i ] ->p_hyper << endl;
	}

	// sort the motifs by ascending order of hypergeometric p-value
	qsort( motifs, num_motifs, sizeof(s_motif*), _compare_p_hyper );

	// mark the motifs that pass FDR <= 0.01
	
	for( i = num_motifs - 1; i >= 0; i -- )
		if( motifs[ i ] ->p_hyper * num_motifs / (i+1.0) <= 0.01 ) // this is the motif that passes the cutoff
		{
			motifs[ i ] ->PFM_optimized = 1;
			break;
		}

	num_motifs = i+1; // i+1 motifs pass this filter
	
	if( num_motifs == initial_num_motifs ) // this was the last round of optimization; so restore the original order of the motifs to ensure compatibility with hmm_scores
		memcpy( motifs, motif_order, sizeof(s_motif*) * num_motifs );

	
	// return the number of motifs that were marked for optimization	
	return num_motifs;
}



/*int _compare_affinities( const void *a, const void *b );
///////////////////////////////////////////////////////////////////////////////////////////
int calculate_optimized_enrichments(
	s_seq *seqs[],
	int num_seqs,
	s_motif *motifs[],
	int num_motifs )
{
	// count the total number of positive sequences
	int K = 0;
	int i, j;
	for( i = 0; i < num_seqs; i ++ )
		if( seqs[ i ] ->positive )
			K ++;

	cout << "Motif" << char(9) << "k" << char(9) << "n" << char(9) << "K" << char(9) << "N" << char(9) << "P-value" << endl;
	for( i = 0; i < num_motifs; i ++ )
	{
		// copy the HMM score of this motif as the active sequence score
		int j;
		for( j = 0; j < num_seqs; j ++ )
			seqs[ j ] ->score = seqs[ j ] ->hmm_scores[ i ];
		// sort the sequences by descending order of the active scores
		qsort( seqs, num_seqs, sizeof(s_seq*), _compare_affinities );

		// count the number of positive sequences among the top 100 that have the highest score for this motif
		int k = 0;
		for( j = 0; j < 100; j ++ )
			if( seqs[ j ] ->positive )
				k ++;
		
		motifs[ i ] ->p_hyper = cumhyper( k, 100, K, num_seqs );
		cout << motifs[ i ] ->name << char(9) << k << char(9) << 100 << char(9) << K
			<< char(9) << num_seqs << char(9) << motifs[ i ] ->p_hyper << endl;
	}

	// store the order of the motifs
	s_motif *motif_order[ MAX_MOTIFS ];
	memcpy( motif_order, motifs, sizeof( s_motif* ) * num_motifs );
	// sort the motifs by ascending order of hypergeometric p-value
	qsort( motifs, num_motifs, sizeof(s_motif*), _compare_p_hyper );

	// mark the motifs that pass FDR <= 0.01
	
	for( i = num_motifs - 1; i >= 0; i -- )
		if( motifs[ i ] ->p_hyper * num_motifs / (i+1.0) <= 0.01 ) // this is the motif that passes the cutoff
		{
			for( j = 0; j <= i; j ++ )
				motifs[ j ] ->PFM_optimized = 1;
				
			break;
		}

	if( i+1 == num_motifs ) // this was the last round of optimization; so restore the original order of the motifs to ensure compatibility with hmm_scores
		memcpy( motifs, motif_order, sizeof( s_motif* ) * num_motifs );
	
	// return the number of motifs that were marked for optimization	
	return i+1;
}
*/