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

///////////////////////////////////////////////////////////////////////////////////////////
int calculate_enrichments(
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
