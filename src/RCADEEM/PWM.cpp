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

int _l = 4; // the length of the recognition sequence for each ZF

//////////////////////////////////////////////////////////
double _get_max_range(
	double *PWM[],
	int PWM_width )
// determines the maximum difference between the max and min values of any column in the PWM
{
	double max_range = 0;	
	// determine the min and max scores in each position, based on which the overall min/max will be deteremined
	int n, x;
	for( x = 0; x < PWM_width; x ++ )
	{
		// get the min and max in position x
		double this_min, this_max;	
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			if( n == 0 ||
				this_min > PWM[ n ][ x ] )
				this_min = PWM[ n ][ x ];
			if( n == 0 ||
				this_max < PWM[ n ][ x ] )
				this_max = PWM[ n ][ x ];
		}
		
		double this_range = this_max - this_min;
		if( x == 0 ||
			max_range < this_range )
			max_range = this_range;
	}
	
	return max_range;
}

//////////////////////////////////////////////////////////
void _convert_PWM_to_likelihood(
	double *PWM[],
	int PWM_width,
	double enrichment_fold )
{	
	double c = _get_max_range( PWM, PWM_width ) / log( enrichment_fold );
	c += 0.00000001; // this pseudocount ensures that there will be no DIV/0
	
	//////////////
	//c = 0.746986509;
	//////////////
	
	// convert to likelihoods
	int x, n;
	for( x = 0; x < PWM_width; x ++ )
	{
		// convert this column values
		
		// pass 1: convert to frequencies
		// initialize the frequencies
		double sum_exp = 0;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			PWM[ n ][ x ] = exp( PWM[ n ][ x ] / c );
			sum_exp += PWM[ n ][ x ];
		}
		// pass 2: convert to log-likelihoods
		// normalize the frequencies, and then covert to log-likelihood
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			PWM[ n ][ x ] /= sum_exp; // normalize the frequencies
			PWM[ n ][ x ] = log10( PWM[ n ][ x ] * NUM_N_LETTERS ); // convert to log-likelihood
		}
	}
}


//////////////////////////////////////////////////////////
void _get_PWM(
	s_motif *motif,
	double *PWM[],
	int *PWM_width,
	double enrichment_fold )
{
	// determine the PWM width based on the number of zinc fingers and the bait length in the model
	*PWM_width = motif ->num_zfs * 3 + _l - 3;
	// set the PWM to zero
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
	{
		// reset all the numbers of this row to zero
		memset( PWM[ n ], 0, sizeof(double) * (*PWM_width) );
	}

	int offset = *PWM_width - _l; // the offset of the zinc finger PWM relative to the overall PWM
	
	// concatenate the PWMs of the ZFs
	int i, x;
	for( i = 0; i < motif ->num_zfs; i ++ )
	{
		for( n = 0; n < NUM_N_LETTERS; n ++ )	// examine each nucleotide
			for( x = 0; x < _l; x ++ ) // examine each bait position
			//for( x = 0; x < 3; x ++ ) // examine only the first three positions
				PWM[ n ][ x + offset ] +=
					motif ->zfs[ i ] ->weight * motif ->zfs[ i ] ->PWM[ n ][ x ];
				
		offset -= 3; // update the offset
	}
	
	// the PWM should be converted to likelihood-based PWM
	if( enrichment_fold > 0 )
		_convert_PWM_to_likelihood( PWM, *PWM_width, enrichment_fold );
}


///////////////////////////////////////////////////////////////////////////////////////////
bool generate_PWMs(
	s_motif *motifs[],
	int num_motifs,
	double enrichment_fold )
{
	// find the maximum number of zinc fingers per motif
	int max_num_zfs = 0;
	int i;
	for( i = 0; i < num_motifs; i ++ )
		if( max_num_zfs < motifs[ i ] ->num_zfs )
			max_num_zfs = motifs[ i ] ->num_zfs;
			
	// initialize the PWM
	double *PWM[ NUM_N_LETTERS ];
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		PWM[ n ] = new double[ max_num_zfs * 3 + _l - 3 ];
		
	// store the PWM for each motif
	for( i = 0; i < num_motifs; i ++ )
	{
		// get the PWM for this motif
		int PWM_width;
		_get_PWM( motifs[ i ], PWM, &PWM_width, enrichment_fold );
		
		// store the PWM and the PFM for this motif
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			motifs[ i ] ->PWM[ n ] = new double[ PWM_width ];
			motifs[ i ] ->PFM[ n ] = new double[ PWM_width ];

			int j;
			for( j = 0; j < PWM_width; j ++ )
			{
				motifs[ i ] ->PWM[ n ][ j ] = PWM[ n ][ j ];
				motifs[ i ] ->PFM[ n ][ j ] = pow( 10, PWM[ n ][ j ] ) / 4.0;
			}
		}
		
		motifs[ i ] ->PWM_width = PWM_width - 1;
		
		// Store the correspondence between each position in the PWM and the position
		// within the "meta" PWM, i.e. the PWM that would be obtained by concatenating of
		// all the motifs of all the ZFs of the protein.
		// Note that the meta-PWM has a reversed order, for simplicity
		
		motifs[ i ] ->metaPFM_pos = new int[ motifs[ i ] ->PWM_width ];
		int j;
		for( j = 0; j < motifs[ i ] ->PWM_width; j ++ )
			motifs[ i ] ->metaPFM_pos[ j ] = ( motifs[ i ] ->first_zf_in_protein - 1 ) * 3 +
				motifs[ i ] ->PWM_width - 1 - j;
		
	}
	
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////
void generate_metaPFM(
	s_motif *motifs[],
	int num_motifs,
	s_motif *metaPFM,
	int metaPFM_width )
// generate a meta-PFM that represents weighted average of all individual PFMs, with the weights corresponding to estimated PFM frequencies
{
	// initiate the meta-PFM object
	_COPY_STR( metaPFM ->protein_name, "meta-PFM" );
	_COPY_STR( metaPFM ->name, "meta-PFM" );
	metaPFM ->PFM_optimized = 1;
	metaPFM ->PWM_width = metaPFM_width;
	
	// initialize the arrays
	int n, x, metax;
	for( n = 0; n < 4; n ++ )
	{
		metaPFM ->opt_PFM[ n ] = new double[ metaPFM_width ];
		memset( metaPFM ->opt_PFM[ n ], 0, sizeof( double ) * metaPFM_width );

		metaPFM ->PFM[ n ] = new double[ metaPFM_width ];
		memset( metaPFM ->PFM[ n ], 0, sizeof( double ) * metaPFM_width );
	}
	
	// calculate the average PFM (the meta-PFM)
	int i, j;
	double total_freq = 0;
	for( i = 0; i < num_motifs; i ++ )
	{
		total_freq += motifs[ i ] ->freq;
		
		for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
			for( n = 0; n < 4; n ++ )
			{
				metax = metaPFM_width - 1 - motifs[ i ] ->metaPFM_pos[ x ]; // the position within the meta-PFM

				metaPFM ->opt_PFM[ n ][ metax ] += motifs[ i ] ->freq * motifs[ i ] ->opt_PFM[ n ][ x ];
				metaPFM ->PFM[ n ][ metax ] += motifs[ i ] ->freq * motifs[ i ] ->PFM[ n ][ x ];
			}
		
		// pad the flanking regions in the meta-PFM with N's
		for( metax = 0; metax < metaPFM_width - 1 - motifs[ i ] ->metaPFM_pos[ 0 ]; metax ++ )
			for( n = 0; n < 4; n ++ )
			{
				metaPFM ->opt_PFM[ n ][ metax ] += motifs[ i ] ->freq * 0.25;
				metaPFM ->PFM[ n ][ metax ] += motifs[ i ] ->freq * 0.25;
			}

		for( metax = metaPFM_width - motifs[ i ] ->metaPFM_pos[ motifs[ i ] ->PWM_width - 1 ]; metax < metaPFM_width; metax ++ )
			for( n = 0; n < 4; n ++ )
		{
				metaPFM ->opt_PFM[ n ][ metax ] += motifs[ i ] ->freq * 0.25;
				metaPFM ->PFM[ n ][ metax ] += motifs[ i ] ->freq * 0.25;
		}
	}
	
	// take the average
	for( metax = 0; metax < metaPFM_width; metax ++ )
		for( n = 0; n < 4; n ++ )
		{
			metaPFM ->opt_PFM[ n ][ metax ] /= total_freq;
			metaPFM ->PFM[ n ][ metax ] /= total_freq;
		}
}
