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

// This script uses numerically stable HMM algorithm described here:
// http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

#define INIFREQ					(0.01/num_motifs/2)
#define INIPRIOR				(0.01/(num_states-1))
#define INIBKG					(0.99)


#define MAX_ITER			10
#define MAX_UNCHANGED		20
#define CONVERGENCE_RATE	0.8

//////////////////////////////////////////////////////////
int _n_index( char nucleotide )
// returns the index of the specified nucleotide
{
	int i;
	for( i = 0; i < 4; i ++ )
		if( nucleotide == __n_letters[ i ] )
			return i;
			
	return -1;
}

//////////////////////////////////////////////////////////
void optimize_PFMs(
	s_seq *seqs[],
	int num_seqs,
	s_motif *motifs[],
	int num_motifs )
{
	// This function uses an HMM to optimize all potential motifs simultaneously.
	// The HMM has the following form:
	// start --> bkg --> motifs
	// the bkg can also transition to itself
	// Each motif is presented using two sets of states, representing the forward and reverse complement matches
	// As much as possible, I will use the same notations as in
	// http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
	// Note that we represent all probabilities as ln(p). This notation is dropped from variables for simplicity


	//////////////////////////////////////////////////////// initialize the transition matrix
	// index 0: bkg
	// index 1: the start of the first motif
	// index 1+w-1: the end of the first motif (w: motif width)
	// index 1+w: the start of the next motif
	// ...

	// First, calculate the number of states:  bkg + 2 * sum( motif_width )
	// Also store the location of motif starts and ends
	int motif_f_start[ MAX_MOTIFS ];
	int motif_f_end[ MAX_MOTIFS ];
	int motif_r_start[ MAX_MOTIFS ];
	int motif_r_end[ MAX_MOTIFS ];
	

	int num_states = 1; //  bkg
	int i, j;
	for( i = 0; i < num_motifs; i ++ )
	{
		motif_f_start[ i ] = num_states;
		motif_f_end[ i ] = num_states + motifs[ i ] ->PWM_width - 1;
		num_states += motifs[ i ] ->PWM_width;

		motif_r_start[ i ] = num_states;
		motif_r_end[ i ] = num_states + motifs[ i ] ->PWM_width - 1;
		num_states += motifs[ i ] ->PWM_width;
	}
	

	// now initialize the matrix of transition probabilities
	double *aij[ MAX_STATES ]; // the ln of transition probabilities; aij[ i ][ j ] : transition from state i to state j
	for( i = 0; i < num_states; i ++ )
	{
		aij[ i ] = new double[ num_states ];
		
		// set the initial values of the transition matrix to zero
		for( j = 0; j < num_states; j ++ )
			aij[ i ][ j ] = 0;
	}

	
	///// initialize the aij values, as well as the prior pointers
	// the prior probability of state i is always equal to the probability of transition from bkg to the motif that state i represents

	double *priors = new double[ 1 + num_motifs ]; // this array will hold the actual priors; index 0: bkg; index 1: motif1; index 2: motif2; etc.
	double *prior[ MAX_STATES ]; // this array will point to the address of the value
	aij[ 0 ][ 0 ] = INIBKG; // the probability of going from bkg to bkg

	prior[ 0 ] = priors + 0; // the address to bkg prior
	*prior[ 0 ] = INIBKG; // the initial value of this prior
	
	for( i = 0; i < num_motifs; i ++ )
	{
		int start = motif_f_start[ i ];
		aij[ 0 ][ start ] = INIFREQ; // the probability of going from bkg to the start of the forward motif
		
		prior[ start ] = priors + i; // the address to motif prior probability
		*prior[ start ] = INIPRIOR; // the initial value of this prior
		
		for( j = 0; j < motifs[ i ] ->PWM_width-1; j ++ )
		{
			aij[ start + j ][ start + j+1 ] = 1; // the probability of going from position j to j+1 of the motif
			prior[ start + j+1 ] = prior[ start + j ]; // the prior probability of any position within a motif is equal to the prior probability of transitioning to that motif
		}
		aij[ start + j ][ 0 ] = 1; // the probability of going from the motif end to bkg

		start = motif_r_start[ i ];
		aij[ 0 ][ start ] = INIFREQ; // the probability of going from bkg to the start of the reverse motif

		prior[ start ] = priors + i; // the address to motif prior probability
		*prior[ start ] = INIPRIOR; // the initial value of this prior

		for( j = 0; j < motifs[ i ] ->PWM_width-1; j ++ )
		{
			aij[ start + j ][ start + j+1 ] = 1; // the probability of going from position j to j+1 of the motif
			prior[ start + j+1 ] = prior[ start + j ]; // the prior probability of any position within a motif is equal to the prior probability of transitioning to that motif
		}
		aij[ start + j ][ 0 ] = 1; // the probability of going from the motif end to bkg
	}
	


	//////////////////////////////////////////////////////// initialize the PFMs and the matrix of the log of emission probabilities
	
	double *bjo[ MAX_STATES ]; // log of emission probabilities; bjo[ j ][ o ]: log of the probability of emission output o from state j
	
	double bkg[ 4 ] = { 0.25, 0.25, 0.25, 0.25 } ; // the log of emission probabilities for background
	bjo[ 0 ] = bkg;

	// copy the initial PFM of each motif, and point the bjo to the corresponding PFM array
	int n, x;
	for( i = 0; i < num_motifs; i ++ )
	{
		// forward PFM
		int start = motif_f_start[ i ];
		for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
		{
			motifs[ i ] ->opt_PFM_f[ x ] = new double[ 4 ];
			for( n = 0; n < 4; n ++ )
				motifs[ i ] ->opt_PFM_f[ x ][ n ] = motifs[ i ] ->PFM[ n ][ x ];
			bjo[ start + x ] = motifs[ i ] ->opt_PFM_f[ x ];
		}

		// reverse PFM
		start = motif_r_start[ i ];
		for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
		{
			motifs[ i ] ->opt_PFM_r[ x ] = new double[ 4 ];
			for( n = 0; n < 4; n ++ )
				motifs[ i ] ->opt_PFM_r[ x ][ n ] = motifs[ i ] ->opt_PFM_f[ motifs[ i ] ->PWM_width-1-x ][ 3-n ];
			bjo[ start + x ] = motifs[ i ] ->opt_PFM_r[ x ];
		}
		
		// initialize the memory for the final PFM, which will be the average of the forward and reverse PFMs
		for( n = 0; n < 4; n ++ )
			motifs[ i ] ->opt_PFM[ n ] = new double[ motifs[ i ] ->PWM_width ];
	}

		
	//////////////////////////////////////////////////////// initialize the alpha and beta matrices
	
	// find the maximum sequence length
	int max_seq_length = -1;
	for( i = 0; i < num_seqs; i ++ )
		if( max_seq_length < seqs[ i ] ->seq_length )
			max_seq_length = seqs[ i ] ->seq_length;
	
	// now initialize memory for alpha, beta, and the scaling factor
	double *alpha[ MAX_SEQ_LENGTH ]; // alpha[ t ][ i ]: log of alpha at time t for state i
	double *beta[ MAX_SEQ_LENGTH ]; // beta[ t ][ i ]: log of beta at time t for state i
	double *scale = new double[ max_seq_length ];
	scale[ 0 ] = 1; // no scaling at t=0;
	
	for( i = 0; i < max_seq_length; i ++ )
	{
		alpha[ i ] = new double[ num_states ];
		beta[ i ] = new double[ num_states ];
	}


	//////////////////////////////////////////////////////// run the expectation maximization algorithm
	double *gamma = new double[ num_states ]; // the probability of being in state i at time t, summed over all times t
	double *gamma0 = new double[ num_states ]; // the probability of being in state i at time 0
	double *g = new double[ num_states ]; // the probability of being in state i at time t
	double *gamma_o[ 4 ]; // the probability of being in state i at time t, summed over all times t that show emission o
	for( i = 0; i < 4; i ++ )
		gamma_o[ i ] = new double[ num_states ];
	double *xi[ MAX_STATES ]; // the probability of being in state i at time t and state j at time t+1, summed over all times t
	double *X[ MAX_STATES ]; // the probability of being in state i at time t and state j at time t+1
	for( i = 0; i < num_states; i ++ )
	{
		xi[ i ] = new double[ num_states ];
		X[ i ] = new double[ num_states ];
	}

	
	int iter, seq, motif;
	for( iter = 0; iter < MAX_ITER; iter ++ ) // perform multiple iterations of the algorithm
	{
		
		// initialize gamma and xi to zero
		for( i = 0; i < num_states; i ++ )
		{
			gamma[ i ] = 0;	
			gamma0[ i ] = 0;	
			
			for( j = 0; j < 4; j ++ )
				gamma_o[ j ][ i ] = 0;
				
			for( j = 0; j < num_states; j ++ )
				xi[ i ][ j ] = 0;
		}

	
		// for each sequence, calculate the alpha and beta separately, and then update gamma, gamma_o, and xi
		for( seq = 0; seq < num_seqs; seq ++ )
		//for( seq = 0; seq < 100; seq ++ )
			if( seqs[ seq ] ->positive )
			{
				cout << "."; cout.flush();
				///////////////////////// forward algorithm - calculate alpha


				////// initialize alpha for t=0
				int t = 0;
				int o = _n_index( seqs[ seq ] ->seq[ t ] ); // the emission output
				for( j = 0; j < num_states; j ++ )
					alpha[ t ][ j ] = *prior[ j ] * ( (o<0)? 1:bjo[ j ][ o ] );

			
				////// update alpha for t=1:T
				double total_scale = 1; // the product of all scaling factors for this sequence
				for( t = 1; t < seqs[ seq ] ->seq_length; t ++ )
				{
					o = _n_index( seqs[ seq ] ->seq[ t ] ); // the emission output
					scale[ t ] = 0;

					int motif_index = 0;
					for( j = 0; j < num_states; j ++ )
					{
						double logalpha = 0;
						if( j == 0 ) // j is the bkg, and can be reached only from itself and from motif ends
						{
							logalpha += alpha[ t-1 ][ 0 ] * aij[ 0 ][ j ];
							for( i = 0; i < num_motifs; i ++ )
							{
								// note that aij[ motif_f_end[ i ] ][ 0 ] = 1
								logalpha += alpha[ t-1 ][ motif_f_end[ i ] ]; // * aij[ motif_f_end[ i ] ][ j ];
								logalpha += alpha[ t-1 ][ motif_r_end[ i ] ]; // * aij[ motif_r_end[ i ] ][ j ];
							}
						}
						else if( motif_index < num_motifs && j == motif_f_start[ motif_index ] ) // j is a motif start, and can only be reached from bkg
							logalpha += alpha[ t-1 ][ 0 ] * aij[ 0 ][ j ];
						else if( motif_index < num_motifs && j == motif_r_start[ motif_index ] ) // j is a motif start, and can only be reached from bkg
						{
							logalpha += alpha[ t-1 ][ 0 ] * aij[ 0 ][ j ];
							motif_index ++;
						}
						else // j is neither the bkg, nor a motif start, and therefore can only be reached from j-1 with p=1
							logalpha += alpha[ t-1 ][ j-1 ]; // since aij[ j-1 ][ j ] = 1, there's no need for multiplication
						
						alpha[ t ][ j ] = logalpha * ( (o<0)? 1:bjo[ j ][ o ] );
						
						scale[ t ] += alpha[ t ][ j ];
					}
					
					for( j = 0; j < num_states; j ++ )
						alpha[ t ][ j ] /= scale[ t ];
						
					total_scale *= scale[ t ];
					
				}
				
				// calculate the probability of observing this sequence given the current model
				double seq_probability = 0;
				for( j = 0; j < num_states; j ++ )
					seq_probability += alpha[ seqs[ seq ] ->seq_length - 1 ][ j ];
				seq_probability *= total_scale;				
				cerr << seq_probability << char(9);


	
				///////////////////////// backward algorithm - calculate beta
			
				////// initialize beta for t=T
				t = seqs[ seq ] ->seq_length - 1;
				for( j = 0; j < num_states; j ++ )
					beta[ t ][ j ] = 1 / scale[ t ];
				
				for( t = seqs[ seq ] ->seq_length - 2; t >= 0; t -- )
				{
					o = _n_index( seqs[ seq ] ->seq[ t+1 ] ); // the emission output

					int motif_index = 0;
					for( i = 0; i < num_states; i ++ )
					{
						double logbeta = 0;
						if( i == 0 ) // i is the bkg, and can only go to itself and to motif starts
						{
							logbeta += beta[ t+1 ][ 0 ] * aij[ i ][ 0 ] * ( (o<0)? 1:bjo[ 0 ][ o ] );
							for( j = 0; j < num_motifs; j ++ )
							{
								logbeta += beta[ t+1 ][ motif_f_start[ j ] ] * aij[ i ][ motif_f_start[ j ] ] * ( (o<0)? 1:bjo[ motif_f_start[ j ] ][ o ] );
								logbeta += beta[ t+1 ][ motif_r_start[ j ] ] * aij[ i ][ motif_r_start[ j ] ] * ( (o<0)? 1:bjo[ motif_r_start[ j ] ][ o ] );
							}
						}
						else if( motif_index < num_motifs && i == motif_f_end[ motif_index ] ) // i is a motif end, and can only go to bkg
							logbeta += beta[ t+1 ][ 0 ] * ( (o<0)? 1:bjo[ 0 ][ o ] ); // * aij[ i ][ 0 ]; // aij[ i ][ 0 ] = 1 for motif ends
						else if( motif_index < num_motifs && i == motif_r_end[ motif_index ] ) // i is a motif end, and can only go to bkg
						{
							logbeta += beta[ t+1 ][ 0 ] * ( (o<0)? 1:bjo[ 0 ][ o ] ); // * aij[ i ][ 0 ]; // aij[ i ][ 0 ] = 1 for motif ends
							motif_index ++;
						}
						else // i is neither the bkg, nor a motif end, and therefore can only go to i+1 with p=1
							logbeta += beta[ t+1 ][ i+1 ] * ( (o<0)? 1:bjo[ i+1 ][ o ] ); // since aij[ i ][ i+1 ] = 1, there's no need for multiplication

						beta[ t ][ i ] = logbeta / scale[ t ];
					}
				}


				///////////////////////// update gamma, gamma_o, and xi
			
				for( t = 0; t < seqs[ seq ] ->seq_length - 1; t ++ )
				{
					//// update gamma and gamma_o
					int o = _n_index( seqs[ seq ] ->seq[ t ] ); // the emission output
					double normalizer = 0;
				
					for( i = 0; i < num_states; i ++ )
					{
						g[ i ] = alpha[ t ][ i ] * beta[ t ][ i ];
						normalizer += g[ i ];
					}
					
					for( i = 0; i < num_states; i ++ )
					{
						gamma[ i ] += g[ i ] / normalizer;

						if( t == 0 )
							gamma0[ i ] += g[ i ] / normalizer;

						if( o < 0 )
							for( n = 0; n < 4; n ++ )
								gamma_o[ n ][ i ] += g[ i ]/normalizer/4;
						else
							gamma_o[ o ][ i ] += g[ i ]/normalizer;
					}
				
					//// update xi
					o = _n_index( seqs[ seq ] ->seq[ t+1 ] ); // the emission output
					normalizer = 0;
				
					for( i = 0; i < num_states; i ++ )
						for( j = 0; j < num_states; j ++ )
						{
							X[ i ][ j ] =
								alpha[ t ][ i ] * aij[ i ][ j ] * ( (o<0)? 1:bjo[ j ][ o ] ) * beta[ t+1 ][ j ];
							normalizer += X[ i ][ j ];
						}

					for( i = 0; i < num_states; i ++ )
						for( j = 0; j < num_states; j ++ )
							xi[ i ][ j ] +=  X[ i ][ j ]/normalizer;
				}
			}
		
		///////////////////////// update the HMM

		// update priors
		double normalizer = 0;
		for( i = 0; i < num_states; i ++ ) // calculate the normalizer, so that the sum of priors is 1
		{
			normalizer += gamma0[ i ];
			*prior[ i ] = 0;
		}

		priors[ 0 ] = gamma0[ i ] / normalizer; // bkg prior

		// for the motifs, the prior of each position is the average of the prior of all positions in forward and reverse PFMs		
		int state_index = 1; // indexes the positions within the motifs
		for( i = 0; i < num_motifs; i ++ )
		{
			for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
			{
				// For every position in a motif, there are two states, corresponding to forward and reverse PFMs
				// therefore, for every position, two states should be added up
				priors[ i ] += gamma0[ state_index ] / normalizer;
				state_index ++;
				priors[ i ] += gamma0[ state_index ] / normalizer;
				state_index ++;
			}
			
			priors[ i ] /= ( motifs[ i ] ->PWM_width * 2 ); // take the average over the motif positions
		}
		
		// update aij
		aij[ 0 ][ 0 ] = xi[ 0 ][ 0 ]/gamma[ 0 ];
		for( i = 0; i < num_motifs; i ++ )
		{
			aij[ 0 ][ motif_f_start[ i ] ] = xi[ 0 ][ motif_f_start[ i ] ]/gamma[ 0 ];
			aij[ 0 ][ motif_r_start[ i ] ] = xi[ 0 ][ motif_r_start[ i ] ]/gamma[ 0 ];
			
			// take the average of the forward and reverse transition probabilities
			aij[ 0 ][ motif_f_start[ i ] ] = aij[ 0 ][ motif_r_start[ i ] ] =
				( aij[ 0 ][ motif_f_start[ i ] ] + aij[ 0 ][ motif_r_start[ i ] ] ) / 2;
		}

		// update bjo
		for( j = 0; j < num_states; j ++ )
			for( n = 0; n < 4; n ++ )
				bjo[ j ][ n ] = ( gamma_o[ n ][ j ]/gamma[ j ] + 0.001 ) / 1.004; // add some uncertainty to avoid zero probabilities
				
		/////////// average the updated forward and reverse PFMs
		for( i = 0; i < num_motifs; i ++ )
			for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
				for( n = 0; n < 4; n ++ )
					motifs[ i ] ->opt_PFM[ n ][ x ] =
					motifs[ i ] ->opt_PFM_f[ x ][ n ] =
					motifs[ i ] ->opt_PFM_r[ motifs[ i ] ->PWM_width-1-x ][ 3-n ] =
						( motifs[ i ] ->opt_PFM_f[ x ][ n ] +
						motifs[ i ] ->opt_PFM_r[ motifs[ i ] ->PWM_width-1-x ][ 3-n ] ) / 2;

		cerr << endl;
	}
}
