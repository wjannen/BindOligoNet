#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "myFold_newpen.h"

#define NEW_RULES 1

#ifndef NEW_RULES 
#error "Choice of terminal mismatch rules unspecified (#define NEW_RULES in myFold_MFT.c)" 
#endif

extern int nucenergy[6];
extern int g0;

void setup(void) {
	int i,j;
	
	/*
	  Allocate memory for the matrices
	*/
	for( i = 0; i < NUM_MATRICES; i++ ){
		the_matrices[i] = (int**) malloc( (len_s+5+DUMMY_LENGTH*2)*sizeof(int*) );
		if(the_matrices[i] == NULL){
			printf("Error allocating memory for matrix %d\n",i);
			exit(1);
		}
		for( j = 0; j <= len_s+DUMMY_LENGTH*2; j++ ){
			the_matrices[i][j] = (int *) malloc( (len_t+5+DUMMY_LENGTH*2)*sizeof(int) );
			if(the_matrices[i][j] == NULL){
				printf("Error allocating memory for matrix %d, row %d\n",i,j);
				exit(1);
			}
		}
	}
	/*
	  Allocate memory for the asymmetric loop matrices
	*/
	lu = (struct asymloopdata**) malloc( (len_s+1+DUMMY_LENGTH*2)*sizeof(struct asymloopdata*) );
	ld = (struct asymloopdata**) malloc( (len_s+1+DUMMY_LENGTH*2)*sizeof(struct asymloopdata*) );
	l1xn = (struct asymloopdata**) malloc( (len_s+1+DUMMY_LENGTH*2)*sizeof(struct asymloopdata*) );
	lnx1 = (struct asymloopdata**) malloc( (len_s+1+DUMMY_LENGTH*2)*sizeof(struct asymloopdata*) );
	for( i = 0; i <= len_s + DUMMY_LENGTH*2; i++){
		lu[i] = (struct asymloopdata*) malloc( (len_t + 1 + DUMMY_LENGTH*2)*sizeof(struct asymloopdata) );
		ld[i] = (struct asymloopdata*) malloc( (len_t + 1 + DUMMY_LENGTH*2)*sizeof(struct asymloopdata) );
		l1xn[i] =  (struct asymloopdata*) malloc( (len_t + 1 + DUMMY_LENGTH*2)*sizeof(struct asymloopdata) );
		lnx1[i] =  (struct asymloopdata*) malloc( (len_t + 1 + DUMMY_LENGTH*2)*sizeof(struct asymloopdata) );
	}
	//assign aliases for quick typing
	pxp = the_matrices[PXP_M];	
	bu = the_matrices[BU_M];
	bd = the_matrices[BD_M];
	bu2 = the_matrices[BU2_M];
	bd2 = the_matrices[BD2_M];
	f = the_matrices[F_M];
	
	//get the nearest neighbor energies from the files
	get_stack_energies();
	get_iloop_energies();
	
	//allocate memory for the pairing matrix
	pairs = (int **)malloc( (len_s+DUMMY_LENGTH)*sizeof(int*) );
	if(pairs == NULL){
		printf( "Failed to allocate memory for pairs\n");
		exit(1);
	}
	for( i = 0; i < len_s + DUMMY_LENGTH; i++ ){
		pairs[i] = (int *)malloc( (len_t+DUMMY_LENGTH)*sizeof(int) );
		
		if(pairs[i] == NULL){
			printf( "Failed to allocate memory for pairs[%d]\n", i);
			exit(1);
		}
		
	}
	
	
	//allocate some initial memory for the structure list
	structures = (struct loop*) malloc( 10*sizeof(struct loop) );
	
	
}

/*
  takes a character string and encodes it into the numerical array
  used in the algorithm. direction is a flag which tells whether the output is
  5' to 3' or 3' to 5'. If direction is 1, the output is 5' to 3'. Otherwise it is 3' to 5'.
  The s should be 5' to 3' and the t should be 3' to 5'.
*/
void encode(int* dest, const char* src, const int length, const int direction){
	int i;
	
	//add the buffers at the beginning and end
	for( i = 0; i < DUMMY_LENGTH; i++ ){
		dest[i] = X; //pad the beginning with 'X's
		dest[i+length+DUMMY_LENGTH] = X; //pad the end with 'X's
	}
	
	//encode the sequences
	// A=0, C=1, G=2, U=3, X=4, *=5 ( X is non-pairing misc. base and * indicates a gap) 
	if(direction == 1){ // the s sequence
		for( i = 0; i < length; i++){
			
			char c = toupper( src[i] ); //convert all characters to upper case
			
			switch(c){
			case 'A':
				dest[i+DUMMY_LENGTH] = A;
				break;
			case 'C':
				dest[i+DUMMY_LENGTH] = C;
				break;
			case 'G':
				dest[i+DUMMY_LENGTH] = G;
				break;
			case 'U':
				dest[i+DUMMY_LENGTH] = U;
				break;
			case 'T':
				dest[i+DUMMY_LENGTH] = U; // rna vs. dna
				break;
			default:
				dest[i+DUMMY_LENGTH] = X; //if not a,c,g or u, make it an 'X'
			}
		}
	}
	
	else{ // or the t sequence
		for( i = 0; i < length; i++){ // t is entered 5' to 3'. Flip it to 3' to 5'
			
			char c = toupper( src[length - i -1] );
			
			switch(c){
			case 'A':
				dest[i+DUMMY_LENGTH] = A;
				break;
			case 'C':
				dest[i+DUMMY_LENGTH] = C;
				break;
			case 'G':
				dest[i+DUMMY_LENGTH] = G;
				break;
			case 'U':
				dest[i+DUMMY_LENGTH] = U;
				break;
			case 'T':
				dest[i+DUMMY_LENGTH] = U; // rna vs. dna
				break;
			default:
				dest[i+DUMMY_LENGTH] = X; //if not a,c,g or u, make it an 'X'
			}
		} 
	}
	
}


/*
  ------------------------ fill_arrays ----------------------
  proceeds through the process of filling the matrices 
*/
void fill_arrays(void){

	int i,j,k;
	
	/*int penalize = 0; //removed 31 Oct 2011 JH: for oligo-oligo binding, use non-MFT version
	if(len_t > PENALTY_THRESHOLD)
		penalize = 1;*/

	/*************************
	 **   INITIALIZIATION   **
	 *************************/

	//initialize all matrices to negative infinity (-9999)
	for( i = 0; i < NUM_MATRICES; i++ ){
		for( j = 0; j <= len_s+DUMMY_LENGTH*2; j++ ){
			for( k = 0; k <= len_t+DUMMY_LENGTH*2; k++ ){
				the_matrices[i][j][k] = NEG_INFINITY;
			}
		}
	}
	
	for( i = 0; i <= len_s+DUMMY_LENGTH*2; i++){
		for( j = 0; j <= len_t+DUMMY_LENGTH*2; j++ ){
			lu[i][j].score = -23423; //NEG_INFINITY*3;
			ld[i][j].score = -28923; //NEG_INFINITY*3;
			
			l1xn[i][j].score = -29387; //NEG_INFINITY*3;
			lnx1[i][j].score = -9823; //NEG_INFINITY*3;
			
			lu[i][j].up = lu[i][j].down = ld[i][j].up = ld[i][j].down = 99;
			l1xn[i][j].up = l1xn[i][j].down = lnx1[i][j].up = lnx1[i][j].down = 99;
		}
	}
	
	//reallocate memory to structures to make it an appropriate size
	structures = (struct loop*) realloc(structures, (len_s+len_t+DUMMY_LENGTH*2)*sizeof(struct loop));
	numstructures = 0;
	
	/***********************
	 **     ALGORITHM     **
	 ***********************/
	
	/*
	  Equivalents to the paper (diagram and supplementary material):
		 au_penalty is the \delta G_{AU/GU} term
		 pxp :  M matrix
		 bu  :  B
		 bu2 :  B2
		 bd  :  b
		 bd2 :  b2
		 lu  :  A
		 ld  :  a
	 */
	//consider the initial conditions
	//when you are computing matrix entries too small for proper bulge loop calculations
	for( i = DUMMY_LENGTH; i <= len_s+DUMMY_LENGTH+1; i++ ){
		for( j = DUMMY_LENGTH; j <= len_t+DUMMY_LENGTH+1; j++ ){
			//printf("(%d,%d) :: %c %c // %d\n",i,j,num_to_base(s[i]),num_to_base(t[j]),match(i+1,j+1));
			int t1, t2, t3, t4, t5, t6; // temporary variables 1-6

			//make sure that matrix only gets filled at valid pairings.
			//if(valid_pair(i,j)) {
				t1 = pxp[i-1][j-1] + match(i,j) + affinePent(j,j); // i, j pair         #aa
				t2 = pxp[i-2][j-2] + i11(i,j) + affinePent(j-1,j); // 1x1 internal loop #bb
				t3 = pxp[i-2][j-3] + i12(i,j) + affinePent(j-2,j); // 1x2 internal loop #cc
				t4 = pxp[i-3][j-2] + i21(i,j) + affinePent(j-1,j); // 2x1 internal loop #dd
				t5 = pxp[i-3][j-3] + i22(i,j) + affinePent(j-2,j); // 2x2 internal loop #ee
				t6 = pxp[i-2][j-1] + matchAlt(i-2,i,j-1,j) - SINGLE_GAP_PENALTY + retnucenergy(j); // bulge loop in s #ff
				pxp[i][j] = MAX(t1, MAX(t2, MAX(t3, MAX(t4, MAX(t5, t6)))));

				t1 = pxp[i-1][j-2] + matchAlt(i-1,i,j-2,j) - SINGLE_GAP_PENALTY + affinePent(j-1,j); // bulge loop in t  #gg
				//t6 = startmatch(i,j) - INITIATION_PENALTY - g0 + retnucenergy(j); // first pair                               #ww
				t6 = startmatch(i,j) - INITIATION_PENALTY - g0 + affinePent(j,j); // first pair                               #ww
				pxp[i][j] = MAX(pxp[i][j], MAX(t1, t6));

			////make sure that matrix only gets filled at valid pairings.
			if(valid_pair(i-1,j-1)) {
				t2 = bu[i-1][j-1] + au_penalty(i,j) + retnucenergy(j); // closing bulge loop in s, i and j binding       #hh
				t3 = bd[i-1][j-1] + au_penalty(i,j) + retnucenergy(j); // bulge loop in t closing, i and j binding       #ii
				t4 = lu[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j); // terminal match after asym internal loop     #jj
				t5 = ld[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j); // terminal match after asym internal loop     #kk
				pxp[i][j] = MAX(pxp[i][j], MAX(t2, MAX(t3, MAX(t4, t5))));
				//pxp[i][j] = MAX(pxp[i][j], MAX(t1, MAX(t2, MAX(t3, MAX(t4, MAX(t5, t6))))));

				t1 = MAX(bu2[i-1][j-1], bd2[i-1][j-1]) + au_penalty(i,j) + retnucenergy(j); // closing bulge loop               #ll
				t2 = MAX(l1xn[i-1][j-1].score, lnx1[i-1][j-1].score) + tmatch2(i,j) + retnucenergy(j); //closing internal loop  #mm
				pxp[i][j] = MAX(pxp[i][j], MAX(t1, t2));
			}

			t1 = pxp[i-2][j] - MULTI_GAP_AFFINE - MULTI_GAP_PENALTY*2 + au_penalty(i-2,j); // closing (<7) s gap     #nn
			t2 = bu[i-1][j] - MULTI_GAP_PENALTY; // increasing length of (<7) s gap                                  #oo
			bu[i][j] = MAX(t1,t2);
			
			//t1 = pxp[i-7][j] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2; // closing (>6) s gap   #pp
			t1 = pxp[i-7][j] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2 + au_penalty(i-7,j); // closing (>6) s gap   #pp
			t2 = bu2[i-1][j] - MULTI_GAP_PENALTY_2; // increasing length of (>6) s gap                               #qq
			bu2[i][j] = MAX(t1, t2);

			t1 = pxp[i][j-2] - MULTI_GAP_AFFINE - MULTI_GAP_PENALTY*2 + au_penalty(i,j-2) + affinePent(j-1,j);  // closing (<7) t gap #rr
			t2 = bd[i][j-1] - MULTI_GAP_PENALTY + retnucenergy(j); // increasing length of (<7) t gap                       #ss
			bd[i][j] = MAX(t1, t2);

			//t1 = pxp[i][j-7] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2 + au_penalty(i,j-7) + 7*mean_pen;// closing (>6) t gap #tt
			t1 = pxp[i][j-7] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2 + au_penalty(i,j-7) + affinePent(j-6,j);// closing (>6) t gap #tt
			t2 = bd2[i][j-1] - MULTI_GAP_PENALTY_2 + retnucenergy(j); // increasing length of (>6) t gap                    #uu
			bd2[i][j] = MAX(t1, t2);
			
			lu_policy(i, j, &lu[i][j]); // secondary structure restriction penalties now embedded in internal loop subs
			ld_policy(i, j, &ld[i][j]); 
			
			l1xn_policy(i, j, &l1xn[i][j]);
			lnx1_policy(i, j, &lnx1[i][j]);
			
			//f[i][j] = pxp[i][j] + finalmatch(i,j);
			f[i][j] = pxp[i-1][j-1] + finalmatch(i,j);
			/*if((pxp[i][j] > 0) && !valid_pair(i-1,j-1))
				printf("(%c,%c)",num_to_base(s[i-1]),num_to_base(t[j-1]));*/
		}
	}	       
	
}

/*
  
  ------------------------------- alignment ---------------------------
  Do the traceback.
  At each entry, figure out which matrix you come from.
  
  pass it the starting location (i,j) and starting matrix, matrix_index
  
  len_align is a global variable determining the current index of the alignment strings
  
*/

int alignment(const int i, const int j, const int matrix_index, const int prefill){
	/*if((matrix_index == PXP_M) || (matrix_index == F_M)) {
		//printf("%d,%d,%d,%c,%c\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH,the_matrices[matrix_index][i][j],num_to_base(s[i-1]),num_to_base(t[j-1]));
		printf("(%d,%d) Total: %d, Match: %d, Removal: %d [%c,%c]\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH,the_matrices[matrix_index][i][j],match(i,j),affinePent(j,j),num_to_base(s[i-1]),num_to_base(t[j-1]));
	}*/
	//printf("matrix_index: %d, i:%d j:%d, G:%d\n", matrix_index,i-DUMMY_LENGTH,j-DUMMY_LENGTH, the_matrices[matrix_index][i][j]);
	/*int mean_pen = 0; //removed 31 Oct 2011 by JH; use non-MFT Bindigo for oligo-oligo binding
	if(len_t > PENALTY_THRESHOLD)
		mean_pen = MEAN_PENALTY;*/

	if( (i==DUMMY_LENGTH) && (j==DUMMY_LENGTH)  ){
		len_align = 0; //if you are here, the traceback is done
	} else if( (i> DUMMY_LENGTH) && (j > DUMMY_LENGTH) && (matrix_index == PXP_M) ){
		if( pxp[i][j] == pxp[i-1][j-1] + match(i,j) + affinePent(j,j)){ // #aa 
			alignment(i-1, j-1, PXP_M, prefill); 
		}
		else if( pxp[i][j] == pxp[i-2][j-2] + i11(i,j) + affinePent(j-1,j)){ // #bb
			
			alignment(i-2, j-2, PXP_M, prefill);
			//printf("1x1: %d,%d\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = num_to_base(s[i-2]);
				aligned_t[len_align] = num_to_base(t[j-2]);
			}
		}
		else if( pxp[i][j] == pxp[i-2][j-3] + i12(i,j) + affinePent(j-2,j)){ // #cc
			alignment(i-2, j-3, PXP_M, prefill);
			//printf("1x2: %d,%d\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = '*';
				aligned_t[len_align] = num_to_base(t[j-3]);
			
				len_align++;

				aligned_s[len_align] = num_to_base(s[i-2]);
				aligned_t[len_align] = num_to_base(t[j-2]); // 09wkj CHANGED FROM j-3 to j-2, swapped *
			}
		}
		else if(  pxp[i][j] == pxp[i-3][j-2] + i21(i,j) + affinePent(j-1,j)){ // #dd
			alignment(i-3, j-2, PXP_M, prefill);
			//printf("2x1: %d,%d\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = num_to_base(s[i-3]);
				aligned_t[len_align] = '*'; 
				
				len_align++;
				aligned_s[len_align] = num_to_base(s[i-2]); // 09wkj CHANGED FROM j-2 to i-2 4/20
				aligned_t[len_align] = num_to_base(t[j-2]); // 09wkj CHANGED FROM j-3 to j-2, swapped *
			}
		}
		else if( pxp[i][j] == pxp[i-3][j-3] + i22(i,j) + affinePent(j-2,j)){ // #ee
			alignment(i-3, j-3, PXP_M, prefill);
			//printf("2x2: %d,%d\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = num_to_base(s[i-3]);
				aligned_t[len_align] = num_to_base(t[j-3]);
				
				len_align++;
				aligned_s[len_align] = num_to_base(s[i-2]);
				aligned_t[len_align] = num_to_base(t[j-2]);
			}
		}
		else if( pxp[i][j] == pxp[i-2][j-1] + matchAlt(i-2,i,j-1,j) - SINGLE_GAP_PENALTY + retnucenergy(j)){ // #ff
			alignment(i-2, j-1, PXP_M, prefill);
			//printf("s bulge: %d,%d\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = num_to_base(s[i-2]);
				aligned_t[len_align] = '*';
			}
			
		}
		else if( pxp[i][j]  == pxp[i-1][j-2] + matchAlt(i-1,i,j-2,j) - SINGLE_GAP_PENALTY + affinePent(j-1,j)){ // #gg
			//printf("t bulge: %d,%d\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
			alignment(i-1, j-2, PXP_M, prefill);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = '*';
				aligned_t[len_align] = num_to_base(t[j-2]);
			}
			
		}
		else if( pxp[i][j] == bu[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #hh
			alignment(i-1, j-1, BU_M, prefill);
		}
		else if( pxp[i][j] == bd[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #ii
			//printf("%d %d ??\n",i,j);
			alignment(i-1, j-1, BD_M, prefill);
		}
		else if( pxp[i][j] == bu2[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #ll (1/2)
			alignment(i-1, j-1, BU2_M, prefill);
		}
		else if( pxp[i][j] == bd2[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #ll (2/2)
			alignment(i-1, j-1, BD2_M, prefill);
		}
		else if( pxp[i][j] == l1xn[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #mm (1/2)
			//printf(" in l1xn \n");
			alignment(i-l1xn[i-1][j-1].up-1, j-l1xn[i-1][j-1].down-1, PXP_M, prefill);
			
			if(prefill)
				bigfill(i,j, l1xn[i-1][j-1].up, l1xn[i-1][j-1].down);
		}
		else if( pxp[i][j] == lnx1[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #mm (2/2)
			//printf(" in lnx1  \n");
			
			alignment(i-lnx1[i-1][j-1].up-1, j-lnx1[i-1][j-1].down-1, PXP_M, prefill);

			if(prefill)
				bigfill(i,j, lnx1[i-1][j-1].up, lnx1[i-1][j-1].down);
		}
		else if( pxp[i][j] == lu[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #jj
			//printf(" in lu \n");
			alignment(i-lu[i-1][j-1].up-1, j-lu[i-1][j-1].down-1, PXP_M, prefill);
			
			if(prefill)
				bigfill(i,j, lu[i-1][j-1].up, lu[i-1][j-1].down);
		}
		else if( pxp[i][j] == ld[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #kk
			//printf("%d,%d in ld  \n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
			
			alignment(i-ld[i-1][j-1].up-1, j-ld[i-1][j-1].down-1, PXP_M, prefill);

			if(prefill)
				bigfill(i,j, ld[i-1][j-1].up, ld[i-1][j-1].down);
		} else { //end of alignment.
			//printf("prefix fill: %d,%d\n",i-DUMMY_LENGTH-1,j-DUMMY_LENGTH-1);
			//len_align = 0;
			if(prefill) {
				fill_prefix(i-1,j-1);
			}
		}
		
		//printf("(m:%d f:%d)[i: %c j:%c]",firstcomp(s[i-1],t[j-1]),prefill,num_to_base(s[i-1]),num_to_base(t[j-1]));
		//printf("(m:%d f:%d)[i: %d j:%d]",firstcomp(s[i-1],t[j-1]),prefill,s[i-1],t[j-1]);


		//printf("paired: %d,%d @ %d\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH,f[i+1][j+1]);
		if(prefill) {
			len_align++;

			aligned_s[len_align] = tolower(num_to_base(s[i-1]));
			aligned_t[len_align] = tolower(num_to_base(t[j-1]));
		}

		pairs[i-DUMMY_LENGTH-1][j-DUMMY_LENGTH-1] = paired;

		if(!prefill) {
			if(f[f_i][f_j] > dgJ[j - DUMMY_LENGTH - 1]) {
				dgJ[j - DUMMY_LENGTH - 1] = f[f_i][f_j];
			}
		}
	} else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BU_M) ){
		
		if( bu[i][j] == pxp[i-2][j] - MULTI_GAP_AFFINE - MULTI_GAP_PENALTY*2 + au_penalty(i-2,j) ){ // #nn
			alignment(i-2, j, PXP_M,prefill);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = num_to_base( s[i-2] );
				aligned_t[len_align] = '*';
				
				len_align++;
				aligned_s[len_align] = num_to_base( s[i-1] );
				aligned_t[len_align] = '*';
			}
		}
		else if( bu[i][j] == bu[i-1][j] - MULTI_GAP_PENALTY){ // #oo
			alignment(i-1, j, BU_M,prefill);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = num_to_base( s[i-1] );
				aligned_t[len_align] = '*';
			}
		}
	}
	else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BD_M) ){
		
		if( bd[i][j] == pxp[i][j-2] - MULTI_GAP_AFFINE - MULTI_GAP_PENALTY*2 + au_penalty(i,j-2) + affinePent(j-1,j) ){ // #rr
			//printf("%d,%d !?\n",i,j);
			alignment(i, j-2, PXP_M, prefill);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = '*';
				aligned_t[len_align] = num_to_base( t[j-2] );
				
				len_align++;
				aligned_s[len_align] = '*';
				aligned_t[len_align] = num_to_base( t[j-1] );
			}
		}
		else if( bd[i][j] == bd[i][j-1] - MULTI_GAP_PENALTY + retnucenergy(j)){ // #ss
			//printf("%d,%d !\n",i,j);
			alignment(i, j-1, BD_M, prefill);
			
			if(prefill) {
				len_align++;
				aligned_s[len_align] = '*';
				aligned_t[len_align] = num_to_base( t[j-1] );
			}
		}
	}
	else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BU2_M) ){
		
		if( bu2[i][j] == pxp[i-7][j] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2 + au_penalty(i-7,j) ){ // #pp
			int ii;
			alignment(i-7, j, PXP_M, prefill);

			//bigfill(i,j,7,0);
			if(prefill) {
				for(ii = 7; ii > 0; --ii){
					len_align++;
					aligned_s[len_align] = num_to_base( s[i-ii] );
					aligned_t[len_align] = '*';
				}
			}
		}
		else if( bu2[i][j] == bu2[i-1][j] - MULTI_GAP_PENALTY_2){ // #qq
			alignment(i-1, j, BU2_M, prefill);

			if(prefill) {
				len_align++;
				aligned_s[len_align] = num_to_base( s[i-1] );
				aligned_t[len_align] = '*';
			}
		}
	} 
	else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BD2_M) ){
		
		if( bd2[i][j] == pxp[i][j-7] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2 + au_penalty(i,j-7) + affinePent(j-6,j)){ // #tt
			int jj;
			alignment(i, j-7, PXP_M, prefill);

			//bigfill(i,j,0,7);			
			if(prefill) {
				for(jj = 7; jj > 0; --jj){
					len_align++;
					aligned_s[len_align] = '*';
					aligned_t[len_align] = num_to_base( t[j-jj] );
				}
			} 
		}
		else if( bd2[i][j] == bd2[i][j-1] - MULTI_GAP_PENALTY_2 + retnucenergy(j)){ // #uu
			alignment(i, j-1, BD2_M, prefill);

			if(prefill) {
				len_align++;
				aligned_s[len_align] = '*';
				aligned_t[len_align] = num_to_base( t[j-1] );
			}
		}
	} 
	
	else if( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH) && (matrix_index == F_M) ){ 
		if( f[i][j] == pxp[i-1][j-1] + finalmatch(i,j) ){// + mean_pen){ // #vv
			alignment(i-1,j-1, PXP_M, prefill);
			
			if(prefill)
				len_align++;
		}
	} else { //end of alignment
		//printf("not paired: %d,%d\n",i-DUMMY_LENGTH-1,j-DUMMY_LENGTH-1);
		
		len_align = 0;
		
		//fill in the beginning of the strings
		if(prefill) { //if we are only interested in o(j), then no prefix fill.
			fill_prefix(i,j);
		}
	}
	
	return 0;
}

int valid_pair(int i,int j) {
	if((s[i] == 0) && (t[j] == 3)) {
		return 1;
	} else if((s[i] == 3) && (t[j] == 0)) {
		return 1;
	} else if((s[i] == 1) && (t[j] == 2)) {
		return 1;
	} else if((s[i] == 2) && (t[j] == 1)) {
		return 1;
	} else if((s[i] == 2) && (t[j] == 3)) {
		return 1;
	} else if((s[i] == 3) && (t[j] == 2)) {
		return 1;
	} else {
		return 0;
	}
}

//In this case, we only wish to perform the traceback and output which sites are paired, rather than print a formatted string.

//int traceback(const int i, const int j, const int matrix_index){
//	
//	//printf("matrix_index: %d, i:%d j:%d, G:%d\n", matrix_index,i-DUMMY_LENGTH,j-DUMMY_LENGTH, the_matrices[matrix_index][i][j]);
//	/*int mean_pen = 0; //removed 31 Oct 2011 by JH; use non-MFT Bindigo for oligo-oligo binding
//	if(len_t > PENALTY_THRESHOLD)
//		mean_pen = MEAN_PENALTY;*/
//	
//	if( (i==DUMMY_LENGTH) && (j==DUMMY_LENGTH)  ){
//		len_align = 0; //if you are here, the traceback is done
//	}
//	else if( (i> DUMMY_LENGTH) && (j > DUMMY_LENGTH) && (matrix_index == PXP_M) ){
//		
//		if( pxp[i][j] == pxp[i-1][j-1] + match(i,j) + retnucenergy(j)){ // #aa
//
//			traceback(i-1, j-1, PXP_M);
//			
//		}
//		else if( pxp[i][j] == pxp[i-2][j-2] + i11(i,j) + affinePent(j-1,j)){ // #bb
//			
//			traceback(i-2, j-2, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base(s[i-2]);
//			aligned_tB[len_align] = num_to_base(t[j-2]);
//		}
//		else if( pxp[i][j] == pxp[i-2][j-3] + i12(i,j) + affinePent(j-2,j)){ // #cc
//			traceback(i-2, j-3, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = '*';
//			aligned_tB[len_align] = num_to_base(t[j-3]);
//			
//			len_align++;
//
//			aligned_sB[len_align] = num_to_base(s[i-2]);
//			aligned_tB[len_align] = num_to_base(t[j-2]); // 09wkj CHANGED FROM j-3 to j-2, swapped *
//		}
//		else if(  pxp[i][j] == pxp[i-3][j-2] + i21(i,j) + affinePent(j-1,j)){ // #dd
//			traceback(i-3, j-2, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base(s[i-3]);
//			aligned_tB[len_align] = '*'; 
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base(s[i-2]); // 09wkj CHANGED FROM j-2 to i-2 4/20
//			aligned_tB[len_align] = num_to_base(t[j-2]); // 09wkj CHANGED FROM j-3 to j-2, swapped *
//		}
//		else if( pxp[i][j] == pxp[i-3][j-3] + i22(i,j) + affinePent(j-2,j)){ // #ee
//			traceback(i-3, j-3, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base(s[i-3]);
//			aligned_tB[len_align] = num_to_base(t[j-3]);
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base(s[i-2]);
//			aligned_tB[len_align] = num_to_base(t[j-2]);
//		}
//		else if( pxp[i][j] == pxp[i-2][j-1] + matchAlt(i-2,i,j-1,j) - SINGLE_GAP_PENALTY + retnucenergy(j)){ // #ff
//			traceback(i-2, j-1, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base(s[i-2]);
//			aligned_tB[len_align] = '*';
//			
//		}
//		else if( pxp[i][j]  == pxp[i-1][j-2] + matchAlt(i-1,i,j-2,j) - SINGLE_GAP_PENALTY + affinePent(j-1,j)){ // #gg
//			traceback(i-1, j-2, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = '*';
//			aligned_tB[len_align] = num_to_base(t[j-2]);
//			
//		}
//		else if( pxp[i][j] == bu[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #hh
//			traceback(i-1, j-1, BU_M);
//		}
//		else if( pxp[i][j] == bd[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #ii
//			traceback(i-1, j-1, BD_M);
//		}
//		else if( pxp[i][j] == bu2[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #ll (1/2)
//			traceback(i-1, j-1, BU2_M);
//		}
//		else if( pxp[i][j] == bd2[i-1][j-1] + au_penalty(i,j) + retnucenergy(j)){ // #ll (2/2)
//			traceback(i-1, j-1, BD2_M);
//		}
//		else if( pxp[i][j] == l1xn[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #mm (1/2)
//			printf(" in l1xn \n");
//			traceback(i-l1xn[i-1][j-1].up-1, j-l1xn[i-1][j-1].down-1, PXP_M);
//			
//			bigfill(i,j, l1xn[i-1][j-1].up, l1xn[i-1][j-1].down);
//		}
//		else if( pxp[i][j] == lnx1[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #mm (2/2)
//			printf(" in lnx1  \n");
//
//			traceback(i-lnx1[i-1][j-1].up-1, j-lnx1[i-1][j-1].down-1, PXP_M);
//			bigfill(i,j, lnx1[i-1][j-1].up, lnx1[i-1][j-1].down);
//		}
//		else if( pxp[i][j] == lu[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #jj
//			printf(" in lu \n");
//			traceback(i-lu[i-1][j-1].up-1, j-lu[i-1][j-1].down-1, PXP_M);
//			
//			bigfill(i,j, lu[i-1][j-1].up, lu[i-1][j-1].down);
//		}
//		else if( pxp[i][j] == ld[i-1][j-1].score + tmatch2(i,j) + retnucenergy(j)){ // #kk
//			printf(" in ld  \n");
//			
//			traceback(i-ld[i-1][j-1].up-1, j-ld[i-1][j-1].down-1, PXP_M);
//			bigfill(i,j, ld[i-1][j-1].up, ld[i-1][j-1].down);
//		} 
//		else { // TODO anything here?
//			fill_prefix(i-1,j-1);
//		}
//		len_align++;
//		aligned_sB[len_align] = tolower(num_to_base(s[i-1])); 
//		aligned_tB[len_align] = tolower(num_to_base(t[j-1]));
//		
//		pairs[i-DUMMY_LENGTH-1][j-DUMMY_LENGTH-1] = paired;
//		
//		// pairs[len_align-1] = '|';
//		// aligned_sB[len_align+1] = num_to_base(s[i]);
//		//aligned_tB[len_align+1] = num_to_base(t[j]);
//		
//	}
//	else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BU_M) ){
//		
//		if( bu[i][j] == pxp[i-2][j] - MULTI_GAP_AFFINE - MULTI_GAP_PENALTY*2 + au_penalty(i-2,j) ){ // #nn
//			traceback(i-2, j, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base( s[i-2] );
//			aligned_tB[len_align] = '*';
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base( s[i-1] );
//			aligned_tB[len_align] = '*';
//		}
//		else if( bu[i][j] == bu[i-1][j] - MULTI_GAP_PENALTY){ // #oo
//			traceback(i-1, j, BU_M);
//			
//			len_align++;
//			aligned_sB[len_align] = num_to_base( s[i-1] );
//			aligned_tB[len_align] = '*';
//		}
//	}
//	else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BD_M) ){
//		
//		if( bd[i][j] == pxp[i][j-2] - MULTI_GAP_AFFINE - MULTI_GAP_PENALTY*2 + au_penalty(i,j-2) + affinePent(j-1,j) ){ // #rr
//			traceback(i, j-2, PXP_M);
//			
//			len_align++;
//			aligned_sB[len_align] = '*';
//			aligned_tB[len_align] = num_to_base( t[j-2] );
//			
//			len_align++;
//			aligned_sB[len_align] = '*';
//			aligned_tB[len_align] = num_to_base( t[j-1] );
//			
//		}
//		else if( bd[i][j] == bd[i][j-1] - MULTI_GAP_PENALTY + retnucenergy(j)){ // #ss
//			traceback(i, j-1, BD_M);
//			
//			len_align++;
//			aligned_sB[len_align] = '*';
//			aligned_tB[len_align] = num_to_base( t[j-1] );
//		}
//	}
//	else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BU2_M) ){
//		
//		if( bu2[i][j] == pxp[i-7][j] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2 + au_penalty(i-7,j) + affinePent(j-6,j) ){ // #pp
//			int ii;
//			traceback(i-7, j, PXP_M);
//
//			//bigfill(i,j,7,0);
//			for(ii = 7; ii > 0; --ii){
//				len_align++;
//				aligned_sB[len_align] = num_to_base( s[i-ii] );
//				aligned_tB[len_align] = '*';
//			}
//		}
//		else if( bu2[i][j] == bu2[i-1][j] - MULTI_GAP_PENALTY_2){ // #qq
//			traceback(i-1, j, BU2_M);
//			len_align++;
//			aligned_sB[len_align] = num_to_base( s[i-1] );
//			aligned_tB[len_align] = '*';
//		}
//	} 
//	else if ( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH ) && (matrix_index == BD2_M) ){
//		
//		if( bd2[i][j] == pxp[i][j-7] - MULTI_GAP_AFFINE_2 - MULTI_GAP_PENALTY_2 + au_penalty(i,j-7) + retnucenergy(j)){ // #tt
//			int jj;
//			traceback(i, j-7, PXP_M);
//
//			//bigfill(i,j,0,7);			
//			for(jj = 7; jj > 0; --jj){
//				len_align++;
//				aligned_sB[len_align] = '*';
//				aligned_tB[len_align] = num_to_base( t[j-jj] );
//			}
//
//		}
//		else if( bd2[i][j] == bd2[i][j-1] - MULTI_GAP_PENALTY_2 + retnucenergy(j)){ // #uu
//			traceback(i, j-1, BD2_M);
//			len_align++;
//			aligned_sB[len_align] = '*';
//			aligned_tB[len_align] = num_to_base( t[j-1] );
//		}
//	} 
//	
//	else if( (i > DUMMY_LENGTH) && (j > DUMMY_LENGTH) && (matrix_index == F_M) ){
//		
//		if( f[i][j] == pxp[i-1][j-1] + finalmatch(i,j) ){// + mean_pen){ // #vv
//			traceback(i-1,j-1, PXP_M);
//			
//			len_align++;
//		}
//	}
//	else  {
//		len_align = 0;
//		
//		//fill in the beginning of the strings
//		fill_prefix(i,j);
//		
//	}
//	
//	return 0;
//}

void fill_prefix(const int i, const int j ){
	int max_prefix = MAX(i,j); //fill in the longest prefix first
	int x; 
    
	// printf( "i: %d, j: %d\n", i, j);
	if( max_prefix == i ){
		for( x = DUMMY_LENGTH; x < i; x++ )
			aligned_s[x] = num_to_base(s[x]);
		
		for( x = DUMMY_LENGTH; x < j; x++ )
			aligned_t[ x + i - j] = num_to_base(t[x]);
		
		len_align = max_prefix - 1;
	} else {
		
		for( x = DUMMY_LENGTH; x < j; x++ )
			aligned_t[x] = num_to_base(t[x]);
		
		for( x = DUMMY_LENGTH; x < i; x++ )
			aligned_s[ x + j - i] = num_to_base(s[x]);
		
		
		len_align = max_prefix - 1;
		
    }
	
}


//convert the sequence code number to the corresponding character
char num_to_base(const int i) {
	
	return code[i];
	
}

//return filename concatenated with the path
void stackstrcat(char * str, char * dest) {
	char * env = getenv("BINPATH");
	char buf[strlen(env) + strlen(str) + 1];

	strcpy(buf,env);

	if(env[strlen(env) - 1] != '/')
		strcat(buf,"/");

	strcat(buf,str); 

	strcpy(dest,buf);
}

/*
  scan in the files, parse the data files and fill the array.
*/
void get_stack_energies() {
	int temp, i, j;
	int count = 6*6*6*6; // stack, tstack, and tstacki all 6*6*6*6 ints

	/*
	  read stack data files
	*/
	FILE *filestack,  *filetstack, *filetstacki;

	char * wd; 
	char * path = malloc(strlen(getenv("BINPATH")) + 20);

	if(wd = getenv("BINPATH")) {
		stackstrcat("stack.bin",path); 
		if(filestack = fopen(path,"r")) {
			temp = fread(&stack_energy, sizeof(int), count, filestack);
			if(temp != count)
				printf("Only wrote %d elements of stack\n",temp);
			fclose(filestack);
		} else {
			printf("Could not open stacking binary (stack.bin) for reading!\n");
			exit(1);
		}
		
		stackstrcat("tstack.bin",path);
		if(filetstack = fopen(path,"r")) {
			temp = fread(&tstack_energy, sizeof(int), 6*6*6*6, filetstack);
			if(temp != count)
				printf("Only wrote %d elements of tstack\n",temp);
			fclose(filetstack);
		} else {
			printf("Could not open terminal stack binary (tstack.bin) for reading!\n");
			exit(1);
		}
		
		stackstrcat("tstacki.bin",path);
		if(filetstacki = fopen(path,"r")) {  
			temp = fread(&tstacki_energy, sizeof(int), 6*6*6*6, filetstacki);
			if(temp != count)
				printf("Only wrote %d elements of tstacki\n",temp);
			fclose(filetstacki);
		} else {
			printf("Could not open internal loop terminal stack binary (tstacki.bin) for reading!\n");
			exit(1);
		}

		
		//TODO: understand this.
		for( i = 1; i < 5; i++ ){
			for( j = 1; j < 5; j++){
				tstack_energy[i][j][5][5] = 0;
				tstack_energy[5][5][i][j] = 0;
			}
		}	
	} else {
		printf("Thermodynamic datatable directory environment variable not set.  Please run export BINPATH=\"/path/to/datatables/\"\n");
		exit(1);
	}

	free(path);
}

/*
    fill the energy arrays for the internal loops by opening the files containing the binary memory images
*/
void get_iloop_energies(){
 
	int temp, count;
	
	FILE* file1x1, *file2x2, *file2x1;
	
	char * wd; 
	char * path = malloc(strlen(getenv("BINPATH")) + 20);
	
	if(wd = getenv("BINPATH")) {
		stackstrcat("int1x1.bin",path);
		if(file1x1 = fopen(path,"r")) {   
			count = 6*6*6*6*6*6;	 //1x1 is 6*6*6*6*6*6 ints
			temp = fread(&iloop11, sizeof(int), count, file1x1);
			if(temp != count)
				printf("Only wrote %d elements of 1x1\n",temp);
			fclose(file1x1);
		} else {
			printf("Could not open internal 1x1 loop binary (int1x1.bin) for reading!\n");
			exit(1);
		}
		
		stackstrcat("int2x2.bin",path);
		if(file2x2 = fopen(path,"r")) {
			count = 6*6*6*6*6*6*6*6; //2x2 is 6*6*6*6*6*6*6*6 ints
			temp = fread(&iloop22, sizeof(int), count, file2x2);
			if(temp != count)
				printf("Only wrote %d elements of 2x2\n",temp);
			fclose(file2x2);
		} else {
			printf("Could not open internal 2x2 loop binary (int2x2.bin) for reading!\n");
			exit(1);
		}
		
		stackstrcat("int2x1.bin",path);
		if(file2x1 = fopen(path,"r")) {
			count = 6*6*6*6*6*6*6;   //2x1 is 6*6*6*6*6*6*6 ints
			temp = fread(&iloop21, sizeof(int), count, file2x1);
			if(temp != count)
				printf("Only wrote %d elements of 2x1\n",temp);
			fclose(file2x1);
		} else {
			printf("Could not open internal 2x1 loop binary (int2x1.bin) for reading!\n");
			exit(1);
		}
	} else {
		printf("Thermodynamic datatable location environment variable not set.  Please cd to this directory and run ./binset, or run export BINPATH=\"/path/to/datatables/\"\n");
		exit(1);
	}
}

// calculate base-dependent secondary structure removal penalty (added JH August 2011)

extern int mftbool;

int affinePent (int a, int b ) {
	int k;
	int p=0;
	for(k=a;k<=b;k++) {
		p += nucenergy[t[k-1]];
	}
	return -p*mftbool;
}

int retnucenergy ( int a ) {
	return -nucenergy[t[a-1]]*mftbool;
}

// returns the match score for s[i-1] and t[j-1]
//	5' s[i-2] s[i-1] 3'
//	3' t[j-2] t[j-1] 5'
int match(const int i, const int j){
	
	//printf( " match(%i,%i)\n %c %c \n %c %c \n %i \n", i,j, num_to_base(s[i-1]), num_to_base(s[i]), num_to_base(t[j-1]), num_to_base(t[j]),  stack_energy[s[i-2]+1][t[j-2]+1][s[i-1]+1][t[j-1]+1]);//stack_energy[s[i-1]][t[j-1]][s[i]][t[j]]);
	return -stack_energy[s[i-2]+1][t[j-2]+1][s[i-1]+1][t[j-1]+1]; //nearest neighbor base-pairing score
	
}

//custom format is
// 	5' s[i] s[i'] 3'
// 	3' t[j] t[j'] 5'
int matchAlt(const int i, const int i_prime, const int j, const int j_prime){
	
	return -stack_energy[s[i-1]+1][t[j-1]+1][s[i_prime-1]+1][t[j_prime-1]+1]; //nearest neighbor base-pairing score
	
}

//returns the terminal mismatch score for s[i] and t[j]
// used for the opening of an internal or bulge loop
int tmatch(const int i, const int j){
	
	return -tstacki_energy[s[i-2]+1][t[j-2]+1][s[i-1]+1][t[j-1]+1];
	
}

int tmatchAlt(const int i, const int i_prime, const int j, const int j_prime){
	
	return -tstacki_energy[s[i-1]+1][t[j-1]+1][s[i_prime-1]+1][t[j_prime-1]+1];
	
}


//returns the terminal mismatch score for s[i] and t[j]
// used for the closing of an internal or bulge loop
int tmatch2(const int i, const int j){
	//   printf( " tmatch2(%d,%d)\n %c %c \n %c %c \n %g \n", i,j, num_to_base(t[j-1]), num_to_base(s[i-2]), num_to_base(s[i-1]), num_to_base(t[j-2]), tstack_energy[t[j-1]][s[i-1]][t[j-2]][s[i-2]] );
	return -tstacki_energy[t[j-1]+1][s[i-1]+1][t[j-2]+1][s[i-2]+1];
	
}

int tmatch2Alt(const int i, const int i_prime, const int j, const int j_prime){
	
	return -tstacki_energy[t[j_prime-1]+1][s[i_prime-1]+1][t[j-1]+1][s[i-1]+1];
	
}

//TODO: check this.
//data arrangement for 1x1 loop tables iloop11[a][b][c][d][e][f]:
//abc
//def
int i11(const int c, const int f){
	
	int a = c-2;
	int d = f-2;
	//the +1 is needed because Mathew's encoding is slightly different
	return -iloop11[s[a-1]+1][s[a]+1][s[a+1]+1][t[d-1]+1][t[d]+1][t[d+1]+1];// + au_penalty(c,f) + au_penalty(a,d);
	
}

//key iloop22[a][b][c][d][j][l][k][m] =
//a j l b
//c k m d
int i22(const int b, const int d){
	
	int a = b - 3;
	int c = d - 3;
	// printf( "(%d,%d) chars (%d %d) = %d\n", a,c,s[a-1], t[c-1], -iloop22[s[a-1]+1][s[a+2]+1][t[c-1]+1][t[c+2]+1][s[a]+1][s[a+1]+1][t[c]+1][t[c+1]+1]);
	return -iloop22[s[a-1]+1][s[a+2]+1][t[c-1]+1][t[c+2]+1][s[a]+1][s[a+1]+1][t[c]+1][t[c+1]+1];// + au_penalty(b,d)+ au_penalty(a+1,c+1);
	
}

//key iloop21[a][b][c][d][e][f][g] =
//a c f
//b d e g
int i12(const int f, const int g){
	
	return -iloop21[s[f-3]+1][t[g-4]+1][s[f-2]+1][t[g-3]+1][t[g-2]+1][s[f-1]+1][t[g-1]+1];
	
}


//21 is a flipped 12
//g e d b
//f c   a
int i21(const int b, const int a){
	
	return -iloop21[ t[a-1]+1 ][ s[b-1]+1 ][ t[a-2]+1 ][ s[b-2]+1 ][ s[b-3]+1 ][ t[a-3]+1 ][ s[b-4]+1 ];
	
}


// return AU_PENALTY if s[i] and t[j] are AU or UA
int au_penalty(const int i, const int j) {
	//TODO: this works if U and *anything else* are (i-1,j-1). A=0 and G=2 and U=3
	//current logic: if (at least one is U) AND (the other isnt also U) ? AU_PENALTY : 0
	return (((s[i-1]-3)*(t[j-1]-3) == 0) && ( s[i-1] + t[j-1] - 6 != 0) ? AU_PENALTY : 0 );

	//                  one is U(3)             given that one is U(3), the other is A(0) or G(2)
	//return (((s[i-1]-3)*(t[j-1]-3) == 0) && ( (s[i-1] + t[j-1] == 3) || (s[i-1] + t[j-1] == 5) ) ? AU_PENALTY : 0 ); // 09wkj ADDED 4/27/10
	
}


//return tstack for initial terminal mismatch
// in data tables, 0=X, 1=A, 2=C, 3=G, 4=T/U, and 5=intermolecular linker, so we add 1
int startmatch(const int i, const int j){
#if NEW_RULES
	return -tstack_energy[t[j-1]+1][s[i-1]+1][t[j-2]+1][s[i-2]+1] + au_penalty(i,j); // 09wkj ADDED 4/28/10
	//return -tstack_energy[s[i-1]+1][t[j-1]+1][s[i-2]+1][t[j-2]+1] + au_penalty(i,j); // 09wkj ADDED 4/27/10 (directions wrong)
	
	// THIS WAS ORIGINAL
	//return -tstack_energy[s[i-2]+1][t[j-2]+1][s[i-1]+1][t[j-1]+1] + au_penalty(i,j);
#else
	//old terminal mismatch rule (sum of 2 dangles)
	return -tstack_energy[s[i-2]+1][5][s[i-1]+1][t[j-1]+1] - tstack_energy[5][t[j-2]+1][s[i-1]+1][t[j-1]+1] + au_penalty(i,j);
#endif
}

//return tstack for final terminal mismatch
// in data tables, 0=X, 1=A, 2=C, 3=G, 4=T/U, and 5=intermolecular linker, so we add 1
int finalmatch(const int i, const int j){
#if NEW_RULES
	return -tstack_energy[s[i-2]+1][t[j-2]+1][s[i-1]+1][t[j-1]+1] + au_penalty(i-1,j-1);	
#else
	//old terminal mismatch rule (sum of 2 dangles)
	return -tstack_energy[s[i-2]+1][t[j-2]+1][5][t[j-1]+1] - tstack_energy[s[i-2]+1][t[j-2]+1][s[i-1]+1][5] + au_penalty(i-1,j-1);
//  return -tstack_energy[s[i-2]+1][t[j-2]+1][0][t[j-1]+1] - tstack_energy[s[i-2]+1][t[j-2]+1][s[i-1]+1][0] + au_penalty(i-1,j-1);
#endif
}

/*
  calculates the score for an asymmetric internal loop ending at (i,j), with the s length of n1 and t length of n2
*/
int asym_iloop(const int i, const int j, const int n1, const int n2){
	
	int sum = n1+n2;
	int initiation, mismatch; //, au_closure;
	
	switch(sum) {
	case(4):
		initiation = -ASYM_ILOOP_4;
		break;
	case(5):
		initiation = -ASYM_ILOOP_5;
		break;
	case(6):
		initiation = -ASYM_ILOOP_6;
		break;
	default:
		initiation = -(ASYM_ILOOP_6 + floor( (ASYM_ILOOP_COEF * RT * log((float) sum/6.0) * 100.0)+.5));
		break;
	}
	
	mismatch = 0;
	
	if( (n1 == 1) || (n2 == 1) ){
		
		//check for UU and GA on the i,j side
		//recall A=0, C=1, G=2, U=3
		if( (s[i-1]  == U) && (t[j-1] == U ) ){
			mismatch -= ASYM_ILOOP_BONUS_UU;
		}
		else if( (s[i-1] + t[j-1] == G) && (s[i-1]*t[j-1] == A) ) {
			mismatch -= ASYM_ILOOP_BONUS_GA;
		}
		
		//check for UU and AG on the i-(n1-1), j-(n2-1) side 
		if( (s[i - (n1 - 1) -1] == U) && (t[j - (n2 - 1) -1] == U) ){
			mismatch -= ASYM_ILOOP_BONUS_UU;
		}
		else if( (s[i - (n1 -1) -1] + t[j - (n2 -1) -1] == G) && (s[i- (n1 - 1) -1]*t[j- (n2 - 1) -1] == A) ){
			mismatch -= ASYM_ILOOP_BONUS_GA;
		}
	}
	
	
	
	// au_closure = ( (s[i-1]-3)*(t[j-1]-3) == 0 ) ? -ASYM_ILOOP_PENALTY_AU - AU_PENALTY : 0;
	// au_closure += ( (s[i-n1-2]-3)*(t[j-n2-2]-3) == 0 ) ? -ASYM_ILOOP_PENALTY_AU - AU_PENALTY: 0;
	
	
	return initiation + mismatch - ASYM_ILOOP_ASYM*abs(n1-n2);// + au_closure;
}

void resetpairs(){
	
	int i,j;
	
	for(i = 0; i < len_s; i++ ){
		for( j = 0; j < len_t; j++ ){
			pairs[i][j] = unpaired;
		}
	}
	
}

/*
  do the fill of the alignment matrices for a large internal loop.
  
  pass it the last two bases of the internal loop and the length of the loop.
  
  it fills everything except the last two matching bases
*/
void bigfill(const int i, const int j, int n1, int n2){
	
	n1--;
	n2--;
	
	int k;
	
	len_align++;
	aligned_s[len_align] = num_to_base(s[i-n1-2]);
	aligned_t[len_align] = num_to_base(t[j-n2-2]);
	
	if(n1 > n2) {
		for(k = 1; k <= n2; k++){
			
			aligned_t[len_align + n1+1 - k] = num_to_base( t[j-k-1] );
			
		}
		for( k = n2+1; k <= n1; k++){
			
			aligned_t[len_align + n1+1 - k] = '*';
			
		}
		for(k = 1; k <= n1; k++){
			
			aligned_s[len_align + n1+1 - k] = num_to_base( s[i-k-1] );
			
		}
	}
	else {
		for(k = 1; k <= n1; k++){
			aligned_s[len_align + n2 +1 -k] = num_to_base( s[i-k-1] );
		}
		for( k = n1+1; k <= n2; k++){
			aligned_s[len_align + n2+1-k] = '*';
		}
		for(k = 1; k <= n2; k++){
			aligned_t[len_align + n2+1- k] = num_to_base( t[j-k-1] );
		}
	}
	len_align += MAX(n1,n2);
}

/*
  Asymmetric internal loop calculations. Corresponds to the A matrix in Fig 2, Bindigo paper
  calculate the new lu[i][j].  pass it lu[i][j] as x.
*/
void lu_policy(const int i, const int j, struct asymloopdata *x){
	
	// 3x2 internal loop
	int score1 = pxp[i-3][j-2] + tmatch(i-2,j-1) + asym_iloop(i,j,3,2) + affinePent(j-1,j);
	
	if( (x->up + x->down ) > 0 ){
		// adding a single base from s to the internal loop
		int score2 = lu[i-1][j].score - 
			asym_iloop(i-1,j,lu[i-1][j].up, lu[i-1][j].down) + 
			asym_iloop(i,j,lu[i-1][j].up+1, lu[i-1][j].down); 
		
		// adding a single base from s to the internal loop
		int score3 = ld[i-1][j].score - 
			asym_iloop(i-1,j,ld[i-1][j].up, ld[i-1][j].down) + 
			asym_iloop(i,j,ld[i-1][j].up+1, ld[i-1][j].down); 

		x->score = MAX( score1, MAX( score2, score3));
		
		if( x->score == score1 ){
			x->up = 3;
			x->down = 2;
		} else if( x->score == score2 ){
			x->up = lu[i-1][j].up+1;
			x->down = lu[i-1][j].down;
		} else if( x->score == score3 ){
			x->up = ld[i-1][j].up+1;
			x->down = ld[i-1][j].down;
		}
		
	} else {
		x->score = score1;
		x->up = 3;
		x->down = 2;
	}
	
}


/*
  calculate the new ld[i][j]. return new ld[i][j]
*/
void ld_policy(const int i, const int j, struct asymloopdata *x){
	
	// 2x3 internal loop
	int score1 = pxp[i-2][j-3] + tmatch(i-1,j-2) + asym_iloop(i,j,2,3) + affinePent(j-2,j);
	
	if( (x->up + x->down) > 0 ){

		// adding a single base from t to the internal loop
		int score2 = lu[i][j-1].score - 
			asym_iloop(i, j-1, lu[i][j-1].up, lu[i][j-1].down) + 
			asym_iloop(i, j, lu[i][j-1].up, lu[i][j-1].down+1) + 
			retnucenergy(j);

		// adding a single base from t to the internal loop
		int score3 = ld[i][j-1].score - 
			asym_iloop(i, j-1, ld[i][j-1].up, ld[i][j-1].down) + 
			asym_iloop(i, j, ld[i][j-1].up, ld[i][j-1].down+1) + 
			retnucenergy(j);

		x->score = MAX( score1, MAX( score2, score3 ) );

		if( x->score == score1 ){
			x->up = 2;
			x->down = 3;
		} else if( x->score == score2 ){
			x->up = lu[i][j-1].up;
			x->down = lu[i][j-1].down+1;
		} else if( x->score == score3 ){
			x->up = ld[i][j-1].up;
			x->down = ld[i][j-1].down+1;
		}

	} else {
		x->score = score1;
		x->up = 2;
		x->down = 3;
	}
}

void l1xn_policy(const int i, const int j, struct asymloopdata *x){
	
	// 1x3 internal loop
	int score1 = pxp[i-1][j-3] + tmatch(i,j-2) + asym_iloop(i,j,1,3) + affinePent(j-2,j);
	
	// adding a base from t to the 1xn internal loop
	int score2 = l1xn[i][j-1].score - 
		asym_iloop(i, j-1, l1xn[i][j-1].up, l1xn[i][j-1].down) + 
		asym_iloop(i, j,   l1xn[i][j-1].up, l1xn[i][j-1].down+1) + 
		retnucenergy(j);

	/*if(i == 26 && j == 29) {
		printf("loop score: %d %d\n",score1,score2);
		printf("sub: %d\n",l1xn[i][j-1].score); 
		printf("#top: %d\n",l1xn[i][j-1].up);
		printf("#bot: %d\n",l1xn[i][j-1].down); 
		printf("first: %d\n",asym_iloop(i, j-1, l1xn[i][j-1].up, l1xn[i][j-1].down)); 
		printf("second: %d\n",asym_iloop(i, j,   l1xn[i][j-1].up, l1xn[i][j-1].down+1)); 

	}*/
	
    x->score = MAX( score1, score2 );

    if( x->score == score1 ){
		x->up = 1;
		x->down = 3;
    } else if( x->score == score2 ){
		x->up = l1xn[i][j-1].up;
		x->down = l1xn[i][j-1].down+1;
    }
}

void lnx1_policy(const int i, const int j, struct asymloopdata *x){
	
	// a 3x1 internal loop
	int score1 = pxp[i-3][j-1] + tmatch(i-2,j) + asym_iloop(i,j,3,1) + retnucenergy(j);
	
	// adding a base from s to the nx1 internal loop
	int score2 = lnx1[i-1][j].score - 
		asym_iloop(i-1, j, lnx1[i-1][j].up,   lnx1[i-1][j].down) + 
		asym_iloop(i,   j, lnx1[i-1][j].up+1, lnx1[i-1][j].down);
	    
    x->score = MAX( score1, score2 );

    if( x->score == score1 ){
		x->up = 3;
		x->down = 1;
    } else if( x->score == score2 ){
		x->up = lnx1[i-1][j].up+1;
		x->down = lnx1[i-1][j].down;
    }
}

void print_matrix(const int matrix_num, const int h, const int w, FILE* output ){
	
	//printf( "printing output for matrix %d\n", matrix_num  );
	
	int i,j;
	fprintf(output, "\n----------------------- %d ---------------------------\n", 
			matrix_num);
	
	if(matrix_num == LU_M){
		for(i = DUMMY_LENGTH; i <= h; i++ ){
			fprintf(output, "\n");
			for( j = DUMMY_LENGTH; j <= w; j++ ){
				
				fprintf(output, "%d,", lu[i][j].score);
				
			}
		}
	}
	else if(matrix_num == LD_M){
		for(i = DUMMY_LENGTH; i <= h; i++ ){
			fprintf(output, "\n");
			for( j = DUMMY_LENGTH; j <= w; j++ ){
				
				fprintf(output, "%d,", ld[i][j].score);
				
			}
		}
	}
	else {
		for(i = DUMMY_LENGTH; i <= h; i++ ){
			fprintf(output, "\n");
			for( j = DUMMY_LENGTH; j <= w; j++ ){
				
				fprintf(output, "%d,", the_matrices[matrix_num][i][j]);
				
			}
		}
	}
	fprintf(output,"\n");
	
	fflush(output);
}

/*
  go through the pair matrix and reconstruct the secondary structure by creating the structure array.
  A nearest neighbor pair is signified by (0,0). A 1x0 is (1,0), etc.
*/
void pairoutput(void){
	
	int i,j;
	
	int started = 0;
	
	//what is each s_i paired to?
	for( i = 0; i < len_s; i++){
		for( j = 0; j < len_t; j++){
			
			if( pairs[i][j] == paired ){
				
				if(!started){
					structures[numstructures].s_start = i;
					structures[numstructures].t_start = j;
					started = 1;
				}
				else {
					structures[numstructures].s_end = i;
					structures[numstructures].t_end = j;
					numstructures++;
					
					structures[numstructures].s_start = i;
					structures[numstructures].t_start = j;
				}
				
			}
			
		}
	}
	
}

/*
   Walks though s vs. aligned_s and t vs. aligned_t
   comparing characters. Returns 1 if there is 
   an inconsistency.
*/
int alignError(){ 


	int o; // original string index
	int a; // aligned string index


	//    S
	a = 0;
	o = 0;
   	while(aligned_s[a] == ' '){
		++a;
	}
	while(num_to_base(s[o]) != 'A' && num_to_base(s[o]) != 'C' && num_to_base(s[o]) != 'G' && num_to_base(s[o]) != 'U'){
		++o;
	}

	while(o < len_s + DUMMY_LENGTH ){
		while(aligned_s[a] == '*'){
			++a;
		}
		if(num_to_base(s[o]) != toupper(aligned_s[a])){
			printf("S -> original[%i]: '%c' != aligned[%i]: '%c'\n", o, num_to_base(s[o]), a, toupper(aligned_s[a])); 
			return 1;
		}
		++a;
		++o;
	}


	//    T
	a = 0;
	o = 0;
	while(aligned_t[a] == ' '){
		++a;
	}
	while(num_to_base(t[o]) != 'A' && num_to_base(t[o]) != 'C' && num_to_base(t[o]) != 'G' && num_to_base(t[o]) != 'U'){
		++o;
	}
	while(o < len_t + DUMMY_LENGTH){
		while(aligned_t[a] == '*'){
			++a;
		}
		if(num_to_base(t[o]) != toupper(aligned_t[a])){
			printf("T -> original[%i]: %c != aligned[%i]: %c\n", o, num_to_base(t[o]), a, toupper(aligned_t[a])); 
			return 1;
		}
		++a;
		++o;
	}
	return 0;
}
