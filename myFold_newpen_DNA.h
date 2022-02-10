/*
  myFold.h
  Nathan Oken Hodas
  Aug 18

  The rewritten version to account for the new PxP, MBU, MBD model.
  Also includes the necessary elements for the loop up and loop down procedures.
*/

/** 
	@file myFold.h, written by Nathan Oken Hodas. 
	-This is the rewritten version to account for the new PxP, MBU, MBD model.
	-It also includes the necessary elements for the loop up and loop down procedures.
	-Added in October 2009 is an average secondary structure penalty, which approximates a 
	penalty for each additional paired (or unpaired) base in the longer string. 
	NOTE: all units are in 1/100ths kcal/mol 
*/


#ifndef MYFOLD_H
#define MYFOLD_H

#define MAX(a,b) ((a) < (b) ? (b) : (a))

/*
  Constants
*/


#define MEAN_PENALTY 30

//PENALTIES FOR SHORT BULGE LOOPS

/** Penalty for a single bulge loop gap */
#define SINGLE_GAP_PENALTY 290 //if you just have 1 gap.  base is "popped out," hence higher.

/** Penalty for n-length bulge loop gaps (n<=6): MULTI_GAP_AFFINE + MULTI_GAP_PENALTY*n */
#define MULTI_GAP_AFFINE 182

/** Penalty for n-length bulge loop gaps (n<=6): MULTI_GAP_AFFINE + MULTI_GAP_PENALTY*n */
#define MULTI_GAP_PENALTY 23


//PENALTIES FOR LONG BULGE LOOPS
//( The transition occurs at 6 )

/** Penalty for n-length bulge loop gaps (n>6): MULTI_GAP_AFFINE_2 + MULTI_GAP_PENALTY_2*(n-5) */
#define MULTI_GAP_AFFINE_2 268

/** Penalty for n-length bulge loop gaps (n>6): MULTI_GAP_AFFINE_2 + MULTI_GAP_PENALTY_2*(n-5) */
#define MULTI_GAP_PENALTY_2 10

//penalties for mismatches (interior loops)

/** (interior loop) affine penalty for mismatch */
#define MISMATCH_AFFINE 9999//160 //affine penalty for mismatch

/** (interior loop) penalty for mismatch */
#define MISMATCH_PENALTY 9999//10

/** The intermolecular initiation energy. It is applied once */
#define INITIATION_PENALTY 100 //g_0 removal penalty is handled separately.

/** The AU and GU helix penalty. Applied at the start and the end of a helix */
#define AU_PENALTY 0 //SL says differently?

//CONSTANTS NEEDED TO CALCULATE ASYMMETRIC INTERNAL LOOP PENALTIES

/** internal loop initiation penalty for n1+n2 = 4 */
#define ASYM_ILOOP_4 310

/** internal loop initiation penalty for n1+n2 = 5 */
#define ASYM_ILOOP_5 350

/** internal loop initiation penalty for n1+n2 = 6 */
#define ASYM_ILOOP_6 390

/** used to extrapolate internal loop initiation penalty for n1+n2>6 (not used) */
#define ASYM_ILOOP_COEF 1.75

/** used to extrapolate internal loop initiation penalty for n1+n2>6. */
/* also anywhere we need a Boltzmann factor */
#define RT 0.00198*310.15

/** the asymmetric penalty (internal loop) */
#define ASYM_ILOOP_ASYM 0 //(nonzero according to SL)

//these lines are shut off, as the parameters don't exist for DNA?
/** special bonuses for opening/closing mismatches (internal loop) */
#define ASYM_ILOOP_BONUS_GA 0

/** special bonuses for opening/closing mismatches (internal loop) */
#define ASYM_ILOOP_BONUS_UU 0

/** AU, GU get different penalties when closing an internal loop */
#define ASYM_ILOOP_PENALTY_AU 0

/** maximum size for one side of an asymmetric internal loop */
#define MAX_ASYM_ILOOP_SIZE 5


//BASE REMOVAL PENALTIES

int nucenergy[6];
int g0;


//MATCHING BONUSES

/** mismatch penalty accounted for by interior loop penalties */
#define MISMATCH -9999



/** 
	indices for the a,a_prime,b,c,d,e matrices. or in the terminology used by the Bindigo paper,
    the M, B, b, B2, b2, F, A, and a matrices of Figure 2.
*/
enum matrixNums { 
	PXP_M=0, ///< the M matrix
	BU_M=1,  ///< the B matrix
	BD_M=2,  ///< the b matrix
	BU2_M=3, ///< the B2 matrix
	BD2_M=4, ///< the b2 matrix
	F_M=5,   ///< the F matrix
	LU_M=6,  ///< the A matrix
	LD_M=7   ///< the a matrix
};

/** the number of matrices stored in the_matrices, which does not include lu and ld 
	because they are of a different data type. */
#define NUM_MATRICES 6

/** negative infinity, used to initialize the dynamic programming tables */
#define NEG_INFINITY -9999

// PARAMETERS FOR SCANNING IN DATA FILES

/** There are this many lines in an MFOLD data file */
#define NUM_LINES 760

/** A slight over estimate of the number of characters in each line */
#define LINE_LENGTH 1100

/** We pad the beginning of the sequences with this many dummy characters. 
	useful for bulge loop calculations at the beginning of the sequence */
#define DUMMY_LENGTH 20
                       
/** If <em>t<\em> is greater than PENALTY_THRESHOLD, then we can expect that it has 
	a secondary structure prior to binding with <em>s<\em>. Binding to <em>t<\em> would 
	therefore require an additional energy of 2.86 per base (in expectation) */
#define PENALTY_THRESHOLD 16



// DECLARE FUNCTIONS

/** 
	Extracts energies from a string that is the contents of an MFOLD data file. Not implemented.
	@param MFOLD_dat_file MFold file to read in
	@param energy matrix to populate.
*/
void parse(char MFOLD_dat_file[], int energy[5][5][5][5]);

/** fill the energy matrices from data files */
void get_stack_energies(); 

/** get the internal loop energies prepared from datreader.cc */
void get_iloop_energies(); 

/** 
	Function to perform a traceback starting at i,j and starting in the indicated matrix
	@param i starting index into the traceback matrix
	@param j starting index into the traceback matrix
	@param matrix_index Matrix to start the traceback. Appropriate values can be selected from \enum matrixNums
*/
int alignment(const int i, const int j, const int matrix_index, const int firstrun);

int alignment2(const int i, const int j, const int matrix_index, FILE* trace);

/** 
	Fills the beginning of the aligned sequence with the appropriate characters.
	@param i index into the aligned_s array
	@param j index into the aligned_t array
*/
void fill_prefix(const int i, const int j); 

/**
   Returns the match score for s[i-1] and t[j-1] (nearest neighbor stacking).
   <tt>
   5'  s[i-2] s[i-1]  3' <br>
   3'  t[j-2] t[j-1]  5' </tt>
	@param i 
	@param j 
 */
int match(const int i, const int j); 

/**
   Same as match, but uses custom sequence positions:
   <tt>
   5'  s[i] s[i']  3' <br>
   3'  t[j] t[j']  5' </tt>

   @param i
   @param i_prime
   @param j
   @param j_prime
 */
int matchAlt(const int i, const int i_prime, const int j, const int j_prime); 

/**
   Returns score of nearest neighbor terminal mismatch. It is used for the opening of an 
   internal or bulge loop.
   @param i
   @param j
 */
int tmatch(const int i, const int j); 

/**
   Returns score of nearest neighbor terminal mismatch (same as tmatch), but uses custom sequence positions.
   @param i
   @param i_prime
   @param j
   @param j_prime
 */
int tmatchAlt(const int i, const int i_prime, const int j, const int j_prime); 

/**
   Returns score of nearest neighbor reversed terminal mismatch (i.e. mismatch, then pairing).
   Used for the closing an internal or bulge loop.
   @param i
   @param j
 */
int tmatch2(const int i, const int j); 

/**
   Returns score of nearest neighbor reversed terminal mismatch (i.e. mismatch, then pairing). The same as tmatch2, but it uses custum sequence positions.
   @param i
   @param i_prime
   @param j
   @param j_prime
 */
int tmatch2Alt(const int i, const int i_prime, const int j, const int j_prime); 

/**
   Returns score for initial terminal mismatch.
   @param i
   @param j
 */
int startmatch(const int i, const int j); 

/**
   Returns score for final terminal mismatch
   @param i
   @param j
*/
int endmatch(const int i, const int j); 

//INTERNAL LOOPS
// TODO: check this vs. 22 and 21
/** data arrangement for 1x1 loop tables. <br> 
    key: <br> <tt> iloop11[a][b][c][d][e][f]: <br> a b c <br> d e f */
int i11(const int a, const int d);

/** data arrangement for 2x2 loop tables. <br> 
    key: <br> <tt> iloop22[a][b][c][d][j][l][k][m] = <br> a j l b <br> c k m d 
*/
int i22(const int a, const int c);

/** data arrangement for 2x1 loop tables. <br> 
    key: <br> <tt> iloop21[a][b][c][d][e][f][g] = <br> a e f b <br> c g   d */
int i21(const int a, const int c);

/** data arrangement for 1x2 loop tables. (1x2 is a flipped 2x1) <br>
    key: <br> <tt> iloop12[a][b][c][d][e][f][g] = <br> a  g  b <br> c e f d */
int i12(const int a, const int c);

/** Calculates the free energy penalty for starting or closing a helix with an AU or GU pair.
	@param i index of the current <em>s<\em> base
	@param j index of the current <em>t<\em> base
	@return AU_PENALTY if needed, otherwise 0
*/
int au_penalty(const int i, const int j);

/** 
	A function to convert an integer to a base according to following assigment: <br>
	'A'=0 'C'=1 'G'=2 'U'=3 'X'=4 '*'=5 '_'=6 <br>
	Note: no bounds checking is done, so numbers not in the range [0,6] will produce 
	unexpected results
	@param i int in the range [0,6] to be converted into a character base
	@return a single character from  {'A' 'C' 'G' 'U' 'X' '*' or '_'}
*/
char num_to_base(const int i);

/** 
	Fill the matrices with their values. The rules correspond to Figure 2 in the Bindigo paper.
*/ 
void fill_arrays(void);

//void print_matrix( int m, int h, int w); //prints a matrix to stdout
//void print_matrix(  int matrix_num, int h, int w );

/** 
	Initialization: creating the matrix, loading the data, etc.
 */
void setup(void);




//GLOBALS -- poor form, but simplest for first attempt

/** length of the <em>s<\em> string */
int len_s;

/** length of the <em>t<\em> string */
int len_t;

/** length of the <em>tt<\em> string */
int len_tt;

/** place-holder for the traceback */
int len_align; 

/** where we enter traceback **/
int f_i;
int f_j;

/** free energy landscape dgJ **/
int * dgJ;

/** The actual matrices. A convenient way of encapsulating the matrices -- an array of matrices.
	The i,j entry of the "a" matrix would be:  the_matrices[A_M][i][j] */
int **the_matrices[NUM_MATRICES];

/** Shorthand for the B matrix stored in the_matrices[BU_M]. */
int **bu;

/** Shorthand for the b matrix stored in the_matrices[BD_M]. */
int **bd; 

/** Shorthand for the B2 matrix stored in the_matrices[BU2_M]. */
int **bu2;

/** Shorthand for the b2 matrix stored in the_matrices[BD2_M]. */
int **bd2;

/** Shorthand for the F matrix stored in the_matrices[F_M]. */
int **f;

/** Shorthand for the M matrix stored in the_matrices[PXP_M]. */
int **pxp;

/** Throughout the program, bases are referenced by indexing into the code array. The codes enum 
	allows the use of the appropriate constant variables throughout the program. */
enum codes {
	A=0, ///< integer value of 'A'
	C=1, ///< integer value of 'C'
	G=2, ///< integer value of 'G'
	U=3, ///< integer value of 'U'
	X=4  ///< integer value of 'X'
};
// A=0, C=1, G=2, U=3, X=4, *=5 (X is non-pairing misc. base and * indicates a gap)

/** The code matrix allows fast decoding. Reverses the directionality of the codes enum. */
const int code[7] = {'A','C','G','U','X','*','_'};

/** An encoded sequence (input) */
int *s;

/** An encoded sequence (input) */
int *t;

/** Full t sequence */
int *tt;

/** A final encoded, aligned sequence (output) */
char *aligned_s;

/** A final encoded, aligned sequence (output) */
char *aligned_t; 

/** Free-energy array for nearest neighbor stacking */
int stack_energy[6][6][6][6];

/** Free-energy array for terminal mismatch nearest neighbor stacking.
	Used for opening and closing internal loops. */
int tstacki_energy[6][6][6][6];

/** Free-energy array for terminal mismatch nearest neighbor stacking.
	Used for initial and final terminal mismatch. */
int tstack_energy[6][6][6][6]; 

/** Free-energy array for 2x2 internal loops. */
int iloop22[6][6][6][6][6][6][6][6];

/** Free-energy array for 2x1 internal loops. */
int iloop21[6][6][6][6][6][6][6];

/** Free-energy array for 1x1 internal loops. */
int iloop11[6][6][6][6][6][6];

/**
   Internal asymmetric loop data. Each matrix entry is a special data type, which 
   includes the score, the number of s bases in the loop, and the number of t bases 
   in the loop.
*/
struct asymloopdata {
	int score; ///< A/a.dG
	int up;   ///< number of s bases in the loop
	int down; ///< number of t bases in the loop
};

//TODO: more specific comment than "loop stuff"
//the arrays to hold the loop stuff
/** The A matrix from Figure 2. */
struct asymloopdata **lu;

/** The a matrix from Figure 2. */
struct asymloopdata **ld;

/** The K matrix from Figure 2. */
struct asymloopdata **lnx1;

/** The k matrix from Figure 2. */
struct asymloopdata **l1xn;

/** returns entropic for generic n1 x n2 internal loop */
int asym_iloop(const int i, const int j, const int n1, const int n2); 

/** a function to perform the asymmetric loop "policy" on lu  */
void lu_policy(const int i, const int j, struct asymloopdata *x);

/** a function to perform the asymmetirc loop "policy" on ld */
void ld_policy(const int i, const int j, struct asymloopdata *x);

/** a function to perform the asymmetric loop "policy" on 1xn loops */
void l1xn_policy(const int i, const int j, struct asymloopdata *x);

/** a function to perform the asymmetric loop "policy" on nx1 loops */
void lnx1_policy(const int i, const int j, struct asymloopdata *x);

/** If the internal loop is one of the special cases, return a very large penalty.
	We don't want generic score overriding the special score. */
int adj_score(struct asymloopdata *loop);

//int getmax(int i, int j);

/**
  Takes a character string and encodes it into the numerical array
  used in the algorithm. 
  @param dest int array where values are placed
  @param src char array where values are polled from
  @param length s and t should be of size length
  @param direction is a flag which tells whether the output is
  5' to 3' or 3' to 5'. If direction is 1, the output is 5' to 3'. Otherwise it is 3' to 5'.
  The s should be 5' to 3' and the t should be 3' to 5'.
 */
void encode(int* dest, const char* src, const int length, const int direction);

/** An extended fill of the alignment matrices. */
void bigfill( const int i, const int j, const int len_up, const int len_down);

/** 
	Writes the specified matrix to the specified file. 
	Note: NO bounds checks are performed. It is up to the user to ensure that height and 
	width do not exceed the dimensions of the specified matrix.
	@param matrix_num Which matrix to print. Options are: PXP_M, BU_M, BD_M, BU2_M, BD2_M, F_M, LU_M, LD_M
	@param h The height (number of lines) of the matrix to print
	@param w The width (number of columns) of the matrix to print
	@param output File to write to
*/
void print_matrix(  const int matrix_num, const int h, const int w, FILE* output );

/*
  pairing info
*/

/** Integer code to indicate whether or not a base is paired. */
enum pairtype {
	unpaired, ///< The base is unpaired
	paired ///< The base is paired
};

/** Matrix to store s-t pairing information */
int **pairs;

/** Sets the entire pairs matrix to unpaired */
void resetpairs();

/**
  This function creates the structures array by going through the pair matrix.
  The structures array can be used to reconstruct the secondary structure.
  A nearest neighbor pair is signified by (0,0). A 1x0 is (1,0), etc.
*/
void pairoutput(void);

/**
   Data structure used to maintain loop information.
 */
struct loop{
	int s_start; ///< index of the start of the loop in the s string
	int s_end;   ///< index of the end of the loop in the s string
    int t_start; ///< index of the start of the loop in the t string
	int t_end;   ///< index of the end of the loop in the t string
};

/**
   A list of the structures in the pairing.
   The structures array can be used to reconstruct the secondary structure.
   A nearest neighbor pair is signified by (0,0). A 1x0 is (1,0), etc.
 */
struct loop *structures; 

/** The number of structures in the pairing */
int numstructures; 

/** The first two or last two bases paired, depending on which way you are doing the traceback */
int ending_s;

/** The first two or last two bases paired, depending on which way you are doing the traceback */
int ending_t; 

/** Checks that the printed characters match the original strings */
int alignError(void);

int finalmatch(const int i, const int j);
int startmatch(const int i, const int j);

#endif 
