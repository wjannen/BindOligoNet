/*
occupational probabilities of u1 on dummy/real splice sites
original BINDIGO code by Nathan O. Hodas, 2004
modified by Julian M. Hess, 2011-2012
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "myFold_MFT-nz.c"

#define RANDOM 2
#define REAL 1
#define DECOY 0

#define BEST_TB 0
#define SING_SING 1
#define PROD_SUM 2
#define TWO_PAIR 3
#define LOCALIZED 4

#define ALIGNMENTS 0
#define DG 1
#define POCCS 2
#define HOTSPOTS 3

#define SCRG 0
#define BEST 1
#define NONE 2

#define SPAN 0
#define GUAC 1

#define THRESH -410 - g0

#ifdef DNA
#define USAGE "Usage: bindoligo-D <options> [s sequence] [t sequence]\n"
#else
#define USAGE "Usage: bindoligo <options> [s sequence] [t sequence]\n"
#endif

//GLOBAL VARIABLES

struct traceback ** tb;

struct traceback {
	double dg; //dG associated with traceback
	int * pocc; //array containing the binary occupancy at each site
	int djt; //size of alignment footprint
	int gubool; //whether pairing occurs at GU/AC
};

struct dgarray * dga;

struct dgarray {
	int * dgs; //array containing all the dG's for given location
	int last; //last spot in array that was filled
};

int mode; //real, decoy, or random
int tb_mode; //BEST_TB, SING_SING, PROD_SUM, TWO_PAIR
int dgout; //how (or whether) to print dG to stdout.
int pairmode; //how to determine what qualifies for a pairing
int singseq_bool; //whether to look at all splice sites or just one
int djtbool; //whether we're outputting compactness
int lineno; //which sequence to look at in single mode
int the_max_i,the_max_j; //maximum value of traceback
int ran; //in the case of the single occupancy, we only need to compute the traceback once.
int omode; //what we're outputting (pocc, dG, hotspots, alignment)

int dynpen = 0; //whether we apply static or dynamic penalty
int mftbool = 1; //whether we calculate dG for secondary structure removal.  1 = yes, 0 = no

int custbool = 0; //whether user has specified custom parameters

#ifdef DNA
int nucenergy[6] = {-2,24,52,4,0,0}; //DNA
#else
int nucenergy[6] = {-23,48,94,16,0,0}; //RNA
#endif


int g0;
int meanpen;

double chem_pot; //chemical potential *divided by kT (RT)*

//FUNCTIONS

void get_t (char * tstring, int lent ) { //t is 5' -> 3'
	t = (int *) realloc(t,(lent + 1 + DUMMY_LENGTH*2)*sizeof(int));
	encode(t,tstring,lent,0); //0 direction corresponds to 5' -> 3'
}

void get_s (char * sstring, int lens ) { //t is 5' -> 3'
	s = (int *) realloc(s,(lens + 1 + DUMMY_LENGTH*2)*sizeof(int));
	encode(s,sstring,lens,1); //1 direction corresponds to 3' -> 5'
}

int filelength(char * file) {
	int i = 0;
	char buf;
	FILE * seqs;

	if(!(seqs = fopen(file, "r"))) {
		printf("File \"%s\" does not exist.\n",file);
		exit(0);
	}

	do {
		buf = fgetc(seqs);
		if(buf == '\n')
			i++;
	} while (buf != EOF);

	fclose(seqs);

	return i;
}

char ** fencode (char * file, int flength) {
	int i,j;
	FILE * seqs;

	int * lengths = malloc(flength*sizeof(int));

	seqs = fopen(file, "r");

	//get the length of each line in the file
	for(i=0;i<flength;i++) {
		j = 1;
		while(fgetc(seqs) != '\n') {
			++j;
		}
		lengths[i] = j + 1;
	}

	rewind(seqs);

	char ** out = (char **) malloc(flength*sizeof(char *));
	for(i=0;i<flength;i++) {
		out[i] = (char *) malloc(lengths[i]*sizeof(char));
	}

	for(i=0;i<flength;i++) {
		fgets(out[i],lengths[i],seqs);
		out[i][strlen(out[i])-1]='\0'; //chomp the newline
	}
	fclose(seqs);

	free(lengths);

	return out;
}

void write_occ_prob ( ) { //the old BINDIGO paper's method, for comparison
	int i,j;
	FILE * dest = fopen("realoutput_mu=2","a");
	for(i=0;i<len_t+1;i++) {
		double scorez=0;
		for(j=0;j<len_s+1;j++) {
			scorez = MAX(scorez,1.0/(1.0+exp( -((double) f[j+DUMMY_LENGTH][i+DUMMY_LENGTH])/(RT*100.0)-chem_pot)));
		}
		fprintf(dest,"%d %g\n",i,scorez);
	}
	fclose(dest);
}

void align_init(int tfill) {
	int j;

	if(tfill) {
		aligned_s = (char*) malloc((len_s+len_t+DUMMY_LENGTH+1)*3*sizeof(char));
		aligned_t = (char*) malloc((len_s+len_t+DUMMY_LENGTH+1)*3*sizeof(char));

		//change garbage characters to spaces
		for(j = 0; j < (DUMMY_LENGTH+len_s+len_t+1)*3; j++){
			aligned_s[j] = ' ';
			aligned_t[j] = ' ';
		}
	}
}

void align_free() {
	free(aligned_s);
	free(aligned_t);
}

/*
double dpf: Calculates the partition function for a given scheme.
Note that this is actually just calculating the denominator in
occupational probability expressions, and may not actually represent
the true partition function.
Potential TODO: deal with PF overflows
*/

double dpf (int mode) {
	int i,j;
	double pf;

	if (mode == SING_SING) {
		pf = 1;
		for(i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
			for(j = DUMMY_LENGTH + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
				if(f[i][j] < THRESH) //variable threshold?
					continue;

				pf += exp(((double) f[i][j])/(RT*100) - chem_pot);
			}
		}
	} else if((mode == PROD_SUM) || (mode == BEST_TB) || (mode == LOCALIZED)) {
		/*In these cases, we do not explicitly calculate a partition function.
		  For the case of PROD_SUM and LOCALIZED, the entire occupancy calculation is done in p_occ
		  For the case of BEST_TB, we are not accounting for multiple states.
		*/
		pf = 1;
	} else if(mode == TWO_PAIR) {
		pf = 1;
		for(i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) { //single occupancy
			for(j = DUMMY_LENGTH + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
				if(f[i][j] < THRESH) //variable threshold?
					continue;

				pf += exp(((double) f[i][j])/(RT*100) - chem_pot);
			}
		}

		/* The two pairing portion of the partition function
			i,j -> i',j'
			ip,jp -> i'',j''
		*/
		int ip,jp;
		for(j = DUMMY_LENGTH + 1; j <= DUMMY_LENGTH + len_t - len_s + 1; j++) { //two pairings
			for(jp = j + len_s + 1; jp <= len_t + DUMMY_LENGTH + 1; jp++) {
				for(i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
					for(ip = DUMMY_LENGTH + 1; ip <= len_s + DUMMY_LENGTH + 1; ip++) {
						if((f[i][j] < THRESH) || (f[ip][jp] < THRESH)) //variable threshold?
							continue;
						pf += (exp(((double) f[i][j] + (double) f[ip][jp])/(RT*100) - 2*chem_pot));
					}
				}
			}
		}
	}

	return pf;
}

/*
Calculates the alignment for a given i,j.
Stores the alignment (i.e. binary occupancy times Boltzmann factor)
in the provided traceback's array.
(When called, this function shoul be passed the valid address
of an allocated strcut traceback.)
*/
int traceback (int i, int j, struct traceback *x, int tfill) {
	int k;

//	int the_min_i = 0;
//	int the_min_j = 0;

	//prep alignment
	len_align = 0;
	align_init(tfill);

	if(tfill) //this is only necessary if we're keeping track of pairs.
		resetpairs();

	//need to specify where we're coming from
	f_i = i;
	f_j = j;

	//execute traceback
	alignment(i,j,F_M,tfill);

	if(tfill) {
		//calculate first column of F_M
		for(k=DUMMY_LENGTH;k<=len_t+DUMMY_LENGTH+1;k++) {
			if(!isspace(aligned_t[k])) {
//				the_min_j = k;
				break;
			}
		}
		for(k=0;k<=len_t+DUMMY_LENGTH+1;k++) {
			if(!isspace(aligned_s[k])) {
//				the_min_i = k;
				break;
			}
		}

		//fill out aligned_t to the end.
		for (k=j; k<=len_t + DUMMY_LENGTH + 1; k++) {
			aligned_t[len_align+k-j]=num_to_base(t[k-1]);
		}

		//fill out aligned_s to the end.
		for (k=i; k<=len_s + DUMMY_LENGTH; k++) {
			aligned_s[len_align+k-i]=num_to_base(s[k-1]);
		}
	}

	//finally, assign values to given structure.
	x->dg = exp(((double) f[i][j])/(RT*100) - chem_pot); //return Boltzmann factor (E_{i,j})
	x->djt = 0; //possible TODO: accurately report spanning region from pairing matrix

	x->gubool = 1; //possible TODO: include interface for pair forcing on command line

	//free alignment arrays (if necessary)
	if((omode != ALIGNMENTS) && tfill)
		align_free();

	return 0;
}

/*
precompute_traceback: iterates traceback starting from each negative f[i][j].
Also used for finding average binding energy at consensus site,
and for binding footprint calculations.
*/
void precompute_traceback(FILE * out, int run) {
	double G = 0;
	int the_max_fij = -9999;

	for(int i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
		for(int j = DUMMY_LENGTH + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
			if(f[i][j] < THRESH)
				continue;

			if(run)
				traceback(i,j,&tb[i][j],omode == HOTSPOTS ? 0 : 1);

			if(dgout != NONE) {
				if((dgout == SCRG)) { // && tb[i][j].gubool) {
					G += tb[i][j].dg;
				} else if((dgout == BEST)) { // && tb[i][j].gubool) {
					if(f[i][j] > the_max_fij)
						the_max_fij = f[i][j];
				}
			}
		}
	}

	if((out != NULL) && (omode != HOTSPOTS)) {
		if((dgout == SCRG) && (G > 0)) {
			fprintf(out,"%g\n",-RT*log(G));
		} else if(dgout == BEST) {
			fprintf(out,"%g\n",(double) -the_max_fij/100.0);
		}
	}
	return;
}

void precompute_dG() { //like precompute_traceback(), but sets all poccs to 1
	int i,j;
	for(i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
		for(j = DUMMY_LENGTH + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
			if(f[i][j] < THRESH)
				continue;

			traceback(i,j,&tb[i][j],1);
			(tb[i][j].pocc)[j] = 1;
		}
	}
}

/*
Actually calculates the occupational probability across entire sequence t for a given mode.  Note that this actually just calculating the numerator in occupational probability expressions, and then dividing by the partition function.
*/
/*
A brief note: should probably pick the best i, otherwise the calculation across an entire t becomes n^3
*/
double p_occ (int jt, int mode, double pf) {
	double occp = 0;

	if (mode == SING_SING) {
		for(int i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
			for(int j = DUMMY_LENGTH + jt + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
				if(f[i][j] < THRESH) //variable threshold?
					continue;

				/*printf("\n(%d,%d)*****\n",i-DUMMY_LENGTH,j-DUMMY_LENGTH);
				for(k = 0; k <= len_t-1; k++) {
					printf("%d",(tb[i][j].pocc)[k]);
				}
				printf("\n*****\n");*/

				occp += ((tb[i][j].pocc)[jt]*tb[i][j].dg);
			}
		}
	} else if(mode == PROD_SUM) {
		double num, denom;
		for(int j = DUMMY_LENGTH + jt + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
			num = 0;
			denom = 1;
			for(int i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
				if(f[i][j] < THRESH)
					continue;

				num += ((tb[i][j].pocc)[jt]*tb[i][j].dg);
				denom += tb[i][j].dg;
			}
			if(denom != 0) //this check is no longer necessary
				occp += (double) (num/denom);
		}
	} else if(mode == TWO_PAIR) {
		/*
		miscellaneous note: should the sums run from DUMMY_LENGTH + 1?
		*/

		double first, second;
		int ip,jp,ipp,jpp;

		//traceback at each i,j has already been precomputed.

		//first factor
		first = 0;
		for(int i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
			for(int j = DUMMY_LENGTH + jt + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
				if(f[i][j] < THRESH)
					continue;

				first += ((tb[i][j].pocc)[jt]*tb[i][j].dg);
			}
		}

		//second/third factors
		second = 0;
		for(int i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++) {
			for(int j = DUMMY_LENGTH + jt + 1; j <= len_t + DUMMY_LENGTH + 1; j++) {
				if(j - DUMMY_LENGTH >= len_t - len_s)
					//continue;

				//ip,jp -> i',j'
				//ipp,jpp -> i'',j'' (ipp for clarity)
				for(jp = j + len_s; jp <= len_t + DUMMY_LENGTH + 1; jp++) {
					for(ip = DUMMY_LENGTH + 1; ip <= len_s + DUMMY_LENGTH + 1; ip++) { //O(n^4) (!)
						if(f[ip][jp] < THRESH)
							continue; //?

						second += ((tb[i][j].pocc)[jt]*tb[i][j].dg*tb[ip][jp].dg);
					}
				}

				for(jpp = DUMMY_LENGTH + jt + 1; jpp <= j - len_s + DUMMY_LENGTH + 1; jpp++) {
					for(ipp = DUMMY_LENGTH + 1; ipp <= len_s + DUMMY_LENGTH + 1; ipp++) { //O(n^4) (!)
						if(f[ipp][jpp] < THRESH)
							continue; //?

						second += ((tb[i][j].pocc)[jt]*tb[i][j].dg*tb[ipp][jpp].dg);
					}
				}
			}
		}

		occp = first + second;
	} else if(mode == BEST_TB) {
		if(!ran) {
			//find the single maximum entry in the scoring matrix.
			the_max_i = 0;
			the_max_j = 0;
			for(int i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH + 1; i++){
				for(int j = DUMMY_LENGTH + 1; j <= len_t + DUMMY_LENGTH + 1; j++){
					if(f[i][j] > f[the_max_i][the_max_j]){
						the_max_i = i;
						the_max_j = j;
					}
				}
			}
			traceback(the_max_i,the_max_j,&tb[the_max_i][the_max_j],1);

			//we only need to run this once.
			ran = 1;
		}

		return ((tb[the_max_i][the_max_j].pocc)[jt]*tb[the_max_i][the_max_j].dg)/(1+tb[the_max_i][the_max_j].dg);
	} else if(mode == LOCALIZED) {
		double num,denom;
		int ip,jp;

		//traceback at each i,j has already been precomputed.

		//i' -> ip, j' -> jp
		for(int i = DUMMY_LENGTH + 1; i <= len_s + DUMMY_LENGTH+1; i++) {
			int bdyr = 0;
			if(jt >= len_t - len_s)
				bdyr = jt + len_s - len_t;

			for(int j = DUMMY_LENGTH + jt + 1; j <= DUMMY_LENGTH + jt + len_s + 1 - bdyr;j++) {
				if(f[i][j] < THRESH)
					continue;

				num = (tb[i][j].pocc)[jt]*tb[i][j].dg;

				denom = 1;

				int bdyl = 0; //prevent jp overflows
				int bdyr = 0;
				if(j - DUMMY_LENGTH - 1 <= len_s)
					bdyl = len_s - j + DUMMY_LENGTH + 1;
				if(j - DUMMY_LENGTH - 1 >= len_t - len_s)
					bdyr = j + len_s - len_t - DUMMY_LENGTH - 1;

				for(ip = DUMMY_LENGTH + 1; ip <= DUMMY_LENGTH + len_s + 1; ip++) {

					for(jp = j - len_s + bdyl; jp <= j + len_s - bdyr; jp++) {
						if(f[ip][jp] < THRESH)
							continue;

						denom += tb[ip][jp].dg;
					}
				}
				occp += (num/denom);
			}
		}
	}

	return occp/pf;
}

char invcase(char base) {
	char out;
	if(isalpha(base)) {
		if(islower(base)) {
			out = toupper(base);
		} else if(isupper(base)) {
			out = tolower(base);
		}
	} else {
		out = base;
	}

	return out;
}

void print_alignment(FILE * dest) {
	int j,k,l;

	int the_min_i = 0;
	int the_min_j = 0;

	the_max_i = 0;
	the_max_j = 0;

	for(j = DUMMY_LENGTH + 1; j <= len_s + DUMMY_LENGTH + 1; j++){
		for(l = DUMMY_LENGTH + 1; l <= len_t + DUMMY_LENGTH + 1; l++){
			if(f[j][l] > f[the_max_i][the_max_j]){
				the_max_i = j;
				the_max_j = l;
			}
		}
	}
	traceback(the_max_i,the_max_j,&tb[the_max_i][the_max_j],1);

	for(k=DUMMY_LENGTH;k<=len_t+DUMMY_LENGTH+1;k++) {
		if(!isspace(aligned_t[k])) {
			the_min_j = k;
			break;
		}
	}
	for(k=0;k<=len_t+DUMMY_LENGTH+1;k++) {
		if(!isspace(aligned_s[k])) {
			the_min_i = k;
			break;
		}
	}

	int tlength = MAX(len_s - the_max_i + DUMMY_LENGTH, len_t - the_max_j + DUMMY_LENGTH ) + len_align + 1;
	int start = -MAX(-the_min_i,-the_min_j);

	fprintf(dest,"\ns:");
	for(j = tlength - 1; j >= start; j--) {
		if(aligned_s[j] != 'X') {
			fprintf(dest,"%c",invcase(aligned_s[j]));
		} else {
			fprintf(dest," ");
		}
	}

	fprintf(dest,"\nt:");
	for(j = tlength - 1; j >= start; j--) {
		if(aligned_t[j] != 'X') {
			fprintf(dest,"%c",invcase(aligned_t[j]));
		} else {
			fprintf(dest," ");
		}
	}

	fprintf(dest,"\n  ");

	//numbering
	for(j = 0; j < tlength; j++) {
		if(j % 10 == 0 ) { fprintf(dest,"|"); } else { fprintf(dest," "); };
	}
}

void struct_alloc() {

	tb = (struct traceback **) malloc((len_s+2+DUMMY_LENGTH)*sizeof(struct traceback *)); //length of s
	for(int j = 0; j <= len_s + DUMMY_LENGTH + 1; j++) { //length of t
		tb[j] = (struct traceback *) malloc((len_t+2+DUMMY_LENGTH)*sizeof(struct traceback));
	}

	//we also have to allocate memory for the arrays within each structure.
	for(int j = 0; j <= len_s + DUMMY_LENGTH + 1; j++) {
		for(int l = 0; l <= len_t + DUMMY_LENGTH + 1; l++) {
			tb[j][l].pocc = (int *) malloc((len_t+2*DUMMY_LENGTH+2)*sizeof(int));
		}
	}

	//allocate memory for array of dg structures // TODO: deprecate this
	/*dga = (struct dgarray *) malloc(len_t*sizeof(struct dgarray));
	//and for the arrays within each structure
	for(i = 0; i <= len_t; i++) {
		dga[i].dgs = (int *) malloc(len_s*len_t*sizeof(int));
		dga[i].last = 0;
	}*/

	dgJ = (int *) calloc(len_t,sizeof(int));
}

void print_hotspots(FILE * out) {
	precompute_traceback(out,1);

	for(int i = 0; i < len_t; i++)
		fprintf(out,"%d: %c, %g\n",i + 1,num_to_base(t[DUMMY_LENGTH + len_t - (i + 1)]),(float) -dgJ[len_t - (i + 1)]/100);
}

void print_poccs(FILE * out) {
	int j;
	ran = 0; //in the case of the single occupancy, we only need to compute the traceback once.

	double prtf = dpf(tb_mode);
	if(tb_mode != BEST_TB)
		precompute_traceback(NULL,1);

	for(j = 1; j <= len_t; j++) {
		fprintf(out,"%c %d %g\n",num_to_base(t[DUMMY_LENGTH + j - 1]),len_t - j + 1,p_occ(j - 1,tb_mode,prtf));
	}
}

void cleanup() {
	int i,j;
	for( i = 0; i < NUM_MATRICES; i++ ){
		for( j = 0; j <= len_s + DUMMY_LENGTH*2; j++){
			free(the_matrices[i][j]);
		}
		free(the_matrices[i]);
	}

	for( i = 0; i <= len_s + DUMMY_LENGTH*2; i++){
		free(lu[i]);
		free(ld[i]);
		free(l1xn[i]);
		free(lnx1[i]);
	}
	free(lu);
	free(ld);
	free(l1xn);
	free(lnx1);


	for( i = 0; i < len_s + DUMMY_LENGTH; i++ ){
		free(pairs[i]);
	}
	free(pairs);

	free(structures);

	for(i=0;i < len_s; i++) {
		for(j=0;j < len_t; j++) {
			free(tb[i][j].pocc);
		}
		free(tb[i]);
	}
	free(tb);

	/*for(i = 0; i < len_t; i++) {
		free(dga[i].dgs);
	}
	free(dga);*/
}

int main(int argc, char **argv){
	int i,l; //dummies

	//file location strings
	char *outputfile = (char *) malloc(200*sizeof(char));
	char *sfile = (char *) malloc(200*sizeof(char));
	char *tfile = (char *) malloc(200*sizeof(char));

	//sequences (if not batchmode)
	char *sa;
	char *ta;

	// deal with arguments
	int stdoutbool; //whether to print to stdout
	int stdinbool; //whether s and t should be read from the command line
	int batchmode; //whether or not we're batching. 0 = no, 1 = yes.
	char custstr[50]; //might want to make this dynamic?

	//checks that options were actually given
	int mbool = 0, obool = 0, fbool = 0, pbool = 0, Pbool = 0, gbool = 0, cbool = 0, sbool = 0, tbool = 0; //, dbool = 0;
	while((l = getopt(argc, argv, "m:o:f:p:g:s:t:c:P:C:")) != -1) { //d: removed since unused
		switch(l) {
			case 'm':
				mbool = 1;
				if(optarg[0] == 's') {
					batchmode = 0;
				} else if(optarg[0] == 'b') {
					batchmode = 1;
				} else {
					printf("Improper mode. (-m)\n");
					printf("%s",USAGE);
					exit(0);
				}
				break;
			case 'o':
				obool = 1;
				switch(optarg[0]) {
					case 'a':
						omode = ALIGNMENTS;
						dgout = NONE;
						break;
					case 'd':
						omode = DG;
						break;
					case 'p':
						omode = POCCS;
						dgout = NONE;
						break;
					case 'l':
						omode = HOTSPOTS;
						break;
					default:
						printf("Improper output type. (-o)\n");
						printf("%s",USAGE);
						exit(0);
						break;
				}
				break;
			case 'f':
				fbool = 1;
				if(optarg) {
					outputfile = optarg;
					stdoutbool = 0;
				} else {
					outputfile = "Magic_Value";
					stdoutbool = 1;
				}
				break;
			case 'p':
				pbool = 1;
				switch(optarg[0]) {
					case 'b':
						tb_mode = BEST_TB;
						break;
					case 's':
						tb_mode = SING_SING;
						break;
					case 'p':
						tb_mode = LOCALIZED;
						break;
					default:
						printf("Improper traceback scheme. (-p)\n");
						printf("%s",USAGE);
						exit(0);
						break;
				}
				break;
			case 'P':
				Pbool = 1;
				switch(optarg[0]) {
					case '5':
						mftbool = 1;
						dynpen = 1;
#ifdef DNA
						g0 = 114; //DNA
#else
						g0 = 160; //RNA
#endif
						break;
					case '2':
						mftbool = 1;
						dynpen = 0;
#ifdef DNA
						meanpen = 21; //DNA
						g0 = 114; //DNA
#else
						meanpen = 34; //RNA
						g0 = 160; //RNA
#endif
						break;
					case '0':
						mftbool = 0;
						g0 = 0;
					case 'c':
						mftbool = 1;
						custbool = 1;
						break;
					case 'C':
						mftbool = 0;
						custbool = 1;
						break;
					default:
						printf("Removal parameters improperly specified.\n");
						printf("%s",USAGE);
						exit(0);
						break;
				}
				break;
			case 'g':
				gbool = 1;
				if(optarg[0] == 'b') {
					dgout = BEST;
				} else if(optarg[0] == 's') {
					dgout = SCRG;
				} else {
					printf("%s",USAGE);
					exit(0);
				}
				break;
//			case 'd':
//				dbool = 1;
//				//will need switch statement here
//				break;
			case 'c':
				cbool = 1;
				if(optarg) {
					chem_pot = -strtod(optarg,NULL); //make sure to specify this as NEGATIVE on command line
				} else {
					printf("Must specify chemical potential.\n");
					printf("%s",USAGE);
					exit(0);
				}
				break;
			case 's':
				sbool = 1;
				strcpy(sfile,optarg);
				break;
			case 't':
				tbool = 1;
				strcpy(tfile,optarg);
				break;
			case 'C':
				strcpy(custstr,optarg);
				break;
			default:
				printf("Improper option (-%c)\n",optopt);
				printf("%s",USAGE);
				exit(0);
				break;
		}
	}
	//will now need to loop over the possible remaining arguments containing s and t
	if(optind != argc) {
		stdinbool = 1;
		if(sbool && tbool) {
			printf("Cannot specify both s and t as arguments and files.\n");
			exit(0);
		}

		if(argv[optind] && argv[optind+1]) {
			batchmode = 0; //specifying both s and t on command line will always override batchmode.
			sa = argv[optind];
			ta = argv[optind+1];
		} else if(argv[optind]) {
			sbool ? (ta = argv[optind]) : (sa = argv[optind]);
		}
	} else {
		stdinbool = 0;
		if(!sbool && !tbool) { //user never gave any s or t
			printf("Sequences must be provided.\n");
			printf("%s",USAGE);
			exit(0);
		}
	}

	//check to see that all the necessary options are present.
	if(!obool) {
		printf("%s",USAGE);
		exit(0);
	} else {
		if(!mbool)
			batchmode = 0;
		if(!fbool) {
			stdoutbool = 1;
			outputfile = "Magic_Value";
		}
		if(!pbool)
			tb_mode = SING_SING;
		if(!Pbool) {
			mftbool = 1;
			dynpen = 1;
			g0 = 160;
		}
		if(!gbool) {
			if((omode == ALIGNMENTS) || (omode == POCCS)) {
				dgout = NONE;
			} else {
				dgout = BEST;
			}
		}
		if(!cbool && (omode == POCCS)) {
			printf("Chemical potential (-c) must be specified in this mode.\n");
			printf("%s",USAGE);
			exit(0);
		} else if(cbool && (omode != POCCS)) {
			chem_pot = 0.0;
		}
	}

	//custom penalties
	if(custbool) {
		char * splitbuf;
		if(mftbool) {
			splitbuf = strtok(custstr,",");
			int j = 0;
			g0 = (int) strtod(splitbuf,NULL)*100;
			while(splitbuf != NULL) {
				if(j > 5) {
					printf("Too many custom parameters specified.\n");
					exit(0);
				}
				nucenergy[j-1] = (int) strtod(splitbuf,NULL)*100;
				splitbuf = strtok(NULL,",");
				j++;
			}
			if(j < 5) {
				printf("Too few custom parameters specified.\n");
				exit(0);
			}
		} else {
			splitbuf = strtok(custstr,",");
			g0 = (int) strtod(splitbuf,NULL)*100;
			splitbuf = strtok(NULL,",");
			meanpen = (int) strtod(splitbuf,NULL)*100;

			splitbuf = strtok(NULL,",");
			if(splitbuf != NULL) {
				printf("Too many custom parameters specified.\n");
				exit(0);
			}
		}
	}

	//static vs. dynamic penalty
	if(dynpen != 1) {
		for(i=0;i<=3;i++) {
			nucenergy[i] = meanpen;
		}
	}

	//encode sequences.
	int flengtht, flengths;
	char ** sseqs;
	char ** tseqs;
	if(batchmode == 1) {
		if(!tbool) { //user must always provide t
			printf("Must specify t as file.\n");
			exit(0);
		}

		if(stdinbool == 0) { //if both batch files are specified as arguments
			flengths = filelength(sfile);
			flengtht = filelength(tfile);

			if(flengths != flengtht) {
				printf("For batch mode, list of s sequences must be same length as list of target sequences. (s list has length %d, t list has length %d)",flengths,flengtht);
				exit(0);
			}

			sseqs = fencode(sfile,flengths);
			tseqs = fencode(tfile,flengtht);
		} else { //they've specified s on the command line but not t
			if(sbool) {
				printf("For single-to-many batch mode, s must be specified on command line, and t provided in file.\n");
				exit(0);
			}

			len_s = strlen(sa);
			s = (int*) malloc( (len_s + 1 + DUMMY_LENGTH*2)*sizeof(int) );
			encode(s, sa, len_s, 1); //1 because we want 3' -> 5' for s

			flengtht = filelength(tfile);
			tseqs = fencode(tfile,flengtht);
		}
	} else if(batchmode == 0) {
		if(stdinbool == 1) { //s and/or t has been specified on the command line
			if(!sbool) {
				len_s = strlen(sa);
				s = (int*) malloc( (len_s + 1 + DUMMY_LENGTH*2)*sizeof(int) );
				encode(s, sa, len_s, 1); //1 because we want 3' -> 5' for s
			} else {
				flengths = filelength(sfile);
				sseqs = fencode(sfile,flengths);
			}

			if(!tbool) {
				len_t = strlen(ta);
				t = (int*) malloc( (len_t + 1 + DUMMY_LENGTH*2)*sizeof(int) );
				encode(t, ta, len_t, 0); //0 because we want 5' -> 3' for t
			} else {
				flengtht = filelength(tfile);
				tseqs = fencode(tfile,flengtht);
			}
		} else { //they've specified both files as single sequences
			flengths = filelength(sfile);
			flengtht = filelength(tfile);

			if(flengths != flengtht) {
				printf("For single mode, list of s sequences must be same length as list of target sequences. (s list has length %d, t list has length %d",flengths,flengtht);
				exit(0);
			}

			sseqs = fencode(sfile,flengths);
			tseqs = fencode(tfile,flengtht);
		}

		if((sbool && flengths != 1) || (tbool && flengtht != 1)) { //TODO: this needs to ensure that sequence files are singleline, but only if applicable (if sbool for s, tbool for t)
			printf("For single mode, each sequence file must contain a single line.\n");
			exit(0);
		}
	}

	aligned_s = (char*) malloc((len_s+len_t+DUMMY_LENGTH+1)*3*sizeof(char));
	aligned_t = (char*) malloc((len_s+len_t+DUMMY_LENGTH+1)*3*sizeof(char));

	struct_alloc();

	/// PRINT OUTPUT
	char farg[512];
	FILE * test;
	FILE * output;
	if(!strcmp(outputfile,"Magic_Value")) {
		sprintf(farg,"rm %s",outputfile);
		test = fopen(outputfile,"r");
		if(test != NULL) { //remove file if it exists.
			if(strcmp(outputfile,"/dev/null") != 0)
				system(farg);
			fclose(test);
		}
	}
	if(stdoutbool == 0) {
		output = fopen(outputfile,"a");
	} else {
		output = stdout;
	}

	if(batchmode == 0) { //if we only want to examine a single sequence ...
		if(sbool) {
			get_s(sseqs[0],strlen(sseqs[0]));
			len_s = strlen(sseqs[0]);
		}
		if(tbool) {
			get_t(tseqs[0],strlen(tseqs[0]));
			len_t = strlen(tseqs[0]);
		}

		setup();
		fill_arrays();

		struct_alloc();

		switch(omode) {
			case ALIGNMENTS:
				print_alignment(output);
				fprintf(output,"\n");
				break;
			case DG:
				ran = 0;
				precompute_traceback(output,0);
				break;
			case POCCS:
				print_poccs(output);
				break;
			case HOTSPOTS:
				print_hotspots(output);
				break;
		}

		cleanup();
	} else { //otherwise, look at multiple sequences
		for(i=0;i<flengtht;i++) {
			if(stdinbool == 0) { //if both batchfiles specified on command line
				len_s = strlen(sseqs[i]);
				get_s(sseqs[i],strlen(sseqs[i]));
			} //else, one s vs. many t
			get_t(tseqs[i],strlen(tseqs[i]));
			len_t = strlen(tseqs[i]);

			setup();
			fill_arrays();

			struct_alloc();

			switch(omode) {
				case ALIGNMENTS:
					print_alignment(output);
					break;
				case DG:
					ran = 0;
					precompute_traceback(output,0);
					break;
				case POCCS:
					print_poccs(output);
					break;
				case HOTSPOTS:
					print_hotspots(output);
					break;
			}
			printf("*\n");
			cleanup();
		}
	}
	fclose(output);

	return 0;
}
