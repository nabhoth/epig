/***************************************************************************
*   Copyright (C) 2006 by martin lukac   				  *
*   lukacm@ece.pdx.edu   						  *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/
//Uncomment to allow multi-threaded computation for smaller circuits 8 < qubits
//#define __SSIZE__
//Uncoment to allow QMDD representation for circuits up to 7 qubits
//#define __QMDD__
//Uncoment to unleash the standard things
//#define __STDR__
//#define __TIMING__

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <pthread.h>
#include <complex>
#include <sstream>

#ifdef __MPI__
   #include "mpi.h"
#endif
#ifdef __CUDA_SDK__
#include <cutil.h>
#include <cutil_inline.h>
#endif
#include "cublas.h"
#include "cuda.h"

#include <gsl_cblas.h>
#include <gsl_complex_math.h>
#include <gsl_blas.h>

using namespace std;

// Thread block size
#define BLOCK_SIZE 16

// max number of individuals
const int MAXNUMBER = 100;
//size of the group for group fitnes evaluation
const int GROUPSIZE = 15;	
// max number generations allowed
const int MAXGEN = 1000000;	
//const float CONSTRAINT = 0.5; // 50% of the wheel has fitness of 1
// max number of gates used
const int MAXNUMOFGATES = 50;	
// max number of gates used
const int MAXNUMOFWIRES = 1000;	
// 2 to the n        er, here set n = 10, so 2^10*2^10. also max number of individual matirx size
const int MAXGATEDEM = 1048576;	
// max number of wires for input of an individual is MAXINDIWIRENUM
const int MAXINDVWIRENUM = 15;	
// max number of segments per individual is MAXSEGNUM
const int MAXSEGNUM = 200;	
// max number of chars in an individual string
const int MAXSTRING = 150;	
// max number of input for any single gate
const int MAXGATEINPUT = 15;	
// max number of measurement gates
const int MAXMEASUREGATE = 50;	
// we are not looking for a single gate solution (before minimization)
const int MINCIRCUITSIZE = 2;	
//maximum  length of a sequence to detect
const int MAXSEQUENCESIZE = 500;
//maximum allowed valuedness for the synthesizer
const int MAXVALUEDNESS = 7;

const int CHARMAX = 126;
const int CHARMIN = 33;
const int REPSIZE = 3;
typedef struct d_list group;

/* global mutex for our program. assignment initializes it */
static pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;


/***********************
* the translation structure - between the GA and the Exahustive searcher
* or other modular extensions
************************/
struct nameT
{
	string normal;
	string nonorm;
};
//structure used in the stringArray
struct tempT
{
// number of gates in this catagory
	int numberGates;	
// for a particular number of input and output, the list of the gates can be used
	char *nameGate[MAXNUMOFGATES];	
};

struct aGate
{
// number of input/output wires
	int numI;
	int *coordI;
	int numO;
	int *coordO;
// cost of the gate
	int Cost;		
// representation
	char *representation;
// LUT representing this gate
	int length_LU;
	int *inLU;	
	int *outLU;	
	char **inBinLU;	
	char **outBinLU;	
//The String Name of this gate
	string my_string;	
	int valuedness;
	string rev_lib_string;
	string rev_lib_char;
	bool negated_vars;
};

struct aConnection
{
	int input_index;
	int output_index;
	//0,1 - logic, -1 - non initiated, -2 - processed
	int current_value;
	bool circuit_input;
	bool circuit_output;
	//connection required and cannot be modified
	bool unmutable;
	bool control;
	//simple index indicating the position of the input of this connection
	int location;
	bool negated;
};

/***********************
* Reversible Circuit
************************/
struct aCircuit
{
	int numIO;
	int numGates;
	int numI;
	int *I;
	int numO;
	int *O;
	int *data_I;
	int *data_O;
	aGate *gates[MAXNUMOFGATES];
	aConnection *connections[MAXNUMOFWIRES];
	double fitness;
	double error;
	double cost;
	bool doubled;
	bool cascade;
};
/************************
*Holds the definition of the result function/behavior
************************/
struct result
{
	float avFitness [100000];
	float avError[100000];
	float avCost[100000];
	int counter;
	int step;
};
// structure representing an individual in the GA population
struct Individual
{
// number of wires of this individual
	int ioNumber;		
// number of segments at the creation, used only for intialization
	int segmentNumber;	
// circuit representation as string
	string my_string;	
// circuit representation as string
	string rev_lib_string;	
//each circuit is now generated with phase
	double phase;	
	double phases[MAXNUMOFGATES];
// fitness value
	float fitness;		
// Cost of the circuit
	float Cost;		
// Error of the circuit
	float Error;		
// for share fitness calcualtion
	int Groups;		
	int Rank;
	//to avoid races copy the gate arrays to every single individual
// all agtes are stored in here
	aGate *gateArray[MAXNUMOFGATES];	
	int valuedness;

};

struct Solution
{
	int ioNumber;
	float error;
	float fitness;
	float cost;
	string my_string;
};

struct d_list
{
	int V;
	struct d_list *next;
};

// compare two gates
int compareColumns(char*, char**, int, int);
double compareColumns(char**, char**, int, int);
double compareColumns(char**, char**, int, int, int);
bool findInArray(int, int*, int);
bool compareRep (char*, char*);	
bool areFromSameGate(aCircuit*, aConnection*, aConnection*);
void copyRep (char*, char*);	
aGate* ccGate (aGate*, aGate*);	
void incRep(char*);
void incRep(char*, int, int);
void initRep(char*, int);
void initRep(char*);
int getGatePos(int, string);
void generateGateIndexes(aGate*, int);
void getGateRepfromStr(int, string, char*);
bool gateInputsSat(aConnection**, int, aGate*, bool);
int getAvailableConn(aConnection**, int);
void generateConnections(aConnection**, int, aGate*, bool, bool);
int getConnByOutput(aConnection**, int, int);
int getConnByInput(aConnection**, int, int);
char getConnValByOutput(aConnection**, int, int);
void copyIntArray(int*, int*, int);
void copyConnections(aConnection**, aConnection**, int);	
// makes a copy oif a quantum gate
void vectorInitZero(int, double *);
void vectorScalarMult(int, double, double *, double *);
aGate tensorProduct (aGate *, aGate *);	
// tensor product of gate A and gate B, returning A
aGate *tensorProduct2 (aGate *, aGate *);	
//generates an Identity Matrix
void destroyGate(aGate *); 
//generates an Identity Matrix
void initGate(aGate *, int, int, int, int); 
void initConnection(aConnection*); 
//generates an Identity Matrix
void initMatrix(aGate *, int); 
//generates empty Matrix
void zeroMatrix(aGate *, int); 
//generates matrix filled with don't cares
void dcMatrix(aGate *, int);
//return the number of proper gates in  the String
int getGatesfromString (string);	
//subtracts two matrices returns the result in the first one
aGate* matrixSubstr(aGate*, aGate*); 
//matrix multiplication using BLAS -> CBLAS
void cblasMatrixProduct(aGate*, aGate*, aGate*);
//matrix multiplication using CUDA BLAS -> CUBLAS
void cublasMatrixProduct (aGate *, aGate *, aGate *, int); 
// matrix product 'matrixProduct3 (aGate *R, aGate *A, aGate *B)' of gate A and gate B into R
void matrixProduct3 (aGate *, aGate *, aGate *);	
// tensor product 'tensorProduct3 (aGate *R, aGate *A, aGate *B)' of gate A and gate B into R
void tensorProduct3 (aGate *, aGate *, aGate *);	
double cuFloatComplex_arg (double);
double cuFloatComplex_abs (double);
double cuFloatComplex_inverse (double);
double cuFloatComplex_logabs (double);
double cuFloatComplex_pow (double, double);
double cuFloatComplex_sqrt (double);
double cuFloatComplex_div (double, double);
double cuFloatComplex_sub (double, double);

short sign(float);
void dec2bin(long, char*, int);
int bin2dec(char*, int);
void cuBlasMMulti(int, aGate*, aGate*, aGate*);
//multiplies a Matrix by a Vector and returns a Vector, Matrix is of size w x w 
void cpuMulv(double*, double*, double*, int); 
void cleanCUDA(double*, double*,double*,double*,double*,double*);
__global__ void matrixMul(double* A, double* B, int wA, int wB, double* C);
__global__ void Muldv(double*, double*, int, double*);
__global__ void vvMul(double*, double*, int, double*, int);
__global__ void Divdv(double*, double*, int, double*);
__global__ void initState(int, double, double*, int);
__global__ void nullState(int, double*);
//requires two host matrices and a host vector, two device matrices and a device vector
void Mulv(double*, double*, int,  double*, double*, double*, double*);
//requires three host matrices and three device matrices pointers 
void Mulm(double*, double*, int, double*, double*, double*, double*);
//matrix multiplication using device matrices pointers
void Mulmnoset(int, double*, double*, double*); 
//matrix multiplication using device matrices pointers
void Mulvnoset(int, double*, double*, double*); 
//matrix multiplication using device matrices pointers
void Mulvvnoset(int, double*, double*, double*); 
void Divvnoset(int, double*, double*, double*);


// holds GA functions
class GA
{
      public:


//The counter of individual in the Threaded mode
	static int populationCounter;	
//the method used to Thread the GA
	static void *startSub (void *arg);	
//measurement switch 0 - turns of measurement, >0 number of measured qubits so as the order of measurement
	static int measurement;	

	void run ();

	  GA ()
	{
	}
	 ~GA ()
	{
#ifdef __CUDA__
		cleanCUDA(d_M1, d_M2, d_M3, d_MV, d_VI, d_VO);
#endif
	};

	 /////////////
	 //GA settings
	 /////////////
	void setpopulationNumber(int);
	int getpopulationNumber();
	 void setsegments(int);
	 int getsegments();
	 void setmincircuitsize(int);
	 int getmincircuitsize();
	 void setalterations(int);
	 int getalterations();
	 void setproba_Mutation(float);
	 float getproba_Mutation();
	 void setproba_Crossover(float);
	 float getproba_Crossover();
	 void setalpha(float);
	 float getalpha();
	 void setbeta(float);
	 float getbeta();
	 void setalpha1(float);
	 float getalpha1();
	 void setbeta1(float);
	 float getbeta1();
	 void setdivider(int);
	 int getdivider();
	 void setphase(int);
	 int getphase();
	 void setdisplaylevel(int);
	 int getdisplaylevel();
	 void setGa(int);
	 int getGa();
	 void setmutation(int);
	 int getmutation();
	 void setcrossover(int);
	 int getcrossover();
	 void setreplicator(int);
	 int getreplicator();
	 void setfitness(int);
	 int getfitness();
	 void setgrouped(int);
	 int getgrouped();
	 void setpareto(int);
	 int getpareto();
	 void setthreshold(int);
	 int getthreshold();
	 void setresultnum(int);
	 int getresultnum();
	 void setmeasurement(int);
	 int getmeasurement();

// initialize the population
	void initializePop (int, string, string);	
// With a 1 point crossover probability recombine 2 parents into 2 new child
	void apply1Crossover (int, int);	
// With a 2 point crossover probability recombine 2 parents into 2 new child
	void apply2Crossover (int, int);	
// A CO operation that preserves the structure a bit more - With a 1 point crossover probability recombine 2 parents into 2 new child
	void apply1CrossoverP (int, int);	
// A CO operation that preserves the structure a bit more - With a 2 point crossover probability recombine 2 parents into 2 new child
	void apply2CrossoverP (int, int);	
	void doMutation(aCircuit*);
	void doCrossover(aCircuit*, aCircuit*);
// Replication of teh selected individuals either SUS or RW
	int applyReplication ();	
	void makeFitness(aCircuit*);
	void evaluateCircuit(aCircuit*, bool);
	void *doMatrixFitness(Individual*,bool);
	void *doMeasureFitness(Individual*,bool);
	void *doMultiMeasureFitness(Individual*,bool);
	void printCircuit(aCircuit*);
void makeThreeWiresSimpleC(aCircuit*);
void makeThreeWiresRestrictedC(aCircuit*);


	void applyRMutation (int, bool);
//A mutation operator that does preserbe the current circuit structure completely
	void applyRMutationP(int, bool);		
	void applyIMutation (int, bool);
//A mutation operator that does preserbe the current circuit structure completely
	void applyIMutationP(int, bool);		
	void evaluatePopulation ();
// somemethod setting final criteria
	void setTerminatingCondition (int);	
// some method returning if the final conditions were reached
	bool terminatingCondition ();	
// calculation of the individual matrix representation for each individual
	void individualMatrix (Individual * i);	
// clear everything in one individual
	void clearAll (int);	
// method minimizing circuit ... not used always
	int minimizeCirc (Individual *);	
	void repairForceCircuit (int);
// method close streams to files
	void closeStream ();	
// writes down the result of the run
	void output (int);	
// writes down the best individual from one generation and/or displays the statistics of the overall run
	void outputBest (bool);	
//copy two indidual by their addresses
	void setIndv (aCircuit *, aCircuit *);	
	void injectNormalSegments (string **);
//translates all parallel blocks to individual NumIO sized circuits
	aGate **getPermutativeSegments ();	
	void calcAverage();

	// size of the population < MAXNUMOFGATES
	int populationNumber;
//the set of probabilistic expectations generated as the result of measurement of the individual states
	double *expectationsAllState[MAXGATEINPUT*MAXVALUEDNESS];   
	// both here bellow fields will
	// be in the future replaced by
	// array of pointers
// population (array of object Individual)
	aCircuit *population[MAXNUMBER];	
// population of individuals used for replication ... non optimal to be improved
	aCircuit *new_population[MAXNUMBER];	
// best Individual from the run
	aCircuit *bestIndv;	
	//method of fitness application
	//0 - normal, 1 - groupped shared
	bool grouped;
	//indicator of pareto rplication
	//0 - no, 1 - yes
	bool pareto;
// looks for the best individual with fitness 1
	bool getFitness ();	
	nameT *normNames[MAXNUMOFGATES];
//the set of meaurement operators defined by parameters from file
	aGate *measurements[MAXMEASUREGATE];	
//the set of meaurement operators defined by parameters from file
	 aGate *measurementsAllState[MAXMEASUREGATE];    
//two fields per input by Max size of the input-output mapping
	double *measureexpected[MAXGATEINPUT*MAXVALUEDNESS]; 
	int measurementQBits[MAXGATEINPUT];
	int measured_desired_output_records;
	int *measured_desired_output_records_idx;
	aGate *permutativeArray[MAXNUMOFGATES];
// all agtes are stored in here
//	aGate **gateArray;
	aGate *gateArray[MAXNUMOFGATES];

// retunrs teh char-representation of the gate at teh int index from gateArray
	int getGate (char*);	
// number of gates to be use < MAXNUMOFGATES
	int numofgates;		
	int numofnormalgates;
// gate/circuit we are looking for
	aGate *finalGate;	
// gate/circuit we are looking for
	Individual *finalIndividual;	
//how many generations does the GA have to do before call a module
	int alterations;	
//display switch
	bool display;
	int progress;
	result finalResult;
	//the level of output details
	//0 - th lowest
	//1 - middle
	//3 - full;
	int displaylevel;

      private:

	ostringstream convert_stream;
// input file stream
	ifstream in_stream;	
// output file stream
	ofstream out_stream;	

	//device matrices for Matrix Matrix multiplication
#ifdef __CUDA__
	double* d_M1;
	double* d_M2;
	double* d_M3;
	//device Matrix and Vector for Matrix Vector multiplication
	double* d_MV;
	double* d_VI;
	double* d_VO;

	double* d_ME_00;
	double* d_ME_01;
	double* d_ME_10;
	double* d_ME_11;

	double* d_VV_M;
	double* d_Value;
#endif
	aGate *measures_o[10], *measures_i[10];
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

	int avrghst[1000];
	int besthst[1000];
	// probabilities of mutation and crossover
	float proba_Mutation;
	float proba_Crossover;
	// scaling for the complex fitness function
	// alpha + beta = 1
	float alpha;
	float beta;
	//scaling fo the pareto replication
	//default alpha1 = 1, beta1 = 1
	float alpha1;
	float beta1;
	//number of segments estimated for solution
	int segments;
	//minimum number of segments - size pressure
	int mincircuitsize;
	// estimated cost divider for the solution circuit
	int divider;
	// type of GA Baldwinian or normal
	//0 - nomral, 1 - baldwinian
	int Ga;
	// use elitism
	//0 - yes, 1 - no
	int elitism;
	//setting different types of mutation
	//0 - normal, 1 - bitwise
	int mutation;
	//setting different types of crossover
	//0 - 1point, 1 - 2point
	int crossover;
	//different types of selection
	//0 - RW, 1- SUS
	int replicator;
	// type of fitness function
	//0 - simple0, 1 - simple1, 2 - complex0, 3 - complex1
	int fitness;
	//Threshold for replication
	//0 - no threshold otherwise must be smaller than 1
	float threshold;
	//more general it is used to allow synthesis from pure gates or rotations
	int phase;
	//number of user specified generations
	int generations;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

	// nullGate - wire
	aGate nullGate;		
	// sorted by input wires in asending order
	tempT *stringArray;	
	// the number of wire for input and output of the result gate
	int resultnum;		
	//indicates if we operate only on variables or also on their complements
	int doubled;		
	//indicates if we do or not cascade the control variable
	int cascade;
	// termination criterian
	int generationCondition;	
	// some boolean condition
	bool condition;		
	float best;
	float (*history)[4];

	void initiateStorage(); 
	int getOutputLUTIndex(aGate*, int);

//matrix manipulators --simple
// matrix product of gate A and gate B
	void matrixProduct2 (aGate *, aGate *);	
// reduce teh total length of the circuit
	void reduceSegNumber (int);	
// reads teh gates from a file and set them in an array
	void initializeGates ();	
// creates string representation for each individual
	void makeIndvWiring (aCircuit*);	
// set an array used only for intialization
	void setStringArray ();	
// method repairing circuit after mutation automatically called from mutation
	void repairCircuit (int);	
#ifdef __CUDA__
	aGate *computeCUDAMatrix (Individual *);
	aGate *computeCUBLASMatrix (Individual *);
	aGate *computeCUDANATIVEMatrix (Individual *);
#else
	aGate *computeMatrix (Individual *);
#endif
	void* computeMeasurement(Individual *);
	//for minimizer of quantum circuits
// compares one quantum gate to all available
	int findGate (aGate *);	
// remove redundantsegments in the circuit
	void removeRedund (Individual *);	
// creates a new segment of gates in string representation
	basic_string < char >synthSegment (int, int);	
// merge two equivalent or equal segments
	void mergeSegs (Individual *, int, int);	
	int getaGate(string);
	void decodeIndvStrToRevLibStr(Individual *);
	void StrToIntArray(string, int*);
	void IntArrayToStr(int* , string);
	void loadGateArray(int, char*);
	int fact(int);

};





