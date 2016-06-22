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

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex>
//#include <bits/cmathcalls.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <pthread.h>
#include <string>
#include <sstream>

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"


const int MAXNUMBER = 500;	// max number of individuals
const int GROUPSIZE = 15;	//size of the group for group fitnes evaluation
const int MAXGEN = 1000000;	// max number generations allowed
//const float CONSTRAINT = 0.5; // 50% of the wheel has fitness of 1
const int MAXNUMOFGATES = 500;	// max number of gates used
const int MAXGATEDEM = 64;	// 2 to the n        er, here set n = 5, so 2^5 = 32. also max number of individual matirx size
const int MAXINDVWIRENUM = 6;	// max number of wires for input of an individual is MAXINDIWIRENUM
const int MAXSEGNUM = 200;	// max number of segments per individual is MAXSEGNUM
const int MAXSTRING = 150;	// max number of chars in an individual string
const int MAXGATEINPUT = 6;	// max number of input for any single gate
const int MAXMEASUREGATE = 12;	// max number of measurement gates
const int MINCIRCUITSIZE = 2;	// we are not looking for a single gate solution (before minimization)
using namespace std;
typedef struct d_list group;
typedef complex<float> cplex;


/***********************
* the translation structure - between the GA and the Exahustive searcher
* or other modular extensions
************************/
struct nameT
{
	string normal;
	string nonorm;
};
// temporary structures used only in incialization
struct tempT
{
	int numberGates;	// number of gates in this catagory
	char nameGate[MAXNUMOFGATES];	// for a particular number of input and output, the list of the gates can be used
};

struct qGate
{
	int numIO;		// number of input/output wires
	int Cost;		// cost of the gate
	char representation, parentA, parentB;	// representation
	  complex < float >gateMatrix1[MAXGATEDEM][MAXGATEDEM];	// matrix representing this gate
	string my_string;	//The String Name of this gate
	int valuedness;

};

struct stringqGate
{
	int wires[MAXGATEINPUT];
	string my_string;	//The String Name of this gate
	char representation;
	int valuedness;
};
struct qState
{
	string initstate;
	string currentstate;
	int currentpos;
	int initpos;
	int valuedness;

};
/************************
*Holds the definition of the result function/behavior
************************/
struct qFunction
{
	complex<float> preBitFunction[MAXGATEINPUT][MAXGATEDEM];
	int measurementQBits[MAXGATEINPUT];
};

struct qDataGate
{
	int numIO;		// number of input/output wires
	int valuedness;
	float trace;
	  complex < float >densityErrorMatrix[MAXGATEDEM][MAXGATEDEM];	// matrix representing this gate
	string my_string;
};

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
	int ioNumber;		// number of wires of this individual
	int segmentNumber;	// number of segments at the creation, used only for intialization
	string my_string;	// circuit representation as string
	  complex < float >phase;	//each circuit is now generated with phase
	  complex < float >phases[MAXNUMOFGATES];
	float fitness;		// fitness value
	float Cost;		// Cost of the circuit
	float Error;		// Error of the circuit
	int Groups;		// for share fitness calcualtion
	int Rank;
	//to avoid races copy the gate arrays to every single individual
	qGate *gateArray[MAXNUMOFGATES];	// all agtes are stored in here
	int valuedness;

};

struct Crc
{
	int ioNumber;		// number of input/output wires
	int Counter;		// counter of all representation so far found
	char representation;	// representation -- initial
	  complex < float >gateMatrix1[MAXGATEDEM][MAXGATEDEM];	// matrix representing this circuit

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

complex < float >getConjugate (complex < float >);
bool compareGate (qGate, qGate);	// compare two gates
qGate copyGate (qGate *);	// makes a copy oif a quantum gate
qGate *cGate (qGate *);
qGate *ccGate(qGate*, qGate*);
void densityErrorMatrix (qGate *, qGate *);
qGate matrixProduct (qGate *, qGate *);	// matrix product of gate A and gate B
complex < float >generatePhase ();
void generatePhases (Individual *);
qGate tensorProduct (qGate *, qGate *);	// tensor product of gate A and gate B
qGate *tensorProduct2 (qGate *, qGate *);	// tensor product of gate A and gate B, returning A
void initMatrix(qGate *, int); //generates an Identity Matrix
int getGatesfromString (string);	//return the number of proper gates in  the String
qGate* matrixSubstr(qGate*, qGate*); //subtracts two matrices returns the result in the first one

class QE
{
      public:
	QE ()
	{
		SolutionCount = 5;
	}
	 ~QE ()
	{
	};

/////////////////////////////////////////////////////////////////////////////
//pointer to the synthetized circuit
	qGate *INDV;
//pointer to teh synthesized circuit properties
	qDataGate *dataINDV;
//field of pointers to all gates read from file
	qGate *gateArray[MAXNUMOFGATES];
//the best gate -- the result we are looking for
	qGate bestOne;
//field with all circuits we are searching for
	Crc *resultCircuits[MAXNUMBER];
//field with names of input gates
	string nameArray[MAXNUMOFGATES];
//field of gates to be used only for specified gates
//and their permutations
	qGate *restricArray[MAXNUMOFGATES];
//field of gates for restricted synthesis
//only pointers to array indexes
	int *restriGate[MAXNUMOFGATES];
//solutions in the form of strings
	string *solutionArray[MAXNUMOFGATES];
//solution 
	Solution bestlocal;
//a character
	char end;
//number fo gates from file
	int numofgates;
//terinating condition
	int generationCondition;
//number of desired segment
	int Segments;
//some global help variable
	int segPointer;
//gates organized b number of wires
	int gateCounter[MAXSEGNUM];
//number of resulting gates
	int population;
//the switch enabling/disabling the use of phase - 

//memory for the Segment manipulator
	qState *segmentstate;

//the maximum parametrized allowed of solutions
	int SolutionCount;
//sometimes we need a general gradient value
	float error;
//streams
	ifstream in_stream;
	ofstream out_stream;

	bool compareIndividual (qGate *);
	group *insertB (group * A, group * m_List);
	void insertA (group * A, group * m_List);
	group *removeFirst (group * m_List);
	group *moveL (group * m_List);
	int sizeL (group * m_List);
	int getGate (char b);
	bool setTerminatingCondition (int type);
	complex < float >assignVal (float x);
	void addGateSet (qGate * array[]);
	void setBestIndividual ();
	void closeStream ();
	void synthSegment (bool begin);
	bool synthPermutations (bool begin);
	bool synthAlterations (bool begin);
	bool synthJumps (bool begin);
	void synthRestric (bool begin);
	group *genPerms (group * head, bool perm);
	void setSegments (int begin);
	void restricSegment ();
	qGate *matrixProduct (qGate * A, qGate * B);
	void matrixProduct2 (qGate * R, qGate * A, qGate * B);
	qGate *densityErrorMatrix (qGate * A, qGate * B);
	bool compareGate (qGate * A);
	void Matrix (bool display);
	void subInit ();
	string **setgateNormArray (qGate **, qGate, int, float);

};

// holds GA functions
class GA
{
      public:


	static int populationCounter;	//The counter of individual in the Threaded mode
	static void *startSub (void *arg);	//the method used to Thread the GA
	static int measurement;	//measurement switch 0 - turns of measurement, >0 number of measured qubits so as the order of measurement

	void run ();

	  GA ()
	{
	}
	 ~GA ()
	{
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

	//methods
	void initializePop (int, string, string);	// initialize the population
	void apply1Crossover (int, int);	// With a 1 point crossover probability recombine 2 parents into 2 new child
	void apply2Crossover (int, int);	// With a 2 point crossover probability recombine 2 parents into 2 new child
	void applyMutation(int, bool);
	void applyReplication ();	// Replication of teh selected individuals either SUS or RW
	//void makeFitness(Individual*, qGate*, float, float);      // according to the input it calculates different fitness values based on the error and the cost 
	void *doFitness();
	void *doMatrixFitness(Individual*,bool);
	void *doMeasureFitness(Individual*,bool);
	void *doMultiMeasureFitness(Individual*,bool);


	void makeFitness (Individual *, bool);
	void makeMeasureFitness(Individual *, bool);
	void applyRMutation (int, bool);
	void evaluatePopulation ();
	void setTerminatingCondition (int);	// somemethod setting final criteria
	bool terminatingCondition ();	// some method returning if the final conditions were reached
	void individualMatrix (Individual * i);	// calculation of the individual matrix representation for each individual
	void clearAll (int);	// clear everything in one individual
	int minimizeCirc (Individual *);	// method minimizing circuit ... not used always
	void repairForceCircuit (int);
	void applyIMutation (int, bool);
	void closeStream ();	// method close streams to files
	void output (int);	// writes down the result of the run
	void outputBest (bool);	// writes down the best individual from one generation and/or displays the statistics of the overall run
	Individual setIndv (Individual);	// copy one Individual and output it
	void setIndv (Individual *, Individual *);	//copy two indidual by their addresses
	void injectNormalCircuits (string **);
	void injectNormalSegments (string **);
	qGate **getPermutativeCircuits ();	//generates all gates from the at the NumIO of the whole circuit
	qGate **getPermutativeSegments ();	//translates all parallel blocks to individual NumIO sized circuits
	void calcAverage();

	// parameters for the GA
	// size of the population < MAXNUMOFGATES
	int populationNumber;
	cplex expectationsAllState[MAXMEASUREGATE*2][MAXGATEINPUT*4];   //the set of probabilistic expectations generated as the result of measurement of ht eindividual states
	// both here bellow fields will 
	// be in the future replaced by
	// array of pointers
	Individual *population[MAXNUMBER];	// population (array of object Individual)
	Individual *new_population[MAXNUMBER];	// population of individuals used for replication ... non optimal to be improved
	Individual *bestIndv;	// best Individual from the run
	qDataGate bestData;
	//method of fitness application
	//0 - normal, 1 - groupped shared
	bool grouped;
	//indicator of pareto rplication
	//0 - no, 1 - yes
	bool pareto;
	bool getFitness ();	// looks for the best individual with fitness 1
	qGate *gateArray[MAXNUMOFGATES];	// all agtes are stored in here
	qGate *permutativeArray[MAXNUMOFGATES];
	nameT *normNames[MAXNUMOFGATES];
	qGate *measurements[MAXMEASUREGATE];	//the set of meaurement operators defined by parameters from file
	 qGate *measurementsAllState[MAXMEASUREGATE];    //the set of meaurement operators defined by parameters from file
	complex<float> measureexpected[MAXGATEINPUT*2][MAXGATEDEM]; //two fields per input by Max size of the input-output mapping
	int measurementQBits[MAXGATEINPUT];

	int getGate (char);	// retunrs teh char-representation of the gate at teh int index from gateArray
	int numofgates;		// number of gates to be use < MAXNUMOFGATES
	int numofnormalgates;
	qGate *finalGate;	// gate/circuit we are looking for
	int alterations;	//how many generations does the GA have to do before call a module
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
	ifstream in_stream;	// input file stream
	ofstream out_stream;	// output file stream

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

	qGate nullGate;		// nullGate - wire
	tempT stringArray[MAXINDVWIRENUM];	// sorted by input wires in asending order
	int resultnum;		// the number of wire for input and output of the result gate
	int generationCondition;	// termination criterian
	bool condition;		// some boolean condition
	float best;
	float (*history)[4]; 
	
	void applySingleObjReplication ();
	void applyParetoReplication ();	// Replication for pareto optimal GA
	void initiateStorage(); //allocates the array of flots for statistics

//matrix manipulators --simple
	void matrixProduct2 (qGate *, qGate *);	// matrix product of gate A and gate B
	//void tensorProduct2 (qGate *, qGate *);	// matrix product of gate A and gate B
	void reduceSegNumber (int);	// reduce teh total length of the circuit
	void initializeGates ();	// reads teh gates from a file and set them in an array
	void makeIndvString (Individual*);	// creates string representation for each individual
	void setStringArray ();	// set an array used only for intialization
	void repairCircuit (int);	// method repairing circuit after mutation automatically called from mutation
	qGate *computeMatrix (Individual *);
	void* computeMeasurement(Individual *);
	//for minimizer of quantum circuits
	int findGate (qGate *);	// compares one quantum gate to all available
	void removeRedund (Individual *);	// remove redundantsegments in the circuit
	basic_string < char >synthSegment (int, int);	// creates a new segment of gates in string representation
	void mergeSegs (Individual *, int, int);	// merge two equivalent or equal segments
	void *generateMeasureOps ();
	void *generateMultiMeasureOps ();
	void generateMeasureExpect ();
	void generateExtendMatrix(char);
	int getqGate(string);
};





