/***************************************************************************
 *   Copyright (C) 2006 by martin lukac   *
 *   lukacm@ece.pdx.edu   *
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


#include "EpiG_mpi.h"


#define WORKTAG 1
#define DIETAG 2
/********************************************
* The one and only application object
********************************************/

/*
struct ThreadInfo{
	bool inUse;
	pthread_t info;
	Individual *ind;
	unsigned int timeout; 
};*/

//ThreadInfo thread_pool[MAXNUMBER];

static void master(void);
static void slave(int);
static float* get_next_work_item(void);
//static void process_results(unit_result_t result);
static float compute();//(unit_of_work_t work);

int main(int argc, char **argv)
{
	int myrank;
	int myteton;

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	/* Find out my identity in the default communicator */
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == 0) {
		cout<< "Starting Master \n"<<endl;
		master();
	} else {
		cout<< "Starting Slave \n"<<endl;
		slave(myrank);
	}

	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}

static void
master(void)
{
	int ntasks, rank;
	float *work;
	float result, *results[100];
	MPI_Status status;
	
	/* Find out how many processes there are in the default
	 * communicator */
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	
	/* Seed the slaves; send one unit of work to each slave. */
	
	for (rank = 1; rank < ntasks; ++rank) {
		
		/* Find the next item of work to do */
		work = get_next_work_item();

		/* Send it to each rank */
		MPI_Send(&work,             /* message buffer */
				1,                 /* one data item */
				MPI_INT,           /* data item is an integer */
				rank,              /* destination process rank */
				WORKTAG,           /* user chosen message tag */
				MPI_COMM_WORLD);   /* default communicator */
	}
	/* Loop over getting new work requests until there is no more work
	 *       to be done */
	work = get_next_work_item();

	if (ntasks > 0)
	while (work != NULL) {
	
	/* Receive results from a slave */
//	cout<<"receiving data from slaves"<<endl;
		MPI_Recv(&result,           /* message buffer */
				1,                 /* one data item */
				MPI_DOUBLE,        /* of type float real */
				MPI_ANY_SOURCE,    /* receive from any sender */
				MPI_ANY_TAG,       /* any type of message */
				MPI_COMM_WORLD,    /* default communicator */
				&status);          /* info about the received message */
		if (result > results[status.MPI_SOURCE][0])
			results[status.MPI_SOURCE][0] = result;

//	cout<<"received data from slaves"<<endl;
//	cout<<"sending data to slaves"<<endl;
		/* Send the slave a new work unit */
		MPI_Send(&work,             /* message buffer */
				1,                 /* one data item */
				MPI_INT,           /* data item is an integer */
				status.MPI_SOURCE, /* to who we just received from */
				WORKTAG,           /* user chosen message tag */
				MPI_COMM_WORLD);   /* default communicator */
		/* Get the next unit of work to be done */
		work = get_next_work_item();
//	cout<<"data to slave "<<status.MPI_SOURCE<<" sent"<<endl;
	}
	/* There's no more work to be done, so receive all the outstanding
	 * 			   *      results from the slaves. */
	for (rank = 1; rank < ntasks; ++rank) {
		MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
				MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		results[status.MPI_SOURCE][0] = result;
	}

	/* Tell all the slaves to exit by sending an empty message with the
	 * 			     *      DIETAG. */
	for (rank = 1; rank < ntasks; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}
}


static void 
slave(int rank)
{
	float work[100];
	float result[20];
	int myrank = rank;
	MPI_Status status;
	GA g;
	QE q;
	srand( time(0));
	int generationCondition = 0;
	pthread_t *childs = new pthread_t[MAXNUMBER];
	int died;
	int datadisplaycounter;
	int pipefd[2];
	float results[MAXNUMBER];
	int bestcounter = 0;
	char str[1024], *cp;
	string inFile;
	//the circuits passed form the GA to various optimization and search modules
	string *transmitted[MAXNUMOFGATES];
	//start the GA
	cout<< "INIT Stage 1 DONE: GA Ready GA rank:>>"<<myrank<<"<<"<<endl;
	g.initializePop(myrank);
	cout<< "INIT Stage 1 DONE: GA Ready GA rank:>>"<<myrank<<"<<"<<endl;
	//evaluate the population of the circuits for the first time
	g.evaluatePopulation();
	cout<< "INIT Stage 2 DONE: First evaluation"<<endl;
	//initiate the BestIndividual by the 0th individual from the population
//	g.setIndv(&g.bestIndv, g.population[0]);
	g.outputBest(true);
	cout<< "INIT Stage 3 DONE: Selected the best individual"<<endl;
	//start to iterate the GA with modules
	int counter = 0;
	g.finalResult.counter = counter;
	//life cycle
	datadisplaycounter = 0;
	generationCondition = 0;

	/* Receive a message from the master */
	MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
			MPI_COMM_WORLD, &status);
	//cout<<"data from master received"<<endl;

	while (!g.terminatingCondition())
	{
	//cout<<"receiving data from master"<<endl;
		
		/* Check the tag of the received message. */
		if (status.MPI_TAG == DIETAG) {
			return;
		}
		datadisplaycounter++;
	//cout<<" current generation is "<<generationCondition<<endl;
		g.applyReplication();
		//cout<<"Mutation Operator"<<endl;
		for (int a = 0;a<g.populationNumber;a++)
		{
			//cout<<"plicator"<<endl;
			g.applyIMutation(a, true);
		}
		bestcounter = 0;
	//cout<<"Startin Threads"<<endl;
		for (int i = 0; i < g.populationNumber; i++)
		{
			g.doMeasureFitness();
		}	
	//cout<<"All threads exited gracefully"<<endl;
		g.evaluatePopulation();
		if (counter == g.alterations)
		{
			//generate output of the best individual
			//add exhaustive search
			//conversion from non normal to normal gates must be here.
			g.injectNormalCircuits(q.setgateNormArray(g.getPermutativeSegments(), g.finalGate, g.bestIndv.segmentNumber, g.bestIndv.Error));
			//exit(0);
			counter = 0;
		}
		
		generationCondition++;
		g.setTerminatingCondition(generationCondition);
		if (g.terminatingCondition())
			g.outputBest(false);
		if (datadisplaycounter >= 10){
			g.calcAverage();
			datadisplaycounter = 0;
		}

		result[0] = g.getfitness();
		/* Send the result back */



		//cout<<"Generation Over: "<<result[0]<<endl;
		//g.outputBest(false);
		counter++;
	}
	//	cout<<"sending data to master"<<endl;
		MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//	cout<<"data to master sent"<<endl;
	
	g.closeStream();
	//---------closing the output stream---------------
}

static float* 
get_next_work_item(void)
{
	float result[20];
	for (int a = 0; a < 20; a++)
	result[a] =  13.4+a;
	  /* Fill in with whatever is relevant to obtain a new unit of work
	   *      suitable to be given to a slave. */
	return result;
}


static void 
process_results(float* result)
{
	  /* Fill in with whatever is relevant to process the results returned
	   *      by the slave */
}
