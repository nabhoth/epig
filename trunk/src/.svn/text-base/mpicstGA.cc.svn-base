/***************************************************************************
 *   Copyright (C) 2005-2008 by martin lukac   *
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
 **************************************************************************/
#include "EpiG.h"

/*
 * global static field for thread counting
 */
int GA::populationCounter;

/* 
 * global mutex for our program. assignment initializes it 
 */
pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

/* 
 * global static switch for measurement/nonmeaurement 
 */
int GA::measurement;
/*******************************************
* finds the fittest individual and if found 
* assigns new best individual 
*******************************************/
void GA::evaluatePopulation()
{
	
	Individual *ind;
	int c = 0;
	for (int t = 1; t < populationNumber; t++){
		if (population[t]->fitness > population[c]->fitness){
			c = t;}
	}
	progress = -1;

	if (bestIndv.fitness >= 1)
	{
		progress = 1;
	}
	else
		if (population[c]->fitness > bestIndv.fitness)
		{
				progress = 1;
		}

	if (population[c]->fitness > bestIndv.fitness){
		//cout<<"Id: "<<fileidentifier<<", Old  best: "<< bestIndv.fitness<< "New best: "<<population[c]->fitness<<endl;
		setIndv(&bestIndv, population[c]);
		outputBest(true);
	}
}

/*******************************************
* static - the function called to invoke the program as a thread
*******************************************/
void* GA::startSub(void *arg)
{
	GA* obj = (GA*)arg;
	obj->doMeasureFitness();
	return 0;
}
/*******************************************
* compute the matrix of a given individual 
* - only the matrix and returns it
*******************************************/
qGate* GA::computeMatrix(Individual *indi)
{

	int phasecounter0 = 0;
	int Time = 0;
	int cost = 0;
	complex<float> temphase;
	int begin = 0;
	int end = indi->my_string.find("p", begin+1);
	int c =0;
	qGate myfinalG, temp0, temp1;
	indi->Cost = 0;
	indi->segmentNumber = 0;
	
	
	while(end > 0) 
	{
		if ((end - begin) > 1)
		{
			c++;
			//get the first  gate of the Segment and copy it
			for (int a = 0;a<= numofgates;a++)
				if (indi->gateArray[a]->representation == indi->my_string.at(begin+1)) 
				{
					temp0= copyGate(indi->gateArray[a]);
					if (phase == 1){
						temphase = indi->phases[phasecounter0++];
						for (int k = 0; k<(int(pow((float)2,(float)temp0.numIO)));k++)
                					for (int l = 0; l<(int(pow((float)2,(float)temp0.numIO)));l++)
                     						temp0.gateMatrix1[k][l] *= temphase;
                     			}
					//continue with the cost
					begin++; cost+= temp0.Cost;
					break;
				}
			//get the next one            
			for (int b = begin+1;b<end;b++) 
			{
				for (int a = 0;a<= numofgates;a++)
					if (indi->gateArray[a]->representation == indi->my_string.at(b)) 
					{
						temp1= copyGate(indi->gateArray[a]);
						temphase = indi->phases[phasecounter0++];
						if (phase == 1){
							for (int k = 0; k<(int(pow((float)2,(float)temp0.numIO)));k++)
								for (int l = 0; l<(int(pow((float)2,(float)temp0.numIO)));l++)
									temp1.gateMatrix1[k][l] *= temphase;
						}
						//continue with the cost
						cost+= temp1.Cost;
						a = numofgates*2;
						//multiply the gates
						//using Kronecker multiplication
						temp0 = tensorProduct(&temp0, &temp1);
					}
			}
			if (Time == 0) 
			{
				//if this is the end of the first segment fill the output value
				myfinalG = copyGate(&temp0);
				Time++;
			}
			else 
			{
				//compile warning is intended
				//multiply using normal Manumofgatestrix product the twou segments
				myfinalG = matrixProduct(&myfinalG, &temp0); 
				Time++;
			}
		}
		//move to the next segment
		begin = indi->my_string.find("p", end);
		end = indi->my_string.find("p", begin+1);
	}

 
	//check if phase is defined and multiply the whole matrix by it
	if (phase == 1){
		if (indi->phase != complex <float> (0,0))
		{
			for (int k = 0; k<(int(pow((float)2,(float)myfinalG.numIO)));k++)
	                	for (int l = 0; l<(int(pow((float)2,(float)myfinalG.numIO)));l++)
        		             myfinalG.gateMatrix1[k][l] *= indi->phase;
		}
	}
	//actualize the segment number
	indi->segmentNumber = c;
	indi->Cost = cost;
	//generate the output resulting gate
	qGate* rGate = cGate(&myfinalG);
	//manual transpose
	for (int k = 0; k<(int(pow((float)2,(float)temp0.numIO)));k++){
		for (int l = 0; l<(int(pow((float)2,(float)temp0.numIO)));l++){
			rGate->gateMatrix1[k][l] = myfinalG.gateMatrix1[l][k];
		}
	}
	return rGate;
}
/*******************************************
* generate the fitness value of the Individual 
* it is set up for
*******************************************/
void* GA::doFitness()
{
	
    	//two dummy gates and one pointer for operations
	qGate* myfinalG;
	Individual *indi;
	int ind;
	int errorcounter = 0;
	int rc = pthread_mutex_lock(&mtx);
	ind = GA::populationCounter;
	GA::populationCounter = (GA::populationCounter+1)%populationNumber;
	rc = pthread_mutex_unlock(&mtx);
	indi = population[ind];
	//null this Individual previous evaluation
	indi->Error = 0;
	indi->fitness = 0;
	
	//threaded code	
	rc = pthread_mutex_lock(&mtx);
		if (rc) {
    		cout<<"Matrix locked: "<<ind<<endl;
		}
		//cout<<"Matrix starting: "<<ind<<endl;
		myfinalG = computeMatrix(indi);
		//cout<<"Matrix done: "<<ind<<endl;
		rc = pthread_mutex_unlock(&mtx);
		//set some more parameters
		myfinalG->numIO = indi->ioNumber;
		//cout<<"From make fitness: "<<indi->segmentNumber <<endl;
		//check if phase is defined and multiply the whole matrix by it
		if (phase == 1){
		if (indi->phase != complex <float> (0,0))
		{
				for (int k = 0; k<(int(pow((float)2,(float)myfinalG->numIO)));k++)
                	for (int l = 0; l<(int(pow((float)2,(float)myfinalG->numIO)));l++)
                     myfinalG->gateMatrix1[k][l] *= indi->phase;
		}
		}
		rc = pthread_mutex_lock(&mtx);
		if (rc) {
    		//cout<<"Result locked: "<<ind<<endl;
		}
        //direct correspondence between number of wires in the goal and current individual
        if (myfinalG->numIO!= finalGate.numIO) 
        {
                indi->fitness = 0;
        }
        else 
        {
			errorcounter = 0;
			
                for (int a = 0; a<(int(pow((float)2,(float)myfinalG->numIO)));a++)
                {
                        for (int b = 0; b<(int(pow((float)2,(float)myfinalG->numIO)));b++)
                        { 
                                //simple comparison -- allowing don't cares as value "10"
                                if (finalGate.gateMatrix1[a][b] != complex<float>(10, 0)
									&& finalGate.gateMatrix1[a][b] != complex<float>(0,0))
								{
									errorcounter++;								
									indi->Error += abs(finalGate.gateMatrix1[a][b]*conj(finalGate.gateMatrix1[b][a]) - myfinalG->gateMatrix1[a][b]*conj(myfinalG->gateMatrix1[b][a]));
								}
                         }
                }
        }

        //normalizing error
	if (errorcounter != 0)
        indi->Error =  indi->Error/errorcounter;
        indi->Cost = (exp(-pow((divider-indi->Cost),2)));
	rc = pthread_mutex_unlock(&mtx);

	//generate fitness value
        if (indi->Error != 0)
		{
			switch(replicator){
					case 0:
						//simple fitness1
						indi->fitness = (1-indi->Error);
						break;
					case 1:
						//simple fitness1
						indi->fitness = (1/(indi->Error+1));
						break;
					case 2:
						//scaled complex fitness1
					    indi->fitness = (alpha*(1-indi->Error) + beta*indi->Cost);
						break;
					case 3:
						//scaled complex fitness2
						indi->fitness = (alpha*(1/(indi->Error+1)) + beta*indi->Cost);
						break;
					case 4:
						indi->fitness = (1-indi->Error);
						break;
			}
		}else 
		{
			indi->Error = 0;
			indi->fitness = 1;	//we found a good individual;
		}
		delete(myfinalG);
		return 0;
}

/*******************************************
* measures the individual and generates the expectations coefficients
* for all possible states
*
* calculates the probabilities for obtaining |0> and |1> for all the measurement operators at a time 
* Right now it is doing measurement in increasing order of all defined bits to measure
*******************************************/
void* GA::doMeasureFitness()
{
	
	//define temporal variables for the results and the manupaltion
	//of the measurement
	int l = (int(pow((float)2,(float)finalGate.numIO)));
	int x = (int(pow((float)2,(float)finalGate.numIO)));;
	int y = (int(pow((float)2,(float)finalGate.numIO)));;
	int p = 0;
	int mes = 0;
	int ind;
	complex<float> inter0[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> inter1[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> logicresult[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> sta[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> expectation0, expectation1, alphas, betas;
	complex<float> resultingTotalError[(int(pow((float)2,(float)finalGate.numIO)))][measurement*2];
	complex<float> resultingTotal[(int(pow((float)2,(float)finalGate.numIO)))][measurement*2];
	complex<float> inter, expect0, expect1;
	qGate *measure0, *measure1,*myfinalG, measure;
	Individual *indi;

	//threaded static counter
	int rc = pthread_mutex_lock(&mtx);
	ind = GA::populationCounter;
	GA::populationCounter = (GA::populationCounter+1)%populationNumber;
	rc = pthread_mutex_unlock(&mtx);
	
	indi = population[ind];
	//null this Individual previous evaluation
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;
	
	//threaded code	
	rc = pthread_mutex_lock(&mtx);
	if (rc) {
    	cout<<"Matrix locked: "<<ind<<endl;
	}
	// compute hte matrix of the circuit
	myfinalG = computeMatrix(indi);
	if (phase == 1){
		if (indi->phase != complex <float> (0,0))
		{
		for (int k = 0; k<(int(pow((float)2,(float)myfinalG->numIO)));k++)
               	  for (int l = 0; l<(int(pow((float)2,(float)myfinalG->numIO)));l++)
                    myfinalG->gateMatrix1[k][l] *= indi->phase;
		}
	}

	//propagate each input through the circuit
	for (int k = 0; k<l; k++){
		
		for (int i = 0; i < l; i++){
			inter0[i] = complex<float>(0,0);
			inter1[i] = complex<float>(0,0);
		}

		//initialize variables
		for (int w = 0; w < l;w++)
			sta[w] = inter0[w] = inter1[w] = logicresult[w] = complex<float>(0,0);
		//set the orthonormal state - k
		sta[k] = complex<float>(1,0);
				
		//propagate the state through the matrix
		for (int i = 0; i < l; i++)
		{    
			for (int j = 0; j < l; j++)
			{
				logicresult[i] += sta[j]*myfinalG->gateMatrix1[i][j];
			}
		}

		//init the measurement operator with respect to the desired output
		//for each measured ouptut state get the expectation values 
		expect0 = expectationsAllState[k][0];
		expect1 = expectationsAllState[k][1];

		//apply the measurement operator for the desired state 
		for (int i = 0; i < l; i++)
		{    
			//rows
			for (int j = 0; j < l; j++)
			{   
				inter0[i] += logicresult[j]*measurementsAllState[2*k]->gateMatrix1[i][j];
				inter1[i] += logicresult[j]*measurementsAllState[2*k+1]->gateMatrix1[i][j];
			}
		}
//p(0) - not really p(0) but rather the desired result ---------------------------------
		expectation0 = 0;
		alphas = complex<float>(0,0);
		//the finnal inner product for p(0)
		for (int j = 0; j < l; j++)
			expectation0 += (getConjugate(logicresult[j])*inter0[j]);
			//state after the measurement for 0
		for (int i = 0; i < l; i++){
			inter0[i] = inter0[i]/sqrt(expectation0);
			alphas +=  inter0[i];
		}
//p(1) ---------------------------------------------------------------------------------
		//vector inner product
		expectation1 = 0;
		betas = complex<float>(0,0);
		//the finnal inner product for p(1)
		for (int i = 0; i < l; i++)
			expectation1 += getConjugate(logicresult[i])*inter1[i];
		//state after the measurement for 1
		for (int i = 0; i < l; i++){
			inter1[i] = inter1[i]/sqrt(expectation1);
			betas +=  inter1[i];
		}
//--------------------------------------------------------------------------------------
		alphas = expectation0;
		betas = expectation1;

//		cout<<"alpha: "<<alphas<<" + beta: "<<betas<<endl;
		//calculate the total state 
		//Total State = M(0)+M(1)State/Measured(0)+Measured(1)
		mes = 0;
		resultingTotal[k][2*mes] = expectation0;
		resultingTotal[k][2*mes+1] = expectation1;
		if (measureexpected[k][2*mes] == complex<float>(0,1) && measureexpected[k][2*mes+1] == complex<float>(0,1)){
			//skip don't cares
			resultingTotalError[k][2*mes] = complex<float>(0,0);
			resultingTotalError[k][2*mes+1] = complex<float>(0,0);
		} else {
			if (real(expectationsAllState[k][mes]) == real(expectationsAllState[k][mes+1])){
				if (real(expectationsAllState[k][mes]) == 1 || real(expectationsAllState[k][mes]) == 0){
					if (real(alphas+betas) == 1 && (real(alphas) == real(betas))){
						resultingTotalError[k][mes*2] = 0;
						resultingTotalError[k][mes*2+1] = 0;
					} else {
							resultingTotalError[k][mes*2] = abs(complex<float>(0.5) - abs(expectation0));
							resultingTotalError[k][mes*2+1] = abs(complex<float>(0.5) - abs(expectation1));
					}
				} else {
					resultingTotalError[k][mes*2] = abs(abs(expect0) - abs(expectation0));
					resultingTotalError[k][mes*2+1] = abs(abs(expect1) - abs(expectation1));
				}
			} else{
				resultingTotalError[k][mes*2] = abs(abs(expect0) - abs(expectation0));
				resultingTotalError[k][mes*2+1] = abs(abs(expect1) - abs(expectation1));
			}

		}
	}

	int m = 0;
	inter = (0,0);
	for (int e = 0; e < l;e++)
	{
//		for (int c = 0; c < measurement; c++){
			if (!(measureexpected[0][e] == complex<float>(0,1) && measureexpected[1][e] == complex<float>(0,1)))
			{
				//get the error over both possible measurements of the state - for 0 and 1 of the first qubit
				inter += (resultingTotalError[e][0] + resultingTotalError[e][1])/complex<float>(1,0);
				//expecting the desired higher probabilities 1
				m++;
			}
//		}
	}
	//indi->Error /= measurement;
	indi->Error = real(inter/complex<float>(m,0));
        indi->Cost = (exp(-pow((divider-indi->Cost),2)));

	//generate fitness value
        if (indi->Error != 0)
	{
		switch(replicator){
			case 0:
				//simple fitness1
				indi->fitness = (1-indi->Error);
				break;
			case 1:
			//simple fitness1
				indi->fitness = (1/(indi->Error+1));
				break;
			case 2:
				//scaled complex fitness1
			    indi->fitness = (alpha*(1-indi->Error) + beta*indi->Cost);
				break;
			case 3:
				//scaled complex fitness2
				indi->fitness = (alpha*(1/(indi->Error+1)) + beta*indi->Cost);
				break;
			case 4:
				indi->fitness = (1-indi->Error);
				break;
		}
	}else 
	{
		indi->Error = 0;
		indi->fitness = 1;//we found a good individual;
	}
	rc = pthread_mutex_unlock(&mtx);
//	cout<<"Fitness: "<<indi->fitness<<endl;
	delete(myfinalG);
}

/*******************************************
 * generate the fitness value of the Individual based on measurement results
 *******************************************/

void GA::makeMeasureFitness(Individual *indi, bool display)
{
	//define temporal variables for the results and the manupaltion
	//of the measurement
	int l = (int(pow((float)2,(float)finalGate.numIO)));
	int x = (int(pow((float)2,(float)finalGate.numIO)));;
	int y = (int(pow((float)2,(float)finalGate.numIO)));;
	int p = 0;
	int m = 0;
	int cost = 0;
	int ind;
	complex<float> inter0[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> inter1[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> logicresult[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> sta[(int(pow((float)2,(float)finalGate.numIO)))];
	complex<float> expectation0, expectation1, alphas, betas;
	complex<float> resultingTotalError[(int(pow((float)2,(float)finalGate.numIO)))][measurement*2];
	complex<float> resultingTotal[(int(pow((float)2,(float)finalGate.numIO)))][measurement*2];
	complex<float> inter, expect0, expect1;
	qGate *measure0, *measure1,*myfinalG, measure;

	//cout<<"got individual: "<<indi<<endl;
	//null this Individual previous evaluation
	indi->Error = 0;
	indi->fitness = 0;
	
	//get the matrix of the circuit
	myfinalG = computeMatrix(indi);
	//set some more parameters
	myfinalG->numIO = indi->ioNumber;
	//check if phase is defined and multiply the whole matrix by it
	if (phase == 1){
		if (indi->phase != complex <float> (0,0))
		{
		for (int k = 0; k<(int(pow((float)2,(float)myfinalG->numIO)));k++)
               	  for (int l = 0; l<(int(pow((float)2,(float)myfinalG->numIO)));l++)
                    myfinalG->gateMatrix1[k][l] *= indi->phase;
		}
	}

	for (int g = 0; g<l; g++)
		for (int h = 0; h<measurement; h++)
	//for every state in the logic state vector
	for (int k = 0; k<l; k++){

		for (int i = 0; i < l; i++){
			inter0[i] = complex<float>(0,0);
			inter1[i] = complex<float>(0,0);
		}

		//set the input state
		for (int w = 0; w < l;w++)
			sta[w] = inter0[w] = inter1[w] = logicresult[w] = complex<float>(0,0);
		sta[k] = complex<float>(1,0);
				
		//propagate the state through the matrix
		for (int i = 0; i < l; i++)
		{    
			for (int j = 0; j < l; j++)
			{
				logicresult[i] += sta[j]*myfinalG->gateMatrix1[i][j];
			}
		}
			
		//init the measurement operator with respect to the desired output
		//for each measured ouptut state get the expectation values 
		expect0 = complex<float>(0,0);
		expect1 = complex<float>(0,0);
		expect0 = expectationsAllState[k][0];
		expect1 = expectationsAllState[k][1];

		//vector product 
		for (int i = 0; i < l; i++)
		{    
			//rows
			for (int j = 0; j < l; j++)
			{   
				inter0[i] += logicresult[j]*measurementsAllState[2*k]->gateMatrix1[i][j];
				inter1[i] += logicresult[j]*measurementsAllState[2*k+1]->gateMatrix1[i][j];
			}
		}
//p(0) - not really p(0) but rather the desired result ------------------------------------------------------------------------------
		expectation0 = 0;
		alphas = complex<float>(0,0);
		//the finnal inner product for p(0)
		for (int j = 0; j < l; j++)
			expectation0 += (getConjugate(logicresult[j])*inter0[j]);
			//state after the measurement for 0
		for (int i = 0; i < l; i++){
			inter0[i] = inter0[i]/sqrt(expectation0);
			alphas +=  inter0[i];
		}
//p(1) ---------------------------------------------------------------------------------
		//vector inner product
		expectation1 = 0;
		betas = complex<float>(0,0);
		//the finnal inner product for p(1)
		for (int i = 0; i < l; i++)
			expectation1 += getConjugate(logicresult[i])*inter1[i];
		//state after the measurement for 1
		for (int i = 0; i < l; i++){
			inter1[i] = inter1[i]/sqrt(expectation1);
			betas +=  inter1[i];
		}
//--------------------------------------------------------------------------------------
		alphas = expectation0;
		betas = expectation1;

//		cout<<"Expected: "<<expect0<<" + "<<expect1<<", Obtained: "<<expectation0<<" + "<<expectation1<<endl;
//		cout<<"alpha: "<<alphas<<" + beta: "<<betas<<endl;
		//calculate the total state 
		//Total State = M(0)+M(1)State/Measured(0)+Measured(1)

		resultingTotal[k][2*m] = expectation0;
		resultingTotal[k][2*m+1] = expectation1;
		if (measureexpected[k][2*m] == complex<float>(0,1) && measureexpected[k][2*m+1] == complex<float>(0,1)){
			//skip don't cares
			resultingTotalError[k][2*m] = complex<float>(0,0);
			resultingTotalError[k][2*m+1] = complex<float>(0,0);
		} else {
			if (real(expectationsAllState[k][m]) == real(expectationsAllState[k][m+1])){
				if (real(expectationsAllState[k][m]) == 1 || real(expectationsAllState[k][m]) == 0){
					if (real(alphas+betas) == 1 && (real(alphas) == real(betas))){
						resultingTotalError[k][m*2] = 0;
						resultingTotalError[k][m*2+1] = 0;
					} else {
							resultingTotalError[k][m*2] = abs(complex<float>(0.5) - abs(expectation0));
							resultingTotalError[k][m*2+1] = abs(complex<float>(0.5) - abs(expectation1));
					}
				} else {
					resultingTotalError[k][m*2] = abs(abs(expect0) - abs(expectation0));
					resultingTotalError[k][m*2+1] = abs(abs(expect1) - abs(expectation1));
				}
			} else{
				resultingTotalError[k][m*2] = abs(abs(expect0) - abs(expectation0));
				resultingTotalError[k][m*2+1] = abs(abs(expect1) - abs(expectation1));
			}

		}
	}

	m = 0;
	inter = (0.0);
	for (int e = 0; e < l;e++)
	{
		if (!(measureexpected[e][0] == complex<float>(0,1) && measureexpected[e][1] == complex<float>(0,1)))
		{	
			inter += (resultingTotalError[e][0] + resultingTotalError[e][1])/complex<float>(1,0);
			m++;
		}  
	}
	indi->Error = real(inter/complex<float>(m,0));
        indi->Cost = (exp(-pow((divider-indi->Cost),2)));

	//generate fitness value
        if (indi->Error != 0)
	{
		switch(replicator){
			case 0:
				//simple fitness1
				indi->fitness = (1-indi->Error);
				break;
			case 1:
			//simple fitness1
				indi->fitness = (1/exp(indi->Error));
				break;
			case 2:
				//scaled complex fitness1
			    indi->fitness = (alpha*(1-indi->Error) + beta*indi->Cost);
				break;
			case 3:
				//scaled complex fitness2
				indi->fitness = (alpha*(1/exp(indi->Error)) + beta*indi->Cost);
				break;
			case 4:
				indi->fitness = (1-indi->Error);
				break;
		}
	}else 
	{
		indi->Error = 0;
		indi->fitness = 1;//we found a good individual;
	}

	/*
	cout<<"Fitness: "<< indi->fitness<<endl;
	cout<<"Error: "<< indi->Error<<endl;
	cout<<"Matrix: " <<endl;
	for (int a = 0; a<(int(pow((float)2,(float)myfinalG->numIO)));a++)
	{
		for (int b = 0; b<(int(pow((float)2,(float)myfinalG->numIO)));b++)
		{
			cout<< myfinalG->gateMatrix1[a][b]<<"";
		}
		cout<<endl;
	}
	cout<<endl;
	cout << " -------- " <<endl;
	*/
         if (display)
        {
		if (displaylevel > 0){
	    		out_stream << " -------- " <<endl;
	    		out_stream <<"Representation String: "<<indi->my_string<<endl;
	    		out_stream << " -------- " <<endl;
			if (phase == 1){
				out_stream <<"Phase for valid gates: "<<endl;
				for(int g = 0; g< getGatesfromString(indi->my_string); g++)
					out_stream <<indi->phases[g]<<", "<<endl;
				out_stream<<endl;
			}
		}else {
			out_stream <<"Representation String: "<<indi->my_string<<endl;
		}
    
		out_stream <<"Fitness: "<< indi->fitness<<endl;
    		out_stream <<"Error: "<< indi->Error<<endl;
    		out_stream <<"Wires: "<< indi->ioNumber<<endl;
		if (displaylevel == 2){
			if (phase == 1)
				out_stream <<"Phase: "<< indi->phase<<endl<<endl;
			out_stream <<"Matrix: " <<endl;
			for (int a = 0; a<(int(pow((float)2,(float)myfinalG->numIO)));a++)
			{
				for (int b = 0; b<(int(pow((float)2,(float)myfinalG->numIO)));b++)
				{
					out_stream<< myfinalG->gateMatrix1[a][b]<<"";
				}
				out_stream<<endl;
			}
			out_stream<<endl;
			out_stream << " -------- " <<endl;
		} else if (displaylevel == 1){
			if (phase == 1)
				out_stream <<"Phase: "<< indi->phase<<endl<<endl;
		}
		if (displaylevel > 0){
			if (measurement != 0){
				out_stream<<"Measured results:  "<<endl;
				for(int u = 0; u < measurement;u++){
					out_stream<<"Measured Qubit is "<<u<<endl;
					m = l;
					inter = (0.0);
					for (int e = 0; e < l;e++)
					{
						if (!(measureexpected[e][2*u] == complex<float>(0,1) && measureexpected[e][2*u+1] == complex<float>(0,1)))
							inter += real(resultingTotalError[e][u*2]+resultingTotalError[e][u*2+1])/2;
						else
							m--;
					}
					inter = inter/complex<float>(m,0);
					out_stream<<"Error 0/1: "<<inter <<endl;
					for (int h = 0; h < (int(pow((float)2,(float)myfinalG->numIO))); h++){
						out_stream<< resultingTotalError[h][u*2]<<"::"<< resultingTotalError[h][u*2+1]<<"   output:  "<<resultingTotal[h][u*2]<<"        "<< resultingTotal[h][u*2+1]<<"   original:  "<<measureexpected[h][u*2]<<"        "<< measureexpected[h][u*2+1]<<endl;
					}
				}
			}
			out_stream << " -------- " <<endl;
		}
        }
	delete(myfinalG);

}

/*******************************************
* the String array represents the gates by the number of wires
* and present a list of gates for each amount of wires
*******************************************/
void GA::setStringArray()
{
	qGate *Q;
	int count[MAXNUMOFGATES];
	
	for (int j = 0; j < MAXGATEINPUT; j++)
	{
		count[j] = 0;
		for (int i = 0; i < numofgates; i++)
		{
			Q = gateArray[i];
			if (Q->numIO == (j+1))
			{
				stringArray[j].nameGate[count[j]] = Q->representation;
				count[j]++;
			}
			stringArray[j].numberGates = count[j];
		}
	}
}

/*******************************************
* this generates the String repressenting every individual.
*it is called upon the structure of the initial population
*******************************************/
void GA::makeIndvString(int i)
{
    int tempNumber,tempNumber2, tempGate, counter, temporyhold;
	char gate;
	Individual *I = population[i];
	counter = 0;
        srand(time(NULL));
	I->my_string = string("");
	for (int j = 0; j < I->segmentNumber; j++) // for each segment of this individual
	{
		if (I->ioNumber > MAXGATEINPUT) // if the input has more than MAXGATEINPUT (10)
			tempGate = MAXGATEINPUT;             // max gate to use
		else
			tempGate = I->ioNumber;
		tempNumber = I->ioNumber; // the number of inputs
		do
		{
			tempNumber2 = rand()%tempGate; // which gate to use
		}while (stringArray[tempNumber2].numberGates == 0);
		
		if (I->ioNumber > 1) // if not only one input
		{
			I->my_string += 'p';
			counter++;
		}
		
		while ((tempNumber-(tempNumber2+1)) > 0) // while more gate could be parallel
		{
			temporyhold = rand()%stringArray[tempNumber2].numberGates;
			gate = stringArray[tempNumber2].nameGate[temporyhold];
			
			I->my_string += gate;
			counter++;
			
			tempNumber = tempNumber - (tempNumber2+1);
			if (tempNumber > tempGate)              // use the smaller one of the two
			{
				temporyhold = rand()%stringArray[tempNumber2].numberGates;
				gate = stringArray[tempNumber2].nameGate[temporyhold];
				do
				{
					tempNumber2 = rand()%(tempGate);
				}while (stringArray[tempNumber2].numberGates == 0);
			}
			else
			{
				do
				{
					tempNumber2 = rand()%(tempNumber);
				}while (stringArray[tempNumber2].numberGates == 0);
			}
		}
		temporyhold = rand()%stringArray[tempNumber2].numberGates;
		gate = stringArray[tempNumber2].nameGate[temporyhold]; // when only can choose one input gate
		
		I->my_string += gate;
		counter++;
		
		if (I->ioNumber > 1)
		{
			I->my_string += 'p';
			counter++;
		}
	}
}


/*******************************************
* this is the main ENd output function
* when the parameter individual = -1;the whole Input/Output Settings, 
* GA settings are dumped into file
* when the parameter individual != -1; The given individual is printed out 
*******************************************/
void GA::output(int individual)
{
	qGate *qG;

	if (individual != -1){
		
      out_stream << "---------------------------" << endl;
      out_stream << " End of run at " <<generationCondition<< endl;
      out_stream << "---------------------------" << endl;
      out_stream << "----- Individual #"<<individual<<"-----" << endl;
 	qG = computeMatrix(population[individual]);
	//get the matrix of the circuit
	out_stream <<"Matrix: " <<endl;
	for (int a = 0; a<(int(pow((float)2,(float)qG->numIO)));a++)
	{
		for (int b = 0; b<(int(pow((float)2,(float)qG->numIO)));b++)
		{
			out_stream<< qG->gateMatrix1[a][b]<<"";
		}
		out_stream<<endl;
	}
	out_stream<<endl;
	out_stream << " -------- " <<endl;

		
	} else {
		
		 out_stream << "There are " << numofgates << " gates from the input file and synthesized.\n";
		 for (int i = 0; i < numofgates; i++)
		 {
			 qG = gateArray[i];
			 out_stream << "The gate number " << i << " is: " << endl;
			 out_stream << "And the gate representation is: " << gateArray[i]->representation << endl;
			 out_stream << "Number of inputs is: " << qG->numIO << endl;
			 out_stream << "Parent A is: " << qG->parentA << endl;
			 out_stream << "Parent B is: " << qG->parentB << endl;
	
			 for (int j = 0; j < (int(pow((float)2,(float)qG->numIO))); j++){
				 for (int k = 0; k < (int(pow((float)2,(float)qG->numIO))); k++){
					 out_stream << qG->gateMatrix1[j][k] << " ";
				 }
				 out_stream << endl;
			 }
		 }
		 out_stream << "The Target gate(s): \n";
		 out_stream << "Number of inputs is: " << finalGate.numIO << endl;
		 for (int m = 0; m < (int(pow((float)2,(float)resultnum))); m++){
			 for (int n = 0; n < (int(pow((float)2,(float)resultnum))); n++){
				 out_stream << finalGate.gateMatrix1[n][m] << " ";
			 }
			 out_stream << endl;
		 }
	 }
	this->display = false;
}


/*******************************************
* outputs the Best Invididual
*******************************************/
void GA::outputBest(bool cond)
{
	//cout<<"Checking for Best individual "<<progress<<"   "<<bestIndv.fitness<<endl;
	qGate *myfinalG;
	if (bestIndv.fitness >= 1 || progress > 0){

		if (displaylevel == 2){
			out_stream << "---------------------------" << endl;
			out_stream << " End of run at " <<generationCondition<< endl;
			out_stream << "---------------------------" << endl;
			out_stream << "----- Best Individual -----" << endl;
		} else if (displaylevel == 1){
			out_stream << "Generation: "<<generationCondition<< endl;
			out_stream << "----- Best Individual -----" << endl;
		} else {
			out_stream << "Generation: "<<generationCondition<< endl;
		}
		//get the matrix of the circuit
		myfinalG = computeMatrix(&bestIndv);
		out_stream <<"Matrix: " <<endl;
		for (int a = 0; a<(int(pow((float)2,(float)myfinalG->numIO)));a++)
		{
			for (int b = 0; b<(int(pow((float)2,(float)myfinalG->numIO)));b++)
			{
				out_stream<< myfinalG->gateMatrix1[a][b]<<"";
			}
			out_stream<<endl;
		}
		out_stream<<endl;
		out_stream<< bestIndv.my_string<<" <- Representation "<<endl;
		out_stream<< bestIndv.fitness<<" <- Fitness "<<endl;
		out_stream<< bestIndv.Error<<" <- Error "<<endl;
		out_stream << " -------- " <<endl;myfinalG = computeMatrix(&bestIndv);
		

		if (displaylevel == 2){
			out_stream << bestIndv.ioNumber<<" <- number of inputs "<<endl;
			out_stream<< bestIndv.my_string<<" <- Representation "<<endl;
			out_stream<< bestIndv.fitness<<" <- Fitness "<<endl;
			out_stream<< bestIndv.Error<<" <- Error "<<endl;
			out_stream<< bestIndv.Cost<<" <- Cost "<<endl;
			out_stream<<"---------------------------" << endl;
			out_stream<<" Data average: " << endl;
			out_stream<<" Fitness " << " Error "<<" Cost "<<endl;
			for (int f = 0; f< finalResult.counter; f++){
				out_stream<<finalResult.avFitness[f]<<" "<<finalResult.avError[f]<<" "<<finalResult.avCost[f]<<endl;
			}
		} else if (displaylevel == 1){
			out_stream<<" Data average: " << endl;
			out_stream<<" Fitness " << " Error "<<" Cost "<<endl;   
			for (int f = 0; f< finalResult.counter; f++){ 
				out_stream<<finalResult.avFitness[f]<<" "<<finalResult.avError[f]<<" "<<finalResult.avCost[f]<<endl;
		       	}
		} else {

		}
	} else {
		if (cond){
			 //get the matrix of the circuit
 			myfinalG = computeMatrix(&bestIndv);
			//get the matrix of the circuit
			out_stream <<"Matrix: " <<endl;
			for (int a = 0; a<(int(pow((float)2,(float)myfinalG->numIO)));a++)
			{
				for (int b = 0; b<(int(pow((float)2,(float)myfinalG->numIO)));b++)
				{
					out_stream<< myfinalG->gateMatrix1[a][b]<<"";
				}
				out_stream<<endl;
			}
			out_stream<<endl;
			out_stream << " -------- " <<endl;myfinalG = computeMatrix(&bestIndv);
 
                        out_stream << bestIndv.ioNumber<<" <- number of inputs "<<endl;
                        out_stream<< bestIndv.my_string<<" <- Representation "<<endl;
                        out_stream<< bestIndv.fitness<<" <- Fitness "<<endl;
                        out_stream<< bestIndv.Error<<" <- Error "<<endl;
                        out_stream<< bestIndv.Cost<<" <- Cost "<<endl;
                        out_stream<<"---------------------------" << endl;
                        out_stream<<" Data average: " << endl;
                        out_stream<<" Fitness " << " Error "<<" Cost "<<endl;
                        for (int f = 0; f< finalResult.counter; f++){
                                out_stream<<finalResult.avFitness[f]<<" "<<finalResult.avError[f]<<" "<<finalResult.avCost[f]<<endl;

			}
		}
	}
}
/*******************************************
* Restrictions:
* 1. no fan out - ex, 5 inputs has to have 5 outputs
* 2. maximum 10 I/O (geometric increase of complexity with more I/Os)
*******************************************/
void GA::initializePop(int rank, string filename)
{
	string inFile = filename;
	string r;
	string outFile;
	r += char(rank+'0');
	string post  = (".out.mach."+r);
	//cout<<"post: "<<post<<" r: "<<r<<endl;
	time_t rawtime;
	string oFile = "";
	struct tm * timeinfo;
	condition= false;
	in_stream.open(inFile.c_str());
	if (in_stream.fail())
	{
		cout << "\"Daa set \" file opening failed.\n";
		exit(1);
	}

	outFile = inFile;
	
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	oFile = asctime (timeinfo);

	unsigned int t = 0;
	while(1){
		//remove spaces
		if (t < oFile.length())
		{
			if(oFile[t] >= '0' && oFile[t] <= 'z'){
				outFile += oFile[t];
 			} else {
				if(oFile[t] == ' '){
				}else {
					break;
				}
			}
		} else break;
		t++;
	}
	
	outFile += post;
	fileidentifier = "";
	for (int h = 0; h < outFile.length();h++)
		 fileidentifier += outFile.at(h);

	out_stream.open(outFile.c_str());
	if (out_stream.fail())
	{
		cout << "Output file opening failed.\n";
        	exit (1);
	}
	cout<<"Generating the array of available gates"<<endl;
	 //initialize the gates
     	initializeGates();
	generateMeasureExpect();
	generateMeasureOps();
	 //initialize the String array representing the gates
	 cout<<"Generating intial individual strings"<<endl;
	setStringArray();
	 cout<<"Generating intial popualtion of circuits"<<endl;
	 //generate the initial population/new_population of individuals
	for (int i = 0; i < populationNumber; i++)
	{
		 population[i] = new Individual;
		 new_population[i] = new Individual;
		//all individuals have the width of the output		 
		 population[i]->ioNumber = finalGate.numIO;



	for (int a = 0; a < numofgates; a++)
	{
		population[i]->gateArray[a] = new qGate;
		//cout<<gateArray[a]<<endl;
		
		population[i]->gateArray[a]->representation = gateArray[a]->representation ; // starting from A
		population[i]->gateArray[a]->numIO = gateArray[a]->numIO ; // starting from A
		for (int j = 0; j < (int(pow((float)2,(float)gateArray[a]->numIO))); j++)
		{
			for (int k = 0; k < (int(pow((float)2,(float)gateArray[a]->numIO))); k++)
			{
				population[i]->gateArray[a]->gateMatrix1[j][k] = complex<float>(gateArray[a]->gateMatrix1[j][k]);
			}
		}
	}

         //population[i]->segmentNumber = rand()%(2*segments)+5; 
         population[i]->segmentNumber = segments; 
		 //cout<<population[i]->segmentNumber<<endl;
		 if (phase == 1){
		 	population[i]->phase  = generatePhase();
		 }
		 //population[i]->segmentNumber = rand()%MAXSEGNUM+1;
		 //generate the String representing this circuit
		 makeIndvString(i);
		 //generate phases for every single gate
		 if (phase == 1){
		 	generatePhases(population[i]);
		 }
		if (measurement > 0)
			makeMeasureFitness(population[i],false);
		else
		 makeFitness(population[i],false);
		 //dump to the output file the intialized individual
		//out_stream<<population[i]->my_string<<endl;
		 
	}

	//init the best individual
	for (int a = 0; a < numofgates; a++)
	{
		bestIndv.gateArray[a] = new qGate;
		bestIndv.gateArray[a]->representation = gateArray[a]->representation; // starting from A
		bestIndv.gateArray[a]->numIO = gateArray[a]->numIO; // starting from A
		for (int j = 0; j < (int(pow((float)2,(float)gateArray[a]->numIO))); j++)
		{
			for (int k = 0; k < (int(pow((float)2,(float)gateArray[a]->numIO))); k++)
			{
				bestIndv.gateArray[a]->gateMatrix1[j][k] = complex<float>(gateArray[a]->gateMatrix1[j][k]);
			}
		}
	}
    out_stream << "\nInitialization done" << endl;
}

/*******************************************
* helper function
*******************************************/
int  GA::getGate(char b)
{
	for (int a = 0;a<numofgates;a++)
		if (gateArray[a]->representation == b) return a;
  	return '0';
}


/*******************************************
* helper function
*******************************************/
int  GA::getqGate(string b)
{
	for (int a = 0;a<numofgates;a++){
		if (b.length() == gateArray[a]->my_string.length()){
			//cout<<a<<"; "<<b<<", "<<gateArray[a]->my_string<<"   "<<(gateArray[a]->my_string).compare(0, b.length(), b) <<endl;
			if ((gateArray[a]->my_string).compare(0, gateArray[a]->my_string.length(), b) >= 0) return a;
		}
	}
  	return -1;
}
/*******************************************
* apply the mutation on the individual of the current population
* take one gate and change it to something else
*******************************************/

void GA::applyRMutation(int indv, bool repair)
{

	Individual *I = population[indv];
	
        srand( time(NULL));
	unsigned int proba, proba2;
	char gate;
	qGate temp;
	proba  = (int)(1-(float(1.0*rand()/((float)RAND_MAX+1.0))));
		
	if (proba <= proba_Mutation)
	{
		proba = rand()%1;
		//change only the phase of the circuit
		if (phase == 1){
		if (proba == 0){
			I->phase = generatePhase();
				return;	
		}
		}
		do {
			proba2 = (rand()%I->my_string.length())-1;
			if (proba2 < 0) proba2 = 0;
				gate =  I->my_string.at(proba2);
		}while(gate != 'p' || proba2>(I->my_string.length() -3));
					
		if (I->my_string.at(proba2+1) == 'p') proba2++;
		int proba3 = I->my_string.find("p", proba2+1);

		I->my_string.erase(proba2, (proba3-proba2+1));

		//tricky part is to keep the width of the circuit constant
		basic_string<char> tempS;
		int V = 0;
		do {
			temp = *gateArray[rand()%numofgates];
			if (temp.numIO <= (I->ioNumber-V) && temp.numIO > 0)
			{
				V += temp.numIO;
				tempS.insert(tempS.length(), 1, temp.representation);
			}
		}while(V != I->ioNumber);  
		tempS.insert(0, "p");
		tempS.insert(tempS.length(), "p");
				

			
		I->my_string.insert(proba2, tempS);
				
		if (repair)
			repairCircuit(indv);			
	}//add here parts moving gate location in the circuit -- this must be done in order to do circuits with controls such 
	 //such as I41_zz do max for circuits with 5 variables.
	//cout <<"mutated string: "<<I->my_string<<endl;
}

/*******************************************
 * apply the mutation on the individual of the current population
 * and for every segment apply the mutation operation
 *******************************************/
void GA::applyIMutation(int indv, bool repair)
{

	Individual *I = population[indv];
	int proba, proba2,temp0, iterations;
	char gate;
	qGate temp;
	complex <float> tempc;
	
        srand( time(NULL));
	temp0 = 0;
	for (unsigned int m = 0; m< I->my_string.size(); m++){
		if (I->my_string.at(m) == 'p')
			temp0++;
	}

	I->segmentNumber = temp0;
	iterations = getGatesfromString(I->my_string);
	//gate wise mutation
	for(int a= 0; a< iterations; a++){

		//generate probabilty for the next index
		proba  = (int)(1-(float(1.0*rand()/((float)RAND_MAX+1.0))));
		
		//try iterator
		if (proba <= proba_Mutation)
			{
				//three parameters to possibly mutate
				proba = rand()%3;
				//change only the phase of the circuit
				switch(proba){
					case 3:
						//modify the global phase
						I->phase = generatePhase();
						break;
					case 2:
						//mutate the n phases of the circuit
						proba2 = rand()%iterations;
						temp0 = 0;
						while(1){
							if (temp0 < I->my_string.length()){
								if (I->my_string[temp0] != 'p')
									proba2--;
							} else break;
							
							if (proba2 == 0 && I->my_string[temp0] != 'p'){
								I->phases[proba2] = generatePhase();
								break;
							}
							temp0++;
						}
						break;
					case 1:
						//mutate the string - single gate
						do {
							proba2 = (rand()%I->my_string.length())-1;
							if (proba2 < 0) proba2 = 0;
							gate =  I->my_string.at(proba2);
						}while(gate != 'p' || proba2>(I->my_string.length() -3));
							
						if (I->my_string.at(proba2+1) == 'p') proba2++;
						int proba3 = I->my_string.find("p", proba2+1);
		
						I->my_string.erase(proba2, (proba3-proba2+1));
		
						//tricky part is to keep the width of the circuit constant
						basic_string<char> tempS;
						int V = 0;
						do {
							temp = *gateArray[rand()%numofgates];
							if (temp.numIO <= (I->ioNumber-V) && temp.numIO > 0)
							{
								V += temp.numIO;
								tempS.insert(tempS.length(), 1, temp.representation);
							}
						}while(V != I->ioNumber);  
						tempS.insert(0, "p");
						tempS.insert(tempS.length(), "p");
						I->my_string.insert(proba2, tempS);
						break;
					}
				}

		}
		if (repair)
			repairCircuit(indv);			
}

/*******************************************
 * check the length
 *******************************************/

void GA::repairCircuit(int a)
{
	
	int counter = 0;
	Individual *I = population[a];
	
//	cout<<I->my_string<<endl;
	if ((int)I->my_string.length() > 0){
	for (int A = 0;A<(int)I->my_string.length();A++)
		if (I->my_string.at(A) == 'p') counter++;
	if (counter < mincircuitsize*2 || counter > segments*4){
		//reinitialize the circuit
		I->segmentNumber = segments;
		makeIndvString(a);
		//cout<<"repaired a circuit "<<a<<endl;
		generatePhases(I);
	}
	} else {
		I->segmentNumber = segments*2;
		makeIndvString(a);
		//cout<<"repaired a circuit "<<a<<endl;
		generatePhases(I);
	}

}

/******************************************* 
 * repair the circuit
 * check the width
 *******************************************/

void GA::repairForceCircuit(int a)
{
	int begin, end, counter = 0;
	Individual *I = population[a];
	qGate temp;
	//one pointer for the beginning of the segment
	begin = 0;
	//onepointer for the current position
	end = 1;
	

	//we put separators on the end and on the beginning
	if (I->my_string.at(begin) != 'p')
		I->my_string.insert(I->my_string.at(begin), 1, 'p');
	if (I->my_string.at(I->my_string.length()-1) != 'p')
		I->my_string.insert(I->my_string.at(I->my_string.length()-1), 1, 'p');
	//detect wrongly placed separators
	//cout<<"repair 0 in: "<<I->my_string<<endl;	
	return;
	if (I->my_string.size() > segments*I->ioNumber)
	{
		end = I->my_string.find('p', (segments*I->ioNumber)/2);
		if (end < I->my_string.size()){
			if (I->my_string.at(end+1) == 'p')
				end++;
			I->my_string.erase(end, I->my_string.length() - end+1);
			return;
		}
		return;		
	} else return;
	counter = I->my_string.length()-2;
	//detect wrongly placed separators
	//cout<<"repair 1 in: "<<I->my_string<<endl;	
	while(true){
		//if we are between 1 and before-last element
		if (counter > 0 && counter <I->my_string.length()-1){
			if (I->my_string.at(counter) == 'p')
			//if this, previous and next are all three 'p' we can erase it
				if (I->my_string.at(counter+1) == 'p')
				{
					if (I->my_string.at(counter-1) == 'p')
						I->my_string.erase(counter, 1);
				}else if (I->my_string.at(counter+1) != 'p')
			//if this is 'p', but previous and next are all gates, add one
						if (I->my_string.at(counter-1) != 'p')
							I->my_string.insert(counter, 1, 'A');
		}else
			break;
		
		counter--;
	}
	//one pointer for the beginning of the segment
	//cout<<"repair 1 out: "<<I->my_string<<endl;	
	begin = 0;
	//onepointer for the current position
	end = 1;
	//counter of the wires in a segment
	counter = 0;

	//check every segment of the string circuit
	while (true)
	{
		
		if (I->my_string.at(begin) != 'p')
			I->my_string.erase(begin, 1);
		else
			if (I->my_string.at(end) == 'p'){
				//if the current end pointer is NOT on gate
				if ((counter != resultnum)){
					//if the resulting wirecount is bigger than desired
					if (counter > resultnum)
						//erase if the segment is too fat
						I->my_string.erase(begin, end-begin+1);
					else {
						//add wires if the given segment is too short on wires
						do {
							temp = *gateArray[0];
		//					if (temp.numIO <= (resultnum - counter))
							counter += temp.numIO;
							I->my_string.insert(I->my_string.at(end-1), 1, 'A');
						}while(counter < I->ioNumber);  
					}
				//if the end is pointing to p and we got correct count -- increment
				}else {
					begin = end +1;
					end = I->my_string.find('p', begin+1);
					counter = 0;
					//cout<<begin<<" "<<end<<endl;	

				}
			//if the current pointer - end - is pointing on a gate
			}else {
				counter+= gateArray[getGate(I->my_string.at(end))]->numIO;
				end++;
			}
	}
		//cout<<"repair 2 out: "<<I->my_string<<endl;	

}

/*******************************************
 * Restrictions:
 * 1. replace gate with same number of I/O
 * 2. replace one gate with n I/O with any combination of smaller gates with
 * same number of I/O
 * 3. conserve strict correspondence of connections
 * void operator=(const class bla & other) { mymember=other.mymember; } 
 *******************************************/

void GA::apply1Crossover(int ind1, int ind2)
{  
     //	apply rules and make crossover
	int proba1, proba2;
	complex <float>temp[MAXNUMOFGATES];
	Individual *I = new_population[ind1];
	Individual *J = new_population[ind2];
	    basic_string<char> temp1;
		basic_string<char> temp2;
	int index1, index2, iterations1, iterations2;

	
        srand( time(NULL));
	//individual 1
		do {
			proba1 = (rand()%(I->my_string.size()));
		}while(proba1 >= I->my_string.length()-1 || I->my_string.at(proba1) != 'p');

		//cout<<proba1;

		if (proba1 <= 0) temp1 = I->my_string.substr(proba1);
			else 
			{				
				if ((I->my_string.at(proba1)) == 'p')
					if ((I->my_string.at(proba1+1)) == 'p') proba1+=1;

				temp1 = I->my_string.substr(proba1);
			}
			
	//individual 2
			
		do {
			proba2 = (rand()%(J->my_string.size()));
		}while(proba2 >= J->my_string.length()-1 || J->my_string.at(proba2) != 'p');

	//cout<<" "<<proba2<<endl;
			if (proba2 <= 0) temp2 = J->my_string.substr(proba2);
			else
			{				
				if ((J->my_string.at(proba2)) == 'p')
					if ((J->my_string.at(proba2+1)) == 'p') proba2+=1;

				temp2 = J->my_string.substr(proba2);
			}

			
	iterations1 = getGatesfromString(I->my_string);
	iterations2 = getGatesfromString(J->my_string);
	
	//get the phase size 
	//we need both indexes in the array of phases
			index1 = 0;
			index2 = 1;
			for (int a = 0;a < I->my_string.length(); a++)
				if( a < proba1){
					if(I->my_string[a] != 'p'){
						index1 +=1;
					}
				}else{
					if(I->my_string[a] != 'p'){
						temp[a-proba1] = I->phases[a];
					}
				}
				
			for (int a = 0;a < proba2; a++)
				if(J->my_string[a] != 'p')
					index2 +=1;

			//	swap the phase sequences
			for (int a = index2;a < (iterations2-index2); a++)
					I->phases[index1+a-index2] = J->phases[a];
			
			for (int a = index1;a < (iterations1-index1); a++)
					J->phases[index2+a-index1] = temp[a-index1];
		
			if (temp1.length() >=3 && temp2.length() >= 3)
			{
				I->my_string.replace(proba1, temp1.length(), temp2, 0, temp2.length());
				J->my_string.replace(proba2, temp2.length(), temp1, 0, temp1.length());
			}

}

/*******************************************
 * two point crossover - no phases swap
 *******************************************/
void GA::apply2Crossover(int ind1, int ind2)
{  
	//apply rules and make crossover
	//first it is used to store segment number
	int proba[4];
	basic_string<char> temp1;
	basic_string<char> temp2;
	Individual *I = new_population[ind1];
	Individual *J = new_population[ind2];
	int a,b,c,d,n,count;
	int pos1, pos2;

        srand( time(NULL));
	if (I->ioNumber == J->ioNumber)
	{
/*		cout<<"Before crossover:"<<endl;
		cout<<"0 individual:"<< I->my_string<<": "<<I->segmentNumber<<endl;
		cout<<"1 individual:"<< J->my_string<<": "<<J->segmentNumber<<endl;
*/		//get string cuting points 
		do {
			//the number of segments to swap on individual 1
			proba[0] = (rand()%(I->segmentNumber/2));
		}while(proba[0] <= 0 );
		a = proba[0]*2;
//		cout<<"1 #segs "<<proba[0]<< "coord: "<<a <<endl;
		do {
			//the location of the cut
			proba[1] = (rand()%(I->segmentNumber));
		}while((proba[1]+proba[0]) >= I->segmentNumber);
		b = proba[1]*2;
//		cout<<"1 loc segs "<<proba[1]<<"coord: "<<b<<endl;
		do {
			//the number of segments to swap on individual 2
			proba[2] = (rand()%(J->segmentNumber/2));
		}while(proba[2] <= 0);
		c = proba[2]*2;
//		cout<<"2 #segs "<<proba[2] <<" coord: "<<c<<endl;
		do {
			//the location of the cut
			proba[3] = (rand()%(J->segmentNumber));
		}while((proba[3]+proba[2]) >= J->segmentNumber);
		d = proba[3]*2;
//		cout<<"2 log segs "<<proba[3] <<" coord: "<<d<<endl;

		n = 0;
		count = 0;
		temp1 = "";
		while(count < (b+a)){

			if (I->my_string.at(n) == 'p')
				count++;
			//the string starts between the indexes of the segments
			if (count >b){
				temp1 += I->my_string.at(n);			
			}
			n++;
		}
		pos1 = I->my_string.find(temp1, 0);

		n = 0;
		count = 0;
		temp2 = "";
		while(count < (c+d)){

			if (J->my_string.at(n) == 'p')
				count++;
			//the string starts between the indexes of the segments
			if (count >d){
				temp2 += J->my_string.at(n);			
			}
			n++;
		}
		pos2 = J->my_string.find(temp2, 0);
/*		cout<<"1 ->"<<temp1 <<endl;
		cout<<"2 ->"<<temp2 <<endl;
*/

		I->my_string.replace(pos1,temp1.length(),temp2,0,temp2.length());
		J->my_string.replace(pos2,temp2.length(),temp1,0,temp1.length());
/*
		cout<<"After crossover:"<<endl;
		cout<<"0 individual:"<< I->my_string<<endl;
		cout<<"1 individual:"<< J->my_string<<endl;
*/

	}
}

/*******************************************
 * replicate the population according to the settings
 *******************************************/
void GA::applyReplication()
{ 

    float fitnessSum = 0;
	float first, second, tempsum, proba;
	int indv1, indv2;

    for (int i = 0; i < populationNumber; i++)
		fitnessSum += population[i]->fitness;

        srand( time(NULL));
   int counter = 0;
    bool condition = false;
    while (!condition)
	{
		if (replicator == 0)
		{
			//roulette wheel
			if (threshold > 0)
				do{
					indv1  = rand()%(populationNumber-1);
					indv2  = rand()%(populationNumber-1);
				} while(population[indv1]->fitness > threshold &&  population[indv2]->fitness > threshold);
			else
				{
					first  = (1-(float(1.0*rand()/((float)RAND_MAX+1.0))))*fitnessSum;
					second  = (1-(float(1.0*rand()/((float)RAND_MAX+1.0))))*fitnessSum;
				}
		}
		else
		{
			//stochastic uniersal sampling
			if (threshold > 0)
				do{
					indv1  = rand()%(populationNumber-1);
					indv2  = rand()%(populationNumber-1);
				} while(indv1 > threshold &&  indv2 > threshold);
			else
				{
					first  = ((1-(float(1.0*rand()/((float)RAND_MAX+1.0))))*fitnessSum)/2;
					second = fitnessSum - first;
				}
		}

		if (threshold == 0)
		{
			tempsum = 0;
			for (int c = 0;c<populationNumber;c++)
			{
				tempsum += population[c]->fitness;
				if (tempsum >= first) 
				{
					indv1 = c; 
					break;
				}
			}
			
			tempsum = 0;
			for (int d = 0;d<populationNumber;d++)
			{
				tempsum += population[d]->fitness;
				if (tempsum >= second) 
				{
					indv2 = d; 
					break;
				}
			}
		}
		
		//cout<<"generating new individual at index: "<<counter<<" from: "<<indv1<<endl;
		if(counter == populationNumber-1)
		{
			setIndv(new_population[counter], population[indv1]);
			counter++;
		} else
		{
			setIndv(new_population[counter], population[indv1]);
			setIndv(new_population[counter+1], population[indv2]);

			proba  = 1-(float(1.0*rand()/((float)RAND_MAX+1.0)));
			if (proba < proba_Crossover && indv1 != indv2)
			{
				 if (crossover == 0)
					 apply1Crossover(counter, counter+1);
				 else apply2Crossover(counter, counter+1);
			}
			
			counter += 2;
		}
		if (counter >= populationNumber-1) condition = true;
	}

    for (int i = 0;i<populationNumber;i++)
	{
	   setIndv(population[i], new_population[i]);
	 
	}
}

/*******************************************
 * Pareto replication with initial settings
 *******************************************/
void GA::applyParetoReplication()
{ 
	float proba;
	int indv1, indv2, tempsum, fitnessSum, first, second;
	int counter = 0;
    bool condition = false;

        srand( time(NULL));
	for (int i = 0; i < populationNumber; i++) population[i]->Rank = -1;
	int temp_rank;
	
	for (int i = 0; i < populationNumber; i++)
    {
		temp_rank = 0;
		for (int j = 0; j < populationNumber; j++)
		{
			if (population[i]->Rank <= -1)
			{
				if ((alpha1*population[i]->Error) <= population[j]->Error && (beta1*population[i]->Cost) <= population[j]->Cost)
					temp_rank++;
			}
		}
		population[i]->Rank = temp_rank;
    }

	fitnessSum = 0;
	for (int i = 0; i < populationNumber; i++)
    {
         fitnessSum += population[i]->Rank;
    }



    while (!condition)
	{
		if (replicator == 0)
		{
			//roulette wheel
			if (threshold > 0)
				do{
					indv1  = rand()%(populationNumber-1);
					indv2  = rand()%(populationNumber-1);
				} while(population[indv1]->fitness > threshold &&  population[indv2]->fitness > threshold);
			else
				{
					first  = (int)(1-(float(1.0*rand()/((float)RAND_MAX+1.0))))*fitnessSum;
					second  = (int)(1-(float(1.0*rand()/((float)RAND_MAX+1.0))))*fitnessSum;
				}
		}
		else
		{
			//stochastic uniersal sampling
			if (threshold > 0)
				do{
					indv1  = rand()%(populationNumber-1);
					indv2  = rand()%(populationNumber-1);
				} while(indv1 > threshold &&  indv2 > threshold);
			else
				{
					first  = (int)((1-(float(1.0*rand()/((float)RAND_MAX+1.0))))*fitnessSum)/2;
					second = fitnessSum - first;
				}
		}

		if (threshold == 0)
		{
			tempsum = 0;
			for (int c = 0;c<populationNumber;c++)
			{
				tempsum += population[c]->Rank;
				if (tempsum >= first) 
				{
					indv1 = c; 
					break;
				}
			}
			
			tempsum = 0;
			for (int d = 0;d<populationNumber;d++)
			{
				tempsum += population[d]->Rank;
				if (tempsum >= second) 
				{
					indv2 = d; 
					break;
				}
			}
		}
		

		if(counter == populationNumber-1){
			  setIndv(new_population[counter], population[indv1]);
			counter++;
		} else
		{

			 setIndv(new_population[counter], population[indv1]);
			 setIndv(new_population[counter+1], population[indv2]);

			proba  = 1-(float(1.0*rand()/((float)RAND_MAX+1.0)));
			if (proba < proba_Crossover && indv1 != indv2)
				if (crossover == 0)
					apply1Crossover(counter, counter+1);
				else apply2Crossover(counter, counter+1);

			counter += 2;
		}
		if (counter >= populationNumber-1) condition = true;
	}
    for (int i = 0;i<populationNumber;i++)
	   setIndv(population[i], new_population[i]);
}

/*******************************************
 * copy one individual to another one - data copy
 *******************************************/
void GA::setIndv(Individual *indv1, Individual *indv2)
{
		indv1->fitness = indv2->fitness;
		indv1->ioNumber = (int)indv2->ioNumber;
		indv1->segmentNumber = (int)indv2->segmentNumber;
		indv1->Cost = (int)indv2->Cost;
		indv1->Error = (float)indv2->Error;
		indv1->Groups = (int)indv2->Groups;			// for share fitness calcualtion
	 	indv1->Rank = (int)indv2->Rank;
		indv1->my_string.erase();
		indv1->my_string = string("");
		indv1->phase = indv2->phase;
		for (int a = 0; a < getGatesfromString(indv2->my_string); a++)
			indv1->phases[a] = indv2->phases[a];
		for (unsigned int a = 0 ; a < indv2->my_string.size(); a++)
			indv1->my_string += indv2->my_string[a];
}

/******************************************* 
 * cuts the amount of segments to 200
 *******************************************/
void GA::reduceSegNumber(int a)
{
	int counter, n;

	counter = 0;
	n = 0;
	
	if (population[a]->ioNumber > 1)
	{
		while(n < (int)population[a]->my_string.length())
		{
			if (population[a]->my_string.at(n) == 'p')
				counter++;
			if (counter >= 200)
			{
				population[a]->my_string.erase(n, population[a]->my_string.length());
				break;
			}
			n++;
		}
	}
}
/*******************************************
 * this works only if the circuits are normalized
 * we generate the Individual from gates and replace them 
 * with real gates.
 * the gates are translated to normal circuit using the 
 * permutativeArray parents
 *******************************************/
void GA::injectNormalCircuits(string **refs)
{
	int counter = 0;
	int count = 5;
	int tcounter = 0;
	
	if (refs == NULL)
		return;
	//cout<<"injecting circuits back to GA"<<endl;
	while(true)
	{
		if (counter >= count-1)
			break;
		else if (refs[counter] == NULL || refs[counter]->size() <= 2)
			break;
		counter++;
	}

	for (int a = counter-1; a >= 0 ; a--)
	{
		//cout<<"* "<<a<<endl;
		//cout<<"* "<<refs[a]<<endl;
		delete(population[populationNumber - a]);
		population[populationNumber - a] = new Individual();
		population[populationNumber - a]->my_string = string("");
		//cout<<"* "<<refs[a]<<"  "<<refs[a]->size()<<endl;
		for (unsigned int b = 0; b< refs[a]->size(); b++)
		{
			tcounter = 0;
			//fetch the tranlsated string from tbe permutative array
			while (permutativeArray[tcounter]->representation !=refs[a]->at(b))
				tcounter++;
			//cout<<"* "<<permutativeArray[tcounter]->my_string<<endl;
			population[populationNumber - a]->my_string += 'p'+permutativeArray[tcounter]->my_string+'p';
		}
		//trivial phase
		population[populationNumber - a]->phase  = complex <float>(0,1);
		population[populationNumber - a]->ioNumber = finalGate.numIO;
		//cout<<"<*> "<<population[populationNumber - a]->my_string<<endl;
		if (measurement > 0)
			makeMeasureFitness(population[populationNumber - a], true);
		else
			makeFitness(population[populationNumber - a], true);
		//cout<<"<> "<<population[populationNumber - a]->my_string<<endl;
	}
}

/*******************************************
 * this works only if the circuits are normalized
 * the gates aretranslated to normal circuit using the 
 * permutativeArray parents
 *******************************************/
void GA::injectNormalSegments(string **refs)
{
	int counter = 0;
	int count = 5;
	int tcounter = 0;
	
	if (refs == NULL)
		return;
	//cout<<"injecting circuits back to GA"<<endl;
	//check if there are valid solutions
	while(true)
	{
		if (counter >= count-1)
			break;
		else if (refs[counter] == NULL || refs[counter]->size() <= 2)
			break;
		counter++;
	}
	//print obtained circuits
	for (int a = 0; a < counter ; a++)
		//cout<<"* "<<refs[a]<<endl;

	//add all solutions into new indivduals and replace them in the GA
	for (int a = counter-1; a >= 0 ; a--)
	{
		//cout<<"* "<<a<<endl;
		//cout<<"* "<<refs[a]<<endl;
		//erase the first n that are the smallest ones.
		delete(population[populationNumber - a]);
		population[populationNumber - a] = new Individual();
		population[populationNumber - a]->my_string = string("");
		//cout<<"* "<<refs[a]<<"  "<<refs[a]->size()<<endl;
		for (unsigned int b = 0; b< refs[a]->size(); b++)
		{
			tcounter = 0;
			//fetch the tranlsated string from tbe permutative array
			while (permutativeArray[tcounter]->representation !=refs[a]->at(b))
				tcounter++;
			//cout<<"* "<<permutativeArray[tcounter]->my_string<<endl;
			//generate teh circuit string in the individual
			population[populationNumber - a]->my_string += 'p'+permutativeArray[tcounter]->my_string+'p';
		}
		population[populationNumber - a]->ioNumber = finalGate.numIO;
		//cout<<"<*> "<<population[populationNumber - a]->my_string<<endl;
		//evaluate the added circuit
		if (measurement > 0)
			makeMeasureFitness(population[populationNumber - a], true);
		else
			makeFitness(population[populationNumber - a], true);
		//cout<<"<> "<<population[populationNumber - a]->my_string<<endl;
	}
}

/*******************************************
 * Returns the aray of gates generated as single gate per parallel segments
 * All gates have thus size CircuitSize and can be simply swapped however the number o fgates grows to maximum
 *******************************************/
qGate** GA::getPermutativeCircuits()
{
	
	qGate *workgate, multigate, *multigate2;
	string names = "";
	int gatecounter = 0;
	int c = 0;
	int rep = 0;
	string s = bestIndv.my_string;
	for (int m = 0; m< MAXNUMOFGATES; m++)
		permutativeArray[m] = NULL;
	//cout<<"starting generating: "<<s<<"  "<<endl;
	for (unsigned int a = 0;a< s.length();a++){
		//cout<<"starting generating "<<s.at(a)<<"  "<<endl;
		rep = 0;
		//we found a real gate
		if (s.at(a) != 'p' && s.at(a) != 'A')
		{
			//cout<<"starting generating 0"<<endl;
			for (int b = 0;b < numofgates;b++){
			//remove the separators and the wire gates
				if (gateArray[b]->representation == s.at(a)){
					//check wether it was already added
					for (unsigned int m = 0; m < names.size(); m++)
						if (names.at(m) == s.at(a))
							rep = 100;
					workgate = gateArray[b];
				}
			}
			//cout<<"starting generating 1"<<endl;
			//if yes do not add this gate again
			if (rep == 0){
				//add the gate
				names += workgate->representation;
				//here comes the transformation for non normal gate to normal ones
				//if the given gate is smaller than the final one
		//cout<<"starting generating 2"<<endl;
				if ( finalGate.numIO > workgate->numIO )
				{
					//generate as many new gates -- normal ones
		//cout<<"starting generating 3"<<endl;
					for (int t = 0; t <= finalGate.numIO - workgate->numIO; t++)
					{
						if (permutativeArray[gatecounter] != NULL)
							delete(permutativeArray[gatecounter]);
						permutativeArray[gatecounter] = new qGate;
						permutativeArray[gatecounter]->my_string = string("");
						//how many W's we have to add to create a shifted one
						for (int j = 0; j <= (finalGate.numIO - workgate->numIO); j++)
						{
							if (j == t)
								permutativeArray[gatecounter]->my_string += workgate->representation;
							else 
								permutativeArray[gatecounter]->my_string += 'A';
						}
						permutativeArray[gatecounter]->representation = char('A'+gatecounter);
						if (t == 0)
							multigate = copyGate(workgate);	
						else 
							multigate = copyGate(gateArray[0]);
						//calculate the matrix of the segment
						c = 0;
						while(true) 
						{
							if (c != t )
								multigate2= gateArray[0];	
							else
								multigate2= workgate;
							//using Kronecker multiplication
							multigate = tensorProduct(&multigate, multigate2);
							c++;
							if (c >=  (finalGate.numIO - workgate->numIO))
								break;
						}
		//cout<<"starting generating 4"<<endl;
						//fill the matrix of the gate
						for (int e = 0; e < (int(pow(float(2),finalGate.numIO))); e++)
							for (int f = 0; f < (int(pow(float(2),finalGate.numIO))); f++)
								permutativeArray[gatecounter]->gateMatrix1[e][f] = multigate.gateMatrix1[e][f];
								permutativeArray[gatecounter]->numIO = finalGate.numIO;
						gatecounter++;
					}
				} else {
						if (permutativeArray[gatecounter] != NULL)
							delete(permutativeArray[gatecounter]);
						permutativeArray[gatecounter] = NULL;
						permutativeArray[gatecounter] = new qGate;
						permutativeArray[gatecounter]->representation = char('A'+gatecounter);
						for (int e = 0; e < (int(pow(float(2),finalGate.numIO))); e++)
							for (int f = 0; f < (int(pow(float(2),finalGate.numIO))); f++)
								permutativeArray[gatecounter]->gateMatrix1[e][f] = workgate->gateMatrix1[e][f];
								permutativeArray[gatecounter]->numIO = finalGate.numIO;
						gatecounter++;
				}
			}
		}
	}
	//cout<<"starting generating 10"<<endl;
	numofnormalgates = gatecounter;
		
	//for (int n = 0; n < numofnormalgates; n++)
		//cout<<permutativeArray[n]->representation<<" : "<<permutativeArray[n]->my_string<<endl;
	//cout<<endl;
	
	return permutativeArray;
}

/*******************************************
* Returns a set of gates representing single segments of the circuit
* this can be used if circuits with big size are generated
* *******************************************/
 
qGate** GA::getPermutativeSegments()
{
	
	qGate *workgate, multigate, *multigate2;
	string names = "";
	int gatecounter = 0;
	string s = bestIndv.my_string;
	bool processing = false;
	for (int m = 0; m< MAXNUMOFGATES; m++){
		permutativeArray[m] = NULL;
	}
	//cout<<"starting generating: "<<s<<"  "<<endl;
	for (unsigned int a = 0;a< s.length();a++){
		//cout<<"starting generating "<<s.at(a)<<"  "<<endl;
		//we found a real gate
		if (s.at(a) == 'p'){
			if (a+1<s.length()){
				if (s.at(a+1) == 'p'){
					//cout<<"Name generated : "<<names<<endl;
					processing = false;
					gatecounter++;
				} else {
					processing = true;
					names = string("");
				}
			}else{ //cout<<"Name generated : "<<names<<endl;
				processing = false;
			}
		}else{
			workgate = gateArray[getGate(s.at(a))];
			names += workgate->representation;
			//do the multiplication when the last gate is registered
			if (a+1<s.length() && s.at(a+1) == 'p'){
				//cout<<"starting generating 0"<<endl;
				//do the final gate multiplication
				for (int n = 0 ; n < names.length(); n++){
					if (n == 0){//set the 0'th gate
						workgate = gateArray[getGate(names.at(n))];
						if (permutativeArray[gatecounter] != NULL)
							delete(permutativeArray[gatecounter]);
						permutativeArray[gatecounter] = new qGate();
						permutativeArray[gatecounter]->representation = char('A'+gatecounter);
						multigate = copyGate(workgate);	
						permutativeArray[gatecounter]->my_string = string(names);
						permutativeArray[gatecounter]->numIO = finalGate.numIO;
						//cout<<"starting generating 10: "<<permutativeArray[gatecounter]->my_string<<endl;
					}else{//add the next gate using the tensor product
						multigate2= gateArray[getGate(names.at(n))];	
						multigate = tensorProduct(&multigate, multigate2);
					}
				}
				for (int e = 0; e < (int(pow(float(2),finalGate.numIO))); e++)
					for (int f = 0; f < (int(pow(float(2),finalGate.numIO))); f++)
						permutativeArray[gatecounter]->gateMatrix1[e][f] = multigate.gateMatrix1[e][f];

				
			}
		}
	}
	gatecounter++; //addlast gate
	//cout<<"starting generating 10"<<endl;
	numofnormalgates = gatecounter;
		
	for (int n = 0; n < numofnormalgates; n++){
		//cout<<permutativeArray[n]->representation<<" : "<<permutativeArray[n]->my_string<<endl;
	}
	//cout<<bestIndv.my_string<<endl;
	//cout<<endl;
	
	return permutativeArray;
}

/*******************************************
 * Returns the terminating condition
 *******************************************/
bool GA::terminatingCondition()
{
     return condition;
}

/*******************************************
 * Sets the terminating condition
 *******************************************/
void GA::setTerminatingCondition(int counter)
{
	generationCondition = counter;
     if (counter >= MAXGEN)
	 {
		 out_stream<<" generation max reached "<<endl<<endl<<endl;
		 condition = true;
	 }
     else  if (getFitness()) 
			{
				out_stream<<" Solution found "<<endl<<endl<<endl;
				condition = true;
			}	
			else condition = false;
}

/*******************************************
 * check wether or not the correct solution was found
 * and if not it refresh the bestIndividual to the newest best value
 *******************************************/
bool GA::getFitness()
{
	
	if (bestIndv.fitness >= 1.0)
	{
		return true;
	}
		else return false;
}

void GA::setpopulationNumber(int populationNumber){
	this->populationNumber = populationNumber;
}
int GA::getpopulationNumber(){
	return this->populationNumber;
}
void GA::setsegments(int segments){
	this->segments = segments;
}
int GA::getsegments(){
	return this->segments;
}
void GA::setmincircuitsize(int mincircuitsize){
	        this->mincircuitsize = mincircuitsize;
}
int GA::getmincircuitsize(){
	        return this->mincircuitsize;
}
void GA::setalterations(int alterations){
	this->alterations = alterations;
}
int GA::getalterations(){
	return this->alterations;
}
void GA::setproba_Mutation(float proba_Mutation){
	this->proba_Mutation = proba_Mutation;
}
float GA::getproba_Mutation(){
	return this->proba_Mutation;
}
void GA::setproba_Crossover(float proba_Crossover){
	this->proba_Crossover = proba_Crossover;
}
float GA::getproba_Crossover(){
	return this->proba_Crossover;
}
void GA::setalpha(float alpha){
	this->alpha = alpha;
}
float GA::getalpha(){
	return this->alpha;
}
void GA::setbeta(float beta){
	this->beta = beta;
}
float GA::getbeta(){
	return this->beta;
}
void GA::setalpha1(float alpha1){
	this->alpha1 = alpha1;
}
float GA::getalpha1(){
	return this->alpha1;
}
void GA::setbeta1(float beta1){
	this->beta1 = beta1;
}
float GA::getbeta1(){
	return this->beta1;
}
void GA::setdivider(int divider){
	this->divider = divider;
}
int GA::getdivider(){
	return this->divider;
}
void GA::setphase(int phase){
	this->phase = phase;
}
int GA::getphase(){
	return this->phase;
}
void GA::setdisplaylevel(int displaylevel){
	this->displaylevel = displaylevel;
}
int GA::getdisplaylevel(){
	return this->displaylevel;
}
void GA::setGa(int Ga){
	this->Ga = Ga;
}
int GA::getGa(){
	return this->Ga;
}
void GA::setmutation(int mutation){
	this->mutation = mutation;
}
int GA::getmutation(){
	return this->mutation;
}
void GA::setcrossover(int crossover){
	this->crossover = crossover;
}
int GA::getcrossover(){
	return this->crossover;
}
void GA::setreplicator(int replicator){
	this->replicator = replicator;
}
int GA::getreplicator(){
	return this->replicator;
}
void GA::setfitness(int fitness){
	this->fitness = fitness;
}
int GA::getfitness(){
	return bestIndv.fitness;
}
void GA::setgrouped(int grouped){
	this->grouped = grouped;
}
int GA::getgrouped(){
	return this->grouped;
}
void GA::setpareto(int pareto){
	this->pareto = pareto;
}
int GA::getpareto(){
	return this->pareto;
}
void GA::setthreshold(int threshold){
	this->threshold = threshold;
}
int GA::getthreshold(){
	return (int)this->threshold;
}
void GA::setresultnum(int resultnum){
	this->resultnum = resultnum;
}
int GA::getresultnum(){
	return this->resultnum;
}
void GA::setmeasurement(int measurement){
	this->measurement = measurement;
}
int GA::getmeasurement(){
	return this->measurement;
}
/*******************************************
 * Initialize the GA gates (input gates, working set, measurements)
 *******************************************/
void GA::initializeGates()
{
	complex <float> I = complex<float>(0., 1.); //  Declare complex float variables
	complex <float> in = complex<float>(0.,0.);
	float x, y;
	int counter = 0;
	char chcounter = char('A');
	string temp;

	//reading the input file 
	cout <<"Reading input parameters"<<endl;
	in_stream>> populationNumber;		//int
	in_stream>> segments;				//int
	in_stream>> mincircuitsize;				//int
	in_stream>> alterations;			//int
	in_stream>> proba_Mutation;		//float
	in_stream>> proba_Crossover;		//float
	in_stream>> alpha;					//float
	in_stream>> beta;					//float
	in_stream>> alpha1;				//float
	in_stream>> beta1;					//float
	in_stream>> divider;				//int
	in_stream>> phase;					//int
	in_stream>> displaylevel;					//int
	//reading file divider
	in_stream >> temp;
	in_stream >> Ga;					//int
	in_stream >> mutation;				//int
	in_stream >> crossover;				//int
	in_stream >> replicator;			//int
	in_stream >> fitness;				//int
	in_stream >> grouped;				//bool
	in_stream >> pareto;				//bool
	in_stream >> threshold;				//float
	in_stream >> resultnum;
	in_stream>> measurement;			//
	if (measurement > 0)
	{
		int i = 0;
		//the input of measured qubits is a sequence of numbers separated by space
		for (int h = 0; h< measurement; h++)
			in_stream >> measurementQBits[i++];
		cout << "Reading Measurements for "<<i<<" qubits"<<endl;
		for (int b = 0; b< measurement*2; b++){
			for (int c = 0; c < (int)(pow((float)2, (float)resultnum)); c++){
				in_stream >> in;
				measureexpected[c][b] = in;
			}
		}
	}

	//reading file divider
	in_stream >> temp;

	cout << " Generating Output " <<endl;
	//generate the header of the output file
	out_stream<< " -------------------------------------------------------- "<<endl;
	out_stream<< "         output file of synthesis: settings               "<<endl;
	out_stream<< " -------------------------------------------------------- "<<endl;
	out_stream<< " Display level:                                           "<<displaylevel<<endl;
	out_stream<< " size of the population:                                  "<<populationNumber<<endl;
	out_stream<< " generation alteration cycle:                             "<<alterations<<endl;
	out_stream<< " number of segments in each circuit:                      "<<segments<<endl;
	out_stream<< " mutation probability is                                  "<< proba_Mutation<<endl;
	out_stream<< " mutation of crossover is                                 "<<proba_Crossover<<endl;
	out_stream<< " factor Alpha is                                          "<<alpha<<endl;
	out_stream<< " factor Beta is                                           "<<beta<<endl;
	out_stream<< " factor Alpha1 is                                         "<<alpha1<<endl;
	out_stream<< " factor Beta1 is                                          "<<beta1<<endl;
	out_stream<< " estimated minimum cost of this circuit                   "<<divider<<endl;
	out_stream<< " phase used                                               "<<phase<<endl;
	out_stream<< " measurement used                                         "<<measurement<<endl;
	out_stream<< " measured bits are                                        ";
	for (int a = 0; a< measurement; a++)
			out_stream<<  measurementQBits[a]<<" ";
	out_stream<< endl;
	for (int b = 0; b< measurement*2; b++){
		for (int c = 0; c < (int)(pow((float)2, (float)resultnum)); c++)
			out_stream << measureexpected[c][b];
		out_stream<<endl;
	} 

	out_stream<< " -------------------------------------- "<<endl;
	out_stream<< " type of GA 0 - normal, 1 - Darwinian                     "<<Ga<<endl;
	out_stream<< " type of mutation 0 - normal, 1 - bitwise                 "<<mutation<<endl;
	out_stream<< " type of crossover 0 - 1point, 1 - 2point                 "<< crossover<<endl;
	out_stream<< " type of replication 0 - RW, 1 - SUS                      "<<replicator<<endl;
	out_stream<< " type of fitness 0 - simplest, 3 - complex                "<<fitness<<endl;
	out_stream<< " type of fitness calculation 0 - individual, 1 - grouped  "<<grouped<<endl;
	out_stream<< " Pareto optimization 0 - no, 1 - yes                      "<< pareto<<endl;
	out_stream<< " Threshold replication 0 - no, other - threshold          "<<threshold<<endl;
	out_stream<< " The number of wires the final circuit has:               "<<resultnum<<endl;
	out_stream<< " -------------------------------------------------------- "<<endl;
	out_stream<< "          output file of synthesis: input gates           "<<endl;
	out_stream<< " -------------------------------------------------------- "<<endl;



	//bedining reading the gates
    	in_stream >> numofgates;
	string s;
	// get first gate
	for (int a = 0; a < numofgates; a++)
	{
		gateArray[a] = new qGate;
		gateArray[a]->representation = char(chcounter++);
		in_stream >> gateArray[a]->numIO;
		in_stream >> gateArray[a]->Cost;
		
		in_stream >>gateArray[a]->my_string;
		for (int j = 0; j < (int(pow((float)2,(float)gateArray[a]->numIO))); j++)
		{
			for (int k = 0; k < (int(pow((float)2,(float)gateArray[a]->numIO))); k++)
			{
				in_stream >> in;
				//in_stream >> y;
				gateArray[a]->gateMatrix1[j][k] = in;//complex <float>(x,y);
			}
		}
		counter++;
	}

	// get the expected result matrix from the file
	in_stream >> resultnum;
	finalGate.numIO = resultnum;
	for (int m = 0; m < (int(pow((float)2,(float)resultnum))); m++)
	{
		for (int n = 0; n < (int(pow((float)2,(float)resultnum))); n++)
		{
			in_stream >> in;
			finalGate.gateMatrix1[m][n] = in;
		}
	}
	
	generateExtendMatrix(chcounter);

	for (int a = 0; a < numofgates; a++)
        {
	out_stream<< gateArray[a]->representation <<endl;
                out_stream<<"  Cost of the gate  "<< gateArray[a]->Cost <<endl;
                out_stream<<"  IO of the gate  "<< gateArray[a]->numIO <<endl;
                out_stream<<"The name of this gate is: "<<gateArray[a]->my_string<<endl;
                for (int j = 0; j < (int(pow((float)2,(float)gateArray[a]->numIO))); j++)
                {
                        for (int k = 0; k < (int(pow((float)2,(float)gateArray[a]->numIO))); k++)
                        {
                                out_stream<< gateArray[a]->gateMatrix1[j][k]<<' ';
                        }
                        out_stream<< endl;
                }
 
        }

	
	out_stream<< "final gate "<<endl;
        for (int m = 0; m < (int(pow((float)2,(float)resultnum))); m++)
        {
                for (int n = 0; n < (int(pow((float)2,(float)resultnum))); n++)
                {
                                out_stream<< finalGate.gateMatrix1[m][n]<<' ';
                }
                out_stream<<endl;
        }


}
/*******************************************
 * Initialize the single qubit measurements
 *******************************************/
void GA::generateMeasureOps()
{
	int l = finalGate.numIO;
	//the counter for measurement matrices is implemented as follows
	//2^(abs) index - numNBGR-1 (to get to the index array values)
	//int k = (int)pow(2.0, (float)(abs(index - (tRule->numNGBR-1))));
	int n,m = 0;
	bool zeros = true;
	int a = 0;
	string srank;
	qGate *measure0, *measure1;
	qGate measureM0, measureM1;
	int p = (int(pow((float)2,(float)l)));

	//generate measurement opertors for every single qubit
	//cout<< "generating measurement circuits: "<< MAXMEASUREGATE;
	while (a < l*2){
		n = 0;
		convert_stream<<m;
		srank = convert_stream.str();

		measurements[a] = new qGate;
		measurements[a]->numIO = l;
		measurements[a]->representation = '0';
		measurements[a]->my_string += "th qubit";
		measurements[a]->my_string = srank+measurements[a]->my_string;
		measurements[a]->my_string = " on the "+ measurements[a]->my_string;
		measurements[a]->my_string = measurements[a]->representation+measurements[a]->my_string;
		measurements[a]->my_string = "measurement for "+measurements[a]->my_string;
		measurements[a+1] = new qGate;
		measurements[a+1]->numIO = l;
		measurements[a+1]->representation = '1';
		measurements[a+1]->my_string += "th qubit";
		measurements[a+1]->my_string = srank+measurements[a+1]->my_string;
		measurements[a+1]->my_string = " on the "+measurements[a+1]->my_string;
		measurements[a+1]->my_string = measurements[a+1]->representation+measurements[a+1]->my_string;
		measurements[a+1]->my_string = "measurement for "+measurements[a+1]->my_string;
		zeros = true;
		for (int x = 0; x < (int)(pow(2.0, (float)l)); x++){
			for (int y = 0; y < (int)(pow(2.0, (float)l)); y++){
				measurements[a]->gateMatrix1[x][y] = complex<float>(0,0);
				measurements[a+1]->gateMatrix1[x][y] = complex<float>(0,0);
				if (x == y){
					if (n <= 0){
						zeros = true;
					} else if (n ==	(int)(pow(2.0, (float)m))){
						zeros = false;
					} 
					if (!zeros){
						measurements[a+1]->gateMatrix1[x][y] = complex<float>(1,0);
						n--;
					}else {
						measurements[a]->gateMatrix1[x][y] = complex<float>(1,0);	
						n++;
					}
				}
			}
		}
		a +=2;
		m++;
	}
	
	convert_stream.flush();
	measure0 = new qGate();
	measure1 = new qGate();
	
	//for every orthonormal state on the input
	cout<< "Generating Measurements For all desired states"<<endl;
	for (int k = 0; k<p; k++){
		initMatrix(&measureM0, resultnum);
		initMatrix(&measureM1, resultnum);
		for (int m = 0; m < measurement;m++)
		{
			measure0 = ccGate(measurements[measurementQBits[m]*2], measure0);
			measure1 = ccGate(measurements[(measurementQBits[m]*2)+1], measure1);
			//check if entanglement is desired
			if (real(measureexpected[k][2*m]) == real(measureexpected[k][2*m+1])){
				if (real(measureexpected[k][2*m]) == 1){
					measureM0 = matrixProduct(&measureM0, measure1);
					measureM1 = matrixProduct(&measureM1, measure0);
				}else if (real(measureexpected[k][2*m]) == 0){
					if (!(imag(measureexpected[k][2*m]) == 1 && imag(measureexpected[k][2*m+1]) == 1)){
						measureM0 = matrixProduct(&measureM0, measure0);
						measureM1 = matrixProduct(&measureM1, measure1);
					}
				}else { 
					if (!(imag(measureexpected[k][2*m]) == 1 && imag(measureexpected[k][2*m+1]) == 1)){
						measureM0 = matrixProduct(&measureM0, measure0);
						measureM1 = matrixProduct(&measureM1, measure1);
					}
				}
			} else{
				if (real(measureexpected[k][2*m]) > real(measureexpected[k][2*m+1])){
					measureM0 = matrixProduct(&measureM0, measure0);
					measureM1 = matrixProduct(&measureM1, measure1);
				} else if (real(measureexpected[k][2*m]) <= real(measureexpected[k][2*m+1])){
					measureM0 = matrixProduct(&measureM0, measure1);
					measureM1 = matrixProduct(&measureM1, measure0);
				}
			}
		}

		convert_stream<<k;
		srank = convert_stream.str();
		measurementsAllState[2*k] = new qGate();
		measurementsAllState[2*k] = ccGate(&measureM0,  measurementsAllState[2*k]);
		measurementsAllState[2*k]->my_string = srank+" Desired input state";
	//	cout<<k<<": measurement operator :"<<measurementsAllState[2*k]<<endl;

		//measurementsAllState[2*k+1] = cGate(&measureM1);
		measurementsAllState[2*k+1] = new qGate();
		measurementsAllState[2*k+1] = ccGate(&measureM1, measurementsAllState[2*k+1]);
		measurementsAllState[2*k+1]->my_string = srank+" Undesired state ";
	//	cout<<k<<": measurement operator :"<<measurementsAllState[2*k+1]<<endl;
        }




	//output the measurement operators
	for (int a = 0; a < p; a++)
	{
	//	cout<<a<<": measurement operator :"<<measurementsAllState[2*a]<<endl;
		out_stream<< "Whole state Measurement Gates:"<<endl;
		out_stream<< measurementsAllState[2*a]->my_string <<endl;
		out_stream<< measurementsAllState[2*a]->representation <<endl;
		out_stream<<"  IO of the gate  "<< measurementsAllState[2*a]->numIO <<endl;
		for (int j = 0; j < (int(pow((float)2,(float)measurementsAllState[2*a]->numIO))); j++)
		{
			for (int k = 0; k < (int(pow((float)2,(float)measurementsAllState[2*a]->numIO))); k++)
			{
				out_stream<< measurementsAllState[2*a]->gateMatrix1[j][k]<<' ';
			}
			out_stream<< endl;
		}

	//	cout<<a<<": measurement operator :"<<measurementsAllState[2*a+1]<<endl;
		out_stream<< measurementsAllState[2*a+1]->my_string <<endl;
		out_stream<< measurementsAllState[2*a+1]->representation <<endl;
		out_stream<<"  IO of the gate  "<< measurementsAllState[2*a+1]->numIO <<endl;
		for (int j = 0; j < (int(pow((float)2,(float)measurementsAllState[2*a+1]->numIO))); j++)
		{
			for (int k = 0; k < (int(pow((float)2,(float)measurementsAllState[2*a+1]->numIO))); k++)
			{
				out_stream<< measurementsAllState[2*a+1]->gateMatrix1[j][k]<<' ';
			}
			out_stream<< endl;
		}

	}
	delete(measure0, measure1);

}
/*******************************************
 * Initialize the expectations probabilities for measured global states
 *******************************************/
void GA::generateMeasureExpect()
{

	int l = finalGate.numIO;
	int p = (int(pow((float)2,(float)l)));
	for (int g = 0; g<p; g++){
		 expectationsAllState[g][0] = cplex(1,0);
                expectationsAllState[g][1] = cplex(1,0);
	}

	for (int k = 0; k<p; k++){
	//	cout << " For input state: "<<k <<endl;
		for (int m = 0; m < measurement;m++)
		{
			//check if entanglement is desired
			if (real(measureexpected[k][2*m]) == real(measureexpected[k][2*m+1])){
				if (real(measureexpected[k][2*m]) != 1 &&  real(measureexpected[k][2*m]) != 0){
					expectationsAllState[k][0] *= measureexpected[k][2*m];
					expectationsAllState[k][1] *= measureexpected[k][2*m+1];
				}
			} else{
				if (real(measureexpected[k][2*m]) > real(measureexpected[k][2*m+1])){
					expectationsAllState[k][0] *= measureexpected[k][2*m];
					expectationsAllState[k][1] *= measureexpected[k][2*m+1];
				} else if (real(measureexpected[k][2*m]) <= real(measureexpected[k][2*m+1])){
					expectationsAllState[k][0] *= measureexpected[k][2*m+1];
					expectationsAllState[k][1] *= measureexpected[k][2*m];
				}
			}
		}
        }
}

/*******************************************
 * Verifies the input gates and generates all missing ones
 * for multiple qubits
 *******************************************/
void GA::generateExtendMatrix(char chcounter)
{

	int counter = 0;
	string name;	
	qGate gate, extgate, *tempgate;
	int topwire = gateArray[numofgates-1]->numIO;
	while(true)
		if (topwire < resultnum)
		{
			topwire++;
			//extgate = new qGate;
			initMatrix(&extgate, 2);
			extgate = *ccGate(gateArray[getqGate("SWAP")], &extgate);
			for (int l = 2; l < topwire; l++){
				extgate = tensorProduct(&extgate, gateArray[0]);		
			}
			counter = 0;
			for (int a = numofgates-1; a >= 0; a--){
				if (gateArray[a]->numIO == topwire-1){
					initMatrix(&gate, gateArray[a]->numIO);
					gate = *ccGate(gateArray[a], &gate);
					name = "";
					for (int m = 0; m < gate.my_string.length();m++)
						name += gate.my_string.at(m);
					gate = tensorProduct(gateArray[0], &gate);
					gate = matrixProduct(&extgate, &gate);
					gate = matrixProduct(&gate, &extgate);
					gate.my_string = "";
					for (int m = 0; m < name.length();m++)
						gate.my_string += name.at(m);
					if (gate.my_string.at(0) == 'C'){
						gate.my_string.insert(1,"I");
					} else if (gate.my_string.at(gate.my_string.length()-1) == 'C'){
                                                gate.my_string.insert(gate.my_string.length()-1,"I");
                                        } else  gate.my_string.insert(1,"I");

					gateArray[numofgates+counter] = cGate(&gate); 
					gateArray[numofgates+counter]->representation = char(chcounter++);
					counter++;

				}
			}
			numofgates += counter;
		} else break;

}
void GA::closeStream()
{
     in_stream.close();
     out_stream.close();
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//methods of the minimizer -- yet all in one class
//no time for improvements
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


/*******************************************
 * compare two gates
 *******************************************/
int GA::findGate(qGate *A)
{
	int a = 0;
	bool matrixTest;
	int result = -1;
	while (a<numofgates)
	{
		out_stream<< a<<" "<<numofgates<<endl;
		matrixTest = true;

		if (A->numIO != gateArray[a]->numIO)
			result = -1;
		else
		{
			for (int i = 0; i < (int(pow((float)2,(float)A->numIO))); i++)
				for (int j = 0; j < (int(pow((float)2,(float)A->numIO))); j++)
                      if (A->gateMatrix1[i][j] != gateArray[a]->gateMatrix1[i][j])
						  matrixTest = false;
			if (matrixTest) 
				result = a;
			else result = -1;

		}

		if (result != -1) break;
		else a++;
	}

	return result;
	
}

/*******************************************
 * the super-method intializing the circuit minimization
 * now the minimzation is done here only virtually
 *******************************************/
int GA::minimizeCirc(Individual *qC)
{
	int begin0, begin1, end0, end1;
	bool equal;
	int a, b;
	int result = 0;

	if (qC->ioNumber > 1)
	{	
		begin0 = 0;
		end0 = qC->my_string.find("p", begin0+1)+1;
		begin1 = end0;
		end1 = qC->my_string.find("p", begin1+1)+1;
		b = 0;

		while (end1 > 0)
		{
			equal = false;
			a = 0;
			if ((qC->my_string.substr(begin0, (end0-begin0))).length() == (qC->my_string.substr(begin1, (end1-begin1))).length() )
			{
				while(a <= (end0-begin0-2))
				{
					if ( (qC->my_string.at(begin0+a)) != 'p')
					{	
						if (gateArray[getGate(qC->my_string.at(begin0+a))]->numIO != gateArray[getGate(qC->my_string.at(begin1+a))]->numIO)
						{
							equal = false;
							break;
						}
						else 
						{
							if (a >= end0-begin0-2) 
							{
								equal = true;
							}
						}
					}
					a++;
				}
			}
			if (!equal)
			{
				for (int a = 1;a<(end0-begin0-1);a++)
				{
					result += gateArray[getGate(qC->my_string.at(begin0+a))]->Cost;
				}
			}
			begin0 = begin1;
			end0 = qC->my_string.find("p", begin0+1)+1;
			begin1 = end0;
			end1 = qC->my_string.find("p", begin1+1)+1;
		}
	}
	for (a = 1;a<(end0-begin0-1);a++)
	{
		result += gateArray[getGate(qC->my_string.at(begin0+a))]->Cost;
	}
	qC->Cost = result;

	return result;
}
/*******************************************
 * removes physically redundant segments in the concenred circuit
 *******************************************/
void GA::removeRedund(Individual *C)
{
	int begin0, begin1, end0, end1;
	bool equal, erased;
	int a, b;

	if (C->ioNumber > 1)
	{	
		begin0 = 0;
		end0 = C->my_string.find("p", begin0+1)+1;
		begin1 = end0;
		end1 = C->my_string.find("p", begin1+1)+1;
		b = 0;

		while (end0<(int)C->my_string.length()-3 || end0 < 0)
		{
			equal = false;
			erased = false;
			a = 0;
			if (C->my_string.substr(begin0, (end0-begin0)) == C->my_string.substr(begin1, (end1-begin1)))
			{
				C->my_string.erase(begin0, (end1-begin0));
				erased = true;
			}
	

			if (erased)
			{
				if (begin0 != 0)
				{
					begin0 = 0;
					end0 = C->my_string.find("p", begin0+1)+1;
					
					begin1 = end0;
					end1 = C->my_string.find("p", begin1+1)+1;
				}
				else 
				{
					begin1 = end0;
					end1 = C->my_string.find("p", begin1+1)+1;
					erased = false;
				}
			}
			else 
			{
				begin0 = begin1;
				end0 = C->my_string.find("p", begin0+1)+1;
				if (end0 < (int)C->my_string.length()-3)
				{
					begin1 = end0;
					end1 = C->my_string.find("p", begin1+1)+1;
				}
			}
		}	
	}
}
/*******************************************
 * Modifies the genome - minimization - lamarckian learning
 *******************************************/
void GA::mergeSegs(Individual *Circ, int first, int second)
{
	basic_string<char> tempStr;
	qGate *temp, *fGate, *sGate;
	int gateCount;

	int endfrsSg = Circ->my_string.find("p", first+1)+1;
	int endscdSg = Circ->my_string.find("p", second+1)+1;
	tempStr.append(1, 'p');

	for (int a = 0;a<(endfrsSg-first-1);a++)
	{
		if (Circ->my_string.at(first+a) != 'p')
		{
			fGate = gateArray[getGate(Circ->my_string.at(first+a))];
			sGate = gateArray[getGate(Circ->my_string.at(second+a))];
			temp = &matrixProduct(fGate, sGate);
			gateCount = findGate(temp);

			if (gateCount != -1)
			{
				tempStr.append(1, gateArray[gateCount]->representation);
			}
			else
			{
				gateArray[numofgates] = new qGate;
				(*gateArray[numofgates]) = copyGate(temp);

				if (gateArray[numofgates-1]->representation == 'o')
					gateArray[numofgates]->representation = gateArray[numofgates-1]->representation+char(2);
				else if (gateArray[numofgates-1]->representation == char(264))
					gateArray[numofgates]->representation = gateArray[numofgates-1]->representation+char(6);
					else gateArray[numofgates]->representation = gateArray[numofgates-1]->representation+char(1);

				gateArray[numofgates]->Cost = fGate->Cost;
				gateArray[numofgates]->parentA = fGate->representation;
				gateArray[numofgates]->parentB = fGate->representation;
				tempStr.append(1, gateArray[numofgates]->representation);
				numofgates++;
			}
		}
	}
	tempStr.append(1, 'p');
	Circ->my_string.erase(first, (endscdSg-first));
	Circ->my_string.insert(first,tempStr);
}
/*******************************************
 * generate the fitness value of the Individual
 *******************************************/
void GA::makeFitness(Individual *indi, bool display)
{
    //dummy variables
	int Time = 0;
	int a, b, c;
	complex <float> temphase;
	int phasemax;
	int	phasecounter0 = 0;
    qGate temp0, temp1, myfinalG;
	
	//null this Individual previous evaluation
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;
	indi->segmentNumber = 0;

	int begin = 0;
	int end = indi->my_string.find("p", begin+1);
	c =0;
	//FIXME: add individual phase for every gate
	// get the number of used gates;
	phasemax = getGatesfromString(indi->my_string);
	while(end > 0) 
	{
		if ((end - begin) > 1)
		{
			c++;
			//get the first  gate of the Segment and copy it
			for (a = 0;a<= numofgates;a++)
				if (gateArray[a]->representation == indi->my_string.at(begin+1)) 
				{
					temp0= copyGate(gateArray[a]);
					if (phase == 1){
						temphase = indi->phases[phasecounter0++];
						for (int k = 0; k<(int(pow((float)2,(float)temp0.numIO)));k++)
                					for (int l = 0; l<(int(pow((float)2,(float)temp0.numIO)));l++)
                     						temp0.gateMatrix1[k][l] *= temphase;
                     			}
					//continue with the cost
					begin++; indi->Cost+= gateArray[a]->Cost;
					break;
				}
			//get the next one            
			for (b = begin+1;b<end;b++) 
			{
				for (a = 0;a<= numofgates;a++)
					if (gateArray[a]->representation == indi->my_string.at(b)) 
					{
						temp1= copyGate(gateArray[a]);
						temphase = indi->phases[phasecounter0++];
						//cout<<"apply the phase2"<<endl;
						if (phase == 1){
							for (int k = 0; k<(int(pow((float)2,(float)temp0.numIO)));k++)
								for (int l = 0; l<(int(pow((float)2,(float)temp0.numIO)));l++)
									temp1.gateMatrix1[k][l] *= temphase;
						}
						//continue with the cost
						indi->Cost+= gateArray[a]->Cost;
						a = numofgates*2;
						//multiply the gates
						//using Kronecker multiplication
						temp0 = tensorProduct(&temp0, &temp1);
					}
			}
			if (Time == 0) 
			{
				//if this is the end of the first segment fill the output value
				myfinalG = copyGate(&temp0);
				Time++;
			}
			else 
			{
				//compile warning is intended
				//multiply using normal Manumofgatestrix product the twou segments
				myfinalG = matrixProduct(&myfinalG, &temp0); 
				Time++;
			}
		}
		//move to the next segment
		begin = indi->my_string.find("p", end);
		end = indi->my_string.find("p", begin+1);
	}
	//set some more parameters
	myfinalG.numIO = indi->ioNumber;
        indi->segmentNumber =c;
	//cout<<"From make fitness: "<<indi->segmentNumber <<endl;
	//check if phase is defined and multiply the whole matrix by it
	if (phase == 1){
		if (indi->phase != complex <float> (0,0))
		{
			for (int k = 0; k<(int(pow((float)2,(float)myfinalG.numIO)));k++)
        	        	for (int l = 0; l<(int(pow((float)2,(float)myfinalG.numIO)));l++)
                     			myfinalG.gateMatrix1[k][l] *= indi->phase;

		}
	}
	temp0= copyGate(&myfinalG);
	for (int k = 0; k<(int(pow((float)2,(float)myfinalG.numIO)));k++)
        	for (int l = 0; l<(int(pow((float)2,(float)myfinalG.numIO)));l++)
			myfinalG.gateMatrix1[k][l] = temp0.gateMatrix1[l][k];
        //direct correspondence between number of wires in the goal and current individual
        if (myfinalG.numIO!= finalGate.numIO) 
        {
                indi->fitness = 0;
        }
        else 
        {
                for (int a = 0; a<(int(pow((float)2,(float)myfinalG.numIO)));a++)
                {
                        for (int b = 0; b<(int(pow((float)2,(float)myfinalG.numIO)));b++)
                        { 
                                //simple comparison -- allowing don't cares as value "10"
                                if (finalGate.gateMatrix1[a][b] != complex<float>(10, 0))
								{
								
									//x = (myfinalG.gateMatrix1[a][b].real()+myfinalG.gateMatrix1[a][b].imag())*(myfinalG.gateMatrix1[a][b].real()-myfinalG.gateMatrix1[a][b].imag());
									//y = (finalGate.gateMatrix1[a][b].real()+finalGate.gateMatrix1[a][b].imag())*(finalGate.gateMatrix1[a][b].real()-finalGate.gateMatrix1[a][b].imag());
                                    
									//evaluation of the error basedon on one to one comaprison
									//indi->Error += abs(finalGate.gateMatrix1[a][b] - myfinalG.gateMatrix1[a][b]);
									
									//evaluation based on the A*transpose(conj(B)) - the density matrix
									//indi->Error += abs(finalGate.gateMatrix1[a][b]*conj(myfinalG.gateMatrix1[b][a]));
									
									//evaluation based only on the A*conj(B) - a simpler version
									//indi->Error += abs(finalGate.gateMatrix1[a][b]*conj(myfinalG.gateMatrix1[a][b]));

									//evaluation based only on the abs(A*conj(A) - B*conj(B)) a simpler version
									indi->Error += abs(finalGate.gateMatrix1[a][b]*conj(finalGate.gateMatrix1[a][b]) - myfinalG.gateMatrix1[a][b]*conj(myfinalG.gateMatrix1[a][b]));
									}
                         }
                }
        }

        //normalizing error
        indi->Error =  indi->Error/pow((float)2,(float)(myfinalG.numIO*2));
        indi->Cost = (exp(-pow((divider-indi->Cost),2)));
		//generate fitness value
        if (indi->Error != 0)
		{
			switch(replicator){
					case 0:
						//simple fitness1
						indi->fitness = (1-indi->Error);
						break;
					case 1:
						//simple fitness1
						indi->fitness = (1/(indi->Error+1));
						break;
					case 2:
						//scaled complex fitness1
					    indi->fitness = (alpha*(1-indi->Error) + beta*indi->Cost);
						break;
					case 3:
						//scaled complex fitness2
						indi->fitness = (alpha*(1/(indi->Error+1)) + beta*indi->Cost);
						break;
					case 4:
						indi->fitness = (1-indi->Error);
						break;
			}
		}else 
		{
			indi->Error = 0;
			indi->fitness = 1;//we found a good individual;
		}
		//dump output for algorithmic debug
         if (display)
        {
            out_stream << " -------- " <<endl;
            out_stream << " -------- " <<endl;
            out_stream << " -------- " <<endl;
            out_stream <<"Representation String: "<<indi->my_string<<endl;
			out_stream <<"Phase for valid gates: "<<end;
			for(int g = 0; g< getGatesfromString(indi->my_string); g++)
				out_stream <<indi->phases[g]<<", "<<end;
			out_stream<<endl;
            out_stream <<"Fitness: "<< indi->fitness<<endl;
            out_stream <<"Error: "<< indi->Error<<endl;
            out_stream <<"Wires: "<< indi->ioNumber<<endl;
            if (phase == 1)
			out_stream <<"Phase: "<< indi->phase<<endl<<endl;
            out_stream <<"Matrix: " <<endl;
        	for (int a = 0; a<(int(pow((float)2,(float)myfinalG.numIO)));a++)
                {
                        for (int b = 0; b<(int(pow((float)2,(float)myfinalG.numIO)));b++)
                        { 
                                out_stream<< myfinalG.gateMatrix1[a][b]<<"";
                         }
                         out_stream<<endl;
                }
                out_stream<<endl;
            out_stream << " -------- " <<endl;
            out_stream << " -------- " <<endl;
            out_stream << " -------- " <<endl;
        }
	delete(&myfinalG, &temp0, &temp1);
}
/*******************************************
 * calculates the average of fitness and error
 *******************************************/
void GA::calcAverage(){
	
	int counter = finalResult.counter;
	float avfitness = 0;
	float averror = 0;
	float avcost = 0;
	for (int t = 0; t < populationNumber; t++){
		avfitness += population[t]->fitness;
		averror += population[t]->Error;
		avcost += population[t]->Cost;
	}
	finalResult.avFitness[counter] = avfitness/populationNumber;
	finalResult.avError[counter] = averror/populationNumber;
	finalResult.avCost[counter] = avcost/populationNumber;
	finalResult.counter = counter + 1;
}


