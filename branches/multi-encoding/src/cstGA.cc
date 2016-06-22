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
#ifdef __MPI__
#include "EpiG_mpi.h"
#else
#include "EpiG.h"
#endif

/*global static field for thread counting */
int GA::populationCounter;

/* global static switch for measurement/nonmeaurement */
int GA::measurement;
/****************************************
 * generates the population order in fitness value:
 * (inc = 0) decreasing order
 * (inc = 1) increasing order
 ****************************************/
void GA::evaluatePopulation() {
	int c = 0;
	float totfit = 0.0;
	float toterr = 0.0;
	float totcost = 0.0;
	for (int t = 0; t < populationNumber; t++) {
		totfit += population[t]->fitness;
		toterr += population[t]->Error;
		totcost += population[t]->Cost;
//	cout<<"currrent : "<<population[t]->fitness<<endl;
		if (population[t]->fitness > population[c]->fitness) {
			c = t;
		}
	}
//	cout<<"currrent best: "<<population[c]->fitness<<endl;
	totfit /= populationNumber;
	totcost /= populationNumber;
	toterr /= populationNumber;
	history[generationCondition][0] = float(totfit);
	history[generationCondition][1] = float(population[c]->fitness);
	history[generationCondition][2] = float(toterr);
	history[generationCondition][3] = float(totcost);
	progress = -1;
	if (bestIndv->fitness >= 1) {
		progress = 1;
	} else if (population[c]->fitness > bestIndv->fitness) {
		cout<<"new best "<<endl;
		progress = 1;
		setIndv(bestIndv, population[c]);
		outputBest(false);
		bestIndv->fitness = 0 + population[c]->fitness;
		cout << "Best new fitness: " << bestIndv->fitness << endl;
	}
	if (Ga > 1)
		sortGates();
}

/****************************************
 * calculate fitness for individual:
 ****************************************/
void GA::makeFitness(Individual *indi) {

	if (indi->Error != 0) {
		switch (replicator) {
		case 0:
			//simple linear fitness
			indi->fitness = (1 - indi->Error);
			break;
		case 1:
			//simple inverse fitness
			indi->fitness = (1 / (indi->Error + 1));
			break;
		case 2:
			//scaled complex fitness1
			indi->fitness = (alpha * (1 - indi->Error) + beta * (exp(-pow(abs(
					divider - indi->Cost), 2))));
			break;
		case 3:
			//scaled complex fitness2
			indi->fitness = (alpha * (1 / (indi->Error + 1)) + beta * (exp(
					-pow(abs(divider - indi->Cost), 2))));
			break;
		case 4:
			indi->fitness = (1 - indi->Error);
			break;
		case 5:
			indi->fitness = exp(-indi->Error);
			break;
		}
	} else {
		indi->Error = 0;
		indi->fitness = 1;//we found a good individual;
	}
}

/****************************************
 * static - the function called to invoke the program as a thread
 ****************************************/
void* GA::startSub(void *arg) {
	GA* obj = (GA*) arg;
	obj->doFitness();
	return 0;
}


/****************************************
 * computes the fitness
 * ****************************************/
void* GA::doFitness() {
	int rc = pthread_mutex_lock(&mtx);
	int ind = GA::populationCounter;
	GA::populationCounter = (GA::populationCounter + 1) % populationNumber;
	rc = pthread_mutex_unlock(&mtx);
//cout<<"ind: "<<ind;
#ifdef __QMDD__

	decodeIndvStrToRevLibStr(population[ind]);
	doQMDDFitness(population[ind], finalIndividual, false);
#else
	if (measurement > 0) {
		if (seq_detect_enabled > 0)
			doMeasureFASeqFitness(population[ind], false);
		else 
			doMeasureFitness(population[ind], false);
	} else {
		doMatrixFitness(population[ind], false);
	}
#endif
//	cout<<"ordering...";
	if(Ga > 1){

		rankInputGates(population[ind]);
	}
//	cout<<"done"<<endl;
}
/****************************************
 * generate the fitness value of the Individual and displays it if desired
 ****************************************/
void* GA::doMatrixFitness(Individual *ind, bool display) {
	//two dummy gates and one pointer for operations
	Individual *indi = ind;
	qGate* myfinalG;
	int rowcount, rowcounter,maxcount;
	int errorcounter = 0;
	float val = indi->valuedness;
	float numIO, error, equal;
	int diagonal_min = 0;
	indi->Error = 0;
	indi->fitness = 0;
	//threaded code
	int rc = pthread_mutex_lock(&mtx);
	if (rc) {
		cout << "Matrix locked: " << ind << endl;
	}

//cout<<"matrix...";
#ifdef __CUBLAS__
	myfinalG = computeCUBLASMatrix(indi);
#else
#ifdef __CUDA__
	myfinalG = computeCUDAMatrix(indi);
#else
	myfinalG = computeMatrix(indi);
#endif
#endif
//cout<<"done"<<endl;
	rc = pthread_mutex_unlock(&mtx);

	//set some more parameters
	myfinalG->numIO = indi->ioNumber;
	numIO = myfinalG->numIO;
	int columns = (int(pow((float) val, (float) numIO)));
	if (phase > 0) {
		maxcount = (int(pow(pow(val, numIO),2)));
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			for (int k = 0; k < maxcount; k++){
					cout<<cuCrealf(myfinalG->gateMatrix1[k])<<","<<cuCimagf(myfinalG->gateMatrix1[k])<<"  ";
					myfinalG->gateMatrix1[k] = cuCmulf(myfinalG->gateMatrix1[k], indi->phase);
			}
		}
	}
	//direct correspondence between number of wires in the goal and current individual
	error = 0;
	errorcounter = 0;
	diagonal_min = 0;
	if (myfinalG->numIO != finalGate->numIO) {
		indi->fitness = 0;
		error = 1;
	} else {
		equal = 0;
		rowcount = (int(pow(val, numIO)));
		if (local_forcing > 0 ){
			for (int a = 0; a < measured_desired_output_records; a++)
				equal += cuCrealf(cuCmulf(myfinalG->gateMatrix1[non_measured_desired_output_records_idx[2*a]+non_measured_desired_output_records_idx[2*a]*columns], cuConjf(myfinalG->gateMatrix1[non_measured_desired_output_records_idx[2*a]+non_measured_desired_output_records_idx[2*a]*columns])));
		}
		if (equal == pow(val, numIO)) {
			error = 20.0;
		} else {
			if(measured_desired_output_records < 1){
				maxcount = (int(pow(pow(val, numIO),2)));
				rowcounter = 0;
				if (force_exact > 0 ){
					for (int a = 0; a < maxcount; a++) {
						//simple comparison -- allowing don't cares as value "-1"
						if (cuCrealf(finalGate->gateMatrix1[a]) != 0 && cuCimagf(finalGate->gateMatrix1[a]) != -1) {
							if (cuCrealf(finalGate->gateMatrix1[rowcounter]) != cuCrealf(myfinalG->gateMatrix1[rowcounter])){
								error += pow (cuCrealf(finalGate->gateMatrix1[rowcounter]) - cuCrealf(myfinalG->gateMatrix1[rowcounter]),2);
							}
							if (cuCimagf(finalGate->gateMatrix1[rowcounter]) != cuCimagf(myfinalG->gateMatrix1[rowcounter])){
								error += pow (cuCimagf(finalGate->gateMatrix1[rowcounter]) - cuCimagf(myfinalG->gateMatrix1[rowcounter]),2);					
							}
							errorcounter++;
						}
					}
				} else {
					for (int a = 0; a < maxcount; a++) {
						//simple comparison -- allowing don't cares as value "-1"
						if (cuCrealf(finalGate->gateMatrix1[a]) != 0 && cuCimagf(finalGate->gateMatrix1[a]) != -1) {
							errorcounter++;
							error += cuCabsf(cuCsubf(cuCmulf(finalGate->gateMatrix1[a], cuConjf(finalGate->gateMatrix1[a])),cuCmulf(myfinalG->gateMatrix1[a], cuConjf(myfinalG->gateMatrix1[a]))));
						}
					}
				}
			} else {
				if (force_exact > 0 ){
					for (int a = 0; a < measured_desired_output_records; a++){
						rowcounter = non_measured_desired_output_records_idx[2*a]+ non_measured_desired_output_records_idx[2*a+1]*columns;
							errorcounter++;
						if (cuCrealf(finalGate->gateMatrix1[rowcounter]) != cuCrealf(myfinalG->gateMatrix1[rowcounter])){
							error += pow (cuCrealf(finalGate->gateMatrix1[rowcounter]) - cuCrealf(myfinalG->gateMatrix1[rowcounter]),2);
						}
						if (cuCimagf(finalGate->gateMatrix1[rowcounter]) != cuCimagf(myfinalG->gateMatrix1[rowcounter])){
							error += pow (cuCimagf(finalGate->gateMatrix1[rowcounter]) - cuCimagf(myfinalG->gateMatrix1[rowcounter]),2);					
						}
					}
				} else {
					for (int a = 0; a < measured_desired_output_records; a++){
						rowcounter = non_measured_desired_output_records_idx[2*a]+ non_measured_desired_output_records_idx[2*a+1]*columns;
						errorcounter++;
						error += cuCabsf(cuCsubf(cuCmulf(finalGate->gateMatrix1[rowcounter], cuConjf(finalGate->gateMatrix1[rowcounter])),cuCmulf(myfinalG->gateMatrix1[rowcounter], cuConjf(myfinalG->gateMatrix1[rowcounter]))));
					}
				}
			}
		}
	}
	//normalizing error
	if (errorcounter != 0)
		indi->Error = error / errorcounter;
	else
		indi->Error = error;
	//generate fitness value
	makeFitness(indi);
	if (display) {
		out_stream << "Matrix: " << endl;
		maxcount = (int(pow(pow(val, myfinalG->numIO),2)));
		rowcount = (int(pow(val, myfinalG->numIO)));
		rowcounter = 0;
		for (int a = 0; a < rowcount; a++) {
			for (int b = 0; b < rowcount; b++) {
				out_stream << "(" << cuCrealf(myfinalG->gateMatrix1[a*rowcount+b])<< "," << cuCimagf(myfinalG->gateMatrix1[a*rowcount+b]) << ")"<< "";

			}
		out_stream << endl;
		}
		out_stream << " -------- " << endl;
	}
//cout<<"cleaning...";
	destroyGate(myfinalG);
	delete (myfinalG);
//cout<<"done"<<endl;
}
/****************************************
 * the String array represents the gates by the number of wires
 * and present a list of gates for each amount of wires
 ****************************************/
void GA::addToStringArray(qGate *Q) {
	int count[MAXNUMOFGATES];
	int j = Q->numIO - 1;
	stringArray[j].numberGates++;
	stringArray[j].nameGate[stringArray[j].numberGates] = new char[4];
	copyRep(Q->representation, stringArray[j].nameGate[stringArray[j].numberGates]);	
}

/****************************************
 * the String array represents the gates by the number of wires
 * and present a list of gates for each amount of wires
 ****************************************/
void GA::setStringArray() {
	qGate *Q;
	int count[MAXNUMOFGATES];
	stringArray = new tempT[MAXGATEINPUT];
	for (int j = 0; j < MAXGATEINPUT; j++) {
		count[j] = 0;
		for (int i = 0; i < numofgates; i++) {
			Q = gateArray[i];
			if (Q->numIO == (j + 1)) {
				stringArray[j].nameGate[count[j]] = new char[4];
				copyRep(Q->representation, stringArray[j].nameGate[count[j]]);	
				//stringArray[j].nameGate[count[j]] = Q->representation;
				count[j]++;
			}
			stringArray[j].numberGates = count[j];
		}
	}
}
/****************************************
 * this generates the String repressenting every individual.
 *it is called upon the structure of the initial population
 ****************************************/
void GA::makeIndvString(Individual *I) {
	int tempNumber, tempNumber2, tempGate, counter, temporyhold, current_wire;
	qGate *qgate;
	char gate[3];
	//Individual *I = population[i];
	counter = 0;
	I->my_string = string("");
	for (int j = 0; j < I->segmentNumber; j++) // for each segment of this individual
	{
		current_wire = 0;
		if (I->ioNumber > MAXGATEINPUT) // if the input has more than MAXGATEINPUT (10)
			tempGate = MAXGATEINPUT; // max gate to use
		else
			tempGate = I->ioNumber;
		tempNumber = I->ioNumber; // the number of inputs
		tempNumber2 = 0;
		if (I->ioNumber > 1) // if not only one input
		{
			I->my_string += 'p';
			counter++;
			while (tempNumber  > 0) // while more gate could be parallel
			{
				do {
					tempNumber2 = rand() % tempGate; // which group of gates to use
				} while (stringArray[tempNumber2].numberGates == 0 || (tempNumber2 + 1) > tempNumber );
				temporyhold = rand() % stringArray[tempNumber2].numberGates;
				qgate = gateArray[getGate(stringArray[tempNumber2].nameGate[temporyhold])];
				if (checkRestrictions(current_wire, qgate)){
					current_wire += tempNumber2+1;
					copyRep(stringArray[tempNumber2].nameGate[temporyhold], gate);
					for(int y = 0; y < REPSIZE; y++)
						I->my_string += gate[y];
					counter++;
					tempNumber = tempNumber - (tempNumber2 + 1);
				}
			}
		} else {
		
			temporyhold = rand() % stringArray[tempNumber2].numberGates;
			copyRep(stringArray[tempNumber2].nameGate[temporyhold], gate);
			for(int y = 0; y < REPSIZE; y++)
				I->my_string += gate[y];
			counter++;
		}

		if (I->ioNumber > 1) {
			I->my_string += 'p';
			counter++;
		}
	}
}
/****************************************
 * this is the main ENd output function
 * when the parameter individual = -1;the whole Input/Output Settings,
 * GA settings are dumped into file
 * when the parameter individual != -1; The given individual is printed out
 ****************************************/
void GA::output(int individual) {
	qGate *qG;
	int columns;

	if (individual != -1) {

		out_stream << "---------------------------" << endl;
		out_stream << " End of run at " << generationCondition << endl;
		out_stream << "---------------------------" << endl;
		out_stream << "----- Individual #" << individual << "-----" << endl;

		if (measurement > 0) {
			//		//this->display = true;
			//		makeMeasureFitness(population[individual], true);
			//	} else {
			//	  makeFitness(population[individual], true);
			//	}
		}

	} else {
		//		if (!(measurement > 0)) {
		out_stream << "There are " << numofgates
		<< " gates from the input file and synthesized.\n";
		float val = (float) gateArray[0]->valuedness;
		out_stream << "The radix is: " << val << "\n";
		for (int i = 0; i < numofgates; i++) {
			qG = gateArray[i];
			out_stream << "The gate number " << i << " is: " << endl;
			out_stream << "And the gate representation is: "
			<< gateArray[i]->representation << endl;
			out_stream << "Number of inputs is: " << qG->numIO << endl;
			out_stream << "Parent A is: " << qG->parentA << endl;
			out_stream << "Parent B is: " << qG->parentB << endl;

			float numIO = (float) gateArray[0]->numIO;
			columns = (int(pow(val, numIO)));
			for (int j = 0; j < columns; j++) {
				for (int k = 0; k < columns; k++) {
					//out_stream << qG->gateMatrix1[j][k] << " ";
					out_stream << "  (" << cuCrealf(qG->gateMatrix1[j+k*columns])<< "," << cuCimagf(qG->gateMatrix1[j+k*columns]) << ")";

				}
				out_stream << endl;
			}
		}
		out_stream << "The Target gate(s): \n";
		out_stream << "Number of inputs is: " << finalGate->numIO << endl;
		columns = (int(pow(val, resultnum)));
		for (int m = 0; m < columns; m++) {
			for (int n = 0; n < columns; n++) {
				//out_stream << finalGate->gateMatrix1[n][m] << " ";
				out_stream << "  (" << cuCrealf(finalGate->gateMatrix1[n+m*columns])<< "," << cuCimagf(finalGate->gateMatrix1[n+m*columns]) << ")";
			}
			out_stream << endl;
		}
	}
	this->display = false;
}
/****************************************
 * outputs the Best Invididual
 ****************************************/
void GA::outputBest(bool end) {
	char gate[REPSIZE];
	if (bestIndv->fitness >= 1 || progress > 0) {
		out_stream << "---------------------------" << endl;
		out_stream << " Generation at " << generationCondition << endl;
		out_stream << "---------------------------" << endl;
		out_stream << "----- Best Individual -----" << endl;
		out_stream << bestIndv->ioNumber << " <- number of inputs " << endl;
		out_stream << bestIndv->my_string << " <- Representation " << endl;
		int location = 0;
//		cout<<bestIndv->my_string<< ',' << bestIndv->my_string.length() << endl;
		while(location < bestIndv->my_string.length()){
			if ( bestIndv->my_string.at(location) == 'p'){
				out_stream << "p";
//				cout<<'p';
				location++;
			} else {
//				cout<<"getting gate at "<<location<<": "<<bestIndv->my_string.at(location)<<", ";
				if (getGateRepfromStr(location, bestIndv->my_string, gate) < 0 ) break;
				cout << gate<<", ";
				out_stream<<" "<<gateArray[getGate(gate)]->my_string<<" ";
				location += 3;
			}
		}
		cout<<bestIndv->my_string<<endl;
		out_stream << bestIndv->my_string << " <- Representation " << endl;
		out_stream << bestIndv->fitness << " <- Fitness " << endl;
		out_stream << bestIndv->Error << " <- Error " << endl;
		out_stream << bestIndv->Cost << " <- Cost " << endl;
		out_stream << bestIndv->valuedness << " <- Radix " << endl;
		if (displaylevel >= 1) {

#ifdef __QMDD__

			decodeIndvStrToRevLibStr(bestIndv);
			doQMDDFitness(bestIndv, finalIndividual, true);
#else
			if (measurement > 0) {
				if (seq_detect_enabled > 0)
					doMeasureFASeqFitness(bestIndv, true);
				else 
					doMeasureFitness(bestIndv, true);
			} else {
				doMatrixFitness(bestIndv, true);
			}

#endif
		}
		out_stream << "---------------------------" << endl;
	}
	if (end && displaylevel > 0) {
		out_stream << "-----Results------" << endl;
		out_stream << "Avg. Fitness, Best Fit.,  Avg. Error,    Avg. Cost  "
		<< endl;
		for (int i = 0; i < generations; i++) {
			for (int j = 0; j < 4; j++) {
				out_stream << history[i][j] << "  ";
			}
			out_stream << endl;
		}
		out_stream << "-----Results------" << endl;

		if (Ga > 1){
			out_stream<<"there are "<<rankV.size()<<" gates"<<endl;
			out_stream<<"Idx Usg Fit"<<endl;
			for (int a = 0; a < rankV.size(); a++){
				out_stream<<' '<<rankV[a].index<<" "<<rankV[a].usage<<" "<<rankV[a].fitness<<" "<<gateArray[rankV[a].index]->representation<<"  "<<gateArray[rankV[a].index]->my_string<<endl;
			}

		}
	}
}
/****************************************
 * Restrictions:
 * 1. no fan out - ex, 5 inputs has to have 5 outputs
 * 2. maximum 10 I/O (geometric increase of complexity with more I/Os)
 ****************************************/
void GA::initializePop(int level, string filename, string out_prefix) {

	string inFile = filename;
	string outFile;
	string post = (".out");
	string oFile = "";
	time_t rawtime;
	struct tm * timeinfo;
	int columns,w;
	condition = false;

	// Open input and output files.
	in_stream.open(inFile.c_str());
	if (in_stream.fail()) {
		cout << "\"Daa set \" file opening failed.\n";
		exit(1);
	}

	outFile = inFile;
/*
int a = 3;
int aa = 3;
	for(int b = 0; b < 9; b++){
		for(int c = 0; c < 9; c++){

			cout<<(c+b*9)<<": element of a: "<<(((c/aa)%a)+(((b/aa)%a)*a))<<" element of b: "<<((c%aa)+((b)%aa)*aa)<<endl;
		}
	}
exit(0);
*/
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	oFile = asctime(timeinfo);

	unsigned int t = 0;
	while (1) {
		//remove spaces
		if (t < oFile.length()) {
			if (oFile[t] >= '0' && oFile[t] <= 'z') {
				outFile += oFile[t];
			} else {
				if (oFile[t] == ' ') {
				} else {
					break;
				}
			}
		} else
			break;
		t++;
	}

	outFile += post;
	outFile = out_prefix + outFile;
	cout << "Output file will be: " << outFile << endl;
	out_stream.open(outFile.c_str());
	if (out_stream.fail()) {
		cout << "Output file opening failed: " << outFile.c_str() << "\n";
		exit(1);
	}
	//initialize the gates

	initializeGates();


	if (measurement > 0){
		cout<<"Observation probabilities Generation\n";
		generateMeasureExpect();
		cout<<" ..... done"<<endl;
	}

	//only single qubit based sequence detection
	for (int r = 0; r < gateArray[0]->valuedness; r++) {
		measures_o[r] = measurements[r];
		measures_i[r] = measurements[(int) gateArray[0]->valuedness + r];
	}

	


	//initialize the String array representing the gates
	setStringArray();
	cout << "Initial population of circuits generation";
	GA::populationCounter = 0;
	//generate the initial population/new_population of individuals
	for (int i = 0; i < populationNumber; i++) {
		population[i] = new Individual;
		new_population[i] = new Individual;
		//all individuals have the width of the output
		population[i]->ioNumber = finalGate->numIO;
		population[i]->valuedness = finalGate->valuedness;

		float val = (float) gateArray[0]->valuedness;
		//population[i]->segmentNumber = rand()%(2*segments)+5;
		population[i]->segmentNumber = segments;
		if (phase > 0) {
			population[i]->phase = generatePhase();
		}
		//population[i]->segmentNumber = rand()%MAXSEGNUM+1;
		//generate the String representing this circuit
		makeIndvString(population[i]);
		//generate phases for every single gate
		if (phase > 0) {
			generatePhases(population[i]);
		}
#ifdef __QMDD__
		decodeIndvStrToRevLibStr(population[i]);
		doQMDDFitness(population[i], finalIndividual, false);
#else
		if (measurement > 0){
			if (seq_detect_enabled > 0)
				doMeasureFASeqFitness(population[i], false);
			else 
				doMeasureFitness(population[i], false);
		}else{
			doMatrixFitness(population[i], false);
		}
		if (Ga > 1)
			rankInputGates(population[i]);
#endif
	}
	cout << " ..... done" << endl;
	//init the best individual
	cout << "Creating best individual storage";
	bestIndv = new Individual();
	bestIndv->segmentNumber = segments;
	bestIndv->ioNumber = finalGate->numIO;
	bestIndv->valuedness = finalGate->valuedness;
	int val = bestIndv->valuedness;
	makeIndvString(bestIndv);
#ifdef __QMDD__
		decodeIndvStrToRevLibStr(bestIndv);
		doQMDDFitness(bestIndv, finalIndividual, false);
#else
	if (measurement > 0){
		if (seq_detect_enabled > 0)
			doMeasureFASeqFitness(bestIndv, false);
		else 
			doMeasureFitness(bestIndv, false);
	}else{
		doMatrixFitness(bestIndv, false);
	}
#endif
	cout << " ..... done" << endl;
	outputBest(false);
	initiateStorage();
	cout << "GA Initialization Successful" << endl<<endl<<endl;
	out_stream << "\nInitialization done" << endl;
}


/****************************************
 * Reads the parameters from the open input streams and generates all gates
 * as well as initiates the measurement operator generation if requested
 ****************************************/
void GA::initializeGates() {
        complex<float> in = complex<float> (0., 0.);

	float x, y;
	int counter = 0;
	int seq_type = 0;
	int dc = 0;
	int ycounter;
	int valuedness = 0;
	int resnum, outnum, a, b, columns;
	bool mem_error = false;
	//char chcounter = char('!');
	char *representation = new char[4];
	initRep(representation);
	std::string temp;
	std::string temp1;
	char buff_temp[100];
	
	char unknown;
	char tempchar;
	char buff[100];
	char binary[100];
	char line[512];

	//reading the input file
	cout << "Reading input parameters";
	in_stream.ignore(256,' ');
	in_stream >> populationNumber; //int
	in_stream.ignore(256,' ');
	in_stream >> segments; //int
	in_stream.ignore(256,' ');
	in_stream >> mincircuitsize; //int
	in_stream.ignore(256,' ');
	in_stream >> segments_multi; //int
	in_stream.ignore(256,' ');
	in_stream >> alterations; //int
	in_stream.ignore(256,' ');
	in_stream >> generations; //int
	in_stream.ignore(256,' ');
	in_stream >> proba_Mutation; //float
	in_stream.ignore(256,' ');
	in_stream >> proba_Crossover; //float
	in_stream.ignore(256,' ');
	in_stream >> alpha; //float
	in_stream.ignore(256,' ');
	in_stream >> beta; //float
	in_stream.ignore(256,' ');
	in_stream >> alpha1; //float
	in_stream.ignore(256,' ');
	in_stream >> beta1; //float
	in_stream.ignore(256,' ');
	in_stream >> divider; //int
	in_stream.ignore(256,' ');
	in_stream >> phase; //int
	in_stream.ignore(256,' ');
	in_stream >> displaylevel; //int
	in_stream.ignore(256,' ');
	in_stream >> Ga; //int
	in_stream.ignore(256,' ');
	in_stream >> elitism; //int
	in_stream.ignore(256,' ');
	in_stream >> mutation; //int
	in_stream.ignore(256,' ');
	in_stream >> crossover; //int
	in_stream.ignore(256,' ');
	in_stream >> replicator; //int
	in_stream.ignore(256,' ');
	in_stream >> fitness; //int
	in_stream.ignore(256,' ');
	in_stream >> grouped; //bool
	in_stream.ignore(256,' ');
	in_stream >> pareto; //bool
	in_stream.ignore(256,' ');
	in_stream >> threshold; //float
	in_stream.ignore(256,' ');
	in_stream >> resultnum;
	in_stream.ignore(256,' ');
	in_stream >> valuedness;
	in_stream.ignore(256,' ');
	in_stream >> force_exact;
	in_stream.ignore(256,' ');
	in_stream >> local_forcing;
	in_stream.ignore(256,' ');
	in_stream >> seq_detect_enabled;
	in_stream.ignore(256,' ');
	in_stream >> measurement; //
	cout<<" ..... done"<<endl;	

	if (measurement > 0) {

		int hh = (int) (pow((float) valuedness,(float) resultnum));
		int g = (int) measurement * valuedness;
		int h = 0;
		in_stream.ignore(256,' ');
		for (h = 0; h < measurement; h++)
			in_stream >> measurementQBits[h];
		cout << "Reading measurements for " << h << " qubits\n";
		in_stream.ignore(256,' ');
		in_stream >> measured_desired_output_records;
		cout << "Reading " << measured_desired_output_records << " output patterns\n";
		if (seq_detect_enabled > 0){
			//read the sequence in as a set of 0 or 1 or MV values with appropriate probabilities
			sequence_desired = new cuComplex[measured_desired_output_records][2];
			for (int m = 0; m < measured_desired_output_records; m++) {
				in_stream >> unknown;
				if (!(unknown == '+' || unknown == '*'))
					a = unknown;

				in_stream >> tempchar;
				seq_type == 0;
				if (tempchar == '+'){
					seq_type == 1; //one or more
					in_stream >> tempchar;
				}
				else if (tempchar == '*'){
					seq_type == 2; //zero or more
					in_stream >> tempchar;
				}
				in_stream >> in;
				sequence_desired[m][0] =  make_cuFloatComplex(int(a-48),0);
				sequence_desired[m][1] =  make_cuFloatComplex(real(in),imag(in));
			}
			in_stream >> temp;
		} else {

			for (int b = 0; b < MAXGATEINPUT*2; b++) {
				cuComplex *cptr  = new cuComplex[hh];
				for (int c = 0; c < hh; c++) {
					cptr[c] = make_cuFloatComplex(0, -1);
				}
				measureexpected[b] = cptr;
			}
			//the input of measured qubits is a sequence of numbers separated by space
			int i = 0;
			measured_desired_output_records_idx = new int[measured_desired_output_records];
			for (int b = 0; b < measured_desired_output_records; b++) {
				cout << " reading ";
				in_stream >> outnum;
				cout << outnum << "  ";
				measured_desired_output_records_idx[b] = outnum;
				in_stream >> tempchar;
				cout << tempchar << "  ";
				measureexpected[outnum] = new cuComplex[measurement*valuedness];
				for (int h = 0; h < measurement; h++) {
					if (h > 0 ){
						in_stream >> tempchar;
					cout << tempchar;
					}
					cout<<endl;
					for (int g = 0; g < valuedness; g++) {
						in_stream >> in;
						cout<<" indexes: "<<outnum<<"  "<<(h * valuedness + g)<<":: ";
						measureexpected[outnum][h * valuedness + g]= make_cuFloatComplex(real(in), imag(in));
						cout <<"("<<cuCrealf(measureexpected[outnum][h*valuedness+g])<<", "<<cuCimagf(measureexpected[outnum][h*valuedness+g]) << ")  ";
					}
					cout<<endl;
				}
				cout<<"done\n";
			cout << endl;
			}
			in_stream >> buff_temp;
        		cout <<" ..... done"<<endl;
		}
	} else {
		int counter = 0;
		do {
			in_stream.getline(line, 512);
		} while (line[0] != '-');
	}

	//reading file divider
	loadGateArray(valuedness, representation);

	int w = (int(pow((float) gateArray[0]->valuedness, (float) resultnum)));
#ifdef __CUBLAS__
        int m_size = w * w * sizeof(cuComplex);
        int v_size = w * sizeof(cuComplex);
        if( cublasAlloc(m_size, sizeof(d_M1[0]), (void**)&d_M1) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_M2[0]), (void**)&d_M2) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_M3[0]), (void**)&d_M3) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if( cublasAlloc(m_size, sizeof(d_M4[0]), (void**)&d_M4) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_M5[0]), (void**)&d_M5) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_M6[0]), (void**)&d_M6) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_MV[0]), (void**)&d_MV) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(v_size, sizeof(d_VI[0]), (void**)&d_VI) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(v_size, sizeof(d_VO[0]), (void**)&d_VO) != CUBLAS_STATUS_SUCCESS) mem_error = true;

        if (cublasAlloc(m_size, sizeof(d_ME_00[0]), (void**)&d_ME_00) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_ME_01[0]), (void**)&d_ME_01) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_ME_10[0]), (void**)&d_ME_10) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(m_size, sizeof(d_ME_11[0]), (void**)&d_ME_11) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(v_size, sizeof(d_VV_M[0]), (void**)&d_VV_M) != CUBLAS_STATUS_SUCCESS) mem_error = true;
        if (cublasAlloc(1, sizeof(cuComplex), (void**)&d_Value) != CUBLAS_STATUS_SUCCESS) mem_error = true;
#else
#ifdef __CUDA__
	cout<<"allocating memory"<<endl;
	int m_size = w * w * sizeof(cuComplex);
	int v_size = w * sizeof(cuComplex);
	if (cudaMalloc((void**)&d_M1, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_M2, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_M3, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_M4, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_M5, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_M6, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_MV, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_VI, v_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_VO, v_size) != cudaSuccess) mem_error = true;

	if (cudaMalloc((void**)&d_ME_00, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_ME_01, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_ME_10, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_ME_11, m_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_VV_M, v_size) != cudaSuccess) mem_error = true;
	if (cudaMalloc((void**)&d_Value, sizeof(cuComplex)) != cudaSuccess) mem_error = true;
#endif
#endif

if(seq_detect_enabled){
#ifdef __CUBLAS__
	if ((cublasSetVector(m_size, sizeof(measures_o[0]->gateMatrix1[0]), measures_o[0]->gateMatrix1, 1, d_ME_00, 1)) != CUBLAS_STATUS_SUCCESS) mem_error = true;
	if ((cublasSetVector(m_size, sizeof(measures_o[1]->gateMatrix1[0]), measures_o[1]->gateMatrix1, 1, d_ME_01, 1)) != CUBLAS_STATUS_SUCCESS) mem_error = true;
	if ((cublasSetVector(m_size, sizeof(measures_i[0]->gateMatrix1[0]), measures_i[0]->gateMatrix1, 1, d_ME_10, 1)) != CUBLAS_STATUS_SUCCESS) mem_error = true;
	if ((cublasSetVector(m_size, sizeof(measures_i[1]->gateMatrix1[0]), measures_i[1]->gateMatrix1, 1, d_ME_11, 1)) != CUBLAS_STATUS_SUCCESS) mem_error = true;
#else
#ifdef __CUDA__
	if ((cudaMemcpy(d_ME_00, measures_o[0]->gateMatrix1, m_size, cudaMemcpyHostToDevice)) != cudaSuccess) mem_error = true;;
	if ((cudaMemcpy(d_ME_01, measures_o[1]->gateMatrix1, m_size, cudaMemcpyHostToDevice)) != cudaSuccess) mem_error = true;;
	if ((cudaMemcpy(d_ME_10, measures_i[0]->gateMatrix1, m_size, cudaMemcpyHostToDevice)) != cudaSuccess) mem_error = true;;
	if ((cudaMemcpy(d_ME_11, measures_i[1]->gateMatrix1, m_size, cudaMemcpyHostToDevice)) != cudaSuccess) mem_error = true;;
#endif
#endif
}
	if (mem_error) {
		cout<<" the circuit you are trying to desing does not fit onto the CUDA hardware,\n either reduce the number of bits or turn of the CUDA acceleration\n";
		exit(0);
	}

	cout<<"Extended gates generation\n";
	generateExtendMatrix(representation);
	cout << " .... done"<<endl;
#ifdef __QMDD__
        cout<<"Generating RevLib Gate Representation:";
	for (int a = 0; a < numofgates; a++) 
		decodeGateToRevLib(gateArray[a]);
        cout<<" ..... done"<<endl;
#endif
	if(measurement > 0)
		generateMeasureOps();
        cout << "Generating Output file"<<endl;


	//generate the header of the output file
	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << "         output file of synthesis: settings               "
	<< endl;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << " Size of the population:                                  "
	<< populationNumber << endl;
	out_stream << " Number of segments in each circuit (approx):             "
	<< segments << endl;
	out_stream << " Minimal Number of segments in each circuit (approx):     "
	<< mincircuitsize << endl;
	out_stream << " Multiplier of the max number of segments:                "
	<< segments_multi << endl;
	out_stream << " Number of GA generations before EX search starts:        "
	<< alterations << endl;
	out_stream << " Total number of GA generations:                          "
	<< generations << endl;
	out_stream << " Mutation probability is:                                 "
	<< proba_Mutation << endl;
	out_stream << " Crossover probability is:                                "
	<< proba_Crossover << endl;
	out_stream << " Factor alpha is (for complex fitness):                   "
	<< alpha << endl;
	out_stream << " Factor beta is (for complex fitness):                    "
	<< beta << endl;
	out_stream << " Factor alpha1 is (for pareto replication):               "
	<< alpha1 << endl;
	out_stream << " Factor beta1 is  (for pareto replication):               "
	<< beta1 << endl;
	out_stream << " Estimated minimum cost of the goal circuit:              "
	<< divider << endl;
	out_stream << " Phase switch (0 - off, 1 - on):                          "
	<< phase << endl;
	out_stream << " Display level:                                           "
	<< displaylevel << endl;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << " Type of GA (0 - normal, 1 - Darwinian):                  "
	<< Ga << endl;
	out_stream << " Use Elitism (0 - no, 1 - yes):                           "
	<< elitism << endl;
	out_stream << " Type of mutation (0 - normal, 1 - bitwise):              "
	<< mutation << endl;
	out_stream << " Type of crossover (0 - 1point, 1 - 2point):              "
	<< crossover << endl;
	out_stream << " Type of replication (0 - RW, 1 - SUS):                   "
	<< replicator << endl;
	out_stream << " Type of fitness (0 - simplest, 3 - complex):             "
	<< fitness << endl;
	out_stream << " Fitness calculation (0 - individual, 1 - grouped):       "
	<< grouped << endl;
	out_stream << " Pareto optimization (0 - no, 1 - yes):                   "
	<< pareto << endl;
	out_stream << " Max-Fitness Threshold for Replication (0 - no, 1< >0):   "
	<< threshold << endl;
	out_stream << " The number of wires the final circuit has:               "
	<< resultnum << endl;
	out_stream << " Valuedness of designed logic:                            "
	<< valuedness << endl;
	out_stream << " Design of Sequential Logic/ Sequence Detector:           "
	<< seq_detect_enabled <<endl;
	out_stream << " Measurement used:                                        "
	<< measurement << endl;



	if (measurement > 0) {
		out_stream
		<< " Measured bits are: ";
		for (int a = 0; a < measurement; a++)
			out_stream << measurementQBits[a] << " ";
		out_stream << endl;

		
		if (seq_detect_enabled > 0){
			for (int c = 0; c < measured_desired_output_records; c++) {
					out_stream << cuCrealf(sequence_desired[c][0])<<": ";
					out_stream << " (" << cuCrealf(sequence_desired[c][1]) << ","<< cuCimagf(sequence_desired[c][1]) << ")";
					out_stream << endl;
				}
				out_stream << endl;
		} else {
			for (int c = 0; c < measured_desired_output_records; c++) {
				out_stream <<c<<": "<<measured_desired_output_records_idx[c]<< ": ";
				for (int b = 0; b < measurement * valuedness; b++){
					out_stream << " (" << cuCrealf(measureexpected[measured_desired_output_records_idx[c]][b]) << ","<< cuCimagf(measureexpected[measured_desired_output_records_idx[c]][b]) << ")";
				}
				out_stream << endl;
			}
		}

	}
	
//	cout << " Writting Output " << endl;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << "            The radix of this logic is:                   "
	<< valuedness << endl;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << "            Number of Input Gates:                        "
	<< numofgates << endl;
	for (int a = 0; a < numofgates; a++) {
		if (gateArray[a]->numIO > 5){
			out_stream<<"Bellow gates are truncated"<<endl;
			out_stream<<"All two or more qubit gates are expanded in the program to the maximum width and are used in the synthesizer"<<endl;
			break;
		}
		out_stream << "Representation: " << gateArray[a]->representation[0]<<gateArray[a]->representation[1]<<gateArray[a]->representation[2]
		<< endl;
		out_stream << " Cost of this gate:  " << gateArray[a]->Cost << endl;
		out_stream << " IO of this gate:  " << gateArray[a]->numIO << endl;
		out_stream << " Real IO of this gate:  " << gateArray[a]->realIO << endl;
		out_stream << " Total number of restrictions of this gate:  " << gateArray[a]->restrictions_number << endl;
		out_stream << " Restrictions of this gate:  ";
		for (int g = 0; g <  gateArray[a]->restrictions_number; g++)
			 out_stream<<"  "<< gateArray[a]->restrictions[g];
		out_stream << endl;
		out_stream << " Allocated wires by this gate:";
		for (int g = 0; g <  gateArray[a]->realIO; g++)
			out_stream <<" "<< gateArray[a]->connections[g];
		out_stream << endl;
		out_stream << " The name of this gate is: " << gateArray[a]->my_string<< endl;
		out_stream << " The RevLib representation of this gate is: " << gateArray[a]->rev_lib_string<< endl;
		columns = (int(pow((float) valuedness,(float) gateArray[a]->numIO)));
		for (int j = 0; j < columns; j++) {
			for (int k = 0; k < columns; k++) {
				out_stream << "  ("
				<< cuCrealf(gateArray[a]->gateMatrix1[j+k*columns]) << ","
				<< cuCimagf(gateArray[a]->gateMatrix1[j+k*columns]) << ")";
			}
			out_stream << endl;
		}
	}
	out_stream << " -------- Generating the desired output evaluation ------ "
	<< endl;
	if (measurement == 0) {

		out_stream << " Desired Gate: " << endl << endl;
		out_stream << " IO of the gate:  " << finalGate->numIO << endl;
		//cout<<endl;
		columns = (int(pow((float) valuedness, (float) resultnum)));
		for (int m = 0; m < columns; m++) {
			for (int n = 0; n < columns; n++) {
				out_stream << "  (" << cuCrealf(finalGate->gateMatrix1[m+n*columns])
				<< "," << cuCimagf(finalGate->gateMatrix1[m+n*columns]) << ")";
			}
			out_stream << endl;
		}
	} else {

		out_stream << " Observables:  " << endl;
		out_stream << " Desirded Qubits of the whole system: "
		<< finalGate->numIO << endl;
		out_stream << " Desirded Valuedness: " << finalGate->valuedness << endl;
		out_stream << endl;

		//output the measurement operators
		float numIO = (float) finalGate->numIO;
		float val = (float) finalGate->valuedness;
		int mes_index = 0;
		int c = measurementQBits[mes_index] * val;
		int val_count = 0;
		for (int a = 0; a < val * measurement; a++) {
			out_stream << "Measurement Gates:" <<endl;
			cout << "Measurement Gates:" << c<<endl;
			if (numIO > 5){
				out_stream<<"Measurement gates truncated due to large size"<<endl;
				break;
			}
			columns = (int(pow(val, (float) measurements[c]->numIO)));
			out_stream << measurements[c]->my_string << endl;
			out_stream << "  IO of the gate  " << measurements[c]->numIO << endl;

			for (int j = 0; j < columns; j++) {
				for (int k = 0; k < columns; k++) {
					out_stream << "  (" << cuCrealf(measurements[c]->gateMatrix1[j+k*columns]) << "," << cuCimagf(measurements[c]->gateMatrix1[j+k*columns]) << ")";
				}
				out_stream << endl;
			}
			c++;
			val_count++;
			if (val_count >= val){
				val_count = 0;
				c = measurementQBits[++mes_index]*val;
			}
		}

	}

	delete [] representation;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << "                   Output of synthesis:                   "
	<< endl;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	cout<<" ..... done"<<endl;
}
/*******************************************
 * Create the measurement oerators for Single Qubits
 * for all orthonormal states in a given valuedness
 *******************************************/
void* GA::generateMeasureOps() {

	float numIO = (float) finalGate->numIO;
	float val = (float) finalGate->valuedness;
	int val_count = 0;
	bool zeros = true;
	char *representation = new char[4];
	int a = 0;
	int m = 0;
	int columns;
	int qubitindex = 0;
	string srank = "";
	int current_step = (int) (pow(val, (float) (qubitindex + 1)));
	int step_length = (int) (pow(val, (float) qubitindex));

	//generate measurement opertors for every single qubit
	cout << "Generating measurement circuits: " << MAXMEASUREGATE << "   "<< numIO<<endl;
	initRep(representation, val_count);
	while (a < resultnum * val) {
		convert_stream.str("");
		convert_stream << qubitindex;
		srank = convert_stream.str();
		cout<<"rank: "<<srank<<", "<<qubitindex<<endl;
		measurements[a] = new qGate;
		initGate(measurements[a], numIO, val, 1);
		copyRep(representation, measurements[a]->representation);
		//measurements[a]->representation = char(val_count);
		measurements[a]->my_string += "th qubit";
		measurements[a]->my_string = srank + measurements[a]->my_string;

		columns = (int) (pow(val, numIO));
		//measurements[a]->gateMatrix1 = new cuComplex[columns*columns];
		for (int x = 0; x < columns; x++) {
			for (int y = 0; y < columns; y++) {
				measurements[a]->gateMatrix1[x+y*columns] = make_cuFloatComplex(0, 0);
				if (x == y) {
					if (step_length > 0) {
						measurements[a]->gateMatrix1[x+y*columns]= make_cuFloatComplex(1, 0);
						step_length--;
					}
					current_step--;
					if (current_step <= 0) {
						current_step = (int) (pow((float) val,
								(float) (qubitindex + 1)));
						step_length = (int) (pow((float) val,
								(float) qubitindex));
					}
				}
			}
		}
		incRep(representation, 0, 10);
		val_count++;
		if (val_count >= val) {
			val_count = 0;
			qubitindex++;
			step_length = (int) (pow((float) val, (float) qubitindex));
			current_step = (int) (pow((float) val, (float) (qubitindex + 1)));
		} else {
			//change the way how it starts for the next observable of this qubit
			step_length = 0;
			//for each observable shift the matrices index
			current_step = val_count * (int) (pow((float) val,
					(float) qubitindex));
		}
		a++;
	}
	cout << " ..... done" << endl;
	delete [] representation;

}
/*******************************************
 * Initialize the expectations probabilities for measured global states
 *******************************************/
void GA::generateMeasureExpect() {

	int l = finalGate->numIO;
	int p = (int(pow((float) 2, (float) l)));
	for (int g = 0; g < p; g++) {
		cuComplex *cptr = new cuComplex[2];
		cptr[0] = make_cuFloatComplex(1, 0);
		cptr[1] = make_cuFloatComplex(1, 0);
		expectationsAllState[g] = cptr;
	}

	for (int k = 0; k < measured_desired_output_records; k++) {
			cout << " For input state: "<<measured_desired_output_records_idx[k] <<endl;
		for (int m = 0; m < measurement; m++) {
			//check if entanglement is desired
			if (cuCrealf(measureexpected[k][2*m]) == cuCrealf(
					measureexpected[k][2*m+1])) {
				if (cuCrealf(measureexpected[k][2*m]) != 1 && cuCrealf(
						measureexpected[k][2*m]) != 0) {
					expectationsAllState[k][0] = cuCmulf(
							expectationsAllState[k][0],
							measureexpected[k][2*m]);
					expectationsAllState[k][1] = cuCmulf(
							expectationsAllState[k][1], measureexpected[k][2*m
							                                               +1]);
				}
			} else {
				if (cuCrealf(measureexpected[k][2*m]) > cuCrealf(
						measureexpected[k][2*m+1])) {
					expectationsAllState[k][0] = cuCmulf(
							expectationsAllState[k][0],
							measureexpected[k][2*m]);
					expectationsAllState[k][1] = cuCmulf(
							expectationsAllState[k][1], measureexpected[k][2*m
							                                               +1]);
				} else if (cuCrealf(measureexpected[k][2*m]) <= cuCrealf(
						measureexpected[k][2*m+1])) {
					expectationsAllState[k][0] = cuCmulf(
							expectationsAllState[k][0], measureexpected[k][2*m
							                                               +1]);
					expectationsAllState[k][1] = cuCmulf(
							expectationsAllState[k][1],
							measureexpected[k][2*m]);
				}
			}
		}
	}
}
/*******************************************
 * Create all mullti state measurements for desired output states
 * This is designed so as we can distinguish  between the most wanted state
 * and all other ones, thus creating E_0 and I-E_0
 *******************************************/
void* GA::generateMultiMeasureOps() {
	int l = finalGate->numIO;
	float numIO = (float) finalGate->numIO;
	//the counter for measurement matrices is implemented as follows
	//2^(abs) index - numNBGR-1 (to get to the index array values)
	//int k = (int)pow(2.0, (float)(abs(index - (tRule->numNGBR-1))));
	int n, m, columns = 0;
	bool zeros = true;
	int val = (float) finalGate->valuedness;
	int a = (int) (resultnum * val);
	string srank;
	qGate *measure, *measure0, *measure1, *measureM0, *measureM1, *temp;
	int p = (int(pow((float) val, (float) l)));

	//generate measurement opertors for every single qubit
	//cout<< "generating measurement circuits: "<< MAXMEASUREGATE;
	cout << "Generating measurements for all desired states" << endl;
	for (int b = 0; b < measured_desired_output_records; b++) {
		n = 0;

		convert_stream.flush();
		measure0 = new qGate();
		initGate(measure0, resultnum, val, 1);
		measure1 = new qGate();
		initGate(measure1, resultnum, val, 1);
		measureM0 = new qGate();
		initGate(measureM0, resultnum, val, 1);
		measureM1 = new qGate();
		initGate(measureM1, resultnum, val, 1);
		initMatrix(measureM1, resultnum);
		temp = new qGate();
		initGate(temp, resultnum, val, 1);
		//for every orthonormal state on the input
		//generate two measurement operators in the form E_0 and I-E_0
		for (int k = 0; k < a; k++) {

			initMatrix(measureM0, resultnum);
			//initMatrix(&measureM1, resultnum);
			for (int m = 0; m < (measurement * val); m++) {
				measure0 = ccGate(measurements[(measurementQBits[m] * val)],
						measure0);
				//check if entanglement is desired
				if (cuCrealf(measureexpected[k][val * m]) == cuCrealf(
						measureexpected[k][2* m + 1])) {
					if (cuCrealf(measureexpected[k][2* m ]) == 1) {
						initMatrix(temp, resultnum);
						cblasMatrixProduct(measureM0, measure1, temp);
						measureM0 = ccGate(temp, measureM0);
					} else if (cuCrealf(measureexpected[k][2* m ]) == 0) {
						if (!(cuCimagf(measureexpected[k][2* m ]) == 1
								&& cuCimagf(measureexpected[k][2* m + 1]) == 1)) {
							initMatrix(temp, resultnum);
							cblasMatrixProduct(measureM0, measure0, temp);
							measureM0 = ccGate(temp, measureM0);
						}
					} else {
						if (!(cuCimagf(measureexpected[k][2* m ]) == 1
								&& cuCimagf(measureexpected[k][2* m + 1]) == 1)) {
							initMatrix(temp, resultnum);
							cblasMatrixProduct(measureM0, measure0, temp);
							measureM0 = ccGate(temp, measureM0);
						}
					}
				} else {
					if (cuCrealf(measureexpected[k][2* m ]) > cuCrealf(
							measureexpected[k][2* m + 1])) {
							initMatrix(temp, resultnum);
							cblasMatrixProduct(measureM0, measure0, temp);
							measureM0 = ccGate(temp, measureM0);
					} else if (cuCrealf(measureexpected[k][2* m ]) <= cuCrealf(
							measureexpected[k][2* m + 1])) {
							initMatrix(temp, resultnum);
							cblasMatrixProduct(measureM0, measure1, temp);
							measureM0 = ccGate(temp, measureM0);
					}
				}
			}

			measureM1 = matrixSubstr(measureM1, measureM0);
			convert_stream << k;
			srank = convert_stream.str();
			measurementsAllState[2* k ] = new qGate();
			initGate(measurementsAllState[2* k ], resultnum, val, 1);
			measurementsAllState[2* k ] = ccGate(measureM0,	measurementsAllState[2* k ]);
			measurementsAllState[2* k ]->my_string = srank
			+ " Desired input state";
			measurementsAllState[2* k + 1] = new qGate();
			initGate(measurementsAllState[2* k +1], resultnum, val, 1);
			measurementsAllState[2* k + 1] = ccGate(measureM1, measurementsAllState[2* k + 1]);
			measurementsAllState[2* k + 1]->my_string = srank
			+ " Undesired state ";
		}
	}

	//output the measurement operators
	for (int a = 0; a < p; a++) {
		out_stream << "Whole state Measurement Gates:" << endl;
		out_stream << measurementsAllState[2* a ]->my_string << endl;
		out_stream << measurementsAllState[2* a ]->representation << endl;
		out_stream << "  IO of the gate  "
		<< measurementsAllState[2* a ]->numIO << endl;
		columns =  (int(pow((float) 2, (float) measurementsAllState[2*a ]->numIO)));
		for (int j = 0; j < columns; j++) {
			for (int k = 0; k < columns; k++) {
				out_stream << "  (" << cuCrealf(measurementsAllState[2* a ]->gateMatrix1[j+k*columns]) << ","<< cuCimagf(measurementsAllState[2* a ]->gateMatrix1[j+k*columns])<< ")";
			}
			out_stream << endl;
		}

		out_stream << measurementsAllState[2* a + 1]->my_string << endl;
		out_stream << measurementsAllState[2* a + 1]->representation << endl;
		out_stream << "  IO of the gate  "
		<< measurementsAllState[2* a + 1]->numIO << endl;
		columns = (int(pow((float) 2, (float) measurementsAllState[2*a + 1]->numIO)));
		for (int j = 0; j < columns; j++) {
			for (int k = 0; k < columns; k++) {
				out_stream << "  (" << cuCrealf(measurementsAllState[2* a + 1]->gateMatrix1[j+k*columns])<< "," << cuCimagf(measurementsAllState[2* a + 1]->gateMatrix1[j+k*columns])<< ")";
			}
			out_stream << endl;
		}

	}
	destroyGate(measure0);
	destroyGate(measure1);
	destroyGate(measureM0);
	destroyGate(measureM1);
	destroyGate(temp);
	delete (measure0, measure1, measureM1, measureM0, temp);

}
/*******************************************
 * Verifies the input gates and generates all missing ones
 * for multiple qubits
 *******************************************/
void GA::generateExtendMatrix(char* representation) {

	int counter = 0;
	int l_count = 0;
	int i = 0;
	bool k = false;
	string name;
	qGate *gate, *extgate, *tempgate, *sgate, *comparegate;
	rGate rgate;
	int controls[MAXINDVWIRENUM];
	int i_placement[MAXINDVWIRENUM];
	int columns, initial, wire_counter;
	int c_count = 0;
	int topwire = 2;
	int function = 0;
	bool found;

	gate = new qGate;
	initGate(gate, resultnum, gateArray[0]->valuedness, 1);
	tempgate = new qGate;
	initGate(tempgate, resultnum, gateArray[0]->valuedness, 1);
	comparegate = new qGate;
	initGate(comparegate, resultnum, gateArray[0]->valuedness, 1);
	extgate = new qGate;
	initGate(extgate, resultnum, gateArray[0]->valuedness, 1);
	sgate = new qGate;
	initGate(sgate, 2, gateArray[0]->valuedness, 1);
	initial = 1;
	initMatrix(sgate, 2);
	columns = (int(pow((float)sgate->valuedness,(float)sgate->numIO)));
	if (gateArray[0]->valuedness == 2){
		sgate->gateMatrix1[1+1*columns] = make_cuFloatComplex(0, 0);
		sgate->gateMatrix1[1+2*columns] = make_cuFloatComplex(1, 0);
		sgate->gateMatrix1[2+1*columns] = make_cuFloatComplex(1, 0);
		sgate->gateMatrix1[2+2*columns] = make_cuFloatComplex(0, 0);
	} else if (gateArray[0]->valuedness == 3){
		zeroMatrix(sgate, 2);
	        sgate->gateMatrix1[0] = make_cuFloatComplex(1, 0);
        	sgate->gateMatrix1[12] = make_cuFloatComplex(1, 0);
	        sgate->gateMatrix1[24] = make_cuFloatComplex(1, 0);
        	sgate->gateMatrix1[28] = make_cuFloatComplex(1, 0);
	        sgate->gateMatrix1[40] = make_cuFloatComplex(1, 0);
        	sgate->gateMatrix1[52] = make_cuFloatComplex(1, 0);
	        sgate->gateMatrix1[56] = make_cuFloatComplex(1, 0);
        	sgate->gateMatrix1[68] = make_cuFloatComplex(1, 0);
	        sgate->gateMatrix1[80] = make_cuFloatComplex(1, 0);
	}
	while (true)
		if (topwire < resultnum) {
			topwire++;
			zeroMatrix(tempgate, extgate->numIO);
			zeroMatrix(comparegate, extgate->numIO);
			counter = 0;
			for (int a = numofgates - 1; a >= 0; a--) {
				if (gateArray[a]->numIO == topwire - 1) {
					c_count = 0;
					wire_counter = 0;
					gate->valuedness = gateArray[a]->valuedness;
					initMatrix(gate, gateArray[a]->numIO);
					gate = ccGate(gateArray[a], gate);
					cout<<"Initial gate: "<<gate->my_string<<endl;
					for (int h = 0 ; h < MAXINDVWIRENUM; h++)
						controls[h] = 0;
					
					for (int m = 0; m < resultnum; m++){
						controls[m] = -10;i_placement[m] = -10;
					}
					name = "";
					function = -10;
					for (int m = 0; m < gate->my_string.length(); m++){
						name += gate->my_string.at(m);
						if (gate->my_string.at(m) == 'C'){
							i_placement[c_count] = m;
							controls[c_count] = wire_counter;
							c_count++;
							wire_counter++;
						} else if (gate->my_string.at(m) == 'I'){
							wire_counter++;
						} else {
							if (function == -10){
								function =  wire_counter;
								wire_counter++;
							}
						}
					}
#ifdef __CUDA__
					Kronm(gateArray[0], gate, (int)(pow((float)(gateArray[0]->valuedness), (float)(gateArray[0]->numIO))), (int)(pow((float)(gate->valuedness), (float)(gate->numIO))), tempgate, d_M1, d_M2, d_M3);
#else
					tensorProduct3(tempgate, gateArray[0] ,gate);
#endif
					
					cout<<"0 gates are equal: "<<compareGate(tempgate, comparegate)<<endl;
					gate = ccGate(tempgate, gate);
					cout<<endl;
					for (int u = 0; u < c_count; u++){
						initMatrix(extgate, resultnum);
						for (int l = 0; l < topwire; l++) {
							if (u == 0 && l == 0 ){
								extgate = ccGate(sgate, extgate);
								l++;
							}else if (l == 0){
								extgate = ccGate(gateArray[0], extgate);
							}else if( u == l){
#ifdef __CUDA__
								Kronm(extgate, sgate, (int)(pow((float)(extgate->valuedness), (float)(extgate->numIO))), (int)(pow((float)(sgate->valuedness), (float)(sgate->numIO))), tempgate, d_M1, d_M2, d_M3);
#else
								tensorProduct3(tempgate, extgate, sgate);
#endif
								cout<<"1 these gates are equal: "<<compareGate(tempgate, comparegate)<<endl;
								extgate = ccGate(tempgate, extgate);
								l++;
							} else {
#ifdef __CUDA__
								Kronm(extgate, gateArray[0], (int)(pow((float)(extgate->valuedness), (float)(extgate->numIO))), (int)(pow((float)(gateArray[0]->valuedness), (float)(gateArray[0]->numIO))),  tempgate, d_M1, d_M2, d_M3);
#else
								tensorProduct3(tempgate, extgate, gateArray[0]);
#endif
								cout<<"2 these gates are equal: "<<compareGate(tempgate, comparegate)<<endl;
								extgate = ccGate(tempgate, extgate);
							}
						}
						extgate->numIO = topwire;
						zeroMatrix(tempgate, extgate->numIO);
						if (extgate->numIO <= BLOCK_MAX_SIZE){
							cblasMatrixProduct(extgate, gate, tempgate);
						} else {
#ifdef __CUBLAS__
						cuBlasMMultiSet((int)(pow((float)(extgate->valuedness), (float)(extgate->numIO))), extgate->gateMatrix1, gate->gateMatrix1, tempgate->gateMatrix1, d_M1, d_M2, d_M3);
#else 
#ifdef __CUDA__
						Mulm(extgate->gateMatrix1, gate->gateMatrix1, (int)(pow((float)(extgate->valuedness), (float)(extgate->numIO))), tempgate->gateMatrix1, d_M1, d_M2, d_M3);
#else
						cblasMatrixProduct(extgate, gate, tempgate);
#endif
#endif
						}

						gate = ccGate(tempgate, gate);

						zeroMatrix(tempgate, extgate->numIO);
						if (extgate->numIO <= BLOCK_MAX_SIZE){
                                                        cblasMatrixProduct(gate, extgate, tempgate);
                                                } else {
#ifdef __CUBLAS__
						cuBlasMMultiSet((int)(pow((float)(extgate->valuedness), (float)(extgate->numIO))), gate->gateMatrix1, extgate->gateMatrix1, tempgate->gateMatrix1, d_M1, d_M2, d_M3);
#else 
#ifdef __CUDA__
						Mulm(gate->gateMatrix1, extgate->gateMatrix1, (int)(pow((float)(extgate->valuedness), (float)(extgate->numIO))), tempgate->gateMatrix1, d_M1, d_M2, d_M3);
#else
						cblasMatrixProduct(gate, extgate, tempgate);
#endif
#endif
						}
						gate = ccGate(tempgate, gate);
						gate->numIO = topwire;
						gate->my_string = "";

						//expand the name by inserting a I
						for (int m = 0; m < name.length(); m++)
							gate->my_string += name.at(m);
						if (i_placement[u] > function)
							gate->my_string.insert(i_placement[u], "I");
						else 
							gate->my_string.insert(i_placement[u]+1, "I");
						//count the wires this gate is located on
						l_count = 0;
						i = 0;
						k = false;
						for (int m = 0; m < gate->my_string.length(); m++){
							if (gate->my_string.at(m) == 'C'){
								k = false;
								gate->connections[i++] = l_count++;
							}
							else if (gate->my_string.at(m) == 'I'){
								l_count++;
								k = false;
							} else if (!k){
								gate->connections[i++] = l_count++;
								k = true;
							}
						}
						
						gate->realIO = i;

						//check if such gate has already been designed
						//possible with MCT like gates
						found = false;
						for (int g = 0 ; g < (numofgates + counter); g++)
							if (gateArray[g]->my_string == gate->my_string){
								found = true;
								break;
							}	
						if (!found){
							gateArray[numofgates + counter] = new qGate;
							rgate.index = (numofgates + counter);
							rankV.push_back(rgate);
							initGate(gateArray[numofgates + counter], topwire, gateArray[0]->valuedness, 1);
							gate->restrictions_number = gateArray[a]->restrictions_number;
							for (int h = 0; h < gate->restrictions_number; h++)
									gate->restrictions[h] = gateArray[a]->restrictions[h];
							gateArray[numofgates + counter] = ccGate(gate,
								gateArray[numofgates + counter]);
							gateArray[numofgates + counter]->Cost = 2;
							copyRep(representation, gateArray[numofgates + counter]->representation);
							incRep(representation);
							counter++;
						}
					}
				}
			}
			numofgates += counter;
		} else
			break;
	destroyGate(extgate);
	destroyGate(tempgate);
	destroyGate(gate);
	destroyGate(sgate);
	delete extgate, tempgate, gate, sgate;
cout<<"done"<<endl;
}
/***********************************
 * reads all the gates from the input file
 * to the storage
 *********************************/
void GA::loadGateArray(int valuedness, char *representation){

	string s;
	int columns,counter,resnum,a,b,c, ycounter,indexa;
	char tempchar;
	string temp;
	char binary[100];
	char buff[100];
	char line[512];
	rGate rgate;
        complex<float> in = complex<float> (0., 0.);
	//bedining reading the gates
	in_stream >> numofgates;
	cout<<"Reading in initial "<<numofgates<<" gates";
	//char chcounter = charcounter;
	// get first gate

	for (int a = 0; a < numofgates; a++) {
		gateArray[a] = new qGate;
		rgate.index = a;
		rankV.push_back(rgate);
		in_stream >> gateArray[a]->numIO;
		gateArray[a]->realIO = gateArray[a]->numIO;
		for (int b = 0; b < gateArray[a]->numIO; b++)
			gateArray[a]->connections[b] = b;
		in_stream >> gateArray[a]->Cost;
		gateArray[a]->valuedness = valuedness;
		initGate(gateArray[a], gateArray[a]->numIO,  gateArray[a]->valuedness, gateArray[a]->Cost);
		copyRep(representation, gateArray[a]->representation);
		incRep(representation);
		indexa = 0;
		do{
			in_stream>>c;
			gateArray[a]->restrictions[indexa++] = c;
		}while (in_stream.peek() != '\n');
		gateArray[a]->restrictions_number = indexa;

		in_stream >> gateArray[a]->my_string;
		columns = (int(pow((float) valuedness,(float) gateArray[a]->numIO)));
		gateArray[a]->gateMatrix1 = new cuComplex[columns*columns];

		for (int j = 0; j < columns; j++) {
			for (int k = 0; k < columns; k++) {
				in_stream >> in;
				gateArray[a]->gateMatrix1[j+k*columns] = make_cuFloatComplex(real(in),imag(in));
			}
		}
		counter++;
	}


	cout<<"... done"<<endl;
	cout<<"Reading in the target gate"<<endl;;
	finalGate = new qGate;
	finalIndividual = new Individual;
	initGate(finalGate, resultnum, valuedness, 1);
	if (measurement == 0) {
		// get the expected result matrix from the file
		in_stream >> resnum;
		in_stream >> measured_desired_output_records;
		if (measured_desired_output_records == 0) {
			columns = (int(pow((float) valuedness, (float) resnum)));
			for (int m = 0; m < columns; m++) {
				for (int n = 0; n < columns; n++) {
					in_stream >> in;
					finalGate->gateMatrix1[m+n*columns] = make_cuFloatComplex(real(in), imag(in));
				}
			}
		} else {
			
#ifdef __QMDD__
				temp = "";
				temp += ".version Target Individual\n.numvars ";
				sprintf(buff, "%d", resultnum);
				counter = 0;
				while(buff[counter] > 0)
					temp += buff[counter++];	
				temp += "\n.variables ";
				for (ycounter = 0; ycounter < resnum; ycounter++){
					sprintf(buff, "%d", ycounter);
					counter = 0;
					temp += 'w';	
					while(buff[counter] > 0)
						temp += buff[counter++];	
					temp += ' ';	
				}
				temp += "\n.inputs ";
				for (ycounter = 0; ycounter < resnum; ycounter++){
					sprintf(buff, "%d", ycounter);
					counter = 0;
					temp += 'w';	
					while(buff[counter] > 0)
						temp += buff[counter++];	
					temp += ' ';	
				}
				temp += "\n.outputs ";
				for (ycounter = 0; ycounter < resnum; ycounter++){
					sprintf(buff, "%d", ycounter);
					counter = 0;
					temp += 'w';	
					while(buff[counter] > 0)
						temp += buff[counter++];	
					temp += ' ';	
				}
				temp += "\n.begin\n";
#endif
				dcMatrix(finalGate, resnum);
				columns = (int(pow((float) valuedness, (float) resnum)));
				non_measured_desired_output_records_idx = new int[2*measured_desired_output_records];
				for (int m = 0; m < measured_desired_output_records; m++) {
					s = "";
					in_stream >> a;
					non_measured_desired_output_records_idx[2*m] = a;
					in_stream >> tempchar;
					in_stream >> b;
					non_measured_desired_output_records_idx[2*m+1] = b;
					dec2bin(b, binary, resnum);
					cout <<binary<<";; "<< a << " " << tempchar << " " << b;
					in_stream >> tempchar;
					in_stream >> in;
					cout << " " << tempchar << " " << in << endl;
					finalGate->gateMatrix1[a+b*columns] = make_cuFloatComplex(real(in),imag(in));
					counter = 0;
#ifdef __QMDD__
					while(binary[counter] > 0)
						temp += binary[counter++];
					temp += "\n";
#endif
				}

				counter = 0;
#ifdef __QMDD__
				temp += ".end\n";
				while (counter < temp.size())
					finalIndividual->rev_lib_string += temp.at(counter++);	
#endif
			}
//		cout << temp<<endl;
	}
	cout<<" ..... done"<<endl;
//	return representation;
}
/***********************************
 * given a three character gate representation
 * returns its index in the storage array
 *********************************/
int GA::getGate(char *b) {
	for (int a = 0; a < numofgates; a++)
		if (compareRep(gateArray[a]->representation, b))
			return a;
	return -1;
}
/****************************************
 * top function to choose mutation
 ***************************************/
void GA::applyMutation() {
	if (mutation == 0)
		for (int a = 0;a<populationNumber;a++)
			applyRMutationP(a);
	else if (mutation == 1)
		for (int a = 0;a<populationNumber;a++)
			applyIMutationP(a);
	else if (mutation == 2)
		for (int a = 0;a<populationNumber;a++)
			applyRMutation(a);
	else if (mutation == 3)
		for (int a = 0;a<populationNumber;a++)
			applyIMutation(a);
}
/****************************************
 * high level circuit optimization
 ***************************************/
void GA::applyOptimization(bool repair) {
	for (int a = 0;a<populationNumber;a++){
		minimizeCirc(population[a]);
		if (numofgates < MAXNUMOFGATES)
			if(Ga > 1)
				restructureCircuit(population[a], true);
			else if (Ga > 0)
				restructureCircuitD(population[a]);
//cout<<"1.5-opt: "<<population[a]->valuedness<<endl;
		if (repair)
			repairCircuit(a);
//cout<<"2-opt: "<<population[a]->valuedness<<endl;
	}
}
/****************************************
 * top function to  choose crossover
 ***************************************/
void GA::applyCrossover(int indv1, int indv2) {
	if (crossover == 0)
		apply1CrossoverP(indv1,indv2);
	else if (crossover == 1)
		apply2CrossoverP(indv1,indv2);
}

/*******************************************
 * givent the string name of a gate returns its index
 *******************************************/
int GA::getqGate(string b) {
	for (int a = 0; a < numofgates; a++) {
		if (b.length() == gateArray[a]->my_string.length()) {
			if ((gateArray[a]->my_string).compare(0,
					gateArray[a]->my_string.length(), b) >= 0)
				return a;
		}
	}
	return -1;
}
/*******************************************
 * helper function - returns the top wires that a gate
 *  resides on  in the circuit
 *******************************************/
int GA::getTopWire(int pos, string str){
	int wires[MAXGATEINPUT];
	char gaterep[REPSIZE];
	qGate * qgate;
	int end = str.find('p',pos);
	int begin = str.rfind('p',pos)+1;
	int counter = 0;
	while (str.at(begin) == 'p')begin++;
	while (str.at(begin) != 'p'){
		getGateRepfromStr(begin, str,gaterep);
		qgate = gateArray[getGate(gaterep)];
		if (begin+3 > pos){
			return (qgate->connections[0]+counter);
		} else {
			counter+= qgate->connections[qgate->realIO-1];	
		} 
		begin += 3;
	}
	return -1;
} 
/***********************************
 * apply the mutation on the individual of the current population
 * take one gate and change it to something else, by preserving
 * the size of the initial gate - simple replacement
 *********************************/
void GA::applyRMutationP(int indv) {

	Individual *I = population[indv];

	unsigned int proba, proba2, ios, pos, top_wire;
	char *gate, gchar;
	qGate *temp;
	gate = new char[3];
	proba = (int) (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))));
	if (proba <= proba_Mutation) {
		proba = rand() % 1;
		//change only the phase of the circuit
		if (phase > 0) {
			if (proba == 0) {
				I->phase = generatePhase();
				return;
			}
		}
		do {
			proba2 = (rand() % (I->my_string.length() - 1));
			if (proba2 < 0)
				proba2 = 0;
			if (proba2 >= I->my_string.length() - 1)
				proba2 = I->my_string.length() - 1;
			gchar = I->my_string.at(proba2);
		} while (gchar == 'p');

		pos = getGatePos(proba2, I->my_string);
		getGateRepfromStr(pos, I->my_string, gate);
		top_wire = getTopWire(pos, I->my_string);
		temp = gateArray[getGate(gate)];
		ios = temp->numIO;
		proba = getRandomGate(ios, true, top_wire,0);
		temp = gateArray[proba];
		I->my_string.replace(pos, 3, temp->representation, 3);
	}
	delete [] gate;
}
/***********************************
 * apply the mutation on the individual of the current population
 * on every gate apply the mutation operation
 * preserve the structure of the circuit
 *********************************/
void GA::applyIMutationP(int indv) {

	Individual *I = population[indv];
	int proba, proba2, temp0, iterations, ios, pos, top_wire;
	char *gate, gchar;
	qGate *temp;
	gate  = new char[3];
	complex<float> tempc;

	temp0 = 0;
	for (unsigned int m = 0; m < I->my_string.size(); m++) {
		if (I->my_string.at(m) == 'p')
			temp0++;
	}

	I->segmentNumber = temp0;
	iterations = getGatesfromString(I->my_string);
	//gate wise mutation
	for (int a = 0; a < iterations; a++) {

		//generate probabilty for the next index
		proba = (int) (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))));

		//try iterator
		if (proba <= proba_Mutation) {
			//three parameters to possibly mutate
			if(phase > 0)
				proba = rand() % 3;
			else
				proba = 0;
			//change only the phase of the circuit
			switch (proba) {
			case 2:
				//modify the global phase
				I->phase = generatePhase();
				break;
			case 1:
				//mutate the n phases of the circuit
				proba2 = rand() % iterations;
				temp0 = 0;
				while (1) {
					if (temp0 < I->my_string.length()) {
						if (I->my_string[temp0] != 'p')
							proba2--;
					} else
						break;

					if (proba2 == 0 && I->my_string[temp0] != 'p') {
						I->phases[proba2] = generatePhase();
						break;
					}
					temp0++;
				}
				break;
			case 0:
				do {
					proba2 = (rand() % (I->my_string.length() - 1));
					if (proba2 < 0)
						proba2 = 0;
					gchar = I->my_string.at(proba2);
				} while (gchar == 'p');

				pos = getGatePos(proba2, I->my_string);
				getGateRepfromStr(pos, I->my_string, gate);
				temp = gateArray[getGate(gate)];
				top_wire = getTopWire(pos, I->my_string);
				ios = temp->numIO;
				proba = getRandomGate(ios, true, top_wire,0);
				temp = gateArray[proba];	
				I->my_string.replace(pos, 3, temp->representation, 3);
				break;
			}
		}

	}
	delete [] gate;
}


/***********************************
 * apply the mutation on the individual of the current population
 * take one segment and change it to something else
 *********************************/
void GA::applyRMutation(int indv) {
	Individual *I = population[indv];
	unsigned int proba, proba2;
	float proba_1;
	char gchar;
	qGate temp;
	proba_1 = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))));
	if (proba_1 <= proba_Mutation) {
		proba = rand() % 1;
		//change only the phase of the circuit
		if (phase > 0) {
			if (proba == 0) {
				I->phase = generatePhase();
				return;
			}
		}
		do {
			proba2 = (rand() % (I->my_string.length() - 1));
			if (proba2 < 0)
				proba2 = 0;
			gchar = I->my_string.at(proba2);
		} while (gchar != 'p' || proba2 > (I->my_string.length() - 3));

		if (I->my_string.at(proba2 + 1) == 'p')
			proba2++;
		int proba3 = I->my_string.find("p", proba2 + 1);

		I->my_string.erase(proba2, (proba3 - proba2 + 1));
		//tricky part is to keep the width of the circuit constant
		basic_string<char> tempS;
		int V = 0;
		do {
			temp = *gateArray[getRandomGate(0, true, V, (I->ioNumber - V))];
			if (temp.numIO <= (I->ioNumber - V)|| ((V+temp.numIO) > I->ioNumber)) {
				V += temp.numIO;
				tempS.insert(tempS.length(), temp.representation, 3);
			}
		} while (V != I->ioNumber);
		tempS.insert(0, "p");
		tempS.insert(tempS.length(), "p");

		I->my_string.insert(proba2, tempS);
	}
}
/***********************************
 *apply the mutation on the individual of the current population
 *on every segment apply the mutation operation
 *********************************/
void GA::applyIMutation(int indv) {
	Individual *I = population[indv];
	int proba2, temp0, iterations, wire_counter, tempg;
	float proba;
	char gate;
	qGate temp;
	complex<float> tempc;
	temp0 = 0;
	for (unsigned int m = 0; m < I->my_string.size(); m++) {
		if (I->my_string.at(m) == 'p')
			temp0++;
	}
	I->segmentNumber = temp0;
	iterations = getGatesfromString(I->my_string);
	//gate wise mutation
	for (int a = 0; a < iterations; a++) {
		//generate probabilty for the next index
		proba =  (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))));
		//try iterator
		if (proba <= proba_Mutation) {
			//mutate the string - single gate
			do {
				proba2 = (rand() % I->my_string.length()) - 1;
				if (proba2 < 0)
					proba2 = 0;
				gate = I->my_string.at(proba2);
			} while (gate != 'p' || proba2 > (I->my_string.length() - 3));

			if (I->my_string.at(proba2 + 1) == 'p')
				proba2++;
			int proba3 = I->my_string.find("p", proba2 + 1);

			I->my_string.erase(proba2, (proba3 - proba2 + 1));
			//tricky part is to keep the width of the circuit constant
			basic_string<char> tempS;
			int V = 0;
			wire_counter = 0;
			do {
				tempg = getRandomGate(0, true, V, (I->ioNumber - V));
				temp = *gateArray[tempg];
cout<<V<<" "<<temp.numIO<<". "<<(I->ioNumber - V)<<" ";
				if (((V+temp.numIO) <= I->ioNumber)) {
				//if (temp.numIO <= (I->ioNumber - V) || ((V+temp.numIO) = I->ioNumber)) {
					V += temp.numIO;
					tempS.insert(tempS.length(), temp.representation, 3);
				}
			} while (V != I->ioNumber);
			tempS.insert(0, "p");
			tempS.insert(tempS.length(), "p");
			I->my_string.insert(proba2, tempS);
		}
	}
}

/**************************
 * randomly seleects a gate with additional constraints
 * wires - specifies precisely over how many wires must be defined
 * restr - specifies if the restrictions should be checked
 * topwire - parameter for the restriction checking
 * limit - is the limit of how many a random gate should be defined at maximum
 * either wires or limit must be equal 0
 ***************************/
int GA::getRandomGate(int wires, bool restr, int topwire, int limit){
	int temp;
	char rg[REPSIZE];
	float totfit = 0;
	float tfit = 0.0;
	if (wires > 0 && limit > 0 ) return -1;
	if (Ga < 1){
		if (wires == 0){
			if (!restr)
				if (limit < 1)
					return rand() % numofgates;
				else do {
					temp = rand() % numofgates;
				} while (gateArray[temp]->numIO > limit);
			else {
				if (limit < 1 )
					do{
						temp = rand() % numofgates;
					} while (!checkRestrictions(topwire, gateArray[temp]));
				else do {
					do{
						temp = rand() % numofgates;
					} while (!checkRestrictions(topwire, gateArray[temp]));
				}while (gateArray[temp]->numIO > limit);
				return temp;
			}
		}
		if (wires != 0){
			if (!restr){
				temp = rand() % stringArray[wires-1].numberGates;	
				return temp;
			} else {
				do {
					temp = rand() % stringArray[wires-1].numberGates;
					temp = getGate(stringArray[wires-1].nameGate[temp]);	
				} while (!checkRestrictions(topwire, gateArray[temp]));
				return temp;
			}
		}
	} else {
		for (int y = 0; y < numofgates; y++)
			totfit += rankV[y].fitness;
//cout<<"get "<<totfit<<", ";
		if (wires == 0){
			if (!restr){
				do {
					tfit = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))))* totfit;
					temp = getBestGateMatch(tfit, topwire, 0, limit, false);
				}while(temp == -1);
				return temp;
			} else {
				do{
					tfit = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))))* totfit;
					temp = getBestGateMatch(tfit, topwire, 0, limit, true);
				} while (temp == -1);// || !checkRestrictions(topwire, gateArray[temp]));
				return temp;
			}
		} else if (wires != 0){
			if (!restr){
				do {
					tfit = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))))* totfit;
					temp = getBestGateMatch(tfit, topwire, wires, 0, false);
				}while(temp == -1);
				return temp;
			} else {
				do {
					tfit = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))))* totfit;
					temp = getBestGateMatch(tfit, topwire, wires, 0, true);
				} while (temp == -1);// || !checkRestrictions(topwire, gateArray[temp]));
				return temp;
			}
		}
	}
	return 0;
}
/**************************
 * check the length of the circuit and
 * either cuts it to maximum allowed
 * either generates a new one if too short
 ***************************/
void GA::repairCircuit(int a) {

	int counter = 0;
	Individual *I = population[a];

	if ((int) I->my_string.length() > 0) {
		for (int A = 0; A < (int) I->my_string.length(); A++)
			if (I->my_string.at(A) == 'p')
				counter++;
		//compare cicuits to the minimum required size and twice the maximum size
		if (counter < (mincircuitsize*2)  || counter > (segments*segments_multi) ) {
			//reinitialize the circuit
			I->segmentNumber = segments;
			makeIndvString(I);
		//	generatePhases(I);
		}
	} else {
		I->segmentNumber = segments;
		makeIndvString(I);
		//generatePhases(I);
	}
}
/************************************************
 *Restrictions:
 *1. replace gate with same number of I/O
 *2. replace one gate with n I/O with any combination of smaller gates with
 *same number of I/O
 *3. conserve strict correspondence of connections
 *void operator=(const class bla & other) { mymember=other.mymember; }
 *************************************************/
void GA::apply1CrossoverP(int ind1, int ind2) {
	//	apply rules and make crossover
	int proba1, proba2, segments_a, segments_b, o_e;
	cuComplex temp[MAXNUMOFGATES];
	Individual *I = new_population[ind1];
	Individual *J = new_population[ind2];
	basic_string<char> temp1;
	basic_string<char> temp2;
	int index1, index2, iterations1, iterations2;

	//individual 1
	do {
		proba1 = (rand() % (I->my_string.size()));
	} while (proba1 >= I->my_string.length() - 1 || I->my_string.at(proba1)
			!= 'p' || proba1 <= 0);

	//cout<<proba1;

	if ((I->my_string.at(proba1)) == 'p')
		if ((I->my_string.at(proba1 + 1)) == 'p')
			proba1 += 1;

	temp1 = I->my_string.substr(proba1);

	segments_a = 0;
	o_e = 0;
	for (int k = 0; k < I->my_string.size(); k++){
		if (k > proba1)
			break;
		if (I->my_string.at(k) == 'p'){
			if (o_e / 2 == 0 && o_e > 0)
				segments_a++;
			o_e++;
		}
	}
	//individual 2
	for (int k = 0; k < J->my_string.size(); k++){
		if (I->my_string.at(k) == 'p'){
			if (o_e / 2 == 0 && o_e > 0)
				segments_b++;
			o_e++;
		}
		if (segments_b == segments_a){
			break;
			proba2 = k;
		}
	}

	if ((J->my_string.at(proba2)) == 'p')
		if ((J->my_string.at(proba2 + 1)) == 'p')
			proba2 += 1;

	temp2 = J->my_string.substr(proba2);

	iterations1 = getGatesfromString(I->my_string);
	iterations2 = getGatesfromString(J->my_string);

	//get the phase size
	//we need both indexes in the array of phases
	index1 = 0;
	index2 = 1;
	for (int a = 0; a < I->my_string.length(); a++)
		if (a < proba1) {
			if (I->my_string[a] != 'p') {
				index1 += 1;
			}
		} else {
			if (phase > 0)
			if (I->my_string[a] != 'p') {
				temp[a - proba1] = make_cuFloatComplex(cuCrealf(I->phases[a]),cuCimagf(I->phases[a]));
			}
		}

	for (int a = 0; a < proba2; a++)
		if (J->my_string[a] != 'p')
			index2 += 1;

	//	swap the phase sequences
	if (phase > 0){
		for (int a = index2; a < (iterations2 - index2); a++)
			I->phases[index1 + a - index2] = J->phases[a];

		for (int a = index1; a < (iterations1 - index1); a++)
			J->phases[index2 + a - index1] = make_cuFloatComplex(cuCrealf(temp[a- index1]), cuCimagf(temp[a - index1]));
	}
	if (temp1.length() >= 3 && temp2.length() >= 3) {
		I->my_string.replace(proba1, temp1.length(), temp2, 0, temp2.length());
		J->my_string.replace(proba2, temp2.length(), temp1, 0, temp1.length());
	}

}

/************************************************
 *Restrictions:
 *1. replace gate with same number of I/O
 *2. replace one gate with n I/O with any combination of smaller gates with
 *same number of I/O
 *3. conserve strict correspondence of connections
 *void operator=(const class bla & other) { mymember=other.mymember; }
 *************************************************/
void GA::apply1Crossover(int ind1, int ind2) {
	//	apply rules and make crossover
	int proba1, proba2;
	cuComplex temp[MAXNUMOFGATES];
	Individual *I = new_population[ind1];
	Individual *J = new_population[ind2];
	basic_string<char> temp1;
	basic_string<char> temp2;
	int index1, index2, iterations1, iterations2;

	//individual 1
	do {
		proba1 = (rand() % (I->my_string.size()));
	} while (proba1 >= I->my_string.length() - 1 || I->my_string.at(proba1)
			!= 'p');

	//cout<<proba1;

	if (proba1 <= 0)
		temp1 = I->my_string.substr(proba1);
	else {
		if ((I->my_string.at(proba1)) == 'p')
			if ((I->my_string.at(proba1 + 1)) == 'p')
				proba1 += 1;

		temp1 = I->my_string.substr(proba1);
	}

	//individual 2

	do {
		proba2 = (rand() % (J->my_string.size()));
	} while (proba2 >= J->my_string.length() - 1 || J->my_string.at(proba2)
			!= 'p');

	//cout<<" "<<proba2<<endl;
	if (proba2 <= 0)
		temp2 = J->my_string.substr(proba2);
	else {
		if ((J->my_string.at(proba2)) == 'p')
			if ((J->my_string.at(proba2 + 1)) == 'p')
				proba2 += 1;

		temp2 = J->my_string.substr(proba2);
	}

	iterations1 = getGatesfromString(I->my_string);
	iterations2 = getGatesfromString(J->my_string);

	//get the phase size
	//we need both indexes in the array of phases
	index1 = 0;
	index2 = 1;
	for (int a = 0; a < I->my_string.length(); a++)
		if (a < proba1) {
			if (I->my_string[a] != 'p') {
				index1 += 1;
			}
		} else {
			if (phase > 0)
			if (I->my_string[a] != 'p') {
				temp[a - proba1] = make_cuFloatComplex(cuCrealf(I->phases[a]),cuCimagf(I->phases[a]));
			}
		}

	for (int a = 0; a < proba2; a++)
		if (J->my_string[a] != 'p')
			index2 += 1;

	//	swap the phase sequences
	if (phase > 0) {
		for (int a = index2; a < (iterations2 - index2); a++)
			I->phases[index1 + a - index2] = J->phases[a];

		for (int a = index1; a < (iterations1 - index1); a++)
			J->phases[index2 + a - index1] = make_cuFloatComplex(cuCrealf(temp[a - index1]), cuCimagf(temp[a - index1]));
	}
	if (temp1.length() >= 3 && temp2.length() >= 3) {
		I->my_string.replace(proba1, temp1.length(), temp2, 0, temp2.length());
		J->my_string.replace(proba2, temp2.length(), temp1, 0, temp1.length());
	}

}
/*******************************************
 * two point crossover but no phases swap yet
 ********************************************/
void GA::apply2CrossoverP(int ind1, int ind2) {
	//apply rules and make crossover
	//first it is used to store segment number
	int proba[4];
	basic_string<char> temp1;
	basic_string<char> temp2;
	Individual *I = new_population[ind1];
	Individual *J = new_population[ind2];
	int a, b, c, d, n, count;
	int pos1, pos2;

	if (J->my_string.length() < I->my_string.length()){
		J = new_population[ind1];
		I = new_population[ind2];
	}
//	cout<<I->my_string<<endl;
//	cout<<J->my_string<<endl;
	if (I->ioNumber == J->ioNumber) {
		 //get string cuting points
		do {
			//the number of segments to swap on individual 1
			proba[0] = (rand() % (I->segmentNumber / 2));
		} while (proba[0] <= 0 && proba[0] > J->segmentNumber / 2);
		a = proba[0];

		do {
			//the location of the cut
			proba[1] = (rand() % (I->segmentNumber));
		} while ((proba[1] + proba[0]) >= I->segmentNumber-1 || (proba[1] + proba[0]) >= J->segmentNumber-1);
		b = proba[1];

		n = 0;
		count = 0;
		temp1 = "";
		while (count < (b + a)) {

			if (I->my_string.at(n) == 'p')
				count++;
			//the string starts between the indexes of the segments
			if (count > b) {
				temp1 += I->my_string.at(n);
			}
			n++;
		}
		pos1 = I->my_string.find(temp1, 0);

		n = 0;
		count = 0;
		temp2 = "";
		while (count < (a + b)) {

			if (J->my_string.at(n) == 'p')
				count++;
			//the string starts between the indexes of the segments
			if (count > b) {
				temp2 += J->my_string.at(n);
			}
			n++;
		}
		pos2 = J->my_string.find(temp2, 0);
		I->my_string.replace(pos1, temp1.length(), temp2, 0, temp2.length());
		J->my_string.replace(pos2, temp2.length(), temp1, 0, temp1.length());
	}
}
/*******************************************
 * two point crossover but no phases swap yet
 ********************************************/
void GA::apply2Crossover(int ind1, int ind2) {
	//apply rules and make crossover
	//first it is used to store segment number
	int proba[4];
	basic_string<char> temp1;
	basic_string<char> temp2;
	Individual *I = new_population[ind1];
	Individual *J = new_population[ind2];
	int a, b, c, d, n, count;
	int pos1, pos2;

	if (I->ioNumber == J->ioNumber) {
		 //get string cuting points
		do {
			//the number of segments to swap on individual 1
			proba[0] = (rand() % (I->segmentNumber / 2));
		} while (proba[0] <= 0);
		a = proba[0];

		do {
			//the location of the cut
			proba[1] = (rand() % (I->segmentNumber));
		} while ((proba[1] + proba[0]) >= I->segmentNumber-1);
		b = proba[1];

		do {
			//the number of segments to swap on individual 2
			proba[2] = (rand() % (J->segmentNumber / 2));
		} while (proba[2] <= 0);
		c = proba[2];

		do {
			//the location of the cut
			proba[3] = (rand() % (J->segmentNumber));
		} while ((proba[3] + proba[2]) >= J->segmentNumber-1);
		d = proba[3];
		n = 0;
		count = 0;
		temp1 = "";
		while (count < (b + a)) {

			if (I->my_string.at(n) == 'p')
				count++;
			//the string starts between the indexes of the segments
			if (count > b) {
				temp1 += I->my_string.at(n);
			}
			n++;
		}
		pos1 = I->my_string.find(temp1, 0);

		n = 0;
		count = 0;
		temp2 = "";
		while (count < (c + d)) {

			if (J->my_string.at(n) == 'p')
				count++;
			//the string starts between the indexes of the segments
			if (count > d) {
				temp2 += J->my_string.at(n);
			}
			n++;
		}
		pos2 = J->my_string.find(temp2, 0);

		I->my_string.replace(pos1, temp1.length(), temp2, 0, temp2.length());
		J->my_string.replace(pos2, temp2.length(), temp1, 0, temp1.length());
		

	}
}
/************************************************
 * replicate the population according to the settings
 ************************************************/
void GA::applyReplication() {
	if (pareto > 0)
		applyParetoReplication();
	else
		applySingleObjReplication();
}

void GA::applySingleObjReplication() {
	float fitnessSum = 0;
	float fitnessavg = 0;
	float first, second, tempsum, proba, best;
	int indv1, indv2;

	for (int i = 0; i < populationNumber; i++)
		fitnessSum += population[i]->fitness;

	fitnessavg = (float)(fitnessSum)/(float)(populationNumber);
	cout<<fitnessavg<<endl;
	int counter = 0;
	bool condition = false;
	while (!condition) {
		if (replicator == 0) {
			//roulette wheel
			//select only individuals with fitness > threshold
			if (threshold > 0)
				do {
					indv1 = rand() % (populationNumber - 1);
					indv2 = rand() % (populationNumber - 1);
				} while (population[indv1]->fitness < threshold
						&& population[indv2]->fitness < threshold);
			else {
				first = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))))
				* fitnessSum;
				second = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))))
				* fitnessSum;
			}
		} else {
			//stochastic universal sampling
			//select only individuals with fitness > threshold
			if (threshold > 0)
				do {
					indv1 = rand() % (populationNumber - 1);
					indv2 = rand() % (populationNumber - 1);
				} while (population[indv1]->fitness < threshold
						&& population[indv2]->fitness < threshold);
			else {
				first = ((1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))))
						* fitnessSum) / 3;
				second = fitnessSum - first;
			}
		}

		if (threshold == 0) {
			tempsum = 0;
			for (int c = 0; c < populationNumber; c++) {
				tempsum += population[c]->fitness;
				if (tempsum >= first) {
					indv1 = c;
					break;
				}
			}

			tempsum = 0;
			for (int d = 0; d < populationNumber; d++) {
				tempsum += population[d]->fitness;
				if (tempsum >= second) {
					indv2 = d;
					break;
				}
			}
		}

		if (counter == populationNumber - 1) {
			setIndv(new_population[counter], population[indv1]);
			counter++;
		} else {
			setIndv(new_population[counter], population[indv1]);
			setIndv(new_population[counter + 1], population[indv2]);

			proba = 1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0)));
			if (proba < proba_Crossover && indv1 != indv2) {
				applyCrossover(counter, counter+1);
			}

			counter += 2;
		}
		if (counter >= populationNumber - 1)
			condition = true;

	}

	if (elitism > 0){

		best = population[0]->fitness;
		indv1 = 0;
		for (int t = 1; t < populationNumber; t++) {
			if (population[t]->fitness > best) {
				best = population[t]->fitness;
				indv1 = t;
			}
		}
		setIndv(population[0], population[indv1]);
		for (int i = 0; i < populationNumber-1; i++) {
			setIndv(population[i+1], new_population[i]);
		}

	} else {
		for (int i = 0; i < populationNumber; i++) {
			setIndv(population[i], new_population[i]);
		
		}
	}

}
/*************************************************
 * Pareto replication with initial settings
 *************************************************/
void GA::applyParetoReplication() {
	float proba;
	int indv1, indv2, tempsum, fitnessSum, first, second;
	int counter = 0;
	bool condition = false;

	for (int i = 0; i < populationNumber; i++)
		population[i]->Rank = -1;
	int temp_rank;

	for (int i = 0; i < populationNumber; i++) {
		temp_rank = 0;
		for (int j = 0; j < populationNumber; j++) {
			if (population[i]->Rank <= -1) {
				if ((alpha1 * population[i]->Error) <= population[j]->Error
						&& (beta1 * population[i]->Cost) <= population[j]->Cost)
					temp_rank++;
			}
		}
		population[i]->Rank = temp_rank;
	}

	fitnessSum = 0;
	for (int i = 0; i < populationNumber; i++) {
		fitnessSum += population[i]->Rank;
	}

	while (!condition) {
		if (replicator == 0) {
			//roulette wheel
			if (threshold > 0)
				do {
					indv1 = rand() % (populationNumber - 1);
					indv2 = rand() % (populationNumber - 1);
				} while (population[indv1]->fitness < threshold
						&& population[indv2]->fitness < threshold);
			else {
				first = (int) (1 - (float(1.0 * rand() / ((float) RAND_MAX
						+ 1.0)))) * fitnessSum;
				second = (int) (1 - (float(1.0 * rand() / ((float) RAND_MAX
						+ 1.0)))) * fitnessSum;
			}
		} else {
			//stochastic uniersal sampling
			if (threshold > 0)
				do {
					indv1 = rand() % (populationNumber - 1);
					indv2 = rand() % (populationNumber - 1);
				} while (population[indv1]->fitness < threshold
						&& population[indv2]->fitness < threshold);
			else {
				first = (int) ((1 - (float(1.0 * rand() / ((float) RAND_MAX
						+ 1.0)))) * fitnessSum) / 2;
				second = fitnessSum - first;
			}
		}

		if (threshold == 0) {
			tempsum = 0;
			for (int c = 0; c < populationNumber; c++) {
				tempsum += population[c]->Rank;
				if (tempsum >= first) {
					indv1 = c;
					break;
				}
			}

			tempsum = 0;
			for (int d = 0; d < populationNumber; d++) {
				tempsum += population[d]->Rank;
				if (tempsum >= second) {
					indv2 = d;
					break;
				}
			}
		}

		if (counter == populationNumber - 1) {
			setIndv(new_population[counter], population[indv1]);
			counter++;
		} else {

			setIndv(new_population[counter], population[indv1]);
			setIndv(new_population[counter + 1], population[indv2]);

			proba = 1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0)));
			if (proba < proba_Crossover && indv1 != indv2)
				applyCrossover(counter, counter+1);
			counter += 2;
		}
		if (counter >= populationNumber - 1)
			condition = true;
	}

	if (elitism > 0){

		best = population[0]->fitness;
		first = 0;
		for (int t = 1; t < populationNumber; t++) {
			if (population[t]->fitness > best) {
				best = population[t]->fitness;
				first = t;
			}
		}
		setIndv(population[0], population[first]);
		for (int i = 0; i < populationNumber-1; i++) {
			setIndv(population[i+1], new_population[i]);
		}

	} else {
		for (int i = 0; i < populationNumber; i++) {
			setIndv(population[i], new_population[i]);
		}
	}

}
/*************************************************
 * copy one individual to another one
 *************************************************/
void GA::setIndv(Individual *indv1, Individual *indv2) {
	indv1->fitness = indv2->fitness;
	indv1->ioNumber = (int) indv2->ioNumber;
	indv1->segmentNumber = (int) indv2->segmentNumber;
	indv1->Cost = (int) indv2->Cost;
	indv1->Error = (float) indv2->Error;
	indv1->Groups = (int) indv2->Groups; // for share fitness calcualtion
	indv1->Rank = (int) indv2->Rank;
	indv1->valuedness = indv2->valuedness;
	//indv1->my_string.erase();
	indv1->my_string = string("");
	indv1->rev_lib_string = string("");
	if (phase > 0){
		indv1->phase = indv2->phase;
		for (int a = 0; a < getGatesfromString(indv2->my_string); a++)
			indv1->phases[a] = indv2->phases[a];
	}
	for (unsigned int a = 0; a < indv2->my_string.size(); a++)
		indv1->my_string += indv2->my_string[a];
	for (unsigned int a = 0; a < indv2->rev_lib_string.size(); a++)
		indv1->rev_lib_string += indv2->rev_lib_string[a];
}
/**************************
 * this works only if the circuits are normalized
 * we generate the Individual from gates and replace them s the
 * n last one in the population.
 * the gates are translated to normal circuit using the
 * permutativeArray parents
 ***************************/
void GA::injectNormalCircuits(string **refs) {
/*	int counter = 0;
	int count = 5;
	int tcounter = 0;

	if (refs == NULL)
		return;
	cout << "injecting circuits back to GA" << endl;
	while (true) {
		if (counter >= count - 1)
			break;
		else if (refs[counter] == NULL || refs[counter]->size() <= 2)
			break;
		counter++;
	}
	for (int a = 0; a < counter; a++)
		cout << "* " << refs[a] << endl;

	for (int a = counter - 1; a >= 0; a--) {
		cout << "* " << a << endl;
		cout << "* " << refs[a] << endl;
		delete (population[populationNumber - a]);
		population[populationNumber - a] = new Individual();
		population[populationNumber - a]->my_string = string("");
		cout << "* " << refs[a] << "  " << refs[a]->size() << endl;
		for (unsigned int b = 0; b < refs[a]->size(); b++) {
			tcounter = 0;
			//fetch the tranlsated string from tbe permutative array
			while (permutativeArray[tcounter]->representation != refs[a]->at(b))
				tcounter++;
			cout << "* " << permutativeArray[tcounter]->my_string << endl;
			population[populationNumber - a]->my_string += 'p'
				+ permutativeArray[tcounter]->my_string + 'p';
		}
		//trivial phase
		population[populationNumber - a]->phase = make_cuFloatComplex(0, 1);
		population[populationNumber - a]->ioNumber = finalGate->numIO;
		cout << "<*> " << population[populationNumber - a]->my_string << endl;
		if (measurement > 0)
			doMeasureFitness(population[populationNumber - a], true);
		else
			doMatrixFitness(population[populationNumber - a], true);
		cout << "<> " << population[populationNumber - a]->my_string << endl;
	}
*/
}
/*******************************************************
 *this works only if the circuits are normalized
 *we generate the Individual from gates and replace them s the
 *n last one in the population.
 *the gates aretranslated to normal circuit using the
 *permutativeArray parents
 *******************************************************/
void GA::injectNormalSegments(string **refs) {
/*	int counter = 0;
	int count = 5;
	int tcounter = 0;

	if (refs == NULL)
		return;
	cout << "injecting circuits back to GA" << endl;
	//check if there are valid solutions
	while (true) {
		if (counter >= count - 1)
			break;
		else if (refs[counter] == NULL || refs[counter]->size() <= 2)
			break;
		counter++;
	}
	//print obtained circuits
	for (int a = 0; a < counter; a++)
		cout << "* " << refs[a] << endl;

	//add all solutions into new indivduals and replace them in the GA
	for (int a = counter - 1; a >= 0; a--) {
		cout << "* " << a << endl;
		cout << "* " << refs[a] << endl;
		//erase the first n that are the smallest ones.
		delete (population[populationNumber - a]);
		population[populationNumber - a] = new Individual();
		population[populationNumber - a]->my_string = string("");
		cout << "* " << refs[a] << "  " << refs[a]->size() << endl;
		for (unsigned int b = 0; b < refs[a]->size(); b++) {
			tcounter = 0;
			//fetch the tranlsated string from tbe permutative array
			while (permutativeArray[tcounter]->representation != refs[a]->at(b))
				tcounter++;
			cout << "* " << permutativeArray[tcounter]->my_string << endl;
			//generate teh circuit string in the individual
			population[populationNumber - a]->my_string += 'p'
				+ permutativeArray[tcounter]->my_string + 'p';
		}
		population[populationNumber - a]->ioNumber = finalGate->numIO;
		cout << "<*> " << population[populationNumber - a]->my_string << endl;
		//evaluate the added circuit
		if (measurement > 0)
			doMeasureFitness(population[populationNumber - a], true);
		else
			doMatrixFitness(population[populationNumber - a], true);
		cout << "<> " << population[populationNumber - a]->my_string << endl;
	}
*/
}
/*************************************************
 * Returns the aray of gates generated as single gate per parallel segments
 * All gates have thus size CircuitSize and can be simply swapped
 *************************************************/
qGate** GA::getPermutativeCircuits() {
/*
	qGate *workgate, *multigate, *multigate2, *tempgate;
	string names = "";
	int gatecounter = 0;
	int c = 0;
	int rep = 0;
	int columns;
	multigate = new qGate;
	initGate(multigate, resultnum, gateArray[0]->valuedness, 1);
	tempgate = new qGate;
	initGate(tempgate, resultnum, gateArray[0]->valuedness, 1);
	string s = bestIndv->my_string;
	for (int m = 0; m < MAXNUMOFGATES; m++)
		permutativeArray[m] = NULL;
	//cout<<"starting generating: "<<s<<"  "<<endl;
	for (unsigned int a = 0; a < s.length(); a++) {
		//cout<<"starting generating "<<s.at(a)<<"  "<<endl;
		rep = 0;
		//we found a real gate
		if (s.at(a) != 'p' && s.at(a) != 'A') {
			//cout<<"starting generating 0"<<endl;
			for (int b = 0; b < numofgates; b++) {
				//remove the separators and the wire gates
				if (gateArray[b]->representation == s.at(a)) {
					//check wether it was already added
					for (unsigned int m = 0; m < names.size(); m++)
						if (names.at(m) == s.at(a))
							rep = 100;
					workgate = gateArray[b];
				}
			}
			//cout<<"starting generating 1"<<endl;
			//if yes do not add this gate again
			if (rep == 0) {
				//add the gate
				names += workgate->representation;
				//here comes the transformation for non normal gate to normal ones
				//if the given gate is smaller than the final one
				//cout<<"starting generating 2"<<endl;
				if (finalGate->numIO > workgate->numIO) {
					//generate as many new gates -- normal ones
					//cout<<"starting generating 3"<<endl;
					for (int t = 0; t <= finalGate->numIO - workgate->numIO; t++) {
						if (permutativeArray[gatecounter] != NULL)
							delete (permutativeArray[gatecounter]);
						permutativeArray[gatecounter] = new qGate;
						initGate(permutativeArray[gatecounter], resultnum, finalGate->valuedness, 1);
						permutativeArray[gatecounter]->my_string = string("");
						//how many W's we have to add to create a shifted one
						for (int j = 0; j <= (finalGate->numIO
								- workgate->numIO); j++) {
							if (j == t)
								permutativeArray[gatecounter]->my_string
								+= workgate->representation;
							else
								permutativeArray[gatecounter]->my_string += 'A';
						}
						permutativeArray[gatecounter]->representation
						= char('A' + gatecounter);
						if (t == 0)
							multigate = ccGate(workgate, multigate);
						else
							multigate = ccGate(gateArray[0], multigate);
						//calculate the matrix of the segment
						c = 0;
						while (true) {
							if (c != t)
								multigate2 = gateArray[0];
							else
								multigate2 = workgate;
							//using Kronecker multiplication
							tensorProduct3(multigate, multigate2, tempgate);
							multigate = ccGate(tempgate, multigate);
							c++;
							if (c >= (finalGate->numIO - workgate->numIO))
								break;
						}
						//cout<<"starting generating 4"<<endl;
						//fill the matrix of the gate
						columns = (int(pow(float(2), finalGate->numIO)));
						for (int e = 0; e < columns; e++)
							for (int f = 0; f < columns; f++)
								permutativeArray[gatecounter]->gateMatrix1[e+f*columns]
								                                              = multigate->gateMatrix1[e+f*columns];
						permutativeArray[gatecounter]->numIO = finalGate->numIO;
						gatecounter++;
					}
				} else {
					if (permutativeArray[gatecounter] != NULL)
						delete (permutativeArray[gatecounter]);
					permutativeArray[gatecounter] = NULL;
					permutativeArray[gatecounter] = new qGate;
					initGate(permutativeArray[gatecounter], finalGate->numIO, finalGate->valuedness, 1);
					permutativeArray[gatecounter]->representation = char('A'
							+ gatecounter);
					columns = (int(pow(float(2), finalGate->numIO)));
					for (int e = 0; e < columns; e++)
						for (int f = 0; f < columns; f++)
							permutativeArray[gatecounter]->gateMatrix1[e+f*columns]
							                                              = workgate->gateMatrix1[e+f*columns];
					gatecounter++;
				}
			}
		}
	}
	//cout<<"starting generating 10"<<endl;
	numofnormalgates = gatecounter;

	for (int n = 0; n < numofnormalgates; n++)
		cout << permutativeArray[n]->representation << " : "
		<< permutativeArray[n]->my_string << endl;
	cout << endl;
	destroyGate(multigate);
	destroyGate(tempgate);
	delete multigate;
	delete tempgate;
	return permutativeArray;
*/
}
/*************************************************
 * Returns a set of gates representing single segments of the circuit
 * this can be used if circuits with big size are generated
 *************************************************/
qGate** GA::getPermutativeSegments() {
/*
	qGate *workgate, *multigate, *multigate2, *tempgate;
	string names = "";
	int gatecounter = 0;
	int columns;
	string s = bestIndv->my_string;
	bool processing = false;
	multigate = new qGate;
	initGate(multigate, resultnum, finalGate->valuedness, 1);
	tempgate = new qGate;
	initGate(tempgate, resultnum, finalGate->valuedness, 1);
	for (int m = 0; m < MAXNUMOFGATES; m++) {
		permutativeArray[m] = NULL;
	}
	cout << "starting generating: " << s << "  " << endl;
	for (unsigned int a = 0; a < s.length(); a++) {
		//cout<<"starting generating "<<s.at(a)<<"  "<<endl;
		//we found a real gate
		if (s.at(a) == 'p') {
			if (a + 1 < s.length()) {
				if (s.at(a + 1) == 'p') {
					cout << "Name generated : " << names << endl;
					processing = false;
					gatecounter++;
				} else {
					processing = true;
					names = string("");
				}
			} else {
				cout << "Name generated : " << names << endl;
				processing = false;
			}
		} else {
			workgate = gateArray[getGate(s.at(a))];
			names += workgate->representation;
			//do the multiplication when the last gate is registered
			if (a + 1 < s.length() && s.at(a + 1) == 'p') {
				cout << "starting generating 0" << endl;
				//do the final gate multiplication
				for (int n = 0; n < names.length(); n++) {
					if (n == 0) {//set the 0'th gate
						workgate = gateArray[getGate(names.at(n))];
						if (permutativeArray[gatecounter] != NULL)
							delete (permutativeArray[gatecounter]);
						permutativeArray[gatecounter] = new qGate();
						permutativeArray[gatecounter]->representation
						= char('A' + gatecounter);
						multigate = ccGate(workgate, multigate);
						permutativeArray[gatecounter]->my_string
						= string(names);
						initGate(permutativeArray[gatecounter], finalGate->numIO, finalGate->valuedness, 1);
						cout << "starting generating 10: "
						<< permutativeArray[gatecounter]->my_string
						<< endl;
					} else {//add the next gate using the tensor product
						multigate2 = gateArray[getGate(names.at(n))];
						initMatrix(tempgate, multigate2->numIO);
						tensorProduct3(multigate, multigate2, tempgate);
						multigate = ccGate(tempgate, multigate);
					}
				}
				columns = (int(pow(float(2), finalGate->numIO)));
				for (int e = 0; e < columns; e++)
					for (int f = 0; f < columns; f++)
						permutativeArray[gatecounter]->gateMatrix1[e+f*columns]
						                                              = multigate->gateMatrix1[e+f*columns];

			}
		}
	}
	gatecounter++; //addlast gate
	//cout<<"starting generating 10"<<endl;
	numofnormalgates = gatecounter;

	for (int n = 0; n < numofnormalgates; n++) {
		cout << permutativeArray[n]->representation << " : "
		<< permutativeArray[n]->my_string << endl;
	}
	cout << bestIndv->my_string << endl;
	cout << endl;
	destroyGate(multigate);
	destroyGate(tempgate);
	delete multigate;
	delete tempgate;

	return permutativeArray;
*/
}
/*************************************************
 *
 *************************************************/
bool GA::terminatingCondition() {
	return condition;
}
/*************************************************
 * the way of stopping things
 *************************************************/
void GA::setTerminatingCondition(int counter) {
	generationCondition = counter;
	if (counter >= MAXGEN || counter >= generations) {
		out_stream << " generation max reached " << endl << endl << endl;
		condition = true;
	} else if (getFitness()) {
		out_stream << " Solution found " << endl << endl << endl;
		condition = true;
	} else
		condition = false;
}
/*************************************************
 * check wether or not the correct solution was found
 * and if not it refresh the bestIndividual to the newest best value
 *************************************************/
bool GA::getFitness() {
	float sigma = 0.0005;
	if ((bestIndv->fitness + sigma) >= 1)
		return true;
	else
		return false;
}

/*******************************************
 * closes the file streams
 *******************************************/
void GA::closeStream() {
	in_stream.close();
	out_stream.close();
}
/*******************************************
 * find if a gate given by A is already stored in gateArray
 *******************************************/
int GA::findGate(qGate *A) {
	int a = 0;
	bool matrixTest;
	int result = -1;
	int columns;
	float val = (float) A->valuedness;
	while (a < numofgates) {
		out_stream << a << " " << numofgates << endl;
		matrixTest = true;
		if (A->numIO != gateArray[a]->numIO)
			result = -1;
		else {
			columns  = (int(pow(val, (float) A->numIO)));
			for (int i = 0; i < columns; i++)
				for (int j = 0; j < columns; j++)
					if (cuCrealf(A->gateMatrix1[i+j*columns]) != cuCrealf(gateArray[a]->gateMatrix1[i+j*columns]) || cuCimagf(A->gateMatrix1[i+j*columns]) != cuCimagf(gateArray[a]->gateMatrix1[i+j*columns])) matrixTest = false;
			if (matrixTest)
				result = a;
			else
				result = -1;
		}
		if (result != -1)
			break;
		else
			a++;
	}
	return result;
}
/*******************************************
 * parse the string circuit into an array where
 * each segment is a sub array with gates being 
 * represented by int
 *******************************************/
void GA::parseCircuitIntoArray(string circuit, int wires, int **intgates){
	int begin, end, wcounter, gcounter, diff, charcounter, wchounter, indwcounter;
	bool condition;
	bool bwires;
	char g[REPSIZE+1];
	int segcount = 0;
	int gate_n;
	begin = 0;
	string s1;
	end = circuit.find("p", begin + 1);
	do {
		diff = end-begin;
		s1 = circuit.substr(begin+1, diff-1);
		wcounter =  0;
		charcounter = 0;
		condition = false;
		bwires = false;
		wchounter = 0;
		gcounter = 0;
		while (wcounter < wires){
			bwires = true;
			getGateRepfromStr(wchounter, s1, g);
			gate_n = getGate(g);
/*			if (gate_n != 0){
				bwires = false;
				break;
			}
*/
//cout<<"parsing circuit: "<<g<<":  "<<gate_n<<" ; ";

			intgates[segcount][gcounter++] = gate_n;
			wcounter += (gateArray[gate_n]->numIO);
			//wcounter += (gateArray[gate_n]->connections[gateArray[gate_n]->realIO-1]+1);
			wchounter +=3;
		}
		for (int j = gcounter; j < wires; j++)
			intgates[segcount][j] = -1;
		segcount++;

/*		wcounter =  0;
		gcounter = 0;
		if (!bwires){
			while (wcounter < wires){
				getGateRepfromStr(charcounter, s1, g);
				gate_n = getGate(g);
				intgates[segcount][gcounter++] = gate_n;
							wcounter += (gateArray[gate_n]->connections[gateArray[gate_n]->realIO-1]+1);
				charcounter +=3;
			}
			for (int j = gcounter; j < wires; j++)
				intgates[segcount][j] = -1;
			segcount++;
		}
*/		if (end < circuit.length()-1){
			begin = end+1;
			end = circuit.find("p", begin + 1);
		} else break;
		if(end < 0) break;
	}while(true);
	intgates[segcount++][0] = -1;
}
/*******************************************
 * Compares two segments represented as arrays of ints
 * each int represent the gate in the gateArray structure
 * It returns two strings - representing two segments
 * after minimization
 *******************************************/
int GA::compareIntSegments(int *intgates1, int *intgates2, char **sseg1, char **sseg2, int **wires1, int **wires2, int wire){
	bool combined = false;
	bool checked = false;	
	int g_counter1 = 0;
	int w_counter1 = 0;
	int g_counter2 = 0;
	int w_counter2 = 0;
	int rgcount1 = 0;
	int wcount1 = 0;
	int rgcount2 = 0;
	int wcount2 = 0;
	int result = 0;
	while(w_counter1 < wire  || w_counter2 < wire){
		while(intgates2[g_counter2] != -1){
			//if the first segment has more gates than the second
			checked = true;
			//if the gates have same top wire
			if(w_counter1 == w_counter2){
				//it the gates have same number of wires merge them
				if(gateArray[intgates1[g_counter1]]->numIO == gateArray[intgates2[g_counter2]]->numIO){
					wires2[w_counter2][0] = 1;
					wires2[w_counter2][1] = w_counter2;
					wires2[w_counter2][2] = w_counter2+gateArray[intgates2[g_counter2]]->numIO;
					wires1[w_counter1][0] = 1;
					wires1[w_counter1][1] = w_counter1;
					wires1[w_counter1][2] = w_counter1+gateArray[intgates1[g_counter1]]->numIO;

					result++;
					for (int j = 0; j < REPSIZE; j++)
						sseg2[wcount2][j] = gateArray[intgates1[g_counter1]]->representation[j];
					sseg2[wcount2][REPSIZE] = ' ';		
					for (int j = 0; j < REPSIZE; j++)
						sseg2[wcount2][j+REPSIZE+1] = gateArray[intgates2[g_counter2]]->representation[j];
					sseg2[wcount2][2*REPSIZE+1] = '\0';
					for (int j = 0; j < gateArray[intgates2[g_counter2]]->numIO;j++){
						sseg1[wcount1][0] = '!';
						sseg1[wcount1][1] = '!';
						sseg1[wcount1][2] = '!';
						sseg1[wcount1++][3] = '\0';
						rgcount1++;
					}
					rgcount2++;
					wcount2 += gateArray[intgates2[g_counter2]]->numIO;
					w_counter1 += gateArray[intgates2[g_counter2]]->numIO;
					w_counter2 += gateArray[intgates2[g_counter2]]->numIO;
					combined = true;
					g_counter1++;
					g_counter2++;
				}  else {
					for (int j = 0; j < REPSIZE; j++)
						sseg1[wcount1][j] = gateArray[intgates1[g_counter1]]->representation[j];
					sseg1[wcount1][REPSIZE] = '\0';

					for (int j = 0; j < REPSIZE; j++)
						sseg2[wcount2][j] = gateArray[intgates2[g_counter2]]->representation[j];
					sseg2[wcount2][REPSIZE] = '\0';
					rgcount1++;
					rgcount2++;
					wcount2 += gateArray[intgates2[g_counter2]]->numIO;
					wcount1 += gateArray[intgates1[g_counter1]]->numIO;
					w_counter1 += gateArray[intgates1[g_counter1]]->numIO;
					w_counter2 += gateArray[intgates2[g_counter2]]->numIO;
					g_counter1++;
					g_counter2++;
				}
			} else break;
			if (g_counter2 >= wire) break;
		}
		//one of the gates was of different size
		//add to both segments as many gates as required to end up with same number of wires
		if (w_counter1 > w_counter2){
			for (int j = 0; j < REPSIZE; j++)
				sseg2[wcount2][j] = gateArray[intgates2[g_counter2]]->representation[j];
			sseg2[wcount2][REPSIZE] = '\0';
			rgcount2++;
			wcount2 += gateArray[intgates2[g_counter2]]->numIO;
			w_counter2 += gateArray[intgates2[g_counter2]]->numIO;
			g_counter2++;
		} else if(w_counter2 > w_counter1){
			for (int j = 0; j < REPSIZE; j++)
				sseg1[wcount1][j] = gateArray[intgates1[g_counter1]]->representation[j];
			sseg1[wcount1][REPSIZE] = '\0';
			rgcount1++;
			wcount1 += gateArray[intgates1[g_counter1]]->numIO;
			w_counter1 += gateArray[intgates1[g_counter1]]->numIO;
			g_counter1++;

		}
		combined = false;
	}
	sseg1[wcount1][0] ='\0';
	sseg2[wcount2][0] ='\0';
//	cout<<"segment comparison done"<<endl;
	return result;
}
/*******************************************
 * used for some type of synthesis process
 * in particular using the GY, RY gates
 * requires that the forst two gates are different than the rest
 * and only the first two gates
 ********************************************/
void GA::restructureCircuitD(Individual *indi) {
	string sresult;
	int result = 0;
	int cost = 0;
	int g_counter, gcount, scount;
	indi->segmentNumber=getSegmentNumber(indi->my_string);
	int *sstring[2*indi->segmentNumber];
	int *lastss[2*indi->segmentNumber];
	int *resultSS[indi->ioNumber+1];


	for (int g = 0; g < (2*indi->segmentNumber); g++){
		sstring[g] = new int[indi->ioNumber+1];
		lastss[g] = new int[indi->ioNumber+1];
	}
	for (int g = 0; g < indi->ioNumber; g++){
		resultSS[g] = new int[2*indi->segmentNumber];
	}
	scount = 0;
	sresult = "";
	g_counter = 0;
	while (g_counter < indi->ioNumber) lastss[g_counter++][0] = '\0';
	parseCircuitIntoArray(indi->my_string, indi->ioNumber, sstring);
/*	cout<<endl;
	for (int h =0; h < indi->ioNumber; h++){
		for (int g =0; g < indi->segmentNumber; g++){
			cout<<sstring[g][h]<<" ";
		}cout<<endl;
	}cout<<endl;

cout<<"circuit Darwinian restructuring start "<<indi->my_string<<endl;
*/	//remove empty columns
	for (int g =0; g < (indi->segmentNumber); g++){
		g_counter = 0;
		for (int h =0; h < indi->ioNumber; h++){
			if (sstring[g][h] == 0)
				g_counter++;
			if (sstring[g][h] == -1 || h == indi->ioNumber-1){
				if (g_counter < indi->ioNumber){
					for (int m = 0; m < indi->ioNumber; m++)
						lastss[scount][m] = sstring[g][m];
					g_counter = 0;
				}
			}
		}	
		//if we did not skip a column of gates
		if (g_counter == 0)
			scount++;
	}
//cout<<"circuit Darwinian restructuring start "<<indi->my_string<<endl;

/*	for (int h =0; h < indi->ioNumber; h++){
		for (int g =0; g < scount; g++){
			cout<<lastss[g][h]<<" ";
		}cout<<endl;
	}cout<<endl;

cout<<"circuit Darwinian restructuring start "<<indi->my_string<<endl;
	for (int g =0; g < scount; g++){
		g_counter = 0;
		for (int h =0; h < indi->ioNumber; h++){
			if (lastss[g][h] == -1) break;
			for(int a = 0; a <  gateArray[lastss[g][h]]->numIO; a++)
				resultSS[a+g_counter][g] = gateArray[lastss[g][h]]->numIO;
			g_counter += gateArray[lastss[g][h]]->numIO;
		}
	}
*/
/*cout<<"circuit Darwinian restructuring start "<<indi->my_string<<endl;
	cout<<endl;
	for (int h =0; h < indi->ioNumber; h++){
		for (int g =0; g < scount; g++){
			cout<<resultSS[h][g]<<" ";
		}cout<<endl;
	}
*/
	//calculate baldwinian cost
	for (int h =0; h < indi->ioNumber; h++)
		for (int g = 1; g < scount; g++){
			if(resultSS[h][g-1] == resultSS[h][g])
				resultSS[h][g-1] = -1;
		}
	for (int h =0; h < indi->ioNumber; h++)
		for (int g = 0; g < scount; g++){
			if (resultSS[h][g] != -1){
				resultSS[h][g] = lastss[g][h];
				if (lastss[g][h] != -1)
					cost+=gateArray[lastss[g][h]]->Cost;
			}
		}
	indi->Cost = cost;
//cout<<"circuit Darwinian restructuring start "<<indi->my_string<<endl;

/*	cout<<endl;
	for (int h =0; h < indi->ioNumber; h++){
		for (int g =0; g < scount; g++){
			cout<<resultSS[h][g]<<" ";
		}cout<<endl;
	}
cout<<"circuit Darwinian restructuring start "<<indi->my_string<<endl;
exit(0);
*/
/*	while(sstring[scount+1][0] != -1){
		//for each segment of the circuit
		for (int g = 0; g < indi->ioNumber; g++) {iresult1[g][0] = -1;iresult2[g][0] = -1;ss1[g][0] = '\0';ss2[g][0] = '\0';}
		iresult1[indi->ioNumber][0] = -1;iresult2[indi->ioNumber][0] = -1;
		result = compareIntSegments(sstring[scount], sstring[scount+1], ss1, ss2, iresult1, iresult2, indi->ioNumber);
cout<<"new seg "<<result<<"   "<<init_time<<" "<<g_counter<<"\n";
		g_counter = 0;
		if (result == 0){
			//if segments are not compatible copy the leftmost to the result string
			sresult += 'p';
			if(init_time){
				//if this is the first segment of the circuit copy it to the resulting circuit
				g_counter = 0;
				cost += analyzeSegment(ss1, indi->ioNumber, rewrite);
				while (g_counter < indi->ioNumber){
					if (ss1[g_counter][0] != '\0'){
						int m = getGate(ss1[g_counter]);
						for (int a = 0; a < REPSIZE; a++){
	                	                       	sresult+= gateArray[m]->representation[a];
						}
					}
					g_counter++;
				}
			} else {
				//copy the last generated segment to the resulting circuit
				g_counter = 0;
				cost += analyzeSegment(lastss, indi->ioNumber, rewrite);
				while (g_counter < indi->ioNumber){
					if (lastss[g_counter][0] != '\0'){
						int m = getGate(lastss[g_counter]);
						for (int a = 0; a < REPSIZE; a++){
                	                	       	sresult+= gateArray[m]->representation[a];
						}
					}
					g_counter++;
				}
			}
			sresult += 'p';
			g_counter = 0;
			//copy the next segment to the last circuit container
			while (g_counter < indi->ioNumber){
				if (ss2[g_counter][0] != '\0'){
					gcount = 0;
					while(ss2[g_counter][gcount] != '\0'){
						lastss[g_counter][gcount] = ss2[g_counter][gcount];
						gcount++;
					}
					lastss[g_counter][gcount] = '\0';
				} else lastss[g_counter][0] = '\0';
				g_counter++;
			}
			//copy the next wire allocation to the last container
			lastss[g_counter][0] = '\0';	
			for (int f = 0; f < indi->ioNumber; f++)
				if(iresult2[f][0] == -1)
					lastiresult[f][0] = iresult2[f][0];
				else for (int u = 0; u < 3; u++)
					lastiresult[f][u] = iresult2[f][u];
			g_counter = 0;
		} else {
			//at least one of the gates between two segments can be combined
			gcount = 0;
			g_counter = 0;
			//chek how many gates are wires
			for (int f = 0; f < indi->ioNumber; f++)
				if (getGate(ss1[f]) == 0)
					g_counter++;
			//if not all gates have been merged and this is first segment of the circuit, copy the left most segment to the resulting string
			if(init_time && g_counter < indi->ioNumber){
				sresult += 'p';
				gcount = 0;
				cost += analyzeSegment(ss1, indi->ioNumber, rewrite);
				while (gcount < indi->ioNumber){
					if (ss1[gcount][0] != '\0'){
						int m = getGate(ss1[gcount]);
						for (int a = 0; a < REPSIZE; a++){
                                	       		sresult+= gateArray[m]->representation[a];
						}
					}
					gcount++;
				}
				sresult += 'p';
				g_counter = 0;
				//copy the next segment to the last circuit container
				while (g_counter < indi->ioNumber){
					if (ss2[g_counter][0] != '\0'){
						gcount = 0;
						while(ss2[g_counter][gcount] != '\0'){
							lastss[g_counter][gcount] = ss2[g_counter][gcount];
							gcount++;
						}
						lastss[g_counter][gcount] = '\0';
					}
					g_counter++;
				}
				lastss[g_counter][0] = '\0';	
				//copy the next wire allocation to the last container
				for (int f = 0; f < indi->ioNumber; f++)
					if(iresult2[f][0] == -1)
						lastiresult[f][0] = iresult2[f][0];
					else for (int u = 0; u < 3; u++)
						lastiresult[f][u] = iresult2[f][u];
			} else {
				//any other cases are parsed here
				if (g_counter > 0){
					//there are some wires 
					temp1 = 0; temp2 = 0;
					wcount = 0;
					for(int y = 0; y < indi->ioNumber; y++){
						condition = true;
						if (iresult1[y][0] < 0){
						 	condition = false;
						}else if (iresult1[y][0] > 0){
							condition = true;
						} 
						//these two gates are mergable
						if(condition){
							temp1 = 0;
							
							while(lastss[y][temp1] != '\0'){
								temp1++;
							}
							temp1 = temp1 - 3;
							if (temp1 < 0) temp1 = 0;
							temp2 = 0;
							while(ss2[y][temp2] != '\0'){
								lastss[y][temp1+temp2] = ss2[y][temp2];
								temp2++;
							}
							lastss[y][temp1+temp2] = '\0';
						} else {
							//if there is no match between the two segments,
							//copy the last merged gate level to the ss1
							// and later copy ss2 the result.
							temp1 = 0;
							if (lastss[y][0] != '\0'){
								while(lastss[y][temp1] != '\0'){
        	                                                        ss1[y][temp1] = lastss[y][temp1];
									temp1++;
                        	                                }
								ss1[y][temp1] = '\0';
							}
							temp1 = 0;
							while(ss2[y][temp1] != '\0'){
                                                                lastss[y][temp1] = ss2[y][temp1];
								temp1++;
                                                        }
							lastss[y][temp1] = '\0';
						}
					}
//cout<<"passed"<<endl;
					//if there is not a complete match
					if (!(g_counter == indi->ioNumber)){
						sresult += 'p';
						gcount = 0;
						cost += analyzeSegment(ss1, indi->ioNumber, rewrite);
//cout<<"passed"<<endl;
						while (gcount < indi->ioNumber){
							if (ss1[gcount][0] != '\0'){
								int m = getGate(ss1[gcount]);
								for (int a = 0; a < REPSIZE; a++){
        	                        		       		sresult+= gateArray[m]->representation[a];
								}
							}
							gcount++;
						}
						sresult += 'p';
					}
        	                        for (int f = 0; f < indi->ioNumber; f++)
                	                        if(iresult2[f][0] == -1)
                        	                        lastiresult[f][0] = iresult2[f][0];
                                	        else for (int u = 0; u < 3; u++)
                                        	        lastiresult[f][u] = iresult2[f][u];
				}
			}
		}
//cout<<"next seg"<<endl;
		scount++;
		if (sstring[scount+1][0] == -1){
			//if this is the last segment
			//finish the circuit string
			cost += analyzeSegment(lastss, indi->ioNumber, rewrite);
			temp1 = 0;
			sresult += 'p';
			while (temp1 < indi->ioNumber){
				if (lastss[temp1][0] != '\0'){
					int n = getGate(lastss[temp1]);
					for (int u = 0; u < REPSIZE; u++){
						sresult+= gateArray[n]->representation[u];;
					}
				}
				temp1++;
			}	
			sresult += 'p';
		}
		init_time = false;
//cout<<"next seg"<<endl;
	}
	sresult += '\0';
	cout<<"circuit restructuring done: "<<sresult<<endl;
//	cout<<"0:"<<sresult<<endl;
//	cout<<"1:"<<sresult<<endl;
	optimized = computeStringCUDAMatrix(sresult,indi->valuedness, indi->ioNumber); 
	orig = computeStringCUDAMatrix(indi->my_string,indi->valuedness, indi->ioNumber); 
//	cout<<"2:"<<sresult<<endl;
	condition = compareGate(orig, optimized);
	cout<<"Equal: "<<condition<<endl;
//	cout<<"3:"<<sresult<<endl;
	cout<<"old"<<endl;
	parseCircuitIntoArray(indi->my_string, indi->ioNumber, sstring);
	cout<<"new"<<endl;
	parseCircuitIntoArray(sresult, indi->ioNumber, sstring);
	int maxcount = (int(pow((float) indi->valuedness, (float) indi->ioNumber)));
	indi->my_string = "";
	for (int m = 0; m < sresult.length(); m++)
		indi->my_string += sresult.at(m);
	indi->Cost = cost;
	cout<<"parsing"<<endl;
//	cout<<"deleting"<<endl;
	for (int g = 0; g < segments; g++)
		delete [] sstring[g];
	cout<<"parsing"<<endl;
//	cout<<"deleting"<<endl;
	for (int g = 0; g < indi->ioNumber+1; g++){
		delete [] iresult1[g];
		delete [] iresult2[g];
		delete [] lastiresult[g];
		delete [] ss1[g];
		delete [] ss2[g];
		delete [] lastss[g];
	}
	cout<<"parsing"<<endl;

	cout<<"circuit restructuring done"<<endl;
*/
}

/*******************************************
 * used for some type of synthesis process
 * in particular using the GY, RY gates
 * requires that the forst two gates are different than the rest
 * and only the first two gates
 ********************************************/
void GA::restructureCircuit(Individual *indi, bool rewrite) {
	string sresult;
	int result = 0;
	int cost = 0;
	int g_counter, gcount, scount, temp1, temp2, wcount;
	bool init_time = true;
	bool condition;

	indi->segmentNumber=getSegmentNumber(indi->my_string);

	int *sstring[2*indi->segmentNumber];
	int *iresult1[indi->ioNumber+1];
	int *iresult2[indi->ioNumber+1];
	int *lastiresult[indi->ioNumber+1];
	char *ss1[indi->ioNumber+1];
	char *ss2[indi->ioNumber+1];
	char *lastss[indi->ioNumber+1];

	qGate *orig;
	qGate *optimized;

	for (int g = 0; g < (2*indi->segmentNumber); g++){
		sstring[g] = new int[indi->ioNumber];
	}
	for (int g = 0; g < indi->ioNumber+1; g++){
		ss1[g] = new char[indi->segmentNumber];
		ss2[g] = new char[indi->segmentNumber];
		lastss[g] = new char[indi->segmentNumber];
		iresult1[g] = new int[3];
		iresult2[g] = new int[3];
		lastiresult[g] = new int[3];
	}
	scount = 0;
	sresult = "";
	g_counter = 0;
cout<<"circuit restructuring start "<<indi->my_string<<endl;
	while (g_counter < indi->ioNumber) lastss[g_counter++][0] = '\0';
	parseCircuitIntoArray(indi->my_string, indi->ioNumber, sstring);

	for (int g =0; g < (indi->segmentNumber); g++){
		for (int h =0; h < indi->ioNumber; h++){
			cout<<sstring[g][h]<<" ";
		}cout<<", ";
	}
	cout<<endl;
	while(sstring[scount+1][0] != -1){
		//for each segment of the circuit
		for (int g = 0; g < indi->ioNumber; g++) {iresult1[g][0] = -1;iresult2[g][0] = -1;ss1[g][0] = '\0';ss2[g][0] = '\0';}
		iresult1[indi->ioNumber][0] = -1;iresult2[indi->ioNumber][0] = -1;
		result = compareIntSegments(sstring[scount], sstring[scount+1], ss1, ss2, iresult1, iresult2, indi->ioNumber);
cout<<"new seg "<<result<<"   "<<init_time<<" "<<g_counter<<"\n";
		g_counter = 0;
		if (result == 0){
			//if segments are not compatible copy the leftmost to the result string
			sresult += 'p';
			if(init_time){
				//if this is the first segment of the circuit copy it to the resulting circuit
				g_counter = 0;
				cost += analyzeSegment(ss1, indi->ioNumber, rewrite);
				while (g_counter < indi->ioNumber){
					if (ss1[g_counter][0] != '\0'){
						int m = getGate(ss1[g_counter]);
						for (int a = 0; a < REPSIZE; a++){
	                	                       	sresult+= gateArray[m]->representation[a];
						}
					}
					g_counter++;
				}
			} else {
				//copy the last generated segment to the resulting circuit
				g_counter = 0;
				cost += analyzeSegment(lastss, indi->ioNumber, rewrite);
				while (g_counter < indi->ioNumber){
					if (lastss[g_counter][0] != '\0'){
						int m = getGate(lastss[g_counter]);
						for (int a = 0; a < REPSIZE; a++){
                	                	       	sresult+= gateArray[m]->representation[a];
						}
					}
					g_counter++;
				}
			}
			sresult += 'p';
			g_counter = 0;
			//copy the next segment to the last circuit container
			while (g_counter < indi->ioNumber){
				if (ss2[g_counter][0] != '\0'){
					gcount = 0;
					while(ss2[g_counter][gcount] != '\0'){
						lastss[g_counter][gcount] = ss2[g_counter][gcount];
						gcount++;
					}
					lastss[g_counter][gcount] = '\0';
				} else lastss[g_counter][0] = '\0';
				g_counter++;
			}
			//copy the next wire allocation to the last container
			lastss[g_counter][0] = '\0';	
			for (int f = 0; f < indi->ioNumber; f++)
				if(iresult2[f][0] == -1)
					lastiresult[f][0] = iresult2[f][0];
				else for (int u = 0; u < 3; u++)
					lastiresult[f][u] = iresult2[f][u];
			g_counter = 0;
		} else {
			//at least one of the gates between two segments can be combined
			gcount = 0;
			g_counter = 0;
			//chek how many gates are wires
			for (int f = 0; f < indi->ioNumber; f++)
				if (getGate(ss1[f]) == 0)
					g_counter++;
			//if not all gates have been merged and this is first segment of the circuit, copy the left most segment to the resulting string
			if(init_time && g_counter < indi->ioNumber){
				sresult += 'p';
				gcount = 0;
				cost += analyzeSegment(ss1, indi->ioNumber, rewrite);
				while (gcount < indi->ioNumber){
					if (ss1[gcount][0] != '\0'){
						int m = getGate(ss1[gcount]);
						for (int a = 0; a < REPSIZE; a++){
                                	       		sresult+= gateArray[m]->representation[a];
						}
					}
					gcount++;
				}
				sresult += 'p';
				g_counter = 0;
				//copy the next segment to the last circuit container
				while (g_counter < indi->ioNumber){
					if (ss2[g_counter][0] != '\0'){
						gcount = 0;
						while(ss2[g_counter][gcount] != '\0'){
							lastss[g_counter][gcount] = ss2[g_counter][gcount];
							gcount++;
						}
						lastss[g_counter][gcount] = '\0';
					}
					g_counter++;
				}
				lastss[g_counter][0] = '\0';	
				//copy the next wire allocation to the last container
				for (int f = 0; f < indi->ioNumber; f++)
					if(iresult2[f][0] == -1)
						lastiresult[f][0] = iresult2[f][0];
					else for (int u = 0; u < 3; u++)
						lastiresult[f][u] = iresult2[f][u];
			} else {
				//any other cases are parsed here
				if (g_counter > 0){
					//there are some wires 
					temp1 = 0; temp2 = 0;
					wcount = 0;
					for(int y = 0; y < indi->ioNumber; y++){
						condition = true;
						if (iresult1[y][0] < 0){
						 	condition = false;
						}else if (iresult1[y][0] > 0){
							condition = true;
						} 
						//these two gates are mergable
						if(condition){
							temp1 = 0;
							
							while(lastss[y][temp1] != '\0'){
								temp1++;
							}
							temp1 = temp1 - 3;
							if (temp1 < 0) temp1 = 0;
							temp2 = 0;
							while(ss2[y][temp2] != '\0'){
								lastss[y][temp1+temp2] = ss2[y][temp2];
								temp2++;
							}
							lastss[y][temp1+temp2] = '\0';
						} else {
							//if there is no match between the two segments,
							//copy the last merged gate level to the ss1
							// and later copy ss2 the result.
							temp1 = 0;
							if (lastss[y][0] != '\0'){
								while(lastss[y][temp1] != '\0'){
        	                                                        ss1[y][temp1] = lastss[y][temp1];
									temp1++;
                        	                                }
								ss1[y][temp1] = '\0';
							}
							temp1 = 0;
							while(ss2[y][temp1] != '\0'){
                                                                lastss[y][temp1] = ss2[y][temp1];
								temp1++;
                                                        }
							lastss[y][temp1] = '\0';
						}
					}
//cout<<"passed"<<endl;
					//if there is not a complete match
					if (!(g_counter == indi->ioNumber)){
						sresult += 'p';
						gcount = 0;
						cost += analyzeSegment(ss1, indi->ioNumber, rewrite);
//cout<<"passed"<<endl;
						while (gcount < indi->ioNumber){
							if (ss1[gcount][0] != '\0'){
								int m = getGate(ss1[gcount]);
								for (int a = 0; a < REPSIZE; a++){
        	                        		       		sresult+= gateArray[m]->representation[a];
								}
							}
							gcount++;
						}
						sresult += 'p';
					}
        	                        for (int f = 0; f < indi->ioNumber; f++)
                	                        if(iresult2[f][0] == -1)
                        	                        lastiresult[f][0] = iresult2[f][0];
                                	        else for (int u = 0; u < 3; u++)
                                        	        lastiresult[f][u] = iresult2[f][u];
				}
			}
		}
//cout<<"next seg"<<endl;
		scount++;
		if (sstring[scount+1][0] == -1){
			//if this is the last segment
			//finish the circuit string
			cost += analyzeSegment(lastss, indi->ioNumber, rewrite);
			temp1 = 0;
			sresult += 'p';
			while (temp1 < indi->ioNumber){
				if (lastss[temp1][0] != '\0'){
					int n = getGate(lastss[temp1]);
					for (int u = 0; u < REPSIZE; u++){
						sresult+= gateArray[n]->representation[u];;
					}
				}
				temp1++;
			}	
			sresult += 'p';
		}
		init_time = false;
//cout<<"next seg"<<endl;
	}
	sresult += '\0';
	cout<<"circuit restructuring done: "<<sresult<<endl;
//	cout<<"0:"<<sresult<<endl;
//	cout<<"1:"<<sresult<<endl;
	optimized = computeStringCUDAMatrix(sresult,indi->valuedness, indi->ioNumber); 
	orig = computeStringCUDAMatrix(indi->my_string,indi->valuedness, indi->ioNumber); 
//	cout<<"2:"<<sresult<<endl;
	condition = compareGate(orig, optimized);
	cout<<"Equal: "<<condition<<endl;
//	cout<<"3:"<<sresult<<endl;
	cout<<"old"<<endl;
	parseCircuitIntoArray(indi->my_string, indi->ioNumber, sstring);
	cout<<"new"<<endl;
	parseCircuitIntoArray(sresult, indi->ioNumber, sstring);
	int maxcount = (int(pow((float) indi->valuedness, (float) indi->ioNumber)));
	indi->my_string = "";
	for (int m = 0; m < sresult.length(); m++)
		indi->my_string += sresult.at(m);
	indi->Cost = cost;
	cout<<"parsing"<<endl;
//	cout<<"deleting"<<endl;
	for (int g = 0; g < segments; g++)
		delete [] sstring[g];
	cout<<"parsing"<<endl;
//	cout<<"deleting"<<endl;
	for (int g = 0; g < indi->ioNumber+1; g++){
		delete [] iresult1[g];
		delete [] iresult2[g];
		delete [] lastiresult[g];
		delete [] ss1[g];
		delete [] ss2[g];
		delete [] lastss[g];
	}
	cout<<"parsing"<<endl;

	cout<<"circuit restructuring done"<<endl;
}

/*************************************************
 * this cuts the amount of segments to 200
 *************************************************/
void GA::reduceSegNumber(int a) {
	int counter, n;

	counter = 0;
	n = 0;

	if (population[a]->ioNumber > 1) {
		while (n < (int) population[a]->my_string.length()) {
			if (population[a]->my_string.at(n) == 'p')
				counter++;
			if (counter >= 200) {
				population[a]->my_string.erase(n,
						population[a]->my_string.length());
				break;
			}
			n++;
		}
	}
}

/*******************************************
 * the circuit minimization
 * eliminates gate segments that are neighboring and are identical
 *******************************************/
int GA::minimizeCirc(Individual *qC) {
	int result = 0;
	int diff0, diff1;
	int begin0, begin1, end0, end1;
	begin0 = 0;
	end0 = qC->my_string.find("p", begin0 + 1);
	do {
		if (end0 < qC->my_string.length()-1){
			begin1 = end0+1;
			end1 = qC->my_string.find("p", begin1 + 1);
		} else break;
		diff0 = end0-begin0;
		diff1 = end1-begin1;
		
		if(diff0 == diff1){
			if (qC->my_string.compare(begin0,diff0, qC->my_string, begin1, diff1) == 0){
				qC->my_string.erase(begin0, diff0+diff1 + 2);
				end0 = qC->my_string.find("p", begin0 + 1);
			}
		}
		if (end1 < qC->my_string.length()-1){
			begin0 = begin1;
			end0 = end1;
		} else break;
	}while(true);
	return result;
}

/*******************************************
 * analyzes a set of gates from a segments 
 * and if necesserry creates a new gate
 * can do - physically modyfies the circuit for lamarckian learning
 *        - simulates baldwinian cost calculation
 *******************************************/
int GA::analyzeSegment(char **segment, int wires, bool create){
	int y, gcount, nscount, scount, state, z, ww, cost, m, wiresifsame;
	char gA[REPSIZE], gB[REPSIZE];
	char newseg[segments];
	int nnewseg[segments];
	qGate *A, *B, *C, *D;
	rGate rgate;
	bool required, initial, equal;
	bool changed[wires];
	string new_my_string;
	string old_my_string;
	C = new qGate;
	D = new qGate;
	initGate(C, resultnum, gateArray[0]->valuedness, 1);
	initGate(D, resultnum, gateArray[0]->valuedness, 1);
	cost = 0;
cout<<"analyzing segment "<<numofgates<<"  "<<create<<endl<<"  "<<segment[0]<<','<<segment[1]<<','<<segment[2]<<','<<segment[3]<<','<<endl;
	if (create){
		//if lamarckian learning is allowed
		for (int x = 0; x < wires; x++){
			y = 0;
			gcount = 0;
			required = false;
			initial = true;
			nscount = 0;
			scount = 0;
			if (segment[x][0] != '\0' && segment[x][0] != ' ') {
				for (int j = 0; j < REPSIZE; j++)
					gA[gcount++] = segment[x][y++];
				gA[gcount] = '\0';
				new_my_string = "";
				changed[x] = false;
				while(true){
					if (segment[x][y] == ' ' && segment[x][y] != '\0') required = true;	
					else required = false;
					state = 0;
					if (required){
						gcount = 0;y++;
						for (int j = 0; j < REPSIZE; j++)
							gB[gcount++] = segment[x][y++];
						gB[gcount] = '\0';
							cout<<"0X "<<gA<<"  "<<gB<<endl;
						if (getGate(gA) == 0){
							cout<<"00"<<endl;
							copyRep(gB, gA);
							changed[x] = true;
							cout<<"01"<<endl;
						}else {
							cout<<"02 "<<gA<<endl;
								nnewseg[scount++] = getGate(gA);
							for (int j = 0; j < REPSIZE; j++)
								newseg[nscount++] = gA[j];
							ww = getGate(gA);
							if (ww > -1)
							for (int i = 0; i < gateArray[ww]->my_string.length(); i++)
								new_my_string += gateArray[ww]->my_string.at(i);
							copyRep(gB, gA);
							changed[x] = true;
							state = 1;
							cout<<"03"<<endl;
						}
					} else {
						cout<<"here"<<endl;
						nnewseg[scount++] = getGate(gA);
						wiresifsame = 1;
						for (int j = 0; j < REPSIZE; j++)
							newseg[nscount++] = gA[j];
						newseg[nscount++] = ' ';
						break;
					}
				}
				newseg[nscount] = '\0';
				nnewseg[scount] = -1;
//				for (int j = 0; j < scount; j++)
//					cout<<nnewseg[j]<<", ";
//				cout<<endl;
//				cout<<"new seg :"<<newseg<<", changed "<<changed[x]<<", scount: "<<scount<<endl;
				if(changed[x]){
					//if a new gate needs to be created
					y = 0;
					m = 0;
					required = false;
					initial = true;
					if (scount == 1){
//			cout<<"wiresifsame "<<wiresifsame<<endl;
						for(int p = 0; p < wiresifsame;p++)
							copyRep(newseg, segment[x+p]);
						for (int i = 0; i < gateArray[nnewseg[0]]->my_string.length(); i++)
							new_my_string += gateArray[nnewseg[0]]->my_string.at(i);
					} else 
					while (m  < scount){
						required = true;
						A = gateArray[nnewseg[m++]];
						
						if (initial){
//				cout<<"0initial new seg :"<<newseg<<endl;
							B = gateArray[nnewseg[m++]];
//				cout<<"0.35initial new seg :"<<newseg<<endl;
							C = ccGate(A,C);
							D = ccGate(A,D);
//				cout<<"0.5initial new seg :"<<newseg<<endl;
							if ((int)(pow((float)(A->valuedness), (float)(A->numIO))) > BLOCK_MAX_SIZE)
								gpuMMmult(A->gateMatrix1,  B->gateMatrix1, A->numIO, A->valuedness, C->gateMatrix1);
							else
								cblasMatrixProduct(A, B, C);
							for (int i = 0; i < A->my_string.length(); i++)
								new_my_string += A->my_string.at(i);	
							for (int i = 0; i < B->my_string.length(); i++)
								new_my_string += B->my_string.at(i);	
							initial = false;
//				cout<<"1initial new seg :"<<newseg<<endl;

						} else {
//				cout<<"2required new seg :"<<newseg<<endl;
							if ((int)(pow((float)(C->valuedness), (float)(C->numIO))) > BLOCK_MAX_SIZE)
								gpuMMmult(C->gateMatrix1, A->gateMatrix1,  A->numIO, A->valuedness, D->gateMatrix1);
							else
								cblasMatrixProduct(C, A, D);

							C = ccGate(D,C);
							for (int i = 0; i < A->my_string.length(); i++)
								new_my_string += A->my_string.at(i);	
						}
//				cout<<"1required new seg :"<<newseg<<endl;
					}
					if (required){
						//check if such gate already exists
						for (int u = 0; u < numofgates; u++){
							equal = false;
							if (gateArray[u]->numIO == C->numIO){
								equal = compareGate(gateArray[u],C);
								if (equal) {
									y = u;
									break;
								}
							}
						}
						if (!equal){
							copyRep(gateArray[numofgates-1]->representation, gB);
							incRep(gB);
							copyRep(gB, C->representation);
							rgate.index = numofgates;
							rankV.push_back(rgate);	
							new_my_string += '\0';
							C->my_string = "";
							cost += C->Cost;
							for (int i = 0; i < new_my_string.length(); i++)
								C->my_string += new_my_string.at(i);	
					 		gateArray[numofgates] = new qGate;
							initGate(gateArray[numofgates], C->numIO, C->valuedness, 1);
							gateArray[numofgates] = ccGate(C, gateArray[numofgates]);
							numofgates++;
							copyRep(gB,newseg);
						} else {
							copyRep(gateArray[y]->representation,newseg);	
							cost += gateArray[y]->Cost;
						}
					}
				} else {
					if (segment[x][0] != '\0' && segment[x][0] != ' ') {
						gcount = 0;
						for (int j = 0; j < REPSIZE; j++)
							gA[gcount++] = segment[x][j];
						gA[gcount] = '\0';
						cost += gateArray[getGate(gA)]->Cost;
					}
				}
				copyRep(newseg, segment[x]);
//				cout<<newseg<<endl;
//				cout<<"done modifying"<<endl;
			}
		}
//cout<<"analyzing segment done"<<endl;
	} else {
		//if only baldwinian cost calculation si allowed
		for (int x = 0; x < wires; x++){
			if (segment[x][0] != '\0' && segment[x][0] != ' ') {
				gcount = 0;
				for (int j = 0; j < REPSIZE; j++)
					gA[gcount++] = segment[x][j];
				gA[gcount] = '\0';
				cost += gateArray[getGate(gA)]->Cost;
			}
		}
	}
	destroyGate(C);
	destroyGate(D);
	delete C;
	C = NULL;
	delete D;
	D = NULL;
//cout<<"finalizing segment"<<endl;//segment[0]<<','<<segment[1]<<','<<segment[2]<<','<<segment[3]<<','<<segment[4]<<endl;
//	cout<<"analysis done"<<endl;
return cost;
}
/*******************************************
 * calculates the average of fitness and error
 *******************************************/
void GA::calcAverage() {

	int counter = finalResult.counter;
	float avfitness = 0;
	float averror = 0;
	float avcost = 0;
	for (int t = 0; t < populationNumber; t++) {
		avfitness += population[t]->fitness;
		averror += population[t]->Error;
		avcost += population[t]->Cost;
	}
	finalResult.avFitness[counter] = avfitness / populationNumber;
	finalResult.avError[counter] = averror / populationNumber;
	finalResult.avCost[counter] = avcost / populationNumber;
	finalResult.counter = counter + 1;
}
/*******************************************
 * recalculates the usage the usage factor and usefulness of used gate
 *******************************************/
void GA::rankInputGates(Individual *indi){
	string s = indi->my_string;
	char *gate_rep = new char[4];
	int begin = 0;
	int csegs = 0;
	int num_g;
	int cost;
	float fit = indi->fitness;
        int end = s.find("p", begin + 1);
        indi->segmentNumber = 0;
        while (end > 0) {
                if ((end - begin) > 1) {
                        csegs++;
                        //get the first  gate of the Segment and copy it
			getGateRepfromStr(begin+1,s, gate_rep);
			num_g = getGate(gate_rep);
			for(int ii=0; ii < numofgates; ii++)
			{
				if (num_g == rankV[ii].index){
					rankV[ii].fitness += fit;
					rankV[ii].fitness /= 2.0;
					rankV[ii].usage++;
					break;
				}
			}

/*			for(int ii=0; ii < numofgates; ii++)
				if (rankV[ii].index == 0){
					rankV[ii].fitness = 0.5;
					break;
				}
  */                      //continue with the cost
                        cost += gateArray[num_g]->Cost;
                        //get the next one
			begin += 4;
                        for (int b = begin; b < end; b+=REPSIZE) {
                        }
                }
                //move to the next segment
                begin = indi->my_string.find("p", end);
                end = indi->my_string.find("p", begin + 1);
        }
	indi->segmentNumber = csegs;
}
/*******************************************
 * Sort the Vector with the gates values
 *******************************************/
void GA::sortGates(){
	int size = numofgates;
	int rank = 0;
	rGate rgate;

	for (int m = 0; m < size; m++){
		rank = m;
		for (int n = m+1; n < size; n++){
			if (rankV[n].fitness > rankV[rank].fitness)
				rank = n;
		}
		if (rank != m){
			rgate.index = rankV[rank].index;
			rgate.usage = rankV[rank].usage;
			rgate.fitness = rankV[rank].fitness;

			rankV[rank].index = rankV[m].index;
			rankV[rank].usage = rankV[m].usage;
			rankV[rank].fitness = rankV[m].fitness;
				
			rankV[m].index = rgate.index;
			rankV[m].usage = rgate.usage;
			rankV[m].fitness = rgate.fitness;
		}
	}
}
/*******************************************
 * Searches for a gate with the fitness closes to the fitdesired
 *******************************************/
int GA::getBestGateMatch(float fitdesired, int topwire, int wires, int limit, bool use_restrictions){
	int s = numofgates;
	int best = 0;
	float fsum = 0;
	bool candidate;
//cout<<"getting gates: "<<fitdesired<<", "<<wires<<", "<<limit<<", ";


	//find gate if no restrictions are required
	if (use_restrictions){
		if (wires == 0){
			if (limit > 0){
				for (int y = s-1; y >=0; y--){
					if (fsum >= fitdesired && gateArray[rankV[y].index]->numIO <= limit)
						return rankV[y].index;		
					if (fsum < fitdesired && gateArray[rankV[y].index]->numIO <= limit)
						best = y;
					fsum += rankV[y].fitness;
				}
//				cout<<best<<endl;
				return rankV[best].index;
			}else 
				for (int y = s-1; y >=0; y--){
					if (fsum >= fitdesired)
						return rankV[y].index;		
					fsum += rankV[y].fitness;
				}
		} else {
			for (int y = s-1; y >=0; y--){
				if (fsum >= fitdesired && gateArray[rankV[y].index]->numIO == wires)
					return rankV[y].index;
				if (gateArray[rankV[y].index]->numIO == wires && fsum < fitdesired)
					best = y;
				fsum += rankV[y].fitness;
			}
			return rankV[best].index;
		}
	} else {
		if (wires == 0){
			if (limit > 0){
				for (int y = s-1; y >=0; y--){
					candidate = true;
					for (int m = 0; m < gateArray[rankV[y].index]->restrictions_number; m++){
         				       if (topwire == gateArray[rankV[y].index]->restrictions[m])
			                        candidate = false;
						break;
				        }
					if (candidate){
						if (fsum >= fitdesired && gateArray[rankV[y].index]->numIO <= limit)
							return rankV[y].index;		
						if (fsum < fitdesired && gateArray[rankV[y].index]->numIO <= limit)
							best = y;
//						if (gateArray[rankV[y].index]->numIO <= limit)
//							best = y;
						fsum += rankV[y].fitness;
					}
				}
				cout<<best<<endl;
				return rankV[best].index;
			}else 
				for (int y = s-1; y >=0; y--){
					candidate = true;
					for (int m = 0; m < gateArray[rankV[y].index]->restrictions_number; m++){
         				       if (topwire == gateArray[rankV[y].index]->restrictions[m])
			                        candidate = false;
						break;
				        }

					if (candidate){
						if (fsum >= fitdesired)
							return rankV[y].index;
						best = y;
						fsum += rankV[y].fitness;
					}
				}
				return rankV[best].index;
		} else {
			for (int y = s-1; y >=0; y--){
				candidate = true;
				for (int m = 0; m < gateArray[rankV[y].index]->restrictions_number; m++){
       				       if (topwire == gateArray[rankV[y].index]->restrictions[m])
		                        candidate = false;
					break;
			        }
				if (candidate){
					if (fsum >= fitdesired && gateArray[rankV[y].index]->numIO == wires)
						return rankV[y].index;
					if (gateArray[rankV[y].index]->numIO == wires && fsum < fitdesired)
						best = y;
					fsum += rankV[y].fitness;
				}
			}
			return rankV[best].index;
		}
	}
	return -1;
}
/*******************************************
 * Initiate the length of the arrays for history records
 *******************************************/
void GA::initiateStorage() {
	cout << "Initiating storage for " << generations << " generations";
	history = new float[generations][4];
	cout<<" ..... done"<<endl;
}

/**********************************
* decodes the name of the gate into a RevLib compatible format
**********************************/
void GA::decodeIndvStrToRevLibStr(Individual *I){
	string rev_lib_circuit;
	int c_count = 0;
	string u;
	int i_gate,i,j,k;
	int wire_counter = 0;
	char buff[10], *gate_rep;
	int y, x, m;
	int indexes[I->ioNumber];
	rev_lib_circuit = ".version Curr Individual\n";
	
	rev_lib_circuit += ".numvars ";
	sprintf(buff, "%d", I->ioNumber);
	x = 0;
	while(buff[x] > 0)
		rev_lib_circuit += buff[x++];	
	rev_lib_circuit += "\n.variables ";
	for (y = 0; y < I->ioNumber; y++){
		sprintf(buff, "%d", y);
		x = 0;
		rev_lib_circuit += 'w';	
		while(buff[x] > 0)
			rev_lib_circuit += buff[x++];	
		rev_lib_circuit += ' ';	
	}
	rev_lib_circuit += "\n.inputs ";
	for (y = 0; y < I->ioNumber; y++){
		sprintf(buff, "%d", y);
		x = 0;
		rev_lib_circuit += 'w';	
		while(buff[x] > 0)
			rev_lib_circuit += buff[x++];	
		rev_lib_circuit += ' ';	
	}
	rev_lib_circuit += "\n.outputs ";
	for (y = 0; y < I->ioNumber; y++){
		sprintf(buff, "%d", y);
		x = 0;
		rev_lib_circuit += 'w';	
		while(buff[x] > 0)
			rev_lib_circuit += buff[x++];	
		rev_lib_circuit += ' ';	
	}
	rev_lib_circuit += "\n.begin\n";
	m = 0;
	gate_rep = new char[3];
	while(true){
		if (m >= I->my_string.size())
			break;
		if (I->my_string.at(m) != 'p'){
			getGateRepfromStr(m, I->my_string, gate_rep);
			i_gate = getGate(gate_rep);	
cout<<"revlib: "<<gate_rep[0]<<gate_rep[1]<<gate_rep[2]<<" <>  "<<gateArray[i_gate]->rev_lib_char<<endl;
			if (i_gate != 0){
				for (i = 0; i < gateArray[i_gate]->rev_lib_char.length(); i++)
					rev_lib_circuit += gateArray[i_gate]->rev_lib_char.at(i);
				rev_lib_circuit += ' ';
				if (wire_counter == 0){
					for (i = 0; i < gateArray[i_gate]->rev_lib_string.length(); i++){
						if (i == 0)
							rev_lib_circuit += 'w';
						rev_lib_circuit += gateArray[i_gate]->rev_lib_string.at(i);
						if (gateArray[i_gate]->rev_lib_string.at(i) == ' ' && i !=  0)
							rev_lib_circuit += 'w';
					}
					rev_lib_circuit += '\n';
					wire_counter += gateArray[i_gate]->numIO;
				}else {
					for (i = 0; i < I->ioNumber; i++)
						indexes[i] = -1;
					StrToIntArray(gateArray[i_gate]->rev_lib_string, indexes);
					for (i = 0; i < I->ioNumber; i++){
						if (indexes[i] == -1)
							break; 					
						indexes[i] += wire_counter;
					}
					y = 0;
					do {
						sprintf(buff, "%d", indexes[y]);
						x = 0;
						rev_lib_circuit += 'w';	
						while(buff[x] > 0)
							rev_lib_circuit += buff[x++];
						rev_lib_circuit += ' ';
					}while(indexes[++y] > -1);
					rev_lib_circuit += '\n';
					wire_counter += gateArray[i_gate]->numIO;
				}
			} else wire_counter++;
			m+=3;
		} else {
			wire_counter = 0;
			m++;
		}
	}
	rev_lib_circuit += ".end";
	I->rev_lib_string = "";
	for (y = 0; y < rev_lib_circuit.length(); y++)
		I->rev_lib_string += rev_lib_circuit.at(y);
	delete [] gate_rep;
	cout<<I->rev_lib_string<<endl;
}

void GA::StrToIntArray(string s, int *data){
	string buff;
	int begin, data_count, y, x;
	begin = 0;
	data_count = 0;
	y = 0;
	do {
		if (s.at(y) == ' '){
			buff = "";
			for (x = begin; x < y; x++)
				buff += s.at(x);
			begin = y;
			data[data_count++] = atoi(buff.c_str()); 
		} else if (y == s.length()-1){
                        buff = "";
                        for (x = begin; x <= y; x++)
                                buff += s.at(x);
                        begin = y;
                        data[data_count++] = atoi(buff.c_str());
                }
	}while(++y <s.length());

}

void GA::IntArrayToStr(int *data, string str){
	char buff[10];
	int y, x;
	y = 0;
	do {
		sprintf(buff, "%d", data[y]);
		x = 0;
		while(buff[x] > 0)
		str += buff[x++];
		str += ' ';

	}while(data[++y] > -1);

	cout<<"in: "<< str.length()<<":out"<<endl;
}
