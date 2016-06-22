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
 * calculate fitness for individual:
 ****************************************/
void GA::makeFitness(aCircuit *indi) {

	if (indi->error != 0) {
		switch (replicator) {
		case 0:
			//simple fitness1
			indi->fitness = (1 - indi->error);
			break;
		case 1:
			//simple fitness1
			indi->fitness = (1 / (indi->error + 1));
			break;
		case 2:
			//scaled complex fitness1
			indi->fitness = (alpha * (1 - indi->error) + beta * (exp(-pow(abs(
					divider - indi->cost), 2))));
			break;
		case 3:
			//scaled complex fitness2
			indi->fitness = (alpha * (1 / (indi->error + 1)) + beta * (exp(
					-pow(abs(divider - indi->cost), 2))));
			break;
		case 4:
			indi->fitness = (1 - indi->error);
			break;
		case 5:
			indi->fitness = exp(-indi->error);
			break;
		}
	} else {
		indi->error = 0;
		indi->fitness = 1;//we found a good individual;
	}

}

/****************************************
 * static - the function called to invoke the program as a thread
 ****************************************/
void* GA::startSub(void *arg) {
	GA* obj = (GA*) arg;
	return 0;
}
void GA::makeThreeWiresRestrictedC(aCircuit *I){
	I->connections[0]->input_index = 0;
	I->connections[0]->output_index = 3;
	I->connections[1]->input_index = 1;
	I->connections[1]->output_index = 4;
	I->connections[2]->input_index = 2;
	I->connections[2]->output_index = 9;
	I->connections[3]->input_index = 5;
	I->connections[3]->output_index = 8;
	I->connections[4]->input_index = 6;
	I->connections[4]->output_index = 14;
	I->connections[5]->input_index = 7;
	I->connections[5]->output_index = 20;
	I->connections[6]->input_index = 10;
	I->connections[6]->output_index = 13;
	I->connections[7]->input_index = 11;
	I->connections[7]->output_index = 19;
	I->connections[8]->input_index = 12;
	I->connections[8]->output_index = 15;
	I->connections[9]->input_index = 16;
	I->connections[9]->output_index = 18;
	I->connections[10]->input_index = 17;
	I->connections[10]->output_index = 25;
	I->connections[11]->input_index = 21;
	I->connections[11]->output_index = 23;
	I->connections[12]->input_index = 22;
	I->connections[12]->output_index = 24;
}
void GA::makeThreeWiresSimpleC(aCircuit *I){
	I->connections[0]->input_index = 0;
	I->connections[0]->output_index = 3;
	I->connections[1]->input_index = 1;
	I->connections[1]->output_index = 4;
	I->connections[2]->input_index = 2;
//	I->connections[2]->output_index = 9;
//	I->connections[3]->input_index = 5;
//	I->connections[3]->output_index = 8;
//	I->connections[6]->input_index = 10;
//	I->connections[6]->output_index = 13;
//	I->connections[9]->input_index = 16;
//	I->connections[9]->output_index = 18;
	I->connections[10]->input_index = 17;
	I->connections[10]->output_index = 25;
	I->connections[11]->input_index = 21;
	I->connections[11]->output_index = 23;
	I->connections[12]->input_index = 22;
	I->connections[12]->output_index = 24;
}
/****************************************
 * this generates the String representing every individual.
 *it is called upon the structure of the initial population
 ****************************************/
void GA::makeIndvWiring(aCircuit *I) {
	int tempNumber, tempNumber2, tempGate, counter, conns;
	int data_output[I->numIO];
	counter = 0;
	for (int m = 0; m < I->numIO; m++){
		data_output[m] = I->data_O[m];
//		cout<<data_output[m]<<"  ";
	}
//	cout<<endl;

	
	for (int a = 0; a < I->numIO; a++){
		I->connections[a]->unmutable = true;
		I->connections[a]->input_index = a;
	}
	if (!doubled){
		for (int a = 0; a < I->numIO; a++){
			I->connections[I->numI-a-1]->unmutable = true;
			I->connections[I->numI-a-1]->output_index =I->data_O[I->numIO-a-1];
		}
	} else {
		for (int a = 0; a < I->numIO-1; a++){
			
			I->connections[I->numI-a-1]->unmutable = true;
			I->connections[I->numI-a-1]->output_index =I->data_O[I->numIO-a-1];
		}
		if (I->data_O[I->numIO-2]-1 == I->data_O[I->numIO-1])
			I->connections[I->numI-I->numIO]->output_index =I->data_O[0];
		else
			I->connections[I->numI-I->numIO-1]->output_index =I->data_O[0];
	}

	for (int j = 0; j < I->numGates; j++) // for each segment of this individual
	{
		while(!gateInputsSat(I->connections, I->numI, I->gates[j], I->doubled)){
			generateConnections(I->connections, I->numI, I->gates[j], I->doubled, I->cascade);
//cout<<"gate: "<<j<<", negated: "<<I->gates[j]->negated_vars<<endl;
        for (int a = 0; a < I->numI; a++)
                cout<<a<<": "<<I->connections[a]->input_index<<", "<<I->connections[a]->output_index<<". "<<I->connections[a]->negated<<endl;

			if (counter > I->numGates)	break;
			counter++;
		}
//	for (int m = 0; m < I->numI; m++)
//		cout<<I->connections[m]->input_index<<": "<<I->connections[m]->output_index<<endl;
	}
//	cout<<"original wireing"<<endl;
//	for (int m = 0; m < I->numI; m++)
//		cout<<I->connections[m]->input_index<<": "<<I->connections[m]->output_index<<endl;
	
	data_output[0] = -1;
	for (int m = 0; m < I->numIO-1; m++){
//	for (int m = 0; m < I->numIO; m++){
		do{
			conns = rand()%I->numIO;
		}while(data_output[conns] == -1);
		for (int k = 0; k < I->numI; k++)
			if(I->connections[k]->output_index == -1){
				I->connections[k]->output_index = data_output[conns];
				data_output[conns] = -1;
				break;
			}
	}
//	cout<<"wireing done"<<endl;
//	for (int m = 0; m < I->numI; m++)
//		cout<<I->connections[m]->input_index<<": "<<I->connections[m]->output_index<<endl;
//	for (int a = 0; a < I->numI; a++)
//		cout<<a<<": "<<I->connections[a]->input_index<<", "<<I->connections[a]->output_index<<". "<<I->connections[a]->negated<<endl;

//exit(0);

}

/****************************************
 * Creates the initial population of individuals
 ****************************************/
void GA::initializePop(int level, string filename, string out_prefix) {

	string inFile = filename;
	string outFile;
	string post = (".out");
	string oFile = "";
	time_t rawtime;
	struct tm * timeinfo;
	int wire_indexes, conn_index, counter;
	condition = false;

	// Open input and output files.
	in_stream.open(inFile.c_str());
	if (in_stream.fail()) {
		cout << "\"Daa set \" file opening failed.\n";
		exit(1);
	}

	outFile = inFile;

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
		cout << "Output file opening failed." << outFile.c_str() << "\n";
		exit(1);
	}
	//initialize the gates
	initializeGates();
	cout << "Intial population of circuits generation\n";
	GA::populationCounter = 0;
	//generate the initial population/new_population of individuals
	for (int i = 0; i < populationNumber; i++) {
		population[i] = new aCircuit;
		new_population[i] = new aCircuit;
		population[i]->numIO = finalGate->numI;
		new_population[i]->numIO = finalGate->numI;
		if (doubled > 0){
			segments *= 2;
			population[i]->doubled = true;
			new_population[i]->doubled = true;
		}
		population[i]->numGates = segments;
		new_population[i]->numGates = segments;
		population[i]->numO = population[i]->numIO+gateArray[0]->numO*segments/2+gateArray[1]->numO*segments/2;
		new_population[i]->numO = new_population[i]->numIO+gateArray[0]->numO*segments/2+gateArray[1]->numO*segments/2;
		population[i]->numI = population[i]->numIO+gateArray[0]->numI*segments/2+gateArray[1]->numI*segments/2;
		new_population[i]->numI = new_population[i]->numIO+gateArray[0]->numI*segments/2+gateArray[1]->numI*segments/2;
		if(cascade > 0){
			population[i]->cascade = true;
			new_population[i]->cascade = true;
		}
		population[i]->I = new int[population[i]->numI];
		new_population[i]->I = new int[new_population[i]->numI];
		population[i]->O = new int[population[i]->numO];
		new_population[i]->O = new int[new_population[i]->numO];
		population[i]->data_I = new int[population[i]->numIO];
		new_population[i]->data_I = new int[new_population[i]->numIO];
		population[i]->data_O = new int[population[i]->numIO];
		new_population[i]->data_O = new int[new_population[i]->numIO];
		for (int a = population[i]->numI-1; a>=0; a--){
			population[i]->connections[a] = new aConnection;
			new_population[i]->connections[a] = new aConnection;
			initConnection(population[i]->connections[a]);
			if (a < population[i]->numIO){
				population[i]->connections[a]->circuit_input = true;
				population[i]->connections[a]->input_index = a;
				population[i]->data_I[a] = a;
			}
		}
		//mark connections belonging to the normal and to the complemented set of gates.
		if (doubled){
			cout<<population[i]->numI<<"  "<<population[i]->numIO<<endl;
			counter = 0;
			for (int a = 0; a < population[i]->numI; a++){
				if (a < population[i]->numIO && (a >= population[i]->numIO/2)){
					population[i]->connections[a]->negated = true;
				} else if (a < (3*segments/2 + population[i]->numIO)){
					if (a%(3*2) > 2){
						population[i]->connections[a]->negated = true;
					}
				} else if (a >= (3*segments/2 + population[i]->numIO)){
					if (counter%(2*2) > 1){
						population[i]->connections[a]->negated = true;
					}
					counter++;
				}
			}
		}
		if (doubled < 1){
			for (int a = 0; a < population[i]->numI; a++){
				if (a > population[i]->numI - population[i]->numIO-1){
					population[i]->connections[a]->circuit_output = true;
					population[i]->data_O[a-population[i]->numI + population[i]->numIO] = (population[i]->numI+population[i]->numO)-(population[i]->numI-a);
				}
			}
		} else {
			counter = population[i]->numIO/2;
			wire_indexes  = population[i]->numIO - 1;
			conn_index = population[i]->numI+population[i]->numO - 1;
			for (int a = population[i]->numI-1; a > 0 ; a--){
				if (counter > 0){
					if (population[i]->connections[a]->output_index == -1 && population[i]->connections[a]->negated)
					{
						cout<<a<<"  "<<conn_index<<"  "<<(wire_indexes)<<"  "<<endl;
						population[i]->connections[a]->circuit_output = true;
	                                        population[i]->data_O[wire_indexes--] = conn_index--;
						counter--;
					}else wire_indexes--;
                                } else break;
			}
/*			for (int a = 0; a < I->numIO; a++){
				I->connections[I->numI-a-1]->unmutable = true;
				I->connections[I->numI-a-1]->output_index =I->data_O[I->numIO-a-1];
			}
*/
			counter = population[i]->numIO/2;
			wire_indexes  = population[i]->numIO - 1;
			for (int a = population[i]->numI-1; a > 0 ; a--){
				if (counter > 0){
					if (population[i]->connections[a]->output_index == -1 && !(population[i]->connections[a]->negated))
					{
						cout<<a<<"  "<<conn_index<<"  "<<(wire_indexes)<<"  "<<endl;
						population[i]->connections[a]->circuit_output = true;
						if (wire_indexes < 0)
							wire_indexes = 0;
	                                        population[i]->data_O[wire_indexes--] = conn_index--;
						counter--;
					}else wire_indexes--;
				} else break;
			}
		}
		wire_indexes = population[i]->numIO;
		conn_index = population[i]->numIO;
		for (int a = 0; a < segments/2; a++){
			population[i]->gates[a] = new aGate;
			new_population[i]->gates[a] = new aGate;
			population[i]->gates[a] = ccGate(gateArray[0], population[i]->gates[a]);
			generateGateIndexes(population[i]->gates[a], wire_indexes);
			for(int j = 0; j < population[i]->gates[a]->numO; j++){
				population[i]->connections[j+conn_index]->input_index = population[i]->gates[a]->coordO[j];
			}
			conn_index += population[i]->gates[a]->numO;
			wire_indexes += (population[i]->gates[a]->numI + population[i]->gates[a]->numO); 
			if (doubled > 0 && (a % 2) != 0){
				population[i]->gates[a]->negated_vars = true;
				new_population[i]->gates[a]->negated_vars = true;
			}
		}
		for (int a = segments/2; a < segments; a++){
			population[i]->gates[a] = new aGate;
			new_population[i]->gates[a] = new aGate;
			population[i]->gates[a] = ccGate(gateArray[1], population[i]->gates[a]);
			generateGateIndexes(population[i]->gates[a], wire_indexes);
			for(int j = 0; j < population[i]->gates[a]->numO; j++){
				population[i]->connections[j+conn_index]->input_index = population[i]->gates[a]->coordO[j];
			}
			conn_index += population[i]->gates[a]->numO;
			wire_indexes += (population[i]->gates[a]->numI + population[i]->gates[a]->numO); 
			if (doubled > 0 && (a % 2) != 0){
				population[i]->gates[a]->negated_vars = true;
				new_population[i]->gates[a]->negated_vars = true;
			}
		}
		makeIndvWiring(population[i]);
		evaluateCircuit(population[i], false);	
		cout<<"gen done"<<endl;
		makeFitness(population[i]);
	}
		
        bestIndv = new aCircuit;
	bestIndv->numIO = finalGate->numI;
	bestIndv->numGates = segments;
	bestIndv->numO = population[0]->numO;
	bestIndv->numI = population[0]->numI;
	bestIndv->I = new int[bestIndv->numI];
	bestIndv->O = new int[bestIndv->numO];
	bestIndv->data_I = new int[bestIndv->numIO];
	bestIndv->data_O = new int[bestIndv->numIO];
	for (int a = bestIndv->numI-1; a>=0; a--){
		bestIndv->connections[a] = new aConnection;
		initConnection(bestIndv->connections[a]);
		if (a < bestIndv->numIO){
			bestIndv->connections[a]->circuit_input = true;
			bestIndv->connections[a]->input_index = a;
			bestIndv->data_I[a] = a;
		}
	}
	for (int a = 0; a < bestIndv->numI; a++){
		if (a > bestIndv->numI - bestIndv->numIO-1){
			bestIndv->connections[a]->circuit_output = true;
			bestIndv->data_O[a-bestIndv->numI + bestIndv->numIO] = (bestIndv->numI+bestIndv->numO)-(bestIndv->numI-a);
		}
	}
	wire_indexes = bestIndv->numIO;
	conn_index = bestIndv->numIO;
	for (int a = 0; a < segments/2; a++){
		bestIndv->gates[a] = new aGate;
		bestIndv->gates[a] = ccGate(gateArray[0], bestIndv->gates[a]);
		generateGateIndexes(bestIndv->gates[a], wire_indexes);
		for(int j = 0; j < bestIndv->gates[a]->numO; j++){
			bestIndv->connections[j+conn_index]->input_index = bestIndv->gates[a]->coordO[j];
		}
		conn_index += bestIndv->gates[a]->numO;
		wire_indexes += (bestIndv->gates[a]->numI + bestIndv->gates[a]->numO);
	}
	for (int a = segments/2; a < segments; a++){
		bestIndv->gates[a] = new aGate;
		bestIndv->gates[a] = ccGate(gateArray[1], bestIndv->gates[a]);
		generateGateIndexes(bestIndv->gates[a], wire_indexes);
		for(int j = 0; j < bestIndv->gates[a]->numO; j++){
			bestIndv->connections[j+conn_index]->input_index = bestIndv->gates[a]->coordO[j];
		}
		conn_index += bestIndv->gates[a]->numO;
		wire_indexes += (bestIndv->gates[a]->numI + bestIndv->gates[a]->numO);
		}
                makeIndvWiring(bestIndv);
                evaluateCircuit(bestIndv, false);
		printCircuit(bestIndv);
                makeFitness(bestIndv);

	initiateStorage();
	cout <<"\n GA Initialization Successful" << endl<<endl<<endl;
	out_stream << "\nInitialization done" << endl;

}


/****************************************
 * Reads the parameters from the open input streams and generates all gates
 * as well as initiates the measurement operator generation if requested
 ****************************************/
void GA::initializeGates() {

	int valuedness = 0;
	int resnum;
	char *representation = new char[4];
	initRep(representation);
	std::string temp;
	char buff_temp[100];
	char binary[100];

	//reading the input file
	cout << "Reading input parameters";
	in_stream >> populationNumber; //int
	in_stream >> segments; //int
	in_stream >> generations; //int
	in_stream >> proba_Mutation; //float
	in_stream >> proba_Crossover; //float
	in_stream >> alpha; //float
	in_stream >> beta; //float
	in_stream >> divider; //int
	in_stream >> displaylevel; //int
	//reading file divider
	in_stream >> temp;
	in_stream >> Ga; //int
	in_stream >> elitism; //int
	in_stream >> mutation; //int
	in_stream >> crossover; //int
	in_stream >> replicator; //int
	in_stream >> fitness; //int
	in_stream >> grouped; //bool
	in_stream >> pareto; //bool
	in_stream >> threshold; //float
	in_stream >> resultnum;
	in_stream >> doubled;
	in_stream >> cascade;
	in_stream >> valuedness;
	//reading file divider
	in_stream >> buff_temp;
	cout<<"... done\n";
	loadGateArray(valuedness, representation);

	cout<<"writting output file ";
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
	out_stream << " Estimated minimum cost of the goal circuit:              "
	<< divider << endl;
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
	out_stream << " Do we use also complemented variables                    "
	<< "\n (doubles the amount of gates):                              "
	<< doubled << endl;
	out_stream << " Do we cascade the control of the gates:                  "
	<< cascade << endl;
	out_stream << " Valuedness of designed logic:                            "
	<< valuedness << endl;


	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << "            The radix of this logic is:                   "
	<< valuedness << endl;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << "            Number of Input Gates:                        "
	<< numofgates << endl;
	for (int a = 0; a < numofgates; a++) {
		out_stream << "Representation: " << gateArray[a]->representation[0]<<gateArray[a]->representation[1]<<gateArray[a]->representation[2]
		<< endl;
		out_stream << " Cost of the gate:  " << gateArray[a]->Cost << endl;
		out_stream << " I/O of the gate:  " << gateArray[a]->numI <<"/"<< gateArray[a]->numO<< endl;
		out_stream << " The name of this gate is: " << gateArray[a]->my_string<< endl;
		for (int j = 0; j < gateArray[a]->length_LU; j++) {
				dec2bin(gateArray[a]->inLU[j], binary,gateArray[a]->numI); 	
				out_stream << binary<<":";
				dec2bin(gateArray[a]->outLU[j], binary,gateArray[a]->numO); 	
				out_stream << binary<<endl;
			}
			out_stream << endl;
	}
	out_stream << " -------- Generating the desired output evaluation ------ "
	<< endl;
	out_stream << " Cost of the gate:  " << finalGate->Cost << endl;
	out_stream << " I/O of the gate:  " << finalGate->numI <<"/"<< finalGate->numO<< endl;
	out_stream << " The name of this gate is: " << finalGate->my_string<< endl;
	for (int j = 0; j < finalGate->length_LU; j++) {
//		cout<<finalGate->inLU[j]<<endl;
		dec2bin(finalGate->inLU[j], binary,finalGate->numI); 	
		out_stream << binary<<":";
		dec2bin(finalGate->outLU[j], binary,finalGate->numO); 	
//		cout<<finalGate->outLU[j]<<endl;
		out_stream << binary<<endl;
	}

	out_stream << endl;


	out_stream << " -------------------------------------------------------- "
	<< endl;
	out_stream << "                   Output of synthesis:                   "
	<< endl;
	out_stream << " -------------------------------------------------------- "
	<< endl;
	cout<<" ..... done"<<endl;
}

/***********************************
 * helper functionnn
 *********************************/
void GA::loadGateArray(int valuedness, char *representation){

	int lines,counter,resnum, inwires, outwires;
	char tempchar;
	char bin[100];
	in_stream >> numofgates;
	cout<<"Reading in initial "<<numofgates<<" gates";;
	for (int a = 0; a < numofgates; a++) {
		gateArray[a] = new aGate;
		in_stream >> gateArray[a]->numI;
		in_stream >> gateArray[a]->numO;
		in_stream >> gateArray[a]->Cost;
		initGate(gateArray[a], gateArray[a]->numI, gateArray[a]->numO, gateArray[a]->Cost, valuedness);
		copyRep(representation, gateArray[a]->representation);
		incRep(representation);
		in_stream >> gateArray[a]->my_string;
		in_stream >> lines;
		gateArray[a]->length_LU = lines;
		for (int j = 0; j < lines; j++) {
			for (int k = 0; k < gateArray[a]->numI; k++) {
				in_stream >> bin[k];
			}
			bin[gateArray[a]->numI] = '\0';
			resnum = bin2dec(bin, gateArray[a]->numI);
			gateArray[a]->inLU[j] = resnum;
			in_stream >> tempchar;
			for (int k = 0; k < gateArray[a]->numO; k++) {
				in_stream >> bin[k];
			}
			bin[gateArray[a]->numO] = '\0';
			int resnum = bin2dec(bin, gateArray[a]->numO);
			gateArray[a]->outLU[j] = resnum;

		}
		counter++;
	}
	cout<<" ..... done"<<endl;
	cout<<"Reading in the desired gates ";

	finalGate = new aGate;
	finalIndividual = new Individual;
	in_stream >> inwires;
	if (doubled > 0)
		finalGate->numI = inwires*2;
	else
		finalGate->numI = inwires;
	in_stream >> outwires;
//	if (doubled > 0)
//		finalGate->numO = outwires*2;
//	else
		finalGate->numO = outwires;
	in_stream >> finalGate->Cost;
	initGate(finalGate, finalGate->numI, finalGate->numO, 1, valuedness);
	in_stream >> lines;
	finalGate->outBinLU = new char*[finalGate->numI];
	finalGate->inBinLU = new char*[finalGate->numI];
	for (int j = 0; j < finalGate->numI; j++){
		finalGate->outBinLU[j] = new char[lines];
		finalGate->inBinLU[j] = new char[lines];
	}
	finalGate->length_LU = lines;
	for (int j = 0; j < lines; j++) {
		for (int k = 0; k < inwires; k++) {
			in_stream >> bin[k];
			finalGate->inBinLU[k][j] = bin[k];
			if (doubled > 0){
				finalGate->inBinLU[k+inwires][j] = (bin[k] == '0') ? '1':'0';
				bin[k+inwires] = (bin[k] == '0') ? '1':'0';
			}
			
		}
		bin[finalGate->numI]= '\0';
//		for (int k = 0; k < finalGate->numI; k++)
//			cout <<"  "<<finalGate->inBinLU[k][j];
//		cout<<endl;
		resnum = bin2dec(bin, finalGate->numI);
//		cout<<"\n  "<<finalGate->numI<< ": "<<resnum;
		finalGate->inLU[j] = resnum;
		in_stream >> tempchar;

		for (int k = 0; k < outwires; k++) {
			in_stream >> bin[k];
			finalGate->outBinLU[k][j] = bin[k];
	//		if (doubled > 0)
	//			finalGate->outBinLU[k+inwires][j] = 1-bin[k];
		}
		bin[finalGate->numO]= '\0';
		resnum = bin2dec(bin, finalGate->numO);
		finalGate->outLU[j] = resnum;
	}
	cout<<" ..... done"<<endl;
}
/***********************************
 * Mutation with structure restrictions
 *********************************/
void GA::doMutation(aCircuit *circuit){
	int con;
	int *columns_match;
	aConnection *acon, *bcon;
	double p;
	int number_m = circuit->numI;
//cout<<"Mutation:";
	if(mutation == 0)
		number_m = 1;
	for (int m = 0; m < number_m; m++){
		p = (float(1.0 * rand() / ((float) RAND_MAX)));
		if (p > proba_Mutation){
			if(mutation == 1){
				do{
					con = rand()%(circuit->numI);
					acon = circuit->connections[con];
					con = rand()%circuit->numI;
					bcon = circuit->connections[con];
				//}while(areFromSameGate(circuit, acon, bcon) || acon->output_index <= bcon->input_index || bcon->unmutable || acon->unmutable);
				}while(areFromSameGate(circuit, acon, bcon) || bcon->output_index > acon->output_index || acon->output_index <= bcon->input_index || bcon->unmutable || acon->unmutable);
//		cout<<acon->input_index<<"-"<<acon->output_index<<"; "<<bcon->input_index<<"-"<<bcon->output_index<<" to ";
				con = acon->output_index;
				acon->output_index = bcon->output_index;
				bcon->output_index = con;
			} else if (mutation == 2){
				con = rand()%(circuit->numI);
                                acon = circuit->connections[con];
				if (con < circuit->numI-1 && circuit->connections[con+1]->input_index == circuit->connections[con]->input_index+1)
					bcon = circuit->connections[con+1];
				else if (con > 0 && circuit->connections[con-1]->input_index == circuit->connections[con]->input_index-1)
					bcon = circuit->connections[con-1];
				con = acon->output_index;
				acon->output_index = bcon->output_index;
				bcon->output_index = con;
			}
//		cout<<acon->input_index<<"-"<<acon->output_index<<"; "<<bcon->input_index<<"-"<<bcon->output_index<<endl;
		}
	}
//cout<<"done"<<endl;
}
/***********************************
 * Crossover with structure restrictions
 *********************************/
void GA::doCrossover(aCircuit *acircuit, aCircuit *bcircuit){
	int acon_counter = 0;
	int bcon_counter = 0;
	int aindexes[acircuit->numI*2], bindexes[bcircuit->numI*2];
	for (int m = 0; m < acircuit->numI; m++){
		if(!acircuit->connections[m]->unmutable)
			aindexes[acon_counter++] = acircuit->connections[m]->output_index;
		if(!bcircuit->connections[m]->unmutable)
			bindexes[bcon_counter++] = bcircuit->connections[m]->output_index;
	}
	for (int m = 0; m < acircuit->numI; m++){
                if(!acircuit->connections[m]->unmutable)
                        acircuit->connections[m]->output_index = bindexes[--bcon_counter];
                if(!bcircuit->connections[m]->unmutable)
                        bcircuit->connections[m]->output_index = aindexes[--acon_counter];
	}
}
/***********************************
 * Evaluates the circuit for given output and also
 * for permutation of the outputs :) !!!
 * Calculates the circuit error
 *********************************/
void GA::evaluateCircuit(aCircuit *circuit, bool disp){
	char binary[100];
	char *bin_table[circuit->numI];
	int *columns_match;
	char *bin_column;
	double error = 0;;
	aGate *gate;
	int con_val, con_idx, resnum, conn_counter, permutations;
	for (int y = 0; y < circuit->numI; y++)
		bin_table[y] = new char[finalGate->length_LU];
	for (int n = 0; n < circuit->numI; n++)
		circuit->connections[n]->current_value = -1;
	for(int k = 0; k < finalGate->length_LU; k++){ 
		dec2bin(finalGate->inLU[k], binary,finalGate->numI); 	
		if (disp) out_stream<<binary<<":";
		for (int m = 0; m < finalGate->numI; m++){
			circuit->connections[m]->current_value = (binary[m] == '0') ? 0:1; 
		}
		for (int m = 0; m < circuit->numGates; m++){
			gate = circuit->gates[m];
			for (int n = 0; n < gate->numI; n++){
				binary[n] = getConnValByOutput(circuit->connections, circuit->numI, gate->coordI[n]);
			}
			binary[gate->numI] = '\0';
		if(disp)	out_stream<<"input to the "<<m<<"th gate "<<binary<<", ";
			resnum = getOutputLUTIndex(gate, bin2dec(binary, gate->numI));
			//resnum = bin2dec(binary, gate->numI);
			dec2bin(gate->outLU[resnum], binary,gate->numO);
		if(disp)	out_stream<<"output of the "<<m<<"th gate "<<binary<<", ";
			for (int n = 0; n < gate->numO; n++){
                                circuit->connections[getConnByInput(circuit->connections, circuit->numI, gate->coordO[n])]->current_value = (binary[n] == '0')?0:1;
                        }
		}
		for (int n = 0; n < circuit->numIO; n++){
			if (circuit->connections[ getConnByOutput(circuit->connections, circuit->numI, circuit->data_O[n])]->current_value == 0)
				binary[n] = '0';
			else binary[n] = '1';
		}
		for (int y = 0; y < circuit->numI; y++)
			bin_table[y][k] = binary[y];
		binary[circuit->numIO] = '\0';
		if (disp) out_stream<<binary<<endl;
		resnum = bin2dec(binary, finalGate->numO);
	}
	if (finalGate->numO != circuit->numO)
		error = compareColumns(bin_table, finalGate->outBinLU,  finalGate->length_LU, circuit->numIO, finalGate->numO);
	else   error = compareColumns(bin_table, finalGate->outBinLU,  finalGate->length_LU, circuit->numIO);
	circuit->error = error;
	circuit->cost = 0;
}
/***********************************
 * helper functionnn
 *********************************/
int GA::fact(int b) {
	if (b < 1 )
		return 0;
	int result = 1;
	for (int a = b; a > 0; a--)
		result *= a;
	return result;
}/***********************************
 * helper functionnn
 *********************************/
int GA::getGate(char *b) {
	for (int a = 0; a < numofgates; a++)
		if (compareRep(gateArray[a]->representation, b))
			return a;
	return '0';
}

/*************************************************
 * helper function
 *************************************************/
bool GA::terminatingCondition() {
        return condition;
}
/*************************************************
 * the way of stopping things
 *************************************************/
void GA::setTerminatingCondition(int counter){
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
                toterr += population[t]->error;
                totcost += population[t]->cost;
                if (population[t]->fitness > population[c]->fitness) {
                        c = t;
                }
        }
        totfit /= populationNumber;
        totcost /= populationNumber;
        toterr /= populationNumber;
        history[generationCondition][0] = float(totfit);
        history[generationCondition][1] = float(population[c]->fitness);
        history[generationCondition][2] = float(toterr);
        history[generationCondition][3] = float(totcost);
        progress = -1;
        cout << "Best current generation fitness: " << population[c]->fitness << endl;
        if (bestIndv->fitness >= 1) {
                progress = 1;
        } else if (population[c]->fitness > bestIndv->fitness) {
                progress = 1;
                cout << "Best current fitness: " << bestIndv->fitness << endl;
                cout << "New best fitness: " << c << " :: " << population[c]->fitness<< endl;
                setIndv(population[c], bestIndv);
		printCircuit(bestIndv);
                bestIndv->fitness = population[c]->fitness;
                cout << "Best new fitness: " << bestIndv->fitness << endl;
        }
//        cout << "Average current fitness: " << totfit << endl;
}
/*************************************************
 * copy one individual to another one
 *************************************************/
void GA::setIndv(aCircuit *indv2, aCircuit *indv1) {
        indv1->fitness = indv2->fitness;
        indv1->error = indv2->error;
        indv1->cost = indv2->cost;
        indv1->numIO = indv2->numIO;
        indv1->numGates = indv2->numGates;
	indv1->numI = indv2->numI;
	indv1->numO = indv2->numO;
	copyIntArray(indv2->I, indv1->I, indv2->numI);	
	copyIntArray(indv2->O, indv1->O, indv2->numO);	
	copyIntArray(indv2->data_I, indv1->data_I, indv2->numIO);	
	copyIntArray(indv2->data_O, indv1->data_O, indv2->numIO);	
	copyConnections(indv2->connections, indv1->connections, indv2->numI);	
}
/*******************************************
 * Initiate the length of the arrays for history records
 *******************************************/
void GA::initiateStorage() {
        cout << "Initiating storage for " << generations << " generations";
        history = new float[generations][4];
        cout<<" ..... done"<<endl;
}

/*******************************************
 * Replication like function
 *******************************************/
int GA::applyReplication(){
        float fitnessSum = 0;
        float fitnessavg = 0;
        float first, second, tempsum, proba, best;
        int indv1, indv2;
        for (int i = 0; i < populationNumber; i++)
                fitnessSum += population[i]->fitness;
        fitnessavg /= populationNumber;
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
                                                * fitnessSum) / 2;
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
                        setIndv(population[indv1], new_population[counter]);
                        counter++;
                } else {
                        setIndv(population[indv1], new_population[counter]);
                        setIndv(population[indv2], new_population[counter + 1]);
			proba = (1 - (float(1.0 * rand() / ((float) RAND_MAX + 1.0))));
			if (proba < proba_Crossover)
				doCrossover(new_population[counter],new_population[counter + 1]);
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
                        setIndv(new_population[i], population[i+1]);
                }

        } else {
                for (int i = 0; i < populationNumber; i++) {
                        setIndv(new_population[i], population[i]);
                }
        }


	return 0;
}

/*******************************************
 * Print out an individual
 *******************************************/
void GA::printCircuit(aCircuit *circuit){

	out_stream<<"-----"<<endl;	
	out_stream<<"fitness: "<<circuit->fitness<<endl;
	out_stream<<"error: "<<circuit->error<<endl;
	out_stream<<"cost: "<<circuit->cost<<endl;
	out_stream<<"Wires: :"<<endl;
	for (int m = 0; m < circuit->numI; m++){
		out_stream<<m<<": "<<circuit->connections[m]->input_index<<"->"<<circuit->connections[m]->output_index<<endl;
	}
	out_stream<<"The functional output:"<<endl;
	evaluateCircuit(circuit, true);
	out_stream<<"-----"<<endl;	
}


/*******************************************
 * Integer encoded-Binary value (input) to int (output LUT index) conversion
 *******************************************/
int GA::getOutputLUTIndex(aGate *gate, int binary){
	for (int m = 0; m < gate->length_LU; m++)
		if (gate->inLU[m] == binary)
			return m;
	return -1;
}
