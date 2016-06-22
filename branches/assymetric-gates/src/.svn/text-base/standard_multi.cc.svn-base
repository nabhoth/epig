#include "EpiG.h"

/****************************************
 * returns the matrix representing the circuit of a given individual
 ****************************************/
#ifdef __STDR__
qGate* GA::computeMatrix(Individual *indi) {

	qGate *temp0, *temp1, *temp2, *myfinalG;
	int phasecounter0 = 0;
	int Time = 0;
	int cost = 0;
	int csegs = 0;
	float numIO, val, maxcount;
	char *gate_rep = new char[4];
	cuComplex temphase;
	temp0 = new qGate;
	temp1 = new qGate;
	temp2 = new qGate;
	myfinalG = new qGate;

	int w = (int(pow((float) indi->valuedness, (float) indi->ioNumber)));
	int N = (int(pow(pow((float) indi->valuedness, (float) indi->ioNumber),2)));
	cuComplex alpha = make_cuFloatComplex(1,0);
	cuComplex beta = make_cuFloatComplex(0,0);

	initGate(temp0, indi->ioNumber, indi->valuedness, 1);	
	initGate(temp1, indi->ioNumber, indi->valuedness, 1);	
	initGate(temp2, indi->ioNumber, indi->valuedness, 1);	
	initGate(myfinalG, indi->ioNumber, indi->valuedness, 1);	
	//int the string counters
	int begin = 0;
	int end = indi->my_string.find("p", begin + 1);
	indi->Cost = 0;
	indi->segmentNumber = 0;

	while (end > 0) {
		if ((end - begin) > 1) {
			csegs++;
			//get the first  gate of the Segment and copy it
				getGateRepfromStr(begin+1,indi->my_string, gate_rep);
			for (int a = 0; a <= numofgates; a++){
				if (compareRep(indi->gateArray[a]->representation, gate_rep)) {
					temp0 = ccGate(indi->gateArray[a], temp0);
					numIO = (float) temp0->numIO;
					val = (float) temp0->valuedness;
					if (phase == 1) {
						temphase = indi->phases[phasecounter0++];
						maxcount = (int(pow(pow(val, numIO),2)));
						for (int k = 0; k < maxcount; k++)
								temp0->gateMatrix1[k] = cuCmulf(
										temp0->gateMatrix1[k], temphase);
					}
					//continue with the cost
					begin++;
					cost += indi->gateArray[a]->Cost;
					break;


				}
			}
			//get the next one
			begin += 3;
			for (int b = begin; b < end; b+=REPSIZE) {
				getGateRepfromStr(b,indi->my_string, gate_rep);
				for (int a = 0; a <= numofgates; a++)
					if (compareRep(indi->gateArray[a]->representation, gate_rep)) {
						temp1 = ccGate(indi->gateArray[a], temp1);
						temphase = indi->phases[phasecounter0++];
						val = (float) temp1->valuedness;
						numIO = (float) temp1->numIO;

						//cout<<"apply the phase2"<<endl;
						if (phase == 1) {
							maxcount = (int(pow(pow(val, numIO),2)));
							for (int k = 0; k < maxcount; k++)
									temp1->gateMatrix1[k] = cuCmulf(
											temp1->gateMatrix1[k], temphase);
						}
						//continue with the cost
						cost += indi->gateArray[a]->Cost;
						a = numofgates * 2;
						tensorProduct3(temp2, temp0, temp1);
						temp0 = ccGate(temp2, temp0);
					}
			}
			if (Time == 0) {
				//if this is the end of the first segment fill the output value
				myfinalG = ccGate(temp0, myfinalG);
				Time++;
			} else {
				//compile warning is intended
				//multiply using normal Manumofgatestrix product the twou segments
				temp2->valuedness = temp0->valuedness;
				zeroMatrix(temp2, temp0->numIO);
				cblasMatrixProduct(temp0, myfinalG, temp2);
				myfinalG = ccGate(temp2, myfinalG);
				Time++;
			}
			begin = indi->my_string.find("p", end);
			if (end != -1 && indi->my_string.at(end) == 'p' && end < indi->my_string.length() - 3)
				end = end - 1;
		} else {
		//move to the next segment
			begin = end+1;
			end = indi->my_string.find("p", begin+1);
			if (end != -1 && indi->my_string.at(end) == 'p' && end < indi->my_string.length() - 3){
				end = end - 1;
				if (indi->my_string.at(end-1) == 'p')
					end = end - 1;
			}
		}
	}

	//check if phase is defined and multiply the whole matrix by it
	if (phase == 1) {
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			val = (float) myfinalG->valuedness;
			numIO = (float) myfinalG->numIO;
			maxcount = (int(pow(pow(val, numIO),2)));
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(
							myfinalG->gateMatrix1[k], indi->phase);
		}
	}
	//actualize the segment number
	indi->segmentNumber = csegs;
	indi->Cost = cost;
	val = ((float) myfinalG->valuedness);
	numIO = ((float) myfinalG->numIO);

	destroyGate(temp0);
	delete temp0;
	destroyGate(temp1);
	delete temp1;
	destroyGate(temp2);
	delete temp2;
	return myfinalG;
}

#else
qGate* GA::computeMatrix(Individual *indi) {



	qGate *temp0, *temp1, *temp2, *myfinalG;
	int phasecounter0 = 0;
	int Time = 0;
	int cost = 0;
	int csegs = 0;
	float numIO, val, maxcount;
	char *gate_rep = new char[4];
	cuComplex temphase;
	temp0 = new qGate;
	temp1 = new qGate;
	temp2 = new qGate;
	myfinalG = new qGate;

	int w = (int(pow((float) indi->valuedness, (float) indi->ioNumber)));
	int N = (int(pow(pow((float) indi->valuedness, (float) indi->ioNumber),2)));
	cuComplex alpha = make_cuFloatComplex(1,0);
	cuComplex beta = make_cuFloatComplex(0,0);

	initGate(temp0, indi->ioNumber, indi->valuedness, 1);	
	initGate(temp1, indi->ioNumber, indi->valuedness, 1);	
	initGate(temp2, indi->ioNumber, indi->valuedness, 1);	
	initGate(myfinalG, indi->ioNumber, indi->valuedness, 1);	
	//int the string counters
	int begin = 0;
	int end = indi->my_string.find("p", begin + 1);
	indi->Cost = 0;
	indi->segmentNumber = 0;

	while (end > 0) {
		if ((end - begin) > 1) {
			csegs++;
			//get the first  gate of the Segment and copy it
				getGateRepfromStr(begin+1,indi->my_string, gate_rep);
			for (int a = 0; a <= numofgates; a++){
				if (compareRep(gateArray[a]->representation, gate_rep)) {
					temp0 = ccGate(gateArray[a], temp0);
					numIO = (float) temp0->numIO;
					val = (float) temp0->valuedness;
					if (phase == 1) {
						temphase = indi->phases[phasecounter0++];
						maxcount = (int(pow(pow(val, numIO),2)));
						for (int k = 0; k < maxcount; k++)
								temp0->gateMatrix1[k] = cuCmulf(
										temp0->gateMatrix1[k], temphase);
					}
					//continue with the cost
					begin++;
					cost += gateArray[a]->Cost;
					break;


				}
			}
			//get the next one
			begin += 3;
			for (int b = begin; b < end; b+=REPSIZE) {
				getGateRepfromStr(b,indi->my_string, gate_rep);
				for (int a = 0; a <= numofgates; a++)
					if (compareRep(gateArray[a]->representation, gate_rep)) {
						temp1 = ccGate(gateArray[a], temp1);
						temphase = indi->phases[phasecounter0++];
						val = (float) temp1->valuedness;
						numIO = (float) temp1->numIO;

						//cout<<"apply the phase2"<<endl;
						if (phase == 1) {
							maxcount = (int(pow(pow(val, numIO),2)));
							for (int k = 0; k < maxcount; k++)
									temp1->gateMatrix1[k] = cuCmulf(
											temp1->gateMatrix1[k], temphase);
						}
						//continue with the cost
						cost += gateArray[a]->Cost;
						a = numofgates * 2;
						tensorProduct3(temp2, temp0, temp1);
						temp0 = ccGate(temp2, temp0);
					}
			}
			if (Time == 0) {
				//if this is the end of the first segment fill the output value
				myfinalG = ccGate(temp0, myfinalG);
				Time++;
			} else {
				//compile warning is intended
				//multiply using normal Manumofgatestrix product the twou segments
				temp2->valuedness = temp0->valuedness;
				zeroMatrix(temp2, temp0->numIO);
				cblasMatrixProduct(temp0, myfinalG, temp2);
				myfinalG = ccGate(temp2, myfinalG);
				Time++;
			}
			begin = indi->my_string.find("p", end);
			if (end != -1 && indi->my_string.at(end) == 'p' && end < indi->my_string.length() - 3)
				end = end - 1;
		} else {
		//move to the next segment
			begin = end+1;
			end = indi->my_string.find("p", begin+1);
			if (end != -1 && indi->my_string.at(end) == 'p' && end < indi->my_string.length() - 3){
				end = end - 1;
				if (indi->my_string.at(end-1) == 'p')
					end = end - 1;
			}
		}
	}

	//check if phase is defined and multiply the whole matrix by it
	if (phase == 1) {
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			val = (float) myfinalG->valuedness;
			numIO = (float) myfinalG->numIO;
			maxcount = (int(pow(pow(val, numIO),2)));
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(
							myfinalG->gateMatrix1[k], indi->phase);
		}
	}
	//actualize the segment number
	indi->segmentNumber = csegs;
	indi->Cost = cost;
	val = ((float) myfinalG->valuedness);
	numIO = ((float) myfinalG->numIO);

	destroyGate(temp0);
	delete temp0;
	destroyGate(temp1);
	delete temp1;
	destroyGate(temp2);
	delete temp2;
	return myfinalG;

}

#endif

/****************************************
 * calculates theprobabilities for obtaining |0> and |1>
 * Right now it is doing measurement in increasing order
 * of all defined bits to measure
 ****************************************/
void* GA::doMeasureFitness(Individual *ind, bool display) {

	//define temporal variables for the results and the manupaltion
	//of the measurement
	qGate *measures[10], *myfinalG;
	int p = 0;
	float val = (float) finalGate->valuedness;
	float numIO = (float) finalGate->numIO;
	float err,maxcount = 0;
	int l = (int(pow(val, numIO)));
	int x = (int(pow(val, numIO)));
	;
	int y = (int(pow(val, numIO)));
	;
	cuComplex inter[(int(pow(val, numIO)))];
	cuComplex logicresult[(int(pow(val, numIO)))];
	cuComplex sta[(int(pow(val, numIO)))];
	cuComplex expectations[(int(pow(val, numIO) * val))], alphas, betas;
	cuComplex resultingTotalError[(int(pow(val, numIO)))][measurement
	                                                      * (int) val];
	cuComplex resultingTotal[(int(pow(val, numIO)))][measurement * (int) val];

	Individual *indi = ind;
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;

	//threaded code
	int rc = pthread_mutex_lock(&mtx);
	if (rc) {
		cout << "Matrix locked: " << ind << endl;
	}
	//myfinalG = computeMatrix(indi);
	//cout << indi->my_string<<endl;
	long time = clock();
	myfinalG = computeMatrix(indi);
	rc = pthread_mutex_unlock(&mtx);
	numIO = (float) myfinalG->numIO;
	//check if phase is defined and multiply the whole matrix by it
	if (phase == 1) {
		maxcount = (int(pow(pow(val, numIO),2)));
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(
							myfinalG->gateMatrix1[k], indi->phase);
		}
	}


	//for single qubit this should always be only 0, i.e. <1
	for (int m = 0; m < measurement; m++) {
		for (int r = 0; r < val; r++) {
			measures[r] = measurements[(measurementQBits[m] * (int) val) + r];
		}
		//for all input states l
		for (int k = 0; k < measured_desired_output_records; k++) {

			for (int w = 0; w < l; w++)
				sta[w] = inter[w] = logicresult[w] = make_cuFloatComplex(0, 0);
			//set the k-th input state - minterm
			sta[measured_desired_output_records_idx[k]] = make_cuFloatComplex(1, 0);
			alphas = make_cuFloatComplex(0, 0);
			betas = make_cuFloatComplex(0, 0);

			//propagate the state through the matrix
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < l; j++) {
					logicresult[i] = cuCaddf(logicresult[i], cuCmulf(sta[j],
							myfinalG->gateMatrix1[i*l+j]));
				}
			}
			//measure for each state of the single given qubit
			for (int r = 0; r < val; r++) {
				for (int w = 0; w < l; w++)
					inter[w] = make_cuFloatComplex(0, 0);
				for (int i = 0; i < l; i++) {
					//rows
					for (int j = 0; j < l; j++) {
						inter[i] = cuCaddf(inter[i], cuCmulf(logicresult[j],
								measures[r]->gateMatrix1[i*l+j]));
					}
				}

/*				cout << measured_desired_output_records_idx[k]<<": ";
				for (int i = 0; i < l; i++) 
				cout <<"O:( "<< cuCrealf( inter[i])<<", "<<cuCimagf( inter[i])<<")";
				cout<<"<>";
*/
				//the finnal inner product for p(r)
				expectations[((m * (int) val) + r)] = make_cuFloatComplex(0, 0);
				for (int j = 0; j < l; j++) {
					expectations[(m * (int) val + r)] = cuCaddf(expectations[(m
							* (int) val + r)], cuCmulf(cuConjf(logicresult[j]),
									inter[j]));
				}

				resultingTotal[measured_desired_output_records_idx[k]][(m * (int) val + r)] = expectations[(m
						* (int) val + r)];
			//for (int i = 0; i < l; i++) 
//				cout <<"R:( "<< cuCrealf( resultingTotal[measured_desired_output_records_idx[k]][(m * (int) val + r)])<<", "<<cuCimagf( resultingTotal[measured_desired_output_records_idx[k]][(m * (int) val + r)])<<"), "<<"E:( "<< cuCrealf( measureexpected[measured_desired_output_records_idx[k]][(m * (int) val + r)])<<", "<<cuCimagf( measureexpected[measured_desired_output_records_idx[k]][(m * (int) val + r)])<<"); ";
//				cout<<"<>";
				
				if (cuCrealf(measureexpected[measured_desired_output_records_idx[k]][(int) val * m + r]) == 0
						&& cuCimagf(measureexpected[measured_desired_output_records_idx[k]][(int) val * m + r]) == 1) {
					//skip don't cares
					resultingTotalError[measured_desired_output_records_idx[k]][m * (int) val + r]
					                       = make_cuFloatComplex(0, 0);
				} else {
					resultingTotalError[measured_desired_output_records_idx[k]][m * (int) val + r]
					                       = cuFloatComplex_sqrt(cuFloatComplex_pow((cuCsubf(resultingTotal[measured_desired_output_records_idx[k]][(m * (int) val) + r],measureexpected[measured_desired_output_records_idx[k]][((int) val * m) + r])), make_cuFloatComplex(2,0)));
				}
				//for (int i = 0; i < l; i++) 
//			cout <<"RE:( "<< cuCrealf( resultingTotalError[measured_desired_output_records_idx[k]][(m * (int) val + r)])<<", "<<cuCimagf( resultingTotalError[measured_desired_output_records_idx[k]][(m * (int) val + r)])<<")";

			}
//			cout<<endl;
		}
	}
//				cout<<endl;

	int m;
	err = 0;
	for (int c = 0; c < measurement; c++) {
		m = l * (int) val;
		inter[0] = make_cuFloatComplex(0, 0);
		for (int e = 0; e < measured_desired_output_records; e++) {
			for (int r = 0; r < val; r++) {
				if (!(cuCrealf(measureexpected[measured_desired_output_records_idx[e]][(int) val * c + r]) == 0
						&& cuCimagf(measureexpected[measured_desired_output_records_idx[e]][(int) val * c + r])	== -1))
					inter[0] = cuCaddf(inter[0], make_cuFloatComplex(cuCrealf(
							resultingTotalError[measured_desired_output_records_idx[e]][c * (int) val + r]), 0));
				else
					m -= (int) val;
			}
		}
		err += cuCrealf(cuFloatComplex_div(inter[0], make_cuFloatComplex(m, 0)));
	}
//	cout<<err<<endl;
	err /= measurement;
	indi->Error = err;

	//generate fitness value
	makeFitness(indi);	
//	cout<<indi->fitness<<endl;
	if (display) {
		out_stream << "Error: " << indi->Error << endl;
		out_stream << "Fitness: " << indi->fitness << endl;
		out_stream << "Valuedness: " << indi->valuedness << endl;
		out_stream << endl;
		out_stream << "Input St. Des_0 Des_1 Out_0 Out_1 Error_0 Error_1 "
		<< endl;

		for (int m = 0; m < measurement; m++) {
			//for all input states l
			for (int k = 0; k < measured_desired_output_records; k++) {
				out_stream << "  " << measured_desired_output_records_idx[k] << "   ";
				for (int r = 0; r < val; r++) {
					out_stream << "  (" << cuCrealf(
							measureexpected[measured_desired_output_records_idx[k]][(int) val * m + r]) << ","
							<< cuCimagf(measureexpected[measured_desired_output_records_idx[k]][(int) val * m + r])
							<< ")";
				}
				//measure for each state of the single given qubit
				for (int r = 0; r < val; r++) {
					out_stream << "  (" << cuCrealf(resultingTotal[measured_desired_output_records_idx[k]][(m
							* (int) val + r)]) << "," << cuCimagf(
									resultingTotal[measured_desired_output_records_idx[k]][(m * (int) val + r)]) << ")";
				}
				for (int r = 0; r < val; r++) {
					out_stream << "  (" << cuCrealf(resultingTotalError[measured_desired_output_records_idx[k]][m
					                                                       * (int) val + r]) << "," << cuCimagf(
					                                                    		   resultingTotalError[measured_desired_output_records_idx[k]][m * (int) val + r]) << ")";
				}
				out_stream << endl;
			}
		}
		out_stream << " -------- " << endl;
	}

	time = clock() -time ;
//	cout<<"CUBLAS fitness done: "<<(float)(time/CLOCKS_PER_SEC)<<", for "<<indi->ioNumber<<" qubits "<<endl;
	destroyGate(myfinalG);
	delete (myfinalG);
}
/****************************************
* Generates a set of states representing a given sequence: the input qubit is the one before the last from the bottom
* the output qubit is the bottom one
* both input and output qubits must be measured after each operation in order to be able to initialize them properly
* for the next step.
* A sequence is detected when for all input elements the output is |0> and for the last input value the output is |1>
 ****************************************/

void* GA::doMeasureFASeqFitness(Individual *ind, bool display) {

	//define temporal variables for the results and the manupaltion
	//of the measurement
	qGate *myfinalG;
	float val = (float) finalGate->valuedness;
	float numIO = (float) finalGate->numIO;
	float err,maxcount,alphas = 0;
	int temp1;
	int l = (int(pow(val, numIO)));
	int x = (int(pow(val, numIO)));
	int y = (int(pow(val, numIO)));
	cuFloatComplex inter[(int(pow(val, numIO)))];
	cuFloatComplex inter_b[(int(pow(val, numIO)))];
	cuFloatComplex logicresult[(int(pow(val, numIO)))];
	cuFloatComplex sta[(int(pow(val, numIO)))];
	cuFloatComplex stb[(int(pow(val, numIO)))];
	cuFloatComplex expectations[(int(pow(val, numIO) * val))];
	cuFloatComplex resultingTotalError[(int(pow(val, numIO)))][measurement * (int) val];
	cuFloatComplex resultingTotal[(int(pow(val, numIO)))][measurement * (int) val];
	cuFloatComplex zero = make_cuFloatComplex(0,0);
	cuFloatComplex one = make_cuFloatComplex(1,0);
	bool set;

	Individual *indi = ind;
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;

	//threaded code
	long time = clock();
	myfinalG = computeMatrix(indi);
	numIO = (float) myfinalG->numIO;
	//check if phase is defined and multiply the whole matrix by it
	if (phase == 1) {
		maxcount = (int(pow(pow(val, numIO),2)));
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(myfinalG->gateMatrix1[k], indi->phase);
		}
	}

if (display){
		out_stream<<indi->my_string<<endl;
}

	//for all input states l initiate the vector at time 0
	vectorInitZero(l, sta);

	for (int k = 0; k < measured_desired_output_records; k++) {


		//empty all intermediary variables
		for (int w = 0; w < l; w++)
			inter[w] = logicresult[w] = stb[w] = make_cuFloatComplex(0,0);
		//set the k-th input from the sequence to desired value and the output to 0
		if (k == 0){
			if (cuCrealf(sequence_desired[k][0]) == 0)
				sta [0] = make_cuFloatComplex(1, 0);
			else
				sta [2] = make_cuFloatComplex(1, 0);
		} else {
			set = false;
			if (cuCrealf(sequence_desired[k][0]) == 0){
				int divider = l/(int(pow(val, numIO-2)));
				for (int w = 0; w < l; w++){
					if (cuCrealf(sta[w]) != 0 || cuCimagf(sta[w]) != 0){
						temp1 = w / divider;
						stb[temp1*divider] = sta[w];
						set = true;
					}
				}
			} else {
				int divider = l/(int(pow(val, numIO-2)));
				for (int w = 0; w < l; w++){
					if (cuCrealf(sta[w]) != 0 || cuCimagf(sta[w]) != 0){
						temp1 = w / divider;
						stb[(temp1*divider)+2] = sta[w];
						set = true;
					}
				}
			}
//	cout<<"Set: "<<set<<endl;
			if (!set)
				if (cuCrealf(sequence_desired[k][0]) == 0){
					stb[0] = make_cuFloatComplex(1, 0);
				} else {
					stb[2] = make_cuFloatComplex(1, 0);
				}
		for (int w = 0; w < l; w++)
			sta[w] = stb[w];
		}

		alphas = 0.0;

		for (int w = 0; w < l; w++)
			inter[w] = inter_b[w] = make_cuFloatComplex(0,0);
if (display){
		cout << "I: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(sta[w])<<","<<cuCimagf(sta[w])<<")";
		cout<<endl;
}
		
		cblas_cgemv(CblasRowMajor, CblasNoTrans, l, l, &one, myfinalG->gateMatrix1, l, sta, 1, &zero, logicresult, 1);

if (display){
		cout << "O: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(logicresult[w])<<","<<cuCimagf(logicresult[w])<<")";
		cout<<endl;
}
		//measure for only the desired state of the single output qubit
		//measure for 0 when non final sequence element
		if (k != measured_desired_output_records-1){
			cblas_cgemv(CblasRowMajor, CblasNoTrans, l, l, &one, measures_o[0]->gateMatrix1, l, logicresult, 1, &zero, inter, 1);
		}else{
			cblas_cgemv(CblasRowMajor, CblasNoTrans, l, l, &one, measures_o[1]->gateMatrix1, l, logicresult, 1, &zero, inter, 1);
		}


if (display){
		cout << "Measure Output: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(inter[w])<<","<<cuCimagf(inter[w])<<")";
		cout<<endl;
}

		//the finnal inner product for p(r) - generate the expectation for the measured output
		cblas_cdotc_sub(l, inter, 1, inter, 1, &resultingTotal[k][0]);
		//skip don't cares
		if (cuCrealf(sequence_desired[k][1]) == 0 && cuCimagf(sequence_desired[k][1]) == -1) {
			resultingTotalError[k][(int) val]     = make_cuFloatComplex(0, 0);
		} else {
			//calculate the error if this is not to be considered as don't care
			resultingTotalError[k][(int) val] = cuFloatComplex_sqrt(cuFloatComplex_pow((cuFloatComplex_sub(resultingTotal[k][0],sequence_desired[k][1])),make_cuFloatComplex(2,0)));
		}
		expectations[0] = make_cuFloatComplex(0,0);
		cblas_zdotc_sub(l, inter, 1, inter, 1, &expectations);
//		cout<<cout<<"("<<cuCrealf(expectations[0])<<","<<cuCimagf(expectations[0])<<")";
		alphas = sqrt(cuCrealf(expectations[0]));
//		cout<<cout<<"("<<(alphas)<<")";
		//alphas = make_cuFloatComplex(real(temp_val), imag(temp_val));

		//if the measurement generated at least some output
		if (!( alphas == 0)){
			//if the output qubit was measured in 0, normalize the vector
			vectorScalarMult(l, cuFloatComplex_div(make_cuFloatComplex(1,0), make_cuFloatComplex(alphas,0)), inter, inter);
		} else {
			//if the measurement was not succesful, use the output of the circuit 
/*			for (int j = 0; j < l; j++) {
				inter[j] = logicresult[j];
			}
*/		}

/*		cout<<cout<<"("<<(alphas)<<")";
                cout<<endl;
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(inter[w])<<","<<cuCimagf(inter[w])<<")";
		cout<<endl;
*/

		//measure the input qubit for 0 - to reset it and allow consequent initialization
		if (k < measured_desired_output_records - 1)
			if (cuCrealf(sequence_desired[k+1][0]) == 0){
				cblas_cgemv(CblasRowMajor, CblasNoTrans, l, l, &one, measures_i[0]->gateMatrix1, l, inter, 1, &zero, inter_b, 1);
			} else {
				cblas_cgemv(CblasRowMajor, CblasNoTrans, l, l, &one, measures_i[1]->gateMatrix1, l, inter, 1, &zero, inter_b, 1);
			}	

if (display){
		cout << "Measure Input: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(inter_b[w])<<","<<cuCimagf(inter_b[w])<<")";
		cout<<endl;
}
		//calculate the expectation value for measuring the input qubit
		expectations[0] = make_cuFloatComplex(0,0);;
		cblas_zdotc_sub(l, inter_b, 1, inter_b, 1, &expectations);
		alphas = sqrt(cuCrealf(expectations[0]));

		//if the input qubit was measured in 0, normalize the vector
		if (!(( alphas) == 0)){

			vectorScalarMult(l, cuFloatComplex_div(make_cuFloatComplex(1,0), make_cuFloatComplex(alphas,0)), inter_b, inter_b);
			//copy the current state vector as the next input
			for (int i = 0; i < l; i++){ 
				sta[i] = make_cuFloatComplex(cuCrealf( inter_b[i]),cuCimagf( inter_b[i]));
			}

		} else {
			//copy the previous intermediary state vector as the next input without normalization
			for (int i = 0; i < l; i++){ 
				sta[i] = make_cuFloatComplex(cuCrealf( inter[i]),cuCimagf( inter[i]));
			}
		}
if (display){
		cout << "At the end: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(sta[w])<<","<<cuCimagf(sta[w])<<")";
		cout<<endl;
}
	}
//	exit(0);

	int m;
	err = 0;
	//calculate the overall error
	for (int c = 0; c < measurement; c++) {
		m = measured_desired_output_records;
		inter[0] = make_cuFloatComplex(0,0);
		for (int e = 0; e < measured_desired_output_records; e++) {
			//inter[0] = cuCmulf(inter[0], resultingTotal[e][0]);
			inter[0] = cuCaddf(inter[0], resultingTotalError[e][(int) val]);
		}
		err += cuCrealf(cuFloatComplex_div(inter[0], make_cuFloatComplex(m, 0)));
	}
	err /= measurement;
	indi->Error = err;

	//generate fitness value
	makeFitness(indi);	
//	cout<<indi->fitness<<endl;
	if (display) {

		out_stream<< "Circuit Matrix:"<<endl;
		for (int i = 0; i < l; i++) {
			for (int j = 0; j < l; j++) {
				out_stream<<"("<<cuCrealf(myfinalG->gateMatrix1[i*l+j])<<","<<cuCimagf(myfinalG->gateMatrix1[i*l+j])<<")";
			}
			out_stream<<endl;
		}
		out_stream<<endl;


		out_stream << "Error: " << indi->Error << endl;
		out_stream << "Fitness: " << indi->fitness << endl;
		out_stream << "Valuedness: " << indi->valuedness << endl;
		out_stream << endl;
		out_stream << "Sequence length: "<<measured_desired_output_records<<endl;;
		out_stream << "Sequence: ";

		for (int m = 0; m < measured_desired_output_records; m++) 
			out_stream<<"  "<<cuCrealf(sequence_desired[m][0])<<"  ";
		
		out_stream << endl<<" Desired: ";
		for (int m = 0; m < measured_desired_output_records; m++) 
			out_stream<<"("<<cuCrealf(sequence_desired[m][1])<<","<<cuCimagf(sequence_desired[m][1])<<")";


		out_stream << endl<<"Observables: ";
		for (int m = 0; m < measured_desired_output_records-1; m++) 
			out_stream<<"  0  ";
		out_stream<<"  1  ";

		out_stream << endl<<"Obtained: ";
		for (int m = 0; m < measured_desired_output_records; m++) 
			out_stream<<"("<<cuCrealf(resultingTotal[m][0])<<","<<cuCimagf(resultingTotal[m][0])<<")";

		out_stream << endl <<" -------- " << endl;
	}

	time = clock() -time ;
	//cout<<"CUBLAS fitness done: "<<(float)(time/CLOCKS_PER_SEC)<<", for "<<indi->ioNumber<<" qubits "<<endl;
	destroyGate(myfinalG);
	delete (myfinalG);
	return 0;
}
/****************************************
 * calculates the probabilities for obtaining two distinct multi qubit states of the form
 * (p_0 + p_1 + .. + p_k)|x_0...x_k> and (r_0 + r_1 + .. + r_k)|y_0....y_k>
 ****************************************/
void* GA::doMultiMeasureFitness(Individual *ind, bool display) {

	//define temporal variables for the results and the manipulation
	//of the measurement
	int l = (int(pow((float) 2, (float) finalGate->numIO)));
	int x = (int(pow((float) 2, (float) finalGate->numIO)));
	int y = (int(pow((float) 2, (float) finalGate->numIO)));
	int p,maxcount,val = 0;
	int mes = 0;
	cuComplex inter0[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex inter1[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex logicresult[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex sta[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex expectation0, expectation1, alphas, betas;
	cuComplex expe_0, expe_1;
	cuComplex
	resultingTotalError[(int(pow((float) 2, (float) finalGate->numIO)))][measurement
	                                                                     * 2];
	cuComplex
	resultingTotal[(int(pow((float) 2, (float) finalGate->numIO)))][measurement
	                                                                * 2];
	cuComplex inter, expect0, expect1;
	qGate *measure0, *measure1, *myfinalG, measure;

	Individual *indi = ind;
	//null this Individual previous evaluation
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;

	//threaded code
	int rc = pthread_mutex_lock(&mtx);
	if (rc) {
		cout << "Matrix locked: " << ind << endl;
	}
	// compute hte matrix of the circuit
	//myfinalG = computeMatrix(indi);
	myfinalG = computeMatrix(indi);
	rc = pthread_mutex_unlock(&mtx);

	if (phase == 1) {
		maxcount = (int(pow(pow(myfinalG->valuedness, myfinalG->numIO),2)));
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(
							myfinalG->gateMatrix1[k], indi->phase);
		}
	}

	//propagate each input through the circuit
	for (int k = 0; k < l; k++) {

		for (int i = 0; i < l; i++) {
			inter0[i] = make_cuFloatComplex(0, 0);
			inter1[i] = make_cuFloatComplex(0, 0);
		}

		//initialize variables
		for (int w = 0; w < l; w++)
			sta[w] = inter0[w] = inter1[w] = logicresult[w]
			                                             = make_cuFloatComplex(0, 0);
		//set the orthonormal state - k
		sta[k] = make_cuFloatComplex(1, 0);

		//propagate the state through the matrix
		for (int i = 0; i < l; i++) {
			for (int j = 0; j < l; j++) {
				logicresult[i] = cuCaddf(logicresult[i], cuCmulf(sta[j],
						myfinalG->gateMatrix1[i*l+j]));
			}
		}

		//init the measurement operator with respect to the desired output
		//for each measured ouptut state get the expectation values
		expect0 = expectationsAllState[k][0];
		expect1 = expectationsAllState[k][1];

		//apply the measurement operator for the desired state
		for (int i = 0; i < l; i++) {
			//rows
			for (int j = 0; j < l; j++) {
				inter0[i] = cuCaddf(inter0[i], cuCmulf(logicresult[j],
						measurementsAllState[2* k ]->gateMatrix1[i*l+j]));
				inter1[i] = cuCaddf(inter1[i], cuCmulf(logicresult[j],
						measurementsAllState[2* k + 1]->gateMatrix1[i*l+j]));
			}
		}
		//p(0) - not really p(0) but rather the desired result ---------------------------------
		expectation0 = make_cuFloatComplex(0, 0);
		alphas = make_cuFloatComplex(0, 0);
		//the finnal inner product for p(0)
		for (int j = 0; j < l; j++)
			expectation0 = cuCaddf(expectation0, cuCmulf(
					cuConjf(logicresult[j]), inter0[j]));
		//state after the measurement for 0
		expe_0
		= make_cuFloatComplex (cuCrealf(expectation0),
				cuCimagf(expectation0));
		expe_0 = cuFloatComplex_sqrt(expe_0);
		expectation0 = make_cuFloatComplex(cuCrealf(expe_0), cuCimagf(expe_0));

		for (int i = 0; i < l; i++) {
			inter0[i] = cuFloatComplex_div(inter0[i], expectation0);
			alphas = cuCaddf(alphas, inter0[i]);
		}
		//p(1) ---------------------------------------------------------------------------------
		//vector inner product
		expectation1 = make_cuFloatComplex(0, 0);
		betas = make_cuFloatComplex(0, 0);
		//the finnal inner product for p(1)
		for (int i = 0; i < l; i++)
			expectation1 = cuCaddf(expectation1, cuCmulf(
					cuConjf(logicresult[i]), inter1[i]));
		//state after the measurement for 1
		expe_1
		= make_cuFloatComplex (cuCrealf(expectation1),
				cuCimagf(expectation1));
		expe_1 = cuFloatComplex_sqrt(expe_1);
		expectation0 = make_cuFloatComplex(cuCrealf(expe_1), cuCimagf(expe_1));
		for (int i = 0; i < l; i++) {
			inter1[i] = cuFloatComplex_div(inter1[i], expectation1);
			betas = cuCaddf(betas, inter1[i]);
		}
		//--------------------------------------------------------------------------------------
		alphas = expectation0;
		betas = expectation1;

		//		cout<<"alpha: "<<alphas<<" + beta: "<<betas<<endl;
		//calculate the total state
		//Total State = M(0)+M(1)State/Measured(0)+Measured(1)
		mes = 0;
		resultingTotal[k][2* mes ] = expectation0;
		resultingTotal[k][2* mes + 1] = expectation1;
		if (cuCimagf(measureexpected[k][2* mes ]) == 1 && cuCimagf(
				measureexpected[k][2* mes + 1]) == 1) {
			//skip don't cares
			resultingTotalError[k][2* mes ] = make_cuFloatComplex(0, 0);
			resultingTotalError[k][2* mes + 1] = make_cuFloatComplex(0, 0);
		} else {
			if (cuCrealf(expectationsAllState[k][mes]) == cuCrealf(
					expectationsAllState[k][mes + 1])) {
				if (cuCrealf(expectationsAllState[k][mes]) == 1 || cuCrealf(
						expectationsAllState[k][mes]) == 0) {
					if ((cuCrealf(cuCaddf(alphas, betas)) == 1) && (cuCrealf(
							alphas) == cuCrealf(betas))) {
						resultingTotalError[k][mes * 2] = make_cuFloatComplex(
								0, 0);
						resultingTotalError[k][mes * 2 + 1]
						                       = make_cuFloatComplex(0, 0);
					} else {
						resultingTotalError[k][mes * 2] = make_cuFloatComplex(
								abs(0.5 - cuFloatComplex_abs(expectation0)), 0);
						resultingTotalError[k][mes * 2 + 1]
						                       = make_cuFloatComplex(abs(0.5 - cuFloatComplex_abs(
						                    		   expectation1)), 0);
					}
				} else {
					resultingTotalError[k][mes * 2] = make_cuFloatComplex(abs(
							cuFloatComplex_abs(expect0) - cuFloatComplex_abs(expectation0)), 0);
					resultingTotalError[k][mes * 2 + 1] = make_cuFloatComplex(
							abs(cuFloatComplex_abs(expect1) - cuFloatComplex_abs(expectation1)), 0);
				}
			} else {
				resultingTotalError[k][mes * 2] = make_cuFloatComplex(abs(
						cuFloatComplex_abs(expect0) - cuFloatComplex_abs(expectation0)), 0);
				resultingTotalError[k][mes * 2 + 1] = make_cuFloatComplex(abs(
						cuFloatComplex_abs(expect1) - cuFloatComplex_abs(expectation1)), 0);
			}

		}
	}

	int m = 0;
	inter = make_cuFloatComplex(0, 0);
	for (int e = 0; e < l; e++) {
		//		for (int c = 0; c < measurement; c++){
		if (!(cuCimagf(measureexpected[0][e]) == 1 && cuCimagf(
				measureexpected[1][e]) == 1)) {
			//get the error over both possible measurements of the state - for 0 and 1 of the first qubit
			inter = cuCaddf(inter, cuFloatComplex_div(cuCaddf(resultingTotalError[e][0],
					resultingTotalError[e][1]), make_cuFloatComplex(1, 0)));
			//expecting the desired higher probabilities 1
			m++;
		}
		//		}
	}
	//indi->Error /= measurement;
	indi->Error = cuCrealf(cuFloatComplex_div(inter, make_cuFloatComplex(m, 0)));
	indi->Cost = (exp(-pow((divider - indi->Cost), 2)));

	//generate fitness value
	makeFitness(indi);	
	delete (myfinalG);
}
