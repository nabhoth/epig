#include "EpiG.h"

/**************************************************************
* Cleans the allocated CUDA ressources
**************************************************************/
void cleanCUDA(cuComplex* d_M1, cuComplex* d_M2, cuComplex* d_M3, cuComplex* d_MV, cuComplex* d_VI, cuComplex* d_VO){
		cudaFree(d_M1);
		cudaFree(d_M2);
		cudaFree(d_M3);
		cudaFree(d_MV);
		cudaFree(d_VI);
		cudaFree(d_VO);
		cudaThreadExit();
}

/**************************************************************
* Square (thus equal size) Matrix multiplication from CUDA Manual
* Host function - 
* it copies the host data to the device storage thus
* requires both the host data as well as the device data pointers
**************************************************************/
void Mulm(const cuComplex* A, const cuComplex* B, int w, cuComplex* C, cuComplex* d_M1, cuComplex* d_M2, cuComplex* d_M3)
{
    int size;
    // Load A and B to the device
    size = w * w * sizeof(cuComplex);
    cudaMemcpy(d_M1, A, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_M2, B, size, cudaMemcpyHostToDevice);
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(w / dimBlock.x, w / dimBlock.y);
    // Launch the device computation
    matrixMul<<<dimGrid, dimBlock>>>(d_M1, d_M2, w, w, d_M3);
    // Read C from the device
    cudaMemcpy(C, d_M3, size, cudaMemcpyDeviceToHost);
    // Free device memory
}

/**************************************************************
* Square (thus equal size) Matrix multiplication using 
* the pointers to the device stored data only
**************************************************************/
void Mulmnoset(int w, cuComplex* d_M1, cuComplex* d_M2, cuComplex* d_M3)
{
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(w / dimBlock.x, w / dimBlock.y);
    // Launch the device computation
    matrixMul<<<dimGrid, dimBlock>>>(d_M1, d_M2, w, w, d_M3);
}

/**************************************************************
* Matrix multiplication from CUDA Manual
* Device function - only for Matrices that have size multiple of 
* Block_SIZE
* Device multiplication function called by Mul()
* Compute C = A * B
* wA is the width of A
* wB is the width of B
**************************************************************/
__global__ void matrixMul(cuComplex* A, cuComplex* B, int wA, int wB, cuComplex* C)
{
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    // Index of the first sub-matrix of A processed by the block
    int aBegin = wA * BLOCK_SIZE * by;
    // Index of the last sub-matrix of A processed by the block
    int aEnd   = aBegin + wA - 1;
    // Step size used to iterate through the sub-matrices of A
    int aStep = BLOCK_SIZE;
    // Index of the first sub-matrix of B processed by the block
    int bBegin = BLOCK_SIZE * bx;
    // Step size used to iterate through the sub-matrices of B
    int bStep = BLOCK_SIZE * wB;
    // The element of the block sub-matrix that is computed
    // by the thread
    cuComplex Csub = make_cuFloatComplex(0, 0);
    // Loop over all the sub-matrices of A and B required to
    // compute the block sub-matrix
    for (int a = aBegin, b = bBegin;
             a <= aEnd;
             a += aStep, b += bStep) {
        // Shared memory for the sub-matrix of A
        __shared__ cuComplex As[BLOCK_SIZE][BLOCK_SIZE];
        // Shared memory for the sub-matrix of B
        __shared__ cuComplex Bs[BLOCK_SIZE][BLOCK_SIZE];
	// Load the matrices from global memory to shared memory;
	// each thread loads one element of each matrix
	As[ty][tx] = A[a + wA * ty + tx];
	Bs[ty][tx] = B[b + wB * ty + tx];
	// Synchronize to make sure the matrices are loaded
	__syncthreads();
	// Multiply the two matrices together;
	// each thread computes one element
	// of the block sub-matrix
	for (int k = 0; k < BLOCK_SIZE; ++k)
		Csub = cuCaddf(Csub, cuCmulf(As[ty][k], Bs[k][tx]));
	// Synchronize to make sure that the preceding
	// computation is done before loading two new
	// sub-matrices of A and B in the next iteration
	__syncthreads();
     }
  // Write the block sub-matrix to global memory;
  // each thread writes one element
  int c = wB * BLOCK_SIZE * by + BLOCK_SIZE * bx;
  C[c + wB * ty + tx] = Csub;
}

/**************************************************************
* Matrix*Vector multiplication Device function, requires both the
* Vector and Matrix size be the multiple of BLOCK_SIZE
* Device multiplication function called by Mul()
* Compute C = A * b
* wA = wB = w is the width,height of A, width of B
**************************************************************/
__global__ void Muldv(cuComplex* A, cuComplex* B, int w, cuComplex* C)
{
    // Index of the first sub-matrix of A processed by the block
    int aBegin = w * threadIdx.x + blockIdx.x * w * blockDim.x;
    // Index of the last sub-matrix of A processed by the block
    int aEnd   = aBegin + w - 1;
    // Step size used to iterate through the sub-matrices of A
    int aStep = 1;
    // Index of the first sub-vector B processed by the block
    int bBegin = 0;
    // Step size used to iterate through the sub-vectors of B
    int bStep = 1;
    // The element of the block sub-vector that is computed
    // by the thread
    cuComplex Csub = make_cuFloatComplex(0, 0);
    // Shared memory for the sub-matrix of A
    //__shared__ cuComplex As;
    // Shared memory for the sub-matrix of B
    //__shared__ cuComplex Bs;
    // Loop over all the sub-matrices of A and subvectors of B required to
    // compute the block sub-matrix
    for (int a = aBegin, b = bBegin;a <= aEnd;a += aStep, b+= bStep) {
	Csub = cuCaddf(Csub, cuCmulf(A[a], B[b]));
	__syncthreads();
     }
  // Write the block sub-matrix to global memory;
  // each thread writes one element
  int c = (blockDim.x * blockIdx.x) + threadIdx.x;
  C[c] = Csub;
}

/**************************************************************
* Square (thus equal size) Matrix multiplication using 
* the pointers to the device stored data 
**************************************************************/
void Mulvnoset(int w, cuComplex* d_MV, cuComplex* d_VI, cuComplex* d_VO)
{
    dim3 dimBlock(BLOCK_SIZE);
    dim3 dimGrid(w / dimBlock.x);
    // Launch the device computation
    Muldv<<<dimGrid, dimBlock>>>(d_MV, d_VI, w, d_VO);
}

/**************************************************************
* Matrix*Vector multiplication
* Host function - copies the host data to the device and performs the
* multiplication
**************************************************************/
void Mulv(cuComplex* A, cuComplex* B, cuComplex* C, int w, cuComplex* d_MV, cuComplex* d_VI, cuComplex* d_VO)
{
//    long time_0 = clock();
    int size;
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock(BLOCK_SIZE);
    dim3 dimGrid(w / dimBlock.x);
    // Launch the device computation
//		long time = clock();
    size = w * w * sizeof(cuComplex);
    (cudaMemcpy(d_MV, A, size, cudaMemcpyHostToDevice));
    size = w * sizeof(cuComplex);
    (cudaMemcpy(d_VI, B, size, cudaMemcpyHostToDevice));
//		cout<<"CUDA Native Matrix : "<<dimBlock.x<<",  "<<dimGrid.x<<", for "<<w<<" matrix size "<<endl;
    Muldv<<<dimBlock, dimGrid>>>(d_MV, d_VI, w, d_VO);
    (cudaMemcpy(C, d_VO, size, cudaMemcpyDeviceToHost));
//		time = clock() -time ;
//		cout<<"CUDA Native Matrix done: "<<time<<",  "<<(double)(time)/(double)(CLOCKS_PER_SEC)<<", for "<<w<<" matrix size "<<endl;

    // Read C from the device
    //cout<<cudaGetErrorString

    // Free device memory
//		time_0 = clock() -time_0 ;
//		cout<<"total CUDA Native Matrix done: "<<time_0<<",  "<<(double)(time_0)/(double)(CLOCKS_PER_SEC)<<", for "<<w<<" matrix size "<<endl;

}

/**************************************************************
* Divide a Vector of cuComplex by a CuComplex using 
* the pointers to the device stored data 
**************************************************************/
void Divvnoset(int w, cuComplex* d_VecI, cuComplex* d_Val, cuComplex* d_VecO)
{
    dim3 dimBlock(1);
    dim3 dimGrid(1);
    Divdv<<<dimGrid, dimBlock>>>(d_VecI, d_Val, w, d_VecO);
}

/**************************************************************
* Device function
* Device multiplication function called by 
* Divvnoset()
* 
* Value*Vector 
* Compute C = A * b
* w is the height of A
**************************************************************/
__global__ void Divdv(cuComplex* A, cuComplex* B, const int w, cuComplex* C)
{
	for (int a = 0 ;a <= w;a++) {
		C[a] = cuCdivf(A[a], B[0]);
	}
	__syncthreads();
}

/**************************************************************
* Vector*Vector (thus equal size) multiplication using 
* the pointers to the device stored data - no host to device copy
**************************************************************/
void Mulvvnoset(int w, cuComplex* d_VecI, cuComplex* d_VecII, cuComplex* d_VecO)
{
    dim3 dimBlock(1);
    dim3 dimGrid(1);
    // Launch the device computation
    vvMul<<<dimGrid, dimBlock>>>(d_VecI, d_VecII, w, d_VecO, 0);
}

/**************************************************************
* conj(Vector)*Vector (thus equal size) multiplication using 
* the pointers to the device stored data - no host to device copy
**************************************************************/
void Mulconjvvnoset(int w, cuComplex* d_VecI, cuComplex* d_VecII,
cuComplex* d_VecO)
{
    dim3 dimBlock(1);
    dim3 dimGrid(1);
    // Launch the device computation
    vvMul<<<dimGrid, dimBlock>>>(d_VecI, d_VecII, w, d_VecO, 1);
}

/**************************************************************
* Device function
* Device multiplication function called by 
* Mulconjvvnoset()
* Mulvvnoset()
* 
* conj(Vector)*Vector (thus equal size) multiplication  for conj_first == 1
* Vector*Vector (thus equal size) multiplication for conj_first == 0
* Compute C = a * b
* w is the height of A, width of B
**************************************************************/
__global__ void vvMul(cuComplex* A, cuComplex* B, const int w, cuComplex* C, int conj_first)
{
	__shared__ cuComplex Csub;
	Csub = make_cuFloatComplex(0, 0);

	if (conj_first < 1){
		for (int a = 0 ;a <= w;a++) {
		Csub = cuCaddf(Csub, cuCmulf(A[a], B[a]));
		__syncthreads();
		}
	} else {
		for (int a = 0 ;a <= w;a++) {
		Csub = cuCaddf(Csub, cuCmulf(cuConjf(A[a]), B[a]));
		__syncthreads();
		}
	}
  C[0] = Csub;
}

/**************************************************************
* Initialize a cuComplex Vector stored on the device to 0...0
**************************************************************/
__global__ void nullState(int w, cuComplex *d_Vec){
	for (int m = 0; m < w; m++)
		d_Vec[m] = make_cuFloatComplex(0,0);
	__syncthreads();
}

/**************************************************************
* Initialize a cuComplex Vector stored on the device to 1...0
**************************************************************/
__global__ void zeroState(int w, cuComplex *d_Vec){
	for (int m = 0; m < w; m++)
		d_Vec[m] = make_cuFloatComplex(0,0);
	d_Vec[0] = make_cuFloatComplex(1,0);
	__syncthreads();
}


/**************************************************************
* Initialize a cuComplex Vector stored on the device to k-th input state
**************************************************************/
__global__ void initState(int k, cuComplex Desired, cuComplex *d_Vec, int length){
	bool set = false;
	if (k == 0){
		if (cuCrealf(Desired) == 0)
			d_Vec[0] = make_cuFloatComplex(1, 0);
		else
			d_Vec[2] = make_cuFloatComplex(1, 0);

	} else {
		for (int w = 0; w < length; w++){
			if (cuCrealf(d_Vec[w]) != 0 || cuCimagf(d_Vec[w]) != 0)
				set = true;
		}
		if (!set)
			if (cuCrealf(Desired) == 0){
				d_Vec[0] = make_cuFloatComplex(1, 0);
			} else {
				d_Vec[2] = make_cuFloatComplex(1, 0);
			}
	}
	__syncthreads();
}
void vectorInitInit(int w, cuComplex *d_Vec){
	dim3 dimBlock(1);
	dim3 dimGrid(1);
	zeroState<<<dimBlock,dimGrid>>>(w, d_Vec);

}
void vectorInitZero(int w, cuComplex *d_Vec){
	dim3 dimBlock(1);
	dim3 dimGrid(1);
	nullState<<<dimBlock,dimGrid>>>(w, d_Vec);

}
void vectorSetSeqNext(int k, int w, cuComplex Desired, cuComplex *d_Vec){
	dim3 dimBlock(1);
	dim3 dimGrid(1);
	initState<<<dimBlock,dimGrid>>>(k, Desired, d_Vec, w);

}

/**************************************************************
* Matrix*Vector multiplication CPU function
**************************************************************/
void cpuMulv(cuComplex* A, cuComplex*B, cuComplex* C, int w){
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < w; j++) {
			C[i] = cuCaddf(C[i], cuCmulf(B[j],A[i*w+j]));
		}
	}
}

/****************************************
 * returns the matrix representing the circuit of a given individual
 ****************************************/
qGate* GA::computeCUDAMatrix(Individual *indi, bool display) {

	qGate *temp0, *temp1, *temp2, *tempprint, *myfinalG;
	int phasecounter0 = 0;
	int Time = 0;
	int cost = 0;
	int csegs = 0;
	float numIO, val, maxcount;
	int rowcount,rowcounter;
	cuComplex temphase;
	temp0 = new qGate;
	initGate(temp0, ga_resultnum, indi->valuedness,1);
	temp1 = new qGate;
	initGate(temp1, ga_resultnum, indi->valuedness,1);
	temp2 = new qGate;
	initGate(temp2, ga_resultnum, indi->valuedness,1);
	tempprint = new qGate;
	initGate(tempprint, ga_resultnum, indi->valuedness,1);
	myfinalG = new qGate;
	initGate(myfinalG, ga_resultnum, indi->valuedness,1);

	int w = (int(pow((float) indi->valuedness, (float) indi->ioNumber)));
	int N = (int(pow(w,2)));
	cuComplex cu_alpha = make_cuFloatComplex(1, 0);
	cuComplex cu_beta = make_cuFloatComplex(0, 0);

	int even_odd = 0;

	//int the string counters
	int begin = 0;
	int end = indi->my_string.find("p", begin + 1);
	indi->Cost = 0;
	indi->segmentNumber = 0;


	/* Allocate device memory for the matrices */

	temp2->valuedness = temp0->valuedness;
	iMatrix(temp2, temp0->numIO);
	(cudaMemcpy(d_M1, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_M2, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));

	while (end > 0) {
		if ((end - begin) > 1) {
			csegs++;
			//get the first  gate of the Segment and copy it
if (display) {
out_stream << "Matrix at "<<csegs<<" : " << endl;
}
			for (int a = 0; a <= numofgates; a++)
				if (gateArray[a]->representation == indi->my_string.at(
				//if (indi->gateArray[a]->representation == indi->my_string.at(
						begin + 1)) {
					temp0 = ccGate(gateArray[a], temp0);
if (display) {
out_stream<<" Gate: "<<temp0->my_string<<", ";
}
					//temp0 = ccGate(indi->gateArray[a], temp0);
					numIO = (float) temp0->numIO;
					val = (float) temp0->valuedness;
					if (ga_phase == 1) {
						temphase = indi->phases[phasecounter0++];
						maxcount = (int(pow(pow(val, numIO),2)));
						for (int k = 0; k < maxcount; k++)
								temp0->gateMatrix1[k] = cuCmulf(
										temp0->gateMatrix1[k], temphase);
					}
					//continue with the cost
					begin++;
					//cost += indi->gateArray[a]->Cost;
					cost += gateArray[a]->Cost;



					break;
				}
			//get the next one
			for (int b = begin + 1; b < end; b++) {
				for (int a = 0; a <= numofgates; a++)
					//if (indi->gateArray[a]->representation
					if (gateArray[a]->representation
							== indi->my_string.at(b)) {
						//temp1 = ccGate(indi->gateArray[a], temp1);
						temp1 = ccGate(gateArray[a], temp1);
if (display) {
out_stream<<" Gate: "<<temp1->my_string<<", ";
}
						temphase = indi->phases[phasecounter0++];
						val = (float) temp1->valuedness;
						numIO = (float) temp1->numIO;

						if (ga_phase == 1) {
							maxcount = (int(pow(pow(val, numIO),2)));
							for (int k = 0; k < maxcount; k++)
									temp1->gateMatrix1[k] = cuCmulf(
											temp1->gateMatrix1[k], temphase);
						}
						//continue with the cost
						cost += gateArray[a]->Cost;
						//cost += indi->gateArray[a]->Cost;
						a = numofgates * 2;
						//multiply the gates
						//using Kronecker multiplication
						tensorProduct3(temp2, temp0, temp1);
						temp0 = ccGate(temp2, temp0);
					}
			}
			
			if (Time == 0) {
				//if this is the end of the first segment fill the output value
				myfinalG = ccGate(temp0, myfinalG);

				/* Initialize the device matrices with the host matrices */
				(cudaMemcpy(d_M1, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
				Time++;

if (display) {
numIO = (float) temp0->numIO;
val = (float) temp0->valuedness;
out_stream << endl;
maxcount = (int(pow(pow(val, numIO),2)));
rowcount = (int(pow(val, numIO)));
rowcounter = 0;
for (int a = 0; a < rowcount; a++) {
for (int b = 0; b < rowcount; b++) {
out_stream << "(" << cuCrealf(temp0->gateMatrix1[a*rowcount+b])<< "," << cuCimagf(temp0->gateMatrix1[a*rowcount+b]) << ")"<< "";
}
out_stream << endl;
}
out_stream << " -------- " << endl;
}
			} else {
if (display) {
numIO = (float) temp0->numIO;
val = (float) temp0->valuedness;
out_stream << "Matrix at "<<csegs<<" : " << endl;
maxcount = (int(pow(pow(val, numIO),2)));
rowcount = (int(pow(val, numIO)));
rowcounter = 0;
for (int a = 0; a < rowcount; a++) {
for (int b = 0; b < rowcount; b++) {
out_stream << "(" << cuCrealf(temp0->gateMatrix1[a*rowcount+b])<< "," << cuCimagf(temp0->gateMatrix1[a*rowcount+b]) << ")"<< "";
}
out_stream << endl;
}
out_stream << " -------- " << endl;
}
				if (even_odd == 0){
					if (Time == 1){
						cudaMemcpy(d_M2, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
						Mulmnoset(w, d_M1, d_M2, d_M3);

if (display) {
cudaMemcpy(tempprint->gateMatrix1, d_M3, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
numIO = (float) tempprint->numIO;
val = (float) tempprint->valuedness;
out_stream << "Result at "<<csegs<<" : " << endl;
maxcount = (int(pow(pow(val, numIO),2)));
rowcount = (int(pow(val, numIO)));
rowcounter = 0;
for (int a = 0; a < rowcount; a++) {
for (int b = 0; b < rowcount; b++) {
out_stream << "(" << cuCrealf(tempprint->gateMatrix1[a*rowcount+b])<< "," << cuCimagf(tempprint->gateMatrix1[a*rowcount+b]) << ")"<< "";
}
out_stream << endl;
}
out_stream << " -------- " << endl;
}
					}else {
						cudaMemcpy(d_M1, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
						Mulmnoset(w, d_M1, d_M3, d_M2);

if (display) {
cudaMemcpy(tempprint->gateMatrix1, d_M2, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
numIO = (float) tempprint->numIO;
val = (float) tempprint->valuedness;
out_stream << "Result at "<<csegs<<" : " << endl;
maxcount = (int(pow(pow(val, numIO),2)));
rowcount = (int(pow(val, numIO)));
rowcounter = 0;
for (int a = 0; a < rowcount; a++) {
for (int b = 0; b < rowcount; b++) {
out_stream << "(" << cuCrealf(tempprint->gateMatrix1[a*rowcount+b])<< "," << cuCimagf(tempprint->gateMatrix1[a*rowcount+b]) << ")"<< "";
}
out_stream << endl;
}
out_stream << " -------- " << endl;
}
					}
					even_odd = 1;
				} else {
					cudaMemcpy(d_M1, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
					Mulmnoset(w, d_M1, d_M2, d_M3);
					even_odd = 0;

if (display) {
cudaMemcpy(tempprint->gateMatrix1, d_M3, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
numIO = (float) tempprint->numIO;
val = (float) tempprint->valuedness;
out_stream << "Result at "<<csegs<<" : " << endl;
maxcount = (int(pow(pow(val, numIO),2)));
rowcount = (int(pow(val, numIO)));
rowcounter = 0;
for (int a = 0; a < rowcount; a++) {
for (int b = 0; b < rowcount; b++) {
out_stream << "(" << cuCrealf(tempprint->gateMatrix1[a*rowcount+b])<< "," << cuCimagf(tempprint->gateMatrix1[a*rowcount+b]) << ")"<< "";
}
out_stream << endl;
}
out_stream << " -------- " << endl;
}
				}
				Time++;
			}
		}
		//move to the next segment
		begin = indi->my_string.find("p", end);
		end = indi->my_string.find("p", begin + 1);
	}


	/* Read the result back */
	if (even_odd == 0)
		cudaMemcpy(myfinalG->gateMatrix1, d_M2, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
	else 
		cudaMemcpy(myfinalG->gateMatrix1, d_M3, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
		


	//check ifga_phase is defined and multiply the whole matrix by it
	if (ga_phase == 1) {
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
	destroyGate(tempprint);
	delete tempprint;

	return myfinalG;
}
/****************************************
* Generates a circuit and executes |phi> = U|000> multiplication. Then it evaluates |phi> using 
* the formula from Coffman  that quantifies the amount of entanglement. 
 ****************************************/
void* GA::doEntanglState(qGate *myfinalG, cuComplex *entang, bool display) {

	float val = (float) finalGate->valuedness;
	float numIO = (float) finalGate->numIO;
	int l = (int(pow(val, numIO)));
	int x = (int(pow(val, numIO)));
	int y = (int(pow(val, numIO)));
	cuComplex logicresult[(int(pow(val, numIO)))];
	short h[2][2][2];

	//a sort of hash code
	h[0][0][0] = 0;
	h[0][0][1] = 1;
	h[0][1][0] = 2;
	h[0][1][1] = 3;
	h[1][0][0] = 4;
	h[1][0][1] = 5;
	h[1][1][0] = 6;
	h[1][1][1] = 7;

	long time = clock();
	//diplay code
	(cudaMemcpy(d_MV, myfinalG->gateMatrix1, l*l*sizeof(cuComplex), cudaMemcpyHostToDevice));
	if (display){
		for (int a = 0; a < l; a++){
			for (int b = 0; b < l; b++){
				cout << "(" << cuCrealf(myfinalG->gateMatrix1[a*l+b])<< "," << cuCimagf(myfinalG->gateMatrix1[a*l+b]) << ")"<< "";
			}
			cout<<endl;
		}
	}

	//execute U|00...0>
	vectorInitInit(l, d_VI);
	Mulvnoset(l, d_MV, d_VI, d_VO);
	cudaMemcpy(logicresult, d_VO, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);

	if (display){
		cout << "O: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(logicresult[w])<<","<<cuCimagf(logicresult[w])<<")";
		cout<<endl;
	}

	cuComplex a000 = logicresult[0];
	cuComplex a001 = logicresult[1];
	cuComplex a010 = logicresult[2];
	cuComplex a011 = logicresult[3];
	cuComplex a100 = logicresult[4];
	cuComplex a101 = logicresult[5];
	cuComplex a110 = logicresult[6];
	cuComplex a111 = logicresult[7];

	entang[0] = cuCaddf(cuCaddf(cuCaddf(cuCmulf(cuCmulf(a000,a000),cuCmulf(a111,a111)),cuCmulf(cuCmulf(a001,a001),cuCmulf(a110,a110))),cuCmulf(cuCmulf(a010,a010),cuCmulf(a101,a101))),cuCmulf(cuCmulf(a100,a100),cuCmulf(a011,a011)));
	entang[1] = cuCaddf(cuCaddf(cuCaddf(cuCaddf(cuCaddf(cuCmulf(cuCmulf(cuCmulf(a000,a111),a011),a100),cuCmulf(cuCmulf(cuCmulf(a000,a111),a101),a010)),cuCmulf(cuCmulf(cuCmulf(a000,a111),a110),a001)),cuCmulf(cuCmulf(cuCmulf(a011,a100),a101),a010)),cuCmulf(cuCmulf(cuCmulf(a011,a100),a110),a001)),cuCmulf(cuCmulf(cuCmulf(a101,a010),a110),a001));
	entang[2] = cuCaddf(cuCmulf(cuCmulf(cuCmulf(a000,a110),a101),a011),cuCmulf(cuCmulf(cuCmulf(a111,a001),a010),a100));


	cuComplex I5 = {0,0};


	for (int i1 = 0;i1<2;i1++)
	for (int j1 = 0;j1<2;j1++)
	for (int k1 = 0;k1<2;k1++)
	for (int i2 = 0;i2<2;i2++)
	for (int j2 = 0;j2<2;j2++)
	for (int k2 = 0;k2<2;k2++)
	for (int i3 = 0;i3<2;i3++)
	for (int j3 = 0;j3<2;j3++)
	for (int k3 = 0;k3<2;k3++){
		I5 = cuCaddf(I5,cuCmulf(logicresult[h[i1][j1][k1]],cuCmulf(logicresult[h[i2][j2][k2]],cuCmulf(logicresult[h[i3][j3][k3]],cuCmulf(cuConjf(logicresult[h[i1][j2][k3]]),cuCmulf(cuConjf(logicresult[h[i2][j3][k1]]), cuConjf(logicresult[h[i3][j1][k2]])))))));	
	}
	entang[6] = I5;

	cuComplex I2 = {0,0};
	cuComplex I3 = {0,0};
	cuComplex I4 = {0,0};

	for (int i1 = 0;i1<2;i1++)
	for (int j1 = 0;j1<2;j1++)
	for (int k1 = 0;k1<2;k1++)
	for (int i2 = 0;i2<2;i2++)
	for (int j2 = 0;j2<2;j2++)
	for (int k2 = 0;k2<2;k2++){
		I2 = cuCaddf(I2,cuCmulf(logicresult[h[i1][j1][k1]],cuCmulf(logicresult[h[i2][j2][k2]],cuCmulf(cuConjf(logicresult[h[i1][j1][k2]]),cuConjf(logicresult[h[i2][j2][k1]])))));
		I3 = cuCaddf(I3,cuCmulf(logicresult[h[i1][j1][k1]],cuCmulf(logicresult[h[i2][j2][k2]],cuCmulf(cuConjf(logicresult[h[i1][j2][k1]]),cuConjf(logicresult[h[i2][j1][k2]])))));
		I4 = cuCaddf(I4,cuCmulf(logicresult[h[i1][j1][k1]],cuCmulf(logicresult[h[i2][j2][k2]],cuCmulf(cuConjf(logicresult[h[i1][j2][k2]]),cuConjf(logicresult[h[i2][j1][k1]])))));
	}
	
	entang[3] = I2;
	entang[4] = I3;
	entang[5] = I4;
	if (display){
		cout << "D4: ";
//			cout<<"("<<cuCrealf(d4)<<","<<cuCimagf(d4)<<")";
		cout<<endl;
	}
	
	return 0;
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
	//of the ga_measurement
	qGate *myfinalG;
	float val = (float) finalGate->valuedness;
	float numIO = (float) finalGate->numIO;
	float err,maxcount = 0;
	int l = (int(pow(val, numIO)));
	int x = (int(pow(val, numIO)));
	int y = (int(pow(val, numIO)));
	cuComplex inter[(int(pow(val, numIO)))];
	cuComplex inter_b[(int(pow(val, numIO)))];
	cuComplex logicresult[(int(pow(val, numIO)))];
	cuComplex stb[(int(pow(val, numIO)))];
	cuComplex expectations[(int(pow(val, numIO) * val))], cu_alphas;
	cuComplex resultingTotalError[(int(pow(val, numIO)))][ga_measurement * (int) val];
	cuComplex resultingTotal[(int(pow(val, numIO)))][ga_measurement * (int) val];
	complex<float> temp_val;

	Individual *indi = ind;
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;

	//threaded code
	long time = clock();
	//accelerated matrix computation
	myfinalG = computeCUDAMatrix(indi, false);
	(cudaMemcpy(d_MV, myfinalG->gateMatrix1, l*l*sizeof(cuComplex), cudaMemcpyHostToDevice));
	numIO = (float) myfinalG->numIO;
	//check ifga_phase is defined and multiply the whole matrix by it
	if (ga_phase == 1) {
		maxcount = (int(pow(pow(val, numIO),2)));
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(
							myfinalG->gateMatrix1[k], indi->phase);
		}
	}

if (display){
		out_stream<<indi->my_string<<endl;
}

	//for all values in the desireed sequence


	(cudaMemcpy(myfinalG->gateMatrix1, d_MV, l*l*sizeof(cuComplex), cudaMemcpyDeviceToHost));

	vectorInitZero(l, d_VI);

	for (int k = 0; k < measured_desired_output_records; k++) {


		//empty all intermediary variables
		for (int w = 0; w < l; w++)
			inter[w] = logicresult[w] = stb[w] = make_cuFloatComplex(0, 0);

		vectorSetSeqNext(k, l, sequence_desired[k][0], d_VI);
if (display){

		cudaMemcpy(stb, d_VI, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);
		cout << "I: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(stb[w])<<","<<cuCimagf(stb[w])<<")";
		cout<<endl;
}
		
		Mulvnoset(l, d_MV, d_VI, d_VO);


if (display){
		cout << "O: ";
		cudaMemcpy(logicresult, d_VO, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(logicresult[w])<<","<<cuCimagf(logicresult[w])<<")";
		cout<<endl;
}
		//measure for only the desired state of the single output qubit
		for (int w = 0; w < l; w++)
			inter[w] = inter_b[w] = make_cuFloatComplex(0, 0);


		//measure for 0 when non final sequence element
		if (k != measured_desired_output_records-1){
			Mulvnoset(l, d_ME_00, d_VO, d_VI);
		}else{
			//measure for 1 when the end of the sequence was reached
			Mulvnoset(l, d_ME_01, d_VO, d_VI);
		}


if (display){
		cout << "Measure Output: ";
		cudaMemcpy(inter, d_VI, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(inter[w])<<","<<cuCimagf(inter[w])<<")";
		cout<<endl;
}

		expectations[0] = make_cuFloatComplex(0, 0);
		//the finnal inner product for p(r) - generate the expectation for the measured output
		Mulconjvvnoset(l,d_VI, d_VI, d_Value);

		//set the probability of obtaining the desired output for this input value
		cudaMemcpy(&resultingTotal[k][0], &d_Value[0], sizeof(cuComplex), cudaMemcpyDeviceToHost);
				
		//skip don't cares
		if (cuCrealf(sequence_desired[k][1]) == 0 && cuCimagf(sequence_desired[k][1]) == -1) {
			resultingTotalError[k][(int) val]     = make_cuFloatComplex(0, 0);
		} else {
			//calculate the error if this is not to be considered as don't care
			resultingTotalError[k][(int) val] = make_cuFloatComplex(sqrt(pow(cuCrealf(cuCsubf(resultingTotal[k][0],sequence_desired[k][1])),(float) 2)), 0);
		}


		temp_val = complex<float> (cuCrealf(expectations[0]), cuCimagf(expectations[0]));
		temp_val = sqrt(temp_val);
		cu_alphas = make_cuFloatComplex(real(temp_val), imag(temp_val));

		//if the ga_measurement generated at least some output
		if (!(cuCrealf( cu_alphas) == 0 && cuCimagf( cu_alphas) == 0)){
			//if the output qubit was measured in 0, normalize the vector
			cudaMemcpy(&d_Value[0], &cu_alphas, sizeof(cuComplex),cudaMemcpyHostToDevice);
			Divvnoset(l, d_VI, d_Value, d_VO);
			cudaMemcpy(d_VI, d_VO, sizeof(cuComplex),cudaMemcpyDeviceToDevice);
		} else {
			//if the ga_measurement was not succesful, use the output of the circuit 
		}



		//measure the input qubit for 0 - to reset it and allow consequent initialization
		if (k < measured_desired_output_records - 1)
			if (cuCrealf(sequence_desired[k+1][0]) == 0){
				Mulvnoset(l, d_ME_10, d_VI, d_VO);
			} else {
				Mulvnoset(l, d_ME_11, d_VI, d_VO);
			}	

if (display){
		cudaMemcpy(inter_b, d_VO, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);
		cout << "Measure Input: ";
		for (int w = 0; w < l; w++)
			cout<<"("<<cuCrealf(inter_b[w])<<","<<cuCimagf(inter_b[w])<<")";
		cout<<endl;
}
		expectations[0] = make_cuFloatComplex(0, 0);
		//calculate the expectation value for measureing the qubit in state |0>
		Mulconjvvnoset(l,d_VO, d_VO, d_Value);
		cudaMemcpy(expectations, d_Value, sizeof(cuComplex), cudaMemcpyDeviceToHost);
		temp_val = complex<float> (cuCrealf(expectations[0]), cuCimagf(expectations[0]));
		temp_val = sqrt(temp_val);
		cu_alphas = make_cuFloatComplex(real(temp_val), imag(temp_val));

		cudaMemcpy(d_VI, d_VO, l*sizeof(cuComplex), cudaMemcpyDeviceToDevice);

		//if the input qubit was measured in 0, normalize the vector
		if (!(cuCrealf( cu_alphas) == 0 && cuCimagf( cu_alphas) == 0)){
			cudaMemcpy(&d_Value[0], &cu_alphas, sizeof(cuComplex),cudaMemcpyHostToDevice);
			Divvnoset(l, d_VI, d_Value, d_VI);

		}

	}

	int m;
	err = 0;
	//calculate the overall error
	for (int c = 0; c < ga_measurement; c++) {
		m = measured_desired_output_records;
		inter[0] = make_cuFloatComplex(0, 0);
		for (int e = 0; e < measured_desired_output_records; e++) {
			//inter[0] = cuCmulf(inter[0], resultingTotal[e][0]);
			inter[0] = cuCaddf(inter[0], make_cuFloatComplex(cuCrealf(resultingTotalError[e][(int) val]), 0));
		}
		err += cuCrealf(cuCdivf(inter[0], make_cuFloatComplex(m, 0)));
	}
	err /= ga_measurement;
	indi->Error = err;

	//generate fitness value
	if (indi->Error != 0) {
		switch (ga_replicator) {
		case 0:
			//simple fitness1
			indi->fitness = (1 - indi->Error);
			break;
		case 1:
			//simple fitness1
			indi->fitness = (1 / (indi->Error + 1));
			break;
		case 2:
			//scaled complex fitness1
			indi->fitness = (ga_alpha * (1 - indi->Error) + ga_beta * (exp(-pow(abs(
					ga_divider - indi->Cost), 2))));
			break;
		case 3:
			//scaled complex fitness2
			indi->fitness = (ga_alpha * (1 / (indi->Error + 1)) + ga_beta * (exp(
					-pow(abs(ga_divider - indi->Cost), 2))));
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
 * calculates theprobabilities for obtaining |0> and |1>
 * Right now it is doing ga_measurement in increasing order
 * of all defined bits to measure
 ****************************************/
void* GA::doMeasureFitness(Individual *ind, bool display) {

	//define temporal variables for the results and the manupaltion
	//of the ga_measurement
	qGate *measures[10], *myfinalG;
	float val = (float) finalGate->valuedness;
	float numIO = (float) finalGate->numIO;
	float err,maxcount = 0;
	int l = (int(pow(val, numIO)));
	int x = (int(pow(val, numIO)));
	int y = (int(pow(val, numIO)));
	cuComplex inter[(int(pow(val, numIO)))];
	cuComplex logicresult[(int(pow(val, numIO)))];
	cuComplex sta[(int(pow(val, numIO)))];
	cuComplex expectations[(int(pow(val, numIO) * val))];
	cuComplex resultingTotalError[(int(pow(val, numIO)))][ga_measurement
	                                                      * (int) val];
	cuComplex resultingTotal[(int(pow(val, numIO)))][ga_measurement * (int) val];

	Individual *indi = ind;
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;

	//myfinalG = computeMatrix(indi);
	//cout << indi->my_string<<endl;
	long time = clock();
	myfinalG = computeCUDAMatrix(indi,false);
	numIO = (float) myfinalG->numIO;
	//check ifga_phase is defined and multiply the whole matrix by it
	if (ga_phase == 1) {
		maxcount = (int(pow(pow(val, numIO),2)));
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(
							myfinalG->gateMatrix1[k], indi->phase);
		}
	}


	//for single qubit this should always be only 0, i.e. <1
	for (int m = 0; m < ga_measurement; m++) {
		for (int r = 0; r < val; r++) {
			measures[r] = measurements[(measurementQBits[m] * (int) val) + r];
		}
		//for all input states l
		for (int k = 0; k < measured_desired_output_records; k++) {

			for (int w = 0; w < l; w++)
				sta[w] = inter[w] = logicresult[w] = make_cuFloatComplex(0, 0);
			//set the k-th input state - minterm
			sta[measured_desired_output_records_idx[k]] = make_cuFloatComplex(1, 0);

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


				//the finnal inner product for p(r)
				expectations[((m * (int) val) + r)] = make_cuFloatComplex(0, 0);
				for (int j = 0; j < l; j++) {
					expectations[(m * (int) val + r)] = cuCaddf(expectations[(m
							* (int) val + r)], cuCmulf(cuConjf(logicresult[j]),
									inter[j]));
				}

				resultingTotal[measured_desired_output_records_idx[k]][(m * (int) val + r)] = expectations[(m
						* (int) val + r)];
				
				if (cuCrealf(measureexpected[measured_desired_output_records_idx[k]][(int) val * m + r]) == 0
						&& cuCimagf(measureexpected[measured_desired_output_records_idx[k]][(int) val * m + r]) == 1) {
					//skip don't cares
					resultingTotalError[measured_desired_output_records_idx[k]][m * (int) val + r]
					                       = make_cuFloatComplex(0, 0);
				} else {
					resultingTotalError[measured_desired_output_records_idx[k]][m * (int) val + r]
					                       = make_cuFloatComplex(sqrt(pow(cuCrealf(cuCsubf(
					                    		   resultingTotal[measured_desired_output_records_idx[k]][(m * (int) val) + r],
					                    		   measureexpected[measured_desired_output_records_idx[k]][((int) val * m) + r])),
					                    		   (float) 2)), 0);
				}
				//for (int i = 0; i < l; i++) 

			}
		}
	}

	int m;
	err = 0;
	for (int c = 0; c < ga_measurement; c++) {
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
		err += cuCrealf(cuCdivf(inter[0], make_cuFloatComplex(m, 0)));
	}
	err /= ga_measurement;
	indi->Error = err;

	//generate fitness value
	if (indi->Error != 0) {
		switch (ga_replicator) {
		case 0:
			//simple fitness1
			indi->fitness = (1 - indi->Error);
			break;
		case 1:
			//simple fitness1
			indi->fitness = (1 / (indi->Error + 1));
			break;
		case 2:
			//scaled complex fitness1
			indi->fitness = (ga_alpha * (1 - indi->Error) + ga_beta * (exp(-pow(abs(
					ga_divider - indi->Cost), 2))));
			break;
		case 3:
			//scaled complex fitness2
			indi->fitness = (ga_alpha * (1 / (indi->Error + 1)) + ga_beta * (exp(
					-pow(abs(ga_divider - indi->Cost), 2))));
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
	if (display) {
		out_stream << "Error: " << indi->Error << endl;
		out_stream << "Fitness: " << indi->fitness << endl;
		out_stream << "Valuedness: " << indi->valuedness << endl;
		out_stream << endl;
		out_stream << "Input St. Des_0 Des_1 Out_0 Out_1 Error_0 Error_1 "
		<< endl;

		for (int m = 0; m < ga_measurement; m++) {
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
	return 0;
}

/****************************************
 * calculates the probabilities for obtaining two distinct multi qubit states of the form
 * (p_0 + p_1 + .. + p_k)|x_0...x_k> and (r_0 + r_1 + .. + r_k)|y_0....y_k>
 ****************************************/
void* GA::doMultiMeasureFitness(Individual *ind, bool display) {

	//define temporal variables for the results and the manipulation
	//of the ga_measurement
	int l = (int(pow((float) 2, (float) finalGate->numIO)));
	int x = (int(pow((float) 2, (float) finalGate->numIO)));
	int y = (int(pow((float) 2, (float) finalGate->numIO)));
	int maxcount = 0;
	int mes = 0;
	cuComplex inter0[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex inter1[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex logicresult[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex sta[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex expectation0, expectation1, ga_alphas, ga_betas;
	complex<float> expe_0, expe_1;
	cuComplex resultingTotalError[(int(pow((float) 2, (float) finalGate->numIO)))][ga_measurement* 2];
	cuComplex resultingTotal[(int(pow((float) 2, (float) finalGate->numIO)))][ga_measurement* 2];

	cuComplex inter, expect0, expect1;
//	qGate *measure0, *measure1;
	qGate *myfinalG, measure;

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
	myfinalG = computeCUDAMatrix(indi,false);
	rc = pthread_mutex_unlock(&mtx);

	if (ga_phase == 1) {
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

		//init the ga_measurement operator with respect to the desired output
		//for each measured ouptut state get the expectation values
		expect0 = expectationsAllState[k][0];
		expect1 = expectationsAllState[k][1];

		//apply the ga_measurement operator for the desired state
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
		ga_alphas = make_cuFloatComplex(0, 0);
		//the finnal inner product for p(0)
		for (int j = 0; j < l; j++)
			expectation0 = cuCaddf(expectation0, cuCmulf(
					cuConjf(logicresult[j]), inter0[j]));
		//state after the ga_measurement for 0
		expe_0
		= complex<float> (cuCrealf(expectation0),
				cuCimagf(expectation0));
		expe_0 = sqrt(expe_0);
		expectation0 = make_cuFloatComplex(real(expe_0), imag(expe_0));

		for (int i = 0; i < l; i++) {
			inter0[i] = cuCdivf(inter0[i], expectation0);
			ga_alphas = cuCaddf(ga_alphas, inter0[i]);
		}
		//p(1) ---------------------------------------------------------------------------------
		//vector inner product
		expectation1 = make_cuFloatComplex(0, 0);
		ga_betas = make_cuFloatComplex(0, 0);
		//the finnal inner product for p(1)
		for (int i = 0; i < l; i++)
			expectation1 = cuCaddf(expectation1, cuCmulf(
					cuConjf(logicresult[i]), inter1[i]));
		//state after the ga_measurement for 1
		expe_1
		= complex<float> (cuCrealf(expectation1),
				cuCimagf(expectation1));
		expe_1 = sqrt(expe_1);
		expectation0 = make_cuFloatComplex(real(expe_1), imag(expe_1));
		for (int i = 0; i < l; i++) {
			inter1[i] = cuCdivf(inter1[i], expectation1);
			ga_betas = cuCaddf(ga_betas, inter1[i]);
		}
		//--------------------------------------------------------------------------------------
		ga_alphas = expectation0;
		ga_betas = expectation1;

		//		cout<<"ga_alpha: "<<ga_alphas<<" + ga_beta: "<<ga_betas<<endl;
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
					if ((cuCrealf(cuCaddf(ga_alphas, ga_betas)) == 1) && (cuCrealf(
							ga_alphas) == cuCrealf(ga_betas))) {
						resultingTotalError[k][mes * 2] = make_cuFloatComplex(
								0, 0);
						resultingTotalError[k][mes * 2 + 1]
						                       = make_cuFloatComplex(0, 0);
					} else {
						resultingTotalError[k][mes * 2] = make_cuFloatComplex(
								abs(0.5 - cuCabsf(expectation0)), 0);
						resultingTotalError[k][mes * 2 + 1]
						                       = make_cuFloatComplex(abs(0.5 - cuCabsf(
						                    		   expectation1)), 0);
					}
				} else {
					resultingTotalError[k][mes * 2] = make_cuFloatComplex(abs(
							cuCabsf(expect0) - cuCabsf(expectation0)), 0);
					resultingTotalError[k][mes * 2 + 1] = make_cuFloatComplex(
							abs(cuCabsf(expect1) - cuCabsf(expectation1)), 0);
				}
			} else {
				resultingTotalError[k][mes * 2] = make_cuFloatComplex(abs(
						cuCabsf(expect0) - cuCabsf(expectation0)), 0);
				resultingTotalError[k][mes * 2 + 1] = make_cuFloatComplex(abs(
						cuCabsf(expect1) - cuCabsf(expectation1)), 0);
			}

		}
	}

	int m = 0;
	inter = make_cuFloatComplex(0, 0);
	for (int e = 0; e < l; e++) {
		//		for (int c = 0; c < ga_measurement; c++){
		if (!(cuCimagf(measureexpected[0][e]) == 1 && cuCimagf(
				measureexpected[1][e]) == 1)) {
			//get the error over both possible ga_measurements of the state - for 0 and 1 of the first qubit
			inter = cuCaddf(inter, cuCdivf(cuCaddf(resultingTotalError[e][0],
					resultingTotalError[e][1]), make_cuFloatComplex(1, 0)));
			//expecting the desired higher probabilities 1
			m++;
		}
		//		}
	}
	//indi->Error /= ga_measurement;
	indi->Error = cuCrealf(cuCdivf(inter, make_cuFloatComplex(m, 0)));
	indi->Cost = (exp(-pow((ga_divider - indi->Cost), 2)));

	//generate fitness value
	if (indi->Error != 0) {
		switch (ga_replicator) {
		case 0:
			//simple fitness1
			indi->fitness = (1 - indi->Error);
			break;
		case 1:
			//simple fitness1
			indi->fitness = (1 / (indi->Error + 1));
			break;
		case 2:
			//scaled complex fitness1
			indi->fitness = (ga_alpha * (1 - indi->Error) + ga_beta * indi->Cost);
			break;
		case 3:
			//scaled complex fitness2
			indi->fitness = (ga_alpha * (1 / (indi->Error + 1)) + ga_beta
					* indi->Cost);
			break;
		case 4:
			indi->fitness = (1 - indi->Error);
			break;
		}
	} else {
		indi->Error = 0;
		indi->fitness = 1;//we found a good individual;
	}
	//	cout<<"Fitness: "<<indi->fitness<<endl;
	destroyGate(myfinalG);
	delete (myfinalG);
	return 0;
}
