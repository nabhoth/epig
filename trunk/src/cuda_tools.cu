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



