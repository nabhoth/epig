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
void Mulm(cuComplex* A, cuComplex* B, int w, cuComplex* C, cuComplex* d_M1, cuComplex* d_M2, cuComplex* d_M3)
{
    int size;
    // Load A and B to the device
    size = w * w * sizeof(cuComplex);
    cudaMemcpy(d_M1, A, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_M2, B, size, cudaMemcpyHostToDevice);
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock(BLOCK_MAX_SIZE, BLOCK_MAX_SIZE);
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
void Mulmnoset(int w, cuComplex* d_1, cuComplex* d_2, cuComplex* d_3)
{
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock(BLOCK_MAX_SIZE, BLOCK_MAX_SIZE);
    dim3 dimGrid(w / dimBlock.x, w / dimBlock.y);
    // Launch the device computation
    matrixMul<<<dimGrid, dimBlock>>>(d_1, d_2, w, w, d_3);
}
#ifdef __CUBLAS__
/**************************************************************
* Interface to cuBlas matrix multiplication
**************************************************************/
void cuBlasMMultiSet(int columns, cuComplex *A, cuComplex *B, cuComplex *C, cuComplex* d_M1, cuComplex* d_M2, cuComplex* d_M3){
	int N = columns*columns;
	cudaMemcpy(d_M1, A, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
	cudaMemcpy(d_M2, B, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
	cublasCgemm ('N','N', columns, columns, columns, make_cuFloatComplex(1, 0),  A, columns, B, columns, make_cuFloatComplex(0, 0), C, columns);
	cudaMemcpy(C, d_M3, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
}
/**************************************************************
* Interface to cuBlas matrix multiplication
**************************************************************/
void cuBlasMMulti(int columns, cuComplex *A, cuComplex *B, cuComplex *C){
	cublasCgemm ('N','N', columns, columns, columns, make_cuFloatComplex(1, 0),  A, columns, B, columns, make_cuFloatComplex(0, 0), C, columns);
}
#endif
/**************************************************************
* Matrix multiplication of matrices of square sizes
* Device function - only for Matrices that have size multiple of 
* Block_SIZE
* Compute C = A KRON B
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
    int aBegin = wA * BLOCK_MAX_SIZE * by;
    // Index of the last sub-matrix of A processed by the block
    int aEnd   = aBegin + wA - 1;
    // Step size used to iterate through the sub-matrices of A
    int aStep = BLOCK_MAX_SIZE;
    // Index of the first sub-matrix of B processed by the block
    int bBegin = BLOCK_MAX_SIZE * bx;
    // Step size used to iterate through the sub-matrices of B
    int bStep = BLOCK_MAX_SIZE * wB;
    // The element of the block sub-matrix that is computed
    // by the thread
    cuComplex Csub = make_cuFloatComplex(0, 0);
    // Loop over all the sub-matrices of A and B required to
    // compute the block sub-matrix
    for (int a = aBegin, b = bBegin;
             a <= aEnd;
             a += aStep, b += bStep) {
        // Shared memory for the sub-matrix of A
        __shared__ cuComplex As[BLOCK_MAX_SIZE][BLOCK_MAX_SIZE];
        // Shared memory for the sub-matrix of B
        __shared__ cuComplex Bs[BLOCK_MAX_SIZE][BLOCK_MAX_SIZE];
	// Load the matrices from global memory to shared memory;
	// each thread loads one element of each matrix
	As[ty][tx] = A[a + wA * ty + tx];
	Bs[ty][tx] = B[b + wB * ty + tx];
	// Synchronize to make sure the matrices are loaded
	__syncthreads();
	// Multiply the two matrices together;
	// each thread computes one element
	// of the block sub-matrix
	for (int k = 0; k < BLOCK_MAX_SIZE; ++k)
		Csub = cuCaddf(Csub, cuCmulf(As[ty][k], Bs[k][tx]));
	// Synchronize to make sure that the preceding
	// computation is done before loading two new
	// sub-matrices of A and B in the next iteration
	__syncthreads();
     }
  // Write the block sub-matrix to global memory;
  // each thread writes one element
  int c = wB * BLOCK_MAX_SIZE * by + BLOCK_MAX_SIZE * bx;
  C[c + wB * ty + tx] = Csub;
}

/**************************************************************
* Square  Kronecker Matrix multiplication
* it copies the host data to the device storage thus
* requires both the host data as well as the device data pointers
**************************************************************/
void Kronm(qGate* A, qGate* B, int wa, int wb, qGate* C, cuComplex* d_M1, cuComplex* d_M2, cuComplex* d_M3)
{
    int size = wa * wa * sizeof(cuComplex);
    // Load A and B to the device
    //cout<<cudaGetErrorString
    (cudaMemcpy(d_M1, A->gateMatrix1, size, cudaMemcpyHostToDevice));
    size = wb * wb * sizeof(cuComplex);
    (cudaMemcpy(d_M2, B->gateMatrix1, size, cudaMemcpyHostToDevice));
    size = wa * wb;
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock;
//cout<<wa<<" "<<wb<<" "<<BLOCK_MAX_SIZE<<endl;
    if ((wb*wa) < BLOCK_MAX_SIZE)
	    dimBlock = dim3((wa*wb), (wa*wb));
    else if (wb < BLOCK_MAX_SIZE)
	    dimBlock = dim3(wb, wb);
    else
	    dimBlock = dim3(BLOCK_MAX_SIZE, BLOCK_MAX_SIZE);
//	dimBlock = dim3(1,1);
    dim3 dimGrid(size / dimBlock.x, size / dimBlock.y);
    // Launch the device computation
//    cout<<"Launching: Grid of "<<dimGrid.x<<"*"<<dimGrid.y<<" each with "<<dimBlock.x<<"*"<<dimBlock.y<<" threads "<<endl;
    matrixKron<<<dimGrid, dimBlock>>>(d_M1, d_M2, wa, wb, d_M3);
    // Read C from the device
    size = size*size*sizeof(cuComplex);
    (cudaMemcpy(C->gateMatrix1, d_M3, size, cudaMemcpyDeviceToHost));
    C->valuedness = A->valuedness;
    C->realIO = (A->realIO + B->realIO);
    C->restrictions_number = (A->restrictions_number + B->restrictions_number);
    C->numIO = (A->numIO + B->numIO);
    C->representation[0] = 'f';
    C->representation[1] = 'f';
    C->representation[2] = 'f';
    C->representation[REPSIZE] = '\0';
    C->my_string = "";
    C->Cost = A->Cost;
    // Free device memory
}

/**************************************************************
* Square Kronecker Matrix multiplication
* it copies the host data to the device storage thus
* requires both the host data as well as the device data pointers
**************************************************************/
void Kronm(cuComplex* A, cuComplex* B, int wa, int wb, cuComplex* C, cuComplex* d_M1, cuComplex* d_M2, cuComplex* d_M3)
{
 //   cout<<wa<<' '<<wb<<' '<<d_M1<<endl;
    int size = wa * wa * sizeof(cuComplex);
    // Load A and B to the device
    cudaMemcpy(d_M1, A, size, cudaMemcpyHostToDevice);
    size = wb * wb * sizeof(cuComplex);
    cudaMemcpy(d_M2, B, size, cudaMemcpyHostToDevice);
    size = wa * wb;
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock;
    if ((wb*wa) < BLOCK_MAX_SIZE)
	    dimBlock = dim3((wa*wb), (wa*wb));
    else if (wb < BLOCK_MAX_SIZE)
	    dimBlock = dim3(wb, wb);
    else
	    dimBlock = dim3(BLOCK_MAX_SIZE, BLOCK_MAX_SIZE);
    dim3 dimGrid(size / dimBlock.x, size / dimBlock.y);
    // Launch the device computation
    //cout<<"Launching: Grid of "<<dimGrid.x<<"*"<<dimGrid.y<<" each with "<<dimBlock.x<<"*"<<dimBlock.y<<" threads "<<endl;
    matrixKron<<<dimGrid, dimBlock>>>(d_M1, d_M2, wa, wb, d_M3);
    // Read C from the device
    size = size*size*sizeof(cuComplex);
    cudaMemcpy(C, d_M3, size, cudaMemcpyDeviceToHost);
    // Free device memory
}

/**************************************************************
* Square (thus equal size) Matrix multiplication using 
* the pointers to the device stored data only
**************************************************************/
void Kronmnoset(int w1, int w2, cuComplex* d_1, cuComplex* d_2, cuComplex* d_3)
{
    int size = w1*w2;
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    //maximum 256 because it is a square of 16x16
    dim3 dimBlock;
    if (w2 < BLOCK_MAX_SIZE)
	    dimBlock = dim3(BLOCK_SIZE, BLOCK_SIZE);
    else
	    dimBlock = dim3(BLOCK_MAX_SIZE, BLOCK_MAX_SIZE);
    dim3 dimGrid(size / dimBlock.x, size / dimBlock.y);
    // Launch the device computation
    matrixKron<<<dimGrid, dimBlock>>>(d_1, d_2, w1, w2, d_3);
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
__global__ void matrixKron(cuComplex* A, cuComplex* B, int wA, int wB, cuComplex* C)
{
    // Block index
//    int bx = blockIdx.x;
//    int by = blockIdx.y;
    // Thread index
//    int tx = threadIdx.x;
//    int ty = threadIdx.y;
//    int bid = bx + by * gridDim.x;
//    int tid = tx + ty * blockDim.x;
    int wc = wA*wB;
    int c = blockIdx.x*blockDim.x+threadIdx.x;
    int b = blockIdx.y*blockDim.y+threadIdx.y;

    int aindex = (((c/wB)%wA)+(((b/wB)%wA)*wA));
    int bindex = ((c%wB)+((b)%wB)*wB);
    int  cindex = (c+b*wc);
    //C[cindex] = make_cuFloatComplex(aindex, bindex);
    //C[0] = make_cuFloatComplex(aindex, bindex);
    C[cindex] = cuCmulf(A[aindex],B[bindex]);
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
    dim3 dimBlock(BLOCK_MAX_SIZE);
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
    // the matrix dimensions are multiples of BLOCK_MAX_SIZE
    dim3 dimBlock(BLOCK_MAX_SIZE);
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
    // Launch the device computation
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
* Initialize a cuComplex Vector stored on the device to k-th input state
* and preserve previous states 
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
* Initialize a cuComplex Vector stored on the device to k-th input state
**************************************************************/
__global__ void initNewState(int k, cuComplex Desired, cuComplex *d_Vec, int length){
		for (int w = 0; w < length; w++)
			if (w == k){
				d_Vec[w] = make_cuFloatComplex(cuCrealf(Desired), cuCimagf(Desired));
			} else {
				d_Vec[w] = make_cuFloatComplex(0, 0);
			}
	__syncthreads();
}
void vectorSetIterState(int k, int w, cuComplex Desired, cuComplex *d_Vec){
	dim3 dimBlock(1);
	dim3 dimGrid(1);
	initNewState<<<dimBlock,dimGrid>>>(k, Desired, d_Vec, w);
}

/**************************************************************
* Matrix*Matrix multiplication GPU wrapper function
**************************************************************/
void GA::gpuMMmult(cuComplex* A, cuComplex* B, int width, int valuedness, cuComplex* C){
        int w = (int(pow((float) valuedness, (float) width)));
        int N = (int(pow(w,2)));
        (cudaMemcpy(d_M1, A, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
        (cudaMemcpy(d_M2, B, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
	Mulmnoset(w, d_M1, d_M2, d_M3);
        (cudaMemcpy(C, d_M3, N*sizeof(cuComplex), cudaMemcpyDeviceToHost));
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

#ifdef __CUBLAS__
/****************************************
 * returns the matrix representing the circuit of a given individual
 * computed using the CUBLAS matrix multiplication 
 ****************************************/
qGate* GA::computeCUBLASMatrix(Individual *indi) {

        qGate *temp0, *temp1, *temp2, *myfinalG;
        int phasecounter0 = 0;
        bool initTime = true;
        int cost = 0;
        int csegs = 0;
        int num_g = 0;
        float numIO, val, maxcount;
	char *gate_rep = new char[4];
        cuComplex temphase;
        temp0 = new qGate;
        initGate(temp0, resultnum, indi->valuedness,1);
        temp1 = new qGate;
        initGate(temp1, resultnum, indi->valuedness,1);
        temp2 = new qGate;
        initGate(temp2, resultnum, indi->valuedness,1);
        myfinalG = new qGate;
        initGate(myfinalG, resultnum, indi->valuedness,1);

        int w = (int(pow((float) indi->valuedness, (float) indi->ioNumber)));
        int N = (int(pow(w,2)));
        cuComplex alpha = make_cuFloatComplex(1, 0);
        cuComplex beta = make_cuFloatComplex(0, 0);

	int even_odd = 0;

        int begin = 0;
        int end = indi->my_string.find("p", begin + 1);
        indi->Cost = 0;
        indi->segmentNumber = 0;
        temp2->valuedness = temp0->valuedness;
        zeroMatrix(temp2, temp0->numIO);
        //cout<<cudaGetErrorString
        (cudaMemcpy(d_M1, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));


//	cout<<"Individual matrix computing:  "<<indi->my_string<<endl;
        while (end > 0) {
                if ((end - begin) > 1) {
                        csegs++;
                        //get the first  gate of the Segment and copy it
			getGateRepfromStr(begin+1,indi->my_string, gate_rep);
			num_g = getGate(gate_rep);
                        temp0 = ccGate(gateArray[num_g], temp0);
                        numIO = (float) temp0->numIO;
                        val = (float) temp0->valuedness;
                        if (phase > 0) {
				temphase = indi->phases[phasecounter0++];
                                maxcount = (int(pow(pow(val, numIO),2)));
                                for (int k = 0; k < maxcount; k++)
                                        temp0->gateMatrix1[k] = cuCmulf(
                                                 temp0->gateMatrix1[k], temphase);
                        }
                        //continue with the cost
                        begin++;
                        cost += gateArray[num_g]->Cost;
                        //get the next one
			begin += 3;
                        for (int b = begin; b < end; b+=REPSIZE) {
				getGateRepfromStr(b,indi->my_string, gate_rep);
				num_g = getGate(gate_rep);
                                temp1 = ccGate(gateArray[num_g], temp1);
                                temphase = indi->phases[phasecounter0++];
                                val = (float) temp1->valuedness;
                                numIO = (float) temp1->numIO;

                                if (phase > 0) {
                                        maxcount = (int(pow(pow(val, numIO),2)));
                                        for (int k = 0; k < maxcount; k++)
                                                       temp1->gateMatrix1[k] = cuCmulf(
                                                                     temp1->gateMatrix1[k], temphase);
                                }
                                //continue with the cost
                                cost += gateArray[num_g]->Cost;
                                //multiply the gates
                                //using Kronecker multiplication
					
				Kronm(temp0, temp1, (int)(pow((float)(temp0->valuedness),(float)(temp0->numIO))), (int)(pow((float)(temp1->valuedness),(float)(temp1->numIO))), temp2, d_M4, d_M5, d_M6);

//                                tensorProduct3(temp2, temp0, temp1);
                                temp0 = ccGate(temp2, temp0);
                        }
                        if (initTime) {
				/* Initialize the device matrices with the host matrices */
				(cudaMemcpy(d_M2, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
				initTime = false;
			} else {
				cudaMemcpy(d_M1, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
				if (even_odd == 0){
					cuBlasMMulti(w, d_M2, d_M1, d_M3);
					even_odd = 1;
				}else {	
					cuBlasMMulti(w, d_M3, d_M1, d_M2);
					even_odd = 0;
				}
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


        //check if phase is defined and multiply the whole matrix by it
        if (phase > 0) {
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
	delete [] gate_rep;
        return myfinalG;
}
#endif
/****************************************
 * returns the matrix representing the circuit of a given individual
 ****************************************/
qGate* GA::computeStringCUDAMatrix(string s, int valuedness, int io) {

	qGate *temp0, *temp1, *temp2, *myfinalG, *temp3;
	bool initTime = true;
	int cost = 0;
	int csegs = 0;
	int num_g;
	float numIO, val, maxcount;
	char *gate_rep = new char[4];
	cuComplex temphase;
	temp0 = new qGate;
	initGate(temp0, resultnum, valuedness,1);
	temp1 = new qGate;
	initGate(temp1, resultnum, valuedness,1);
	temp2 = new qGate;
	initGate(temp2, resultnum, valuedness,1);
	myfinalG = new qGate;
	initGate(myfinalG, resultnum, valuedness,1);

	int w = (int(pow((float) valuedness, (float) io)));
	int N = (int(pow(w,2)));
	cuComplex alpha = make_cuFloatComplex(1, 0);
	cuComplex beta = make_cuFloatComplex(0, 0);

	int even_odd = 0;

	//int the string counters
	int begin = 0;
	int end = s.find("p", begin + 1);


	/* Allocate device memory for the matrices */
	temp2->valuedness = temp0->valuedness;
	zeroMatrix(temp2, temp0->numIO);
	/* Initialize the device matrices with the host matrices */
	//cout<<cudaGetErrorString
	(cudaMemcpy(d_M1, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_M2, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_M3, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));

	while (end > 0) {
		if ((end - begin) > 1) {
			csegs++;
			//get the first  gate of the Segment and copy it
			getGateRepfromStr(begin+1,s, gate_rep);
			num_g = getGate(gate_rep);
			temp0 = ccGate(gateArray[num_g], temp0);
//			printqGate(temp0);
			numIO = (float) temp0->numIO;
			val = (float) temp0->valuedness;
			//continue with the cost
			begin++;
			cost += temp0->Cost;
			//get the next one
			begin += 3;
			for (int b = begin; b < end; b+=REPSIZE) {
				getGateRepfromStr(b,s, gate_rep);
				num_g = getGate(gate_rep);
				temp1 = ccGate(gateArray[num_g], temp1);
				temp3 = gateArray[num_g];
				val = (float) valuedness;
				numIO = (float) temp1->numIO;

				if (phase > 0) {
					maxcount = (int(pow(pow(val, numIO),2)));
					for (int k = 0; k < maxcount; k++)
						temp1->gateMatrix1[k] = cuCmulf(
							temp1->gateMatrix1[k], temphase);
				}
				//continue with the cost
				cost += gateArray[num_g]->Cost;
				//multiply the gates using Kronecker multiplication
//				printqGate(temp3);
				Kronm(temp0, temp3, (int)(pow((float)(temp0->valuedness),(float)(temp0->numIO))), (int)(pow((float)(temp3->valuedness),(float)(temp3->numIO))), temp2, d_M4, d_M5, d_M6);
				temp0 = ccGate(temp2, temp0);
			}
//			printqGate(temp0);
			if (initTime) {
				//if this is the end of the first segment fill the output value
				//cout<<cudaGetErrorString
				(cudaMemcpy(d_M2, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
				initTime=false;
			} else {
				cudaMemcpy(d_M1, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
				if (even_odd == 0){
					Mulmnoset(w, d_M2, d_M1, d_M3);
					even_odd = 1;
				} else {
					Mulmnoset(w, d_M3, d_M1, d_M2);
					even_odd = 0;
				}
			}
		}
		//move to the next segment
		begin = s.find("p", end);
		end = s.find("p", begin + 1);
	}

	/* Read the result back */
	if (even_odd == 0)
		cudaMemcpy(myfinalG->gateMatrix1, d_M2, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
	else 
		cudaMemcpy(myfinalG->gateMatrix1, d_M3, N*sizeof(cuComplex), cudaMemcpyDeviceToHost);
		


	//actualize the segment number
	val = ((float) myfinalG->valuedness);
	numIO = ((float) myfinalG->numIO);

	destroyGate(temp0);
	delete temp0;
	destroyGate(temp1);
	delete temp1;
	destroyGate(temp2);
	delete temp2;
	delete [] gate_rep;

	return myfinalG;
}

/****************************************
 * returns the matrix representing the circuit of a given individual
 ****************************************/
qGate* GA::computeCUDAMatrix(Individual *indi) {

	qGate *temp0, *temp1, *temp2, *myfinalG, *temp3;
	int phasecounter0 = 0;
	bool initTime = true;
	int cost = 0;
	int csegs = 0;
	int num_g;
	float numIO, val, maxcount;
	char *gate_rep = new char[4];
	cuComplex temphase;

	temp0 = new qGate;
	initGate(temp0, resultnum, indi->valuedness,1);
	temp1 = new qGate;
	initGate(temp1, resultnum, indi->valuedness,1);
	temp2 = new qGate;
	initGate(temp2, resultnum, indi->valuedness,1);
	myfinalG = new qGate;
	initGate(myfinalG, resultnum, indi->valuedness,1);

	int w = (int(pow((float) indi->valuedness, (float) indi->ioNumber)));
	int N = (int(pow(w,2)));
	cuComplex alpha = make_cuFloatComplex(1, 0);
	cuComplex beta = make_cuFloatComplex(0, 0);

	int even_odd = 0;

	//int the string counters
	int begin = 0;
	int end = indi->my_string.find("p", begin + 1);
	indi->Cost = 0;
	indi->segmentNumber = 0;

	/* Allocate device memory for the matrices */
	temp2->valuedness = temp0->valuedness;
	zeroMatrix(temp2, temp0->numIO);
	/* Initialize the device matrices with the host matrices */
	//cout<<cudaGetErrorString
	(cudaMemcpy(d_M1, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_M2, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_M3, temp2->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));

	while (end > 0) {
		if ((end - begin) > 1) {
			csegs++;
			//get the first  gate of the Segment and copy it
			getGateRepfromStr(begin+1,indi->my_string, gate_rep);
			num_g = getGate(gate_rep);
			temp0 = ccGate(gateArray[num_g], temp0);
			numIO = (float) temp0->numIO;
			val = (float) temp0->valuedness;
			if (phase > 0) {
				temphase = indi->phases[phasecounter0++];
				maxcount = (int(pow(pow(val, numIO),2)));
				for (int k = 0; k < maxcount; k++)
					temp0->gateMatrix1[k] = cuCmulf(
						temp0->gateMatrix1[k], temphase);
			}
			//continue with the cost
			begin++;
			cost += temp0->Cost;
			//get the next one
			begin += 3;
			for (int b = begin; b < end; b+=REPSIZE) {
				getGateRepfromStr(b,indi->my_string, gate_rep);
				num_g = getGate(gate_rep);
				temp1 = ccGate(gateArray[num_g], temp1);
				temp3 = gateArray[num_g];
				if (phase > 0)
					temphase = indi->phases[phasecounter0++];
				val = (float) temp1->valuedness;
				numIO = (float) temp1->numIO;

				if (phase > 0) {
					maxcount = (int(pow(pow(val, numIO),2)));
					for (int k = 0; k < maxcount; k++)
						temp1->gateMatrix1[k] = cuCmulf(
							temp1->gateMatrix1[k], temphase);
				}
				//continue with the cost
				cost += gateArray[num_g]->Cost;
				//multiply the gates using Kronecker multiplication
				Kronm(temp0, temp3, (int)(pow((float)(temp0->valuedness),(float)(temp0->numIO))), (int)(pow((float)(temp3->valuedness),(float)(temp3->numIO))), temp2, d_M4, d_M5, d_M6);
//				tensorProduct3(temp2, temp0, temp3);
				temp0 = ccGate(temp2, temp0);
			}
			if (initTime) {
				//if this is the end of the first segment fill the output value
				//cout<<cudaGetErrorString
				(cudaMemcpy(d_M2, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice));
				initTime=false;
			} else {
				cudaMemcpy(d_M1, temp0->gateMatrix1, N*sizeof(cuComplex), cudaMemcpyHostToDevice);
				if (even_odd == 0){
					Mulmnoset(w, d_M2, d_M1, d_M3);
					even_odd = 1;
				} else {
					Mulmnoset(w, d_M3, d_M1, d_M2);
					even_odd = 0;
				}
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
		


	//check if phase is defined and multiply the whole matrix by it
	if (phase > 0) {
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
	delete [] gate_rep;

/*		cout << "Matrix: w: "<<w<<endl;;
		for (int y = 0; y < w; y++){
			for (int x = 0; x < w; x++){
				cout<<"("<<cuCrealf(myfinalG->gateMatrix1[x+w*y])<<","<<cuCimagf(myfinalG->gateMatrix1[x+w*y])<<")";
			}
			cout<<endl;
		}
		cout<<endl;

exit(0);
*/
	return myfinalG;
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
	float err,maxcount = 0;
	int l = (int(pow(val, numIO)));
	int x = (int(pow(val, numIO)));
	int y = (int(pow(val, numIO)));
	cuComplex inter[(int(pow(val, numIO)))];
	cuComplex inter_b[(int(pow(val, numIO)))];
	cuComplex logicresult[(int(pow(val, numIO)))];
	cuComplex stb[(int(pow(val, numIO)))];
	cuComplex expectations[(int(pow(val, numIO) * val))], alphas;
	cuComplex resultingTotalError[(int(pow(val, numIO)))][measurement * (int) val];
	cuComplex resultingTotal[(int(pow(val, numIO)))][measurement * (int) val];
	complex<float> temp_val;

	Individual *indi = ind;
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;

	//threaded code
	long time = clock();
	//accelerated matrix computation
#ifdef __CUBLAS__
	myfinalG = computeCUBLASMatrix(indi);
#else
	myfinalG = computeCUDAMatrix(indi);
#endif
	(cudaMemcpy(d_MV, myfinalG->gateMatrix1, l*l*sizeof(cuComplex), cudaMemcpyHostToDevice));
	numIO = (float) myfinalG->numIO;
	//check if phase is defined and multiply the whole matrix by it
	if (phase > 0) {
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
		alphas = make_cuFloatComplex(real(temp_val), imag(temp_val));

		//if the measurement generated at least some output
		if (!(cuCrealf( alphas) == 0 && cuCimagf( alphas) == 0)){
			//if the output qubit was measured in 0, normalize the vector
			cudaMemcpy(&d_Value[0], &alphas, sizeof(cuComplex),cudaMemcpyHostToDevice);
			Divvnoset(l, d_VI, d_Value, d_VO);
			cudaMemcpy(d_VI, d_VO, sizeof(cuComplex),cudaMemcpyDeviceToDevice);
		} else {
			//if the measurement was not succesful, use the output of the circuit 
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
		alphas = make_cuFloatComplex(real(temp_val), imag(temp_val));

		cudaMemcpy(d_VI, d_VO, l*sizeof(cuComplex), cudaMemcpyDeviceToDevice);

		//if the input qubit was measured in 0, normalize the vector
		if (!(cuCrealf( alphas) == 0 && cuCimagf( alphas) == 0)){
			cudaMemcpy(&d_Value[0], &alphas, sizeof(cuComplex),cudaMemcpyHostToDevice);
			Divvnoset(l, d_VI, d_Value, d_VI);

		}

	}

	int m;
	err = 0;
	//calculate the overall error
	for (int c = 0; c < measurement; c++) {
		m = measured_desired_output_records;
		inter[0] = make_cuFloatComplex(0, 0);
		for (int e = 0; e < measured_desired_output_records; e++) {
			//inter[0] = cuCmulf(inter[0], resultingTotal[e][0]);
			inter[0] = cuCaddf(inter[0], make_cuFloatComplex(cuCrealf(resultingTotalError[e][(int) val]), 0));
		}
		err += cuCrealf(cuCdivf(inter[0], make_cuFloatComplex(m, 0)));
	}
	err /= measurement;
	indi->Error = err;

	//generate fitness value
	makeFitness(indi);
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
 * Right now it is doing measurement in increasing order
 * of all defined bits to measure
 ****************************************/
void* GA::doMeasureFitness(Individual *ind, bool display) {

	//define temporal variables for the results and the manupaltion
	//of the measurement
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
	cuComplex resultingTotalError[(int(pow(val, numIO)))][measurement
	                                                      * (int) val];
	cuComplex resultingTotalAllVals;
	cuComplex resultingTotal[(int(pow(val, numIO)))][measurement * (int) val];
	double penalty = 0.0;

	Individual *indi = ind;
	indi->Error = 0;
	indi->fitness = 0;
	indi->Cost = 0;

//	cout << indi->my_string<<endl;
	long time = clock();
#ifdef __CUBLAS__
	myfinalG = computeCUBLASMatrix(indi);
#else
	myfinalG = computeCUDAMatrix(indi);
#endif
	(cudaMemcpy(d_MV, myfinalG->gateMatrix1, l*l*sizeof(cuComplex), cudaMemcpyHostToDevice));
	numIO = (float) myfinalG->numIO;
	//check if phase is defined and multiply the whole matrix by it
	if (phase > 0) {
		maxcount = (int(pow(pow(val, numIO),2)));
		if (cuCrealf(indi->phase) != 0 && cuCimagf(indi->phase) != 0) {
			for (int k = 0; k < maxcount; k++)
					myfinalG->gateMatrix1[k] = cuCmulf(
							myfinalG->gateMatrix1[k], indi->phase);
		}
	}

//cout<<measurement << ' ' << val << ' ' << measured_desired_output_records<<endl;
	//for single qubit this should always be only 0, i.e. <1
	for (int m = 0; m < measurement; m++) {
		for (int r = 0; r < val; r++) {
			measures[r] = measurements[(measurementQBits[m] * (int) val) + r];
		}
		//for all input states l
		vectorInitZero(l, d_VI);
		for (int k = 0; k < measured_desired_output_records; k++) {

			//empty all intermediary variables
			for (int w = 0; w < l; w++)
				sta[w] = inter[w] = logicresult[w] = make_cuFloatComplex(0, 0);
			//set the k-th input state - minterm
			sta[measured_desired_output_records_idx[k]] = make_cuFloatComplex(1, 0);
			vectorSetIterState(measured_desired_output_records_idx[k], l, sta[measured_desired_output_records_idx[k]], d_VI);
			//cudaMemcpy(d_VI, sta, l*sizeof(cuComplex), cudaMemcpyHostToDevice);

//			cout << "Input: ";
//			for (int w = 0; w < l; w++)
//				cout<<"("<<cuCrealf(sta[w])<<","<<cuCimagf(sta[w])<<")";
//			cout<<endl;
//			cudaMemcpy(sta, d_VI, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);
//
//			cout << "Input: ";
//			for (int w = 0; w < l; w++)
//				cout<<"("<<cuCrealf(sta[w])<<","<<cuCimagf(sta[w])<<")";
//			cout<<endl<<endl;


			//propagate the state through the matrix
			Mulvnoset(l, d_MV, d_VI, d_VO);
			
//			for (int i = 0; i < l; i++) {
//				for (int j = 0; j < l; j++) {
//					logicresult[i] = cuCaddf(logicresult[i], cuCmulf(sta[j],
//							myfinalG->gateMatrix1[i*l+j]));
//				}
//			}
//
//			cout << "Output: ";
//			for (int w = 0; w < l; w++)
//				cout<<"("<<cuCrealf(logicresult[w])<<","<<cuCimagf(logicresult[w])<<")";
//			cout<<endl<<endl;

			cudaMemcpy(logicresult, d_VO, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);
//			cout << "Output: ";
//			for (int w = 0; w < l; w++)
//				cout<<"("<<cuCrealf(logicresult[w])<<","<<cuCimagf(logicresult[w])<<")";
//			cout<<endl;


			
			//measure for each state of the single given qubit
			resultingTotalAllVals = make_cuFloatComplex(0, 0);
			for (int r = 0; r < val; r++) {
				for (int w = 0; w < l; w++)
					inter[w] = make_cuFloatComplex(0, 0);

				(cudaMemcpy(d_M1, measures[r]->gateMatrix1, l*l*sizeof(cuComplex), cudaMemcpyHostToDevice));
				cudaMemcpy(d_VV_M, logicresult, l*sizeof(cuComplex), cudaMemcpyHostToDevice);
				Mulvnoset(l, d_M1, d_VV_M, d_VO);
				cudaMemcpy(inter, d_VO, l*sizeof(cuComplex), cudaMemcpyDeviceToHost);

/*				for (int i = 0; i < l; i++) {
					//rows
					for (int j = 0; j < l; j++) {
						inter[i] = cuCaddf(inter[i], cuCmulf(logicresult[j],
								measures[r]->gateMatrix1[i*l+j]));
					}
				}
*/
/*			cout << "O: ";
			for (int w = 0; w < l; w++)
				cout<<"("<<cuCrealf(inter[w])<<","<<cuCimagf(inter[w])<<")";
			cout<<endl;

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
				
//				cout <<"indexes :" <<measured_desired_output_records_idx[k]<<"  "<<(int)((m* (int) val ) + r) <<":: "<< cuCrealf(measureexpected[measured_desired_output_records_idx[k]][(m* (int) val ) + r]) << ':';
				resultingTotalAllVals =  cuCaddf(resultingTotalAllVals,expectations[(m * (int) val + r)]);
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
//				for (int i = 0; i < l; i++) 
//					cout<<k<< ' ' << m << ' ' <<cuCrealf(resultingTotalError[measured_desired_output_records_idx[k]][m * (int) val + r])<< ' ';

			}
//			cout<<endl;
			if (cuCrealf(resultingTotalAllVals) == 0)
				penalty++;	
		}
	}

//cout<<"next"<<endl;
	int m;
	err = 0;
	for (int c = 0; c < measurement; c++) {
		m = l * (int) val;
		inter[0] = make_cuFloatComplex(0, 0);
		for (int e = 0; e < measured_desired_output_records; e++) {
			for (int r = 0; r < val; r++) {
				if (!(cuCrealf(measureexpected[measured_desired_output_records_idx[e]][(int) val * c + r]) == 0
						&& cuCimagf(measureexpected[measured_desired_output_records_idx[e]][(int) val * c + r])	== -1)){
					inter[0] = cuCaddf(inter[0], make_cuFloatComplex(cuCrealf(
							resultingTotalError[measured_desired_output_records_idx[e]][c * (int) val + r]), 0));
				} else {
					m -= (int) val;
				}
//				cout<<e<<' '<< c << ' ' <<cuCrealf(resultingTotalError[measured_desired_output_records_idx[e]][c * (int) val + r])<< ' ';
			}
//			cout<<endl ;
		}
//		cout<<m<< ' ' <<cuCrealf(inter[0])<< ' ' <<cuCimagf(inter[0])<<endl;
		err += cuCrealf(cuCdivf(inter[0], make_cuFloatComplex(m, 0)));
	}
	//cout<<cuCrealf(expectations[(m * (int) val + r)])<< ' ' <<cuCimagf(expectations[(m * (int) val + r)])<< endl;
//	cout<<err<< endl;
//exit(0);

/*		cout<< "Circuit Matrix:"<<endl;
		for (int i = 0; i < l; i++) {
			for (int j = 0; j < l; j++) {
				cout<<"("<<cuCrealf(myfinalG->gateMatrix1[i*l+j])<<","<<cuCimagf(myfinalG->gateMatrix1[i*l+j])<<")";
			}
			cout<<endl;
		}
		cout<<endl;
		cout<<endl;
*/
	err /= measurement;
	indi->Error = err+(penalty/measured_desired_output_records);

	//generate fitness value
	makeFitness(indi);
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
//	cout<<"CUBLAS fitness done: "<<(float)(time/CLOCKS_PER_SEC)<<", for "<<indi->ioNumber<<" qubits, fitness "<<indi->fitness<< ' ' <<indi->Error<<endl;
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
	int maxcount = 0;
	int mes = 0;
	cuComplex inter0[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex inter1[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex logicresult[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex sta[(int(pow((float) 2, (float) finalGate->numIO)))];
	cuComplex expectation0, expectation1, alphas, betas;
	complex<float> expe_0, expe_1;
	cuComplex
	resultingTotalError[(int(pow((float) 2, (float) finalGate->numIO)))][measurement
	                                                                     * 2];
//	cuComplex
//	resultingTotal[(int(pow((float) 2, (float) finalGate->numIO)))][measurement
//	                                                                * 2];
	cuComplex inter, expect0, expect1;
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
	myfinalG = computeCUDAMatrix(indi);
	rc = pthread_mutex_unlock(&mtx);

	if (phase > 0) {
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
		= complex<float> (cuCrealf(expectation0),
				cuCimagf(expectation0));
		expe_0 = sqrt(expe_0);
		expectation0 = make_cuFloatComplex(real(expe_0), imag(expe_0));

		for (int i = 0; i < l; i++) {
			inter0[i] = cuCdivf(inter0[i], expectation0);
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
		= complex<float> (cuCrealf(expectation1),
				cuCimagf(expectation1));
		expe_1 = sqrt(expe_1);
		expectation0 = make_cuFloatComplex(real(expe_1), imag(expe_1));
		for (int i = 0; i < l; i++) {
			inter1[i] = cuCdivf(inter1[i], expectation1);
			betas = cuCaddf(betas, inter1[i]);
		}
		//--------------------------------------------------------------------------------------
		alphas = expectation0;
		betas = expectation1;

		//		cout<<"alpha: "<<alphas<<" + beta: "<<betas<<endl;
		//calculate the total state
		//Total State = M(0)+M(1)State/Measured(0)+Measured(1)
		mes = 0;
//		resultingTotal[k][2* mes ] = expectation0;
//		resultingTotal[k][2* mes + 1] = expectation1;
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
		//		for (int c = 0; c < measurement; c++){
		if (!(cuCimagf(measureexpected[0][e]) == 1 && cuCimagf(
				measureexpected[1][e]) == 1)) {
			//get the error over both possible measurements of the state - for 0 and 1 of the first qubit
			inter = cuCaddf(inter, cuCdivf(cuCaddf(resultingTotalError[e][0],
					resultingTotalError[e][1]), make_cuFloatComplex(1, 0)));
			//expecting the desired higher probabilities 1
			m++;
		}
		//		}
	}
	//indi->Error /= measurement;
	indi->Error = cuCrealf(cuCdivf(inter, make_cuFloatComplex(m, 0)));
	indi->Cost = (exp(-pow((divider - indi->Cost), 2)));

	//generate fitness value
	makeFitness(indi);
	//	cout<<"Fitness: "<<indi->fitness<<endl;
	destroyGate(myfinalG);
	delete (myfinalG);
	return 0;
}
