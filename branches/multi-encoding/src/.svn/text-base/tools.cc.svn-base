
/****************************************
 * These are various functions for cstGA
 * **************************************/


#include "EpiG.h"

/****************************************
 * These are various functions for cstGA
 * taken from the GSL library and changed for our purpose and data type 
 * **************************************/

double cuFloatComplex_arg (cuComplex z)
{                               /* return arg(z),  -pi < arg(z) <= +pi */
  double x = cuCrealf (z);
  double y = cuCimagf (z);

  if (x == 0.0 && y == 0.0)
    {
      return 0;
    }

  return atan2 (y, x);
}


double cuFloatComplex_abs (cuComplex z)
{                               /* return |z| */
  return hypot (cuCrealf (z), cuCimagf (z));
}


cuComplex cuFloatComplex_inverse (cuComplex a)
{                               /* z=1/a */
  double s = 1.0 / cuFloatComplex_abs (a);

  cuComplex z;
  z = make_cuFloatComplex ((cuCrealf (a) * s) * s, -(cuCimagf (a) * s) * s);
  return z;
}

double cuFloatComplex_logabs (cuComplex z)
{                               /* return log|z| */
  double xabs = fabs (cuCrealf (z));
  double yabs = fabs (cuCimagf (z));
  double max, u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }

  /* Handle underflow when u is close to 0 */

  return log (max) + 0.5 * log1p (u * u);
}

cuComplex cuFloatComplex_sub (cuComplex a, cuComplex b)
{                               /* z=a-b */
  double ar = cuCrealf (a), ai = cuCimagf (a);
  double br = cuCrealf (b), bi = cuCimagf (b);

  cuComplex z;
  z = make_cuFloatComplex (ar - br, ai - bi);
  return z;
}

cuComplex cuFloatComplex_div (cuComplex a, cuComplex b)
{                               /* z=a/b */
  double ar = cuCrealf (a), ai = cuCimagf (a);
  double br = cuCrealf (b), bi = cuCimagf (b);

  double s = 1.0 / cuFloatComplex_abs (b);

  double sbr = s * br;
  double sbi = s * bi;

  double zr = (ar * sbr + ai * sbi) * s;
  double zi = (ai * sbr - ar * sbi) * s;

  cuComplex z;
  z = make_cuFloatComplex (zr, zi);
  return z;
}



cuComplex cuFloatComplex_pow (cuComplex a, cuComplex b)
{                               /* z=a^b */
  cuComplex z;

  if (cuCrealf (a) == 0 && cuCimagf (a) == 0.0)
    {
      if (cuCrealf (b) == 0 && cuCimagf (b) == 0.0)
        {
          z = make_cuFloatComplex (1.0, 0.0);
        }
      else 
        {
          z = make_cuFloatComplex (0.0, 0.0);
        }
    }
  else if (cuCrealf (b) == 1.0 && cuCimagf (b) == 0.0) 
    {
      return a;
    }
  else if (cuCrealf (b) == -1.0 && cuCimagf (b) == 0.0) 
    {
      return cuFloatComplex_inverse (a);
    }
  else
    {
      double logr = cuFloatComplex_logabs (a);
      double theta = cuFloatComplex_arg (a);

      double br = cuCrealf (b), bi = cuCimagf (b);

      double rho = exp (logr * br - bi * theta);
      double beta = theta * br + bi * logr;

      z = make_cuFloatComplex (rho * cos (beta), rho * sin (beta));
    }

  return z;
}

cuComplex cuFloatComplex_sqrt (cuComplex a)
{                               /* z=sqrt(a) */
  cuComplex z;

  if (cuCrealf (a) == 0.0 && cuCimagf (a) == 0.0)
    {
      z = make_cuFloatComplex (0, 0);
    }
  else
    {
      double x = fabs (cuCrealf (a));
      double y = fabs (cuCimagf (a));
      double w;

      if (x >= y)
        {
          double t = y / x;
          w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + t * t)));
        }
      else
        {
          double t = x / y;
          w = sqrt (y) * sqrt (0.5 * (t + sqrt (1.0 + t * t)));
        }

      if (cuCrealf (a) >= 0.0)
        {
          double ai = cuCimagf (a);
          z = make_cuFloatComplex ( w, ai / (2.0 * w));
        }
      else
        {
          double ai = cuCimagf (a);
          double vi = (ai >= 0) ? w : -w;
          z = make_cuFloatComplex (ai / (2.0 * vi), vi);
        }
    }

  return z;
}
/*******************************************
 * extracts the number of gates
 * from this segment 
 *******************************************/
int getGatesNumber(string b){
	int counter = 0;
	int index = 0;
	while (b.at(index++) == 'p');
	while (b.at(index) != 'p'){
		counter++;
		index += 3;
	}
	return counter;
}
/****************************************
 * exatracts characters representing a gate 
 * from the current position of the 
 * string indv 
 * **************************************/
int getGatesfromString(string indv) {
	int chars = 0;
	for (int a = 0; a < indv.length(); a++)
		if (indv[a] != 'p')
			chars++;
	return chars;
}
/****************************************
 * Compares two gates
 * **************************************/
bool compareGate(qGate *A, qGate *B) {
	int maxcount = (int(pow((float) A->valuedness, (float) A->numIO)));
//	if (A->numIO != B->numIO)
//		return false;
//	if (A->valuedness != B->valuedness)
//		return false;
//	else {
		for (int i = 0; i < maxcount; i++)
			for (int j = 0; j < maxcount; j++){
				if(cuCabsf(A->gateMatrix1[i*maxcount+j]) != cuCabsf(B->gateMatrix1[i*maxcount+j]))
					return false;
				if ((cuCrealf(A->gateMatrix1[i*maxcount+j]) != cuCrealf(B->gateMatrix1[i*maxcount+j])) || (cuCimagf(A->gateMatrix1[i*maxcount+j]) != cuCimagf(B->gateMatrix1[i*maxcount+j]))){
					return false;
				}
			}

//	}
	return true;

}
/****************************************
 * Copies the content from one gate to another pointer
 * **************************************/
qGate* ccGate(qGate *Circ, qGate *Dest) {
	int val = Circ->valuedness;
	int maxcount = (int(pow(pow((float) val, (float) Circ->numIO),2)));
	Dest->numIO = Circ->numIO;
	Dest->realIO = Circ->realIO;
	Dest->restrictions_number = Circ->restrictions_number;
	for (int i = 0; i < Circ->restrictions_number; i++)
			Dest->restrictions[i] = Circ->restrictions[i];
	for (int j = 0; j < Circ->realIO; j++)
		Dest->connections[j] = Circ->connections[j];
	Dest->valuedness = Circ->valuedness;
	for (int r = 0; r < maxcount; r++){
			Dest->gateMatrix1[r] = Circ->gateMatrix1[r];
	}

	Dest->Cost = Circ->Cost;
	copyRep(Circ->representation,Dest->representation);
	Dest->my_string = "";
	for (int a = 0; a < Circ->my_string.length(); a++) {
		Dest->my_string += Circ->my_string.at(a);
	}

	return Dest;
}

/****************************************
 * Copies the content from one gate to another temporary one
 * **************************************/
qGate copyGate(qGate *Circ) {
	qGate Circf;
	int val = Circ->valuedness;
	int maxcount = (int(pow(pow((float) val, (float) Circ->numIO),2)));
	for (int i = 0; i < maxcount; i++)
			Circf.gateMatrix1[i] = Circ->gateMatrix1[i];

	Circf.numIO = Circ->numIO;
	Circf.realIO = Circ->realIO;
	Circf.restrictions_number = Circ->restrictions_number;
	for (int i = 0; i < Circ->restrictions_number; i++)
			Circf.restrictions[i] = Circ->restrictions[i];
	for (int i = 0; i < Circ->realIO; i++)
		Circf.connections[i] = Circ->connections[i];
	Circf.valuedness = Circ->valuedness;
	Circf.Cost = Circ->Cost;
	copyRep(Circ->representation, Circf.representation);
	//Circf.representation = Circ->representation;
	Circf.my_string = "";
	for (int a = 0; a < Circ->my_string.length(); a++) {
		Circf.my_string += Circ->my_string.at(a);
	}
	return Circf;
}
/****************************************
 * In development - calculates the error for a density matrix
 * **************************************/

void densityErrorMatrix(qGate *A, qGate *B, qGate *C) {
	if (A->numIO != B->numIO)
		return;
	C->numIO = A->numIO;
	cblasMatrixProduct(A, B, C);
}

/****************************************
 * Matrix Multiplication using CBLAS libs
 * **************************************/
void cblasMatrixProduct(qGate *A, qGate *B, qGate *C) {
	int val = A->valuedness;
	int rows = int(pow((float) val, (float) A->numIO));
	cuComplex alpha = make_cuFloatComplex(1, 0);
	cuComplex beta = make_cuFloatComplex(0, 0);
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows, &alpha, A->gateMatrix1, rows, B->gateMatrix1, rows, &beta, C->gateMatrix1, rows);
	C->numIO = A->numIO;
}


/***************************************
 *Genereates a random phase for a circuit based on simple rules
 ****************************************/
cuComplex generatePhase() {
	//calculate random phase of a circuit
	float nom = rand() % 90;
	float den = rand() % 100;
	if (den == 0)
		den = 1;
	float sign = rand() % 1;
	if (sign > 0)
		sign = 1;
	else
		sign = -1;
	float pi = 4* atan (1.0);
	//calculate complex exponent
	float r = exp(0.0);
	return make_cuFloatComplex(r * cos(sign * pi * nom / den), r * sin(sign * pi 	* nom / den));
}
/***************************************
 *Genereates a random phase for a whole circuit (for every unitary gate) based on simple rules
 ****************************************/
void generatePhases(Individual *I) {
	//calculate phase of the each circuit
	int phasemax = getGatesfromString(I->my_string);
	for (int a = 0; a < phasemax; a++) {
		float nom = rand() % 90;
		float den = rand() % 100;
		if (den == 0)
			den = 1;
		float sign = rand() % 1;
		if (sign > 0)
			sign = 1;
		else
			sign = -1;
		float pi = 4* atan (1.0);
		//calculate complex exponent
		float r = exp(0.0);
		I->phases[a] = make_cuFloatComplex(r * cos(sign * pi * nom / den), r * sin(sign * pi * nom / den));

	}
}

/****************************************
 * Fills the content of a matrix size resultnum to I
 * **************************************/
void zeroMatrix(qGate *A, int resultnum) {
	A -> numIO = resultnum;
	A -> Cost = 1;
	int val = A->valuedness;
	int maxcount = (int(pow(pow((float) val, (float) A->numIO),2)));
	for (int i = 0; i < maxcount; i++) {
		A->gateMatrix1[i] = make_cuFloatComplex(0, 0);
	}
}
/****************************************
 * Fills the content of a matrix size resultnum to -i 
 * **************************************/

void dcMatrix(qGate *A, int resultnum) {
	A -> Cost = 1;
	A -> numIO = resultnum;
	int val = A->valuedness;
	int maxcount = (int(pow(pow((float) val, (float) A->numIO),2)));
	for (int i = 0; i < maxcount; i++) {
		A->gateMatrix1[i] = make_cuFloatComplex(0, -1);

	}
}
/****************************************
 * Init a Vector to all 0
 * **************************************/
#ifdef __CUDA__
#else
void vectorInitZero(int w, cuComplex *d_Vec){
	for (int m = 0; m < w; m++)
                d_Vec[m] = make_cuFloatComplex(0, 0);
}
#endif
/****************************************
 * Dot Multiply two vectors
 * **************************************/
void vectorScalarMult(int w, cuComplex scalar, cuComplex *a_Vec, cuComplex *b_Vec){
	cuComplex result = make_cuFloatComplex(0, 0);
        for (int m = 0; m < w; m++)
                b_Vec[m] =  cuCmulf(a_Vec[m], scalar);
}
/****************************************
 * Sets the initialized gate to be the Identity gate
 * **************************************/
void initMatrix(qGate *A, int resultnum) {
	A -> numIO = resultnum;
	A -> Cost = 1;
	A -> restrictions_number = 0;
	int val = A->valuedness;
	int rows = (int(pow((float) val, (float) A->numIO)));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {
			(i == j) ? A->gateMatrix1[i*rows+j] = make_cuFloatComplex(1, 0): A->gateMatrix1[i*rows+j] = make_cuFloatComplex(0, 0);
		}
	}
}

/****************************************
 * Creates the gate and set the appropriate parameters
 * **************************************/
void initGate(qGate *A, int resultnum, int valuedness, int cost) {
	A -> numIO = resultnum;
	A -> realIO = resultnum;
	int val = A->valuedness = valuedness;
	A -> Cost = cost;
	int maxcount = (int)(pow(pow((float) valuedness, (float) resultnum),2.0));
	A->gateMatrix1 = new cuComplex[maxcount];
	A->representation = new char[REPSIZE+1];
	A->parentA = new char[REPSIZE+1];
	A->parentB = new char[REPSIZE+1];
}

/****************************************
 * Deletes the gate matrix
 * **************************************/
void destroyGate(qGate *A){
	delete [] A->gateMatrix1;
	delete [] A->representation;
	delete [] A->parentA;
	delete [] A->parentB;
	A->gateMatrix1 = NULL;
}

/**********************************
*RowMajor Kronecker Product
 **********************************/
void tensorProduct2(qGate *R, qGate *A, qGate *B) {
	int val = A->valuedness;
	R->valuedness = val;
	int index;
	int k = 0;
	int maxcount_a = (int(pow(pow((float) val, (float) A->numIO),2)));
	int maxcount_b = (int(pow(pow((float) val, (float) B->numIO),2)));
	int dim_a = (int(pow((float) val, (float) A->numIO)));
	int dim_b = (int(pow((float) val, (float) B->numIO)));
	int new_row_dim = (int(dim_a*dim_b));
	int p_a,d_a,p_b,d_b;
	for (int i = 0; i < maxcount_a; i++) {
		p_a = i % dim_a;
		d_a = i / dim_a;
		k = 0;
		for (int m = 0; m < maxcount_b / dim_b; m++)
		{
			for (int j = 0; j < dim_b; j++) {
				p_b = k / dim_b;
				index = (p_a*dim_b)+(d_a*dim_a*maxcount_b)+(p_b*new_row_dim)+j;
				R->gateMatrix1[index] = cuCmulf(A->gateMatrix1[i], B->gateMatrix1[k]);
				k++;
			}
		}
	}
	R->realIO = (A->realIO + B->realIO);
	R->restrictions_number = (A->restrictions_number + B->restrictions_number);
	R->numIO = (A->numIO + B->numIO);
}
/**********************************
*Column Major Kronecker Product
 **********************************/
void tensorProduct3(qGate *R, qGate *A, qGate *B) {
	int val = A->valuedness;
	R->valuedness = val;
	int index;
	int k = 0;
	int maxcount_a = (int(pow(pow((float) val, (float) A->numIO),2)));
	int maxcount_b = (int(pow(pow((float) val, (float) B->numIO),2)));
	int dim_a = (int(pow((float) val, (float) A->numIO)));
	int dim_b = (int(pow((float) val, (float) B->numIO)));
	int new_column_dim = (int(dim_a*dim_b));
	int p_a,d_a,p_b,d_b;
	for (int i = 0; i < maxcount_a; i++) {
		p_a = i % dim_a;
		d_a = i / dim_a;
		k = 0;
		for (int m = 0; m < maxcount_b / dim_b; m++)
		{
			for (int j = 0; j < dim_b; j++) {
				p_b = k / dim_b;
				index = (p_a*dim_b)+(d_a*dim_b*new_column_dim)+(p_b*new_column_dim)+j;
				R->gateMatrix1[index] = cuCmulf(A->gateMatrix1[i], B->gateMatrix1[k]);
				k++;
			}
		}
	}
	R->realIO = (A->realIO + B->realIO);
	R->restrictions_number = (A->restrictions_number + B->restrictions_number);
	R->numIO = (A->numIO + B->numIO);
	R->representation[0] = 'f';
	R->representation[1] = 'f';
	R->representation[2] = 'f';
	R->representation[REPSIZE] = '\0';
	R->my_string = "";
	R->Cost = 1;
}


/**********************************
 * substract two matrices, returns first argument
 **********************************/
qGate* matrixSubstr(qGate *A, qGate *B) {
	int val = A->valuedness;
	int maxcount = (int(pow(pow((float) val, (float) A->numIO),2)));
	for (int i = 0; i < maxcount; i++) {
		A->gateMatrix1[i] = cuCsubf(A->gateMatrix1[i], B->gateMatrix1[i]);
	}
	return A;
}


short sign(float c){
	if (c >= 0 )
		return 1;
	else return -1;
}
/**********************************
 * Compare two matrices
 **********************************/
bool matrixCompare(int w, int val, cuComplex* A, cuComplex* B){
	int maxcount = (int(pow(pow((float) val, (float) w),2)));
	for (int i = 0; i < maxcount; i++) {
		if ((cuCrealf(A[i]) != cuCrealf(B[i])) || (cuCrealf(A[i]) != cuCrealf(B[i])))
			return false;
	}
	return true;
}
/**********************************
 * Analyzes a generated merged sub-segment and creates a new matrix
 **********************************/
void segAnaAndMerge(char *subsegment){
	int c = 0;
	char g[REPSIZE];
	while (subsegment[c] != '!'){
		for(int m = 0; m < REPSIZE;m++){
			g[m] = subsegment[m+c];
		}
		g[REPSIZE] = '\0';

	}

}
/**********************************
* accepts a decimal integer and returns a binary coded string
* only for unsigned values
 **********************************/

void dec2bin(long decimal, char *binary, int radix)
{
	int k = 0, n = 0;
	int neg_flag = 0;
	int remain;
	int old_decimal; // for test
	char temp[80];
	// take care of negative input

	if (decimal < 0)
	{
		decimal = -decimal;
		neg_flag = 1;
	}
	do
	{
		old_decimal = decimal; // for test
		remain = decimal % 2;
		// whittle down the decimal number
		decimal = decimal / 2;
		// this is a test to show the action
		//printf("%d/2 = %d remainder = %d\n", old_decimal, decimal, remain);
		// converts digit 0 or 1 to character '0' or '1'
		temp[k++] = remain + '0';
	} while (decimal > 0);

//	if (neg_flag)
//		temp[k++] = '-'; // add - sign
//	else
//		temp[k++] = ' '; // space

	// reverse the spelling
	if (k < radix)
		while (n < (radix - k))
			binary[n++] = '0';
	while (k >= 0)
		binary[n++] = temp[--k];

	binary[n-1] = 0; // end with NULL
}

/**********************************
* copy the representation of one gate to another
 **********************************/
void copyRep(char *src, char *dst){
	for (int a = 0; a < REPSIZE; a++){
		dst[a] = src[a];
	}
	dst[REPSIZE] = '\0';
}
/**********************************
* copy the representation of one gate to another string
 **********************************/
void copyRepUnchecked(string dst, char *src, int indx){
	for (int a = 0; a < REPSIZE; a++){
		dst[indx+a]= src[a];
	}
}
/**********************************
* compare the representation of one gate to another
 **********************************/
bool compareRep(char *src, char *dst){
	for (int a = 0; a < REPSIZE; a++){
		if (dst[a] != src[a]) return false;
	}
	return true;
}
/**********************************
* increment the representation with the provided seed
 **********************************/
void initRep(char *dst){
	for (int a = 0; a < REPSIZE; a++){
		dst[a] = char(CHARMIN);
	}
	dst[REPSIZE] = '\0';
}/**********************************
* increment the representation with the provided seed
 **********************************/
void initRep(char *dst, int seed){
	for (int a = 0; a < REPSIZE; a++){
		dst[a] = char(seed);
	}
	dst[REPSIZE] = '\0';
}

/**********************************
* increment the representation within the allowed
* range and skip the 'p' character
 **********************************/
void incRep(char *dest){
	if (dest[0] < CHARMAX){
		dest[0] = char(dest[0]+1);
		if (dest[0] == 'p') dest[0] = char(dest[0]+1);
	} else {
		dest[0] = CHARMIN;
		if (dest[1] < CHARMAX){
			dest[1] = char(dest[1]+1);
			if (dest[1] == 'p') dest[1] = char(dest[1]+1);
		} else {
			dest[1] = CHARMIN;
			if (dest[2] < CHARMAX){
				dest[2] = char(dest[2]+1);
				if (dest[2] == 'p') dest[2] = char(dest[2]+1);
			} else {
				cout<<"Cannot incerement anymore the gate counter, Exiting!!!"<<endl;
			}
		}
	}
	dest[REPSIZE] = '\0';
}
/**********************************
* increment the representation within the allowed
* range and skip the 'p' character parametrized
 **********************************/
void incRep(char *dest, int min, int max){
	if (dest[0] < max){
		dest[0] = char(dest[0]+1);
		if (dest[0] == 'p') dest[0] = char(dest[0]+1);
	} else {
		dest[0] = min;
		if (dest[1] < max){
			dest[1] = char(dest[1]+1);
			if (dest[1] == 'p') dest[1] = char(dest[1]+1);
		} else {
			dest[1] = min;
			if (dest[2] < max){
				dest[2] = char(dest[2]+1);
				if (dest[2] == 'p') dest[2] = char(dest[2]+1);
			} else {
				cout<<"Cannot incerement anymore the gate counter, Exiting!!!"<<endl;
			}
		}
	}
	dest[REPSIZE] = '\0';
}
/**********************************
* checks if the restrictions of this gate match
* the provided wire number
 **********************************/
bool checkRestrictions(int wire, qGate *g){
	if (g->restrictions[0] == -1)
		return true;
	for (int m = 0; m < g->restrictions_number; m++){
		if (wire == g->restrictions[m])
			return true;		
	}
	return false;
}
/**********************************
* return the beginning index of a gate in a genome
 **********************************/
int getGatePos(int pos, string s){
	int end = s.find('p',pos);
	int begin = s.rfind('p',pos)+1;
	int real_pos = pos-begin;
	int pos_div = real_pos/3;
	int pos_mod = real_pos%3;
	return pos-pos_mod;
}

/**********************************
* return the pointer to the gate encoding
 **********************************/
int getGateRepfromStr(int pos, string str, char *representation){
	int lpos = pos;
	while (str.at(lpos) == 'p') {lpos++; if (lpos >= str.length()) return -1;}
	if (lpos >= (str.length()-1)) return -1;
	representation[0] = str.at(lpos);
	representation[1] = str.at(lpos+1);
	representation[2] = str.at(lpos+2);
	representation[REPSIZE] = '\0';
	return 0;	
}

bool isIn(int val, int arr[], int l){
	for (int m = 0; m < l; m++)
		if (val == arr[m])
			return true;
	return false;
}

int isInAt(int val, int arr[], int l){
	for (int m = 0; m < l; m++)
		if (val == arr[m])
			return m;
	return -1;
}

void printqGate(qGate *gate){
	cout<<" IO: "<<gate->numIO<<endl;
	int n = (int)(pow((float)(gate->valuedness),(float)(gate->numIO)));
	for(int j = 0; j < n; j ++){
		for (int i = 0; i < n; i ++){
			cout<<cuCrealf(gate->gateMatrix1[i+j*n])<<','<<cuCimagf(gate->gateMatrix1[i+j*n])<<' ';
		}
		cout<<endl;
	}
	cout<<"----"<<endl;
}

int getSegmentNumber(string circuit){
	int count = 0;
	for (int a = 0; a < circuit.length(); a++)
		if (circuit.at(a) == 'p') count++;
	return count/2;

}
