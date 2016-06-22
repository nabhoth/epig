
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
/****************************************
 * These are various functions for cstGA
 * returns the complex conjugate of a complex number 
 * **************************************/
int getGatesfromString(string indv) {
	int chars = 0;
	for (int a = 0; a < indv.length(); a++)
		if (indv[a] >= '0' && indv[a] <= 'z' && indv[a] != 'p')
			chars++;
	return chars;
}
/****************************************
 * These are various functions for cstGA
 * **************************************/
bool compareGate(qGate A, qGate B) {
	int val = A.valuedness;
	int maxcount = (int(pow(pow((float) val, (float) A.numIO),2)));
	if (A.numIO != B.numIO)
		return false;
	else {
		for (int i = 0; i < maxcount; i++)
			if (cuCrealf(A.gateMatrix1[i]) != cuCrealf(B.gateMatrix1[i]) || cuCimagf(A.gateMatrix1[i]) != cuCimagf(B.gateMatrix1[i])) return false;

	}
	return true;

}
/****************************************
 * Copies the content from one gate to another pointer
 * **************************************/
qGate* ccGate(qGate *Circ, qGate *Dest) {
	int val = Circ->valuedness;
	int maxcount = (int(pow(pow((float) val, (float) Circ->numIO),2)));
	Dest->numIO = Circ->numIO;
	Dest->valuedness = Circ->valuedness;
	for (int i = 0; i < maxcount; i++){
			Dest->gateMatrix1[i] = Circ->gateMatrix1[i];
	}

	Dest->Cost = Circ->Cost;
	Dest->representation = Circ->representation;
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
	Circf.valuedness = Circ->valuedness;
	Circf.Cost = Circ->Cost;
	Circf.representation = Circ->representation;
	Circf.my_string = "";
	for (int a = 0; a < Circ->my_string.length(); a++) {
		Circf.my_string += Circ->my_string.at(a);
	}

	return Circf;
}
/****************************************
 * creates a new gate by copying it self
 * **************************************/

qGate* cGate(qGate *Circ) {
	qGate *Circf = new qGate();
	int val = Circ->valuedness;
	int maxcount = (int(pow(pow((float) val, (float) Circ->numIO),2)));
	Circf->numIO = Circ->numIO;
	Circf->valuedness = Circ->valuedness;
	for (int i = 0; i < maxcount; i++)
			Circf->gateMatrix1[i] = Circ->gateMatrix1[i];

	Circf->Cost = Circ->Cost;
	Circf->representation = Circ->representation;
	Circf->my_string = "";
	for (int a = 0; a < Circ->my_string.length(); a++) {
		Circf->my_string += Circ->my_string.at(a);
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
void vectorInitZero(int w, cuComplex *d_Vec){
	for (int m = 0; m < w; m++)
                d_Vec[m] = make_cuFloatComplex(0, 0);
}
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
	int val = A->valuedness = valuedness;
	A -> Cost = cost;
	int maxcount = (int)(pow(pow((float) valuedness, (float) resultnum),2.0));
	A->gateMatrix1 = new cuComplex[maxcount];
}

/****************************************
 * Deletes the gate matrix
 * **************************************/
void destroyGate(qGate *A){
	delete [] A->gateMatrix1;
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
	R->numIO = (A->numIO + B->numIO);
	R->representation = 'f';
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
	for (int a = 0; a < 2; a++){
		dst[a] = src[a];
	}
}
/**********************************
* compare the representation of one gate to another
 **********************************/
bool compareRep(char *src, char *dst){
	for (int a = 0; a < 2; a++){
		if (dst[a] != src[a]) return false;
	}
	return true;
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
}
