
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
		if (indv[a] != 'p')
			chars++;
	return chars;
}
/****************************************
 * Matrix Multiplication using CBLAS libs
 * **************************************/
void cblasMatrixProduct(aGate *A, aGate *B, aGate *C) {
	int val = A->valuedness;
//	int rows = int(pow((float) val, (float) A->numIO));
//	cuComplex alpha = make_cuFloatComplex(1, 0);
//	cuComplex beta = make_cuFloatComplex(0, 0);
//	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows, &alpha, A->gateMatrix1, rows, B->gateMatrix1, rows, &beta, C->gateMatrix1, rows);
//	C->numIO = A->numIO;
}


#ifdef __CUDA__
#else
void vectorInitZero(int w, double *d_Vec){
	for (int m = 0; m < w; m++)
                d_Vec[m] = 0.0;
}
#endif
/****************************************
 * Copies the content from one gate to another pointer
 * **************************************/
aGate* ccGate(aGate *Circ, aGate *Dest) {
	int val = Circ->valuedness;
        Dest->numI = Circ->numI;
	Dest->coordI = new int[Dest->numI];
        Dest->numO = Circ->numO;
	Dest->coordO = new int[Dest->numO];
        Dest->Cost = Circ->Cost;
        Dest->length_LU = Circ->length_LU;
        Dest->valuedness = Circ->valuedness;
	Dest->inLU = new int[Dest->length_LU];
	Dest->outLU = new int[Dest->length_LU];
        Dest->representation = new char[REPSIZE];
	copyRep(Circ->representation, Dest->representation);
        for (int i = 0; i < Dest->length_LU; i++){
		Dest->inLU[i] = Circ->inLU[i];
		Dest->outLU[i] = Circ->outLU[i];
        }
        Dest->my_string = "";
        for (int a = 0; a < Circ->my_string.length(); a++) {
                Dest->my_string += Circ->my_string.at(a);
        }

        return Dest;
}
/****************************************
 * Generate individual indexes for the gate inputs and outputs
 * **************************************/
void generateGateIndexes(aGate *gate, int current_index){
	for(int m = 0; m < gate->numI; m++){
		gate->coordI[m] = current_index++;
	}
	for(int n = 0; n < gate->numO; n++){
		gate->coordO[n] = current_index++;
	}
}

/****************************************
 * Returns true if all the inputs to a gate have been properly conencted
 * **************************************/
bool gateInputsSat(aConnection **connections, int cl, aGate *gate, bool doubled){
	int countr = 0;
	int con_in_idx = 0;
	int n = 0;
	while(con_in_idx <= gate->coordO[gate->numO-1]){ 
		if (n >= cl) break;
		con_in_idx = connections[n]->input_index;
		for (int m = 0; m < gate->numI; m++){
			if (connections[n]->output_index == gate->coordI[m]){
				countr++;
			}
		}
		if (countr >= gate->numI)
			return true;
		n++;
	}
	return false;
}

/****************************************
 * Returns all non connected connections
 * **************************************/
int getAvailableConn(aConnection **connections, int count){
	int res = 0;
	for (int a = 0; a < count; a++)
		if (connections[a]->output_index == -1)
			res++;
	return res;
}

/****************************************
 * connect the given gates by connections
 * **************************************/
void generateConnections(aConnection **connections, int count, aGate *gate, bool doubled, bool cascade){
	int counter = 0;
	int con_idx = 0;
	int con_idx_max = 0;
	bool c_required;
	int counter_required;
	int n = 0;
	int c = 0;
	int *required  = new int[gate->numI];
	int *available_cox = new int[MAXNUMOFWIRES];
//	cout<<"oth gates connex: "<<gate->coordI[0]<<endl;
	while (connections[con_idx]->input_index < gate->coordI[0])
	{
		if (connections[con_idx]->output_index == -1){
			available_cox[con_idx_max++] = con_idx;
//			cout<<"Avail: "<<con_idx<<" ";
		}
		con_idx++;
	}
	n = 0;
	counter_required = 0;
	for (int y = 0; y < gate->numI; y++){
		c_required = true;
		for (int x = 0; x < count; x++)
			if (gate->coordI[y] == connections[x]->output_index)
				c_required = false;
		if (c_required){
			counter_required++;
			required[n++] = gate->coordI[y];
		}
	}
	n = 0;
	if (doubled){
		while(counter < counter_required){ 
			do {
				c = rand()%con_idx_max;
				n = rand()%counter_required;
			}while (required[n] == -1 || available_cox[c] == -1 || (gate->negated_vars != connections[available_cox[c]]->negated));
			connections[available_cox[c]]->output_index = required[n];
			required[n] = -1;
			available_cox[c] = -1;
			counter++;		
		}
	} else {
		while(counter < counter_required){ 
			do {
				c = rand()%con_idx_max;
				n = rand()%counter_required;
			}while (required[n] == -1 || available_cox[c] == -1);
			connections[available_cox[c]]->output_index = required[n];
			required[n] = -1;
			available_cox[c] = -1;
			counter++;		
		}
	}
}

/****************************************
 * determines if two connections are from the same gate
 * **************************************/
bool  areFromSameGate(aCircuit *circuit, aConnection *acon, aConnection *bcon){
	int in_a, in_b, out_a, out_b;
	for (int m = 0; m < circuit->numGates; m++){
		for(int n = 0; n < circuit->gates[m]->numI; n++){
			if (circuit->gates[m]->coordI[n] == acon->output_index) 
				out_a = m;
			else if (circuit->gates[m]->coordI[n] == bcon->output_index)
				out_b = m;
		}
		for(int n = 0; n < circuit->gates[m]->numO; n++){
			if (circuit->gates[m]->coordO[n] == acon->input_index)
				in_a = m;
			else if (circuit->gates[m]->coordO[n] == bcon->input_index)
				in_b = m;
		}
		if (in_a == out_b || in_b == out_a) return true;
	}
	
	return false;
}
/****************************************
 * returns a connection id from the connection input index
 * **************************************/
int compareColumns(char *output, char **desired, int length, int width){
	int matches, column;
	column = -1;
	for (int j = 0; j < width; j++){
		matches = 0;
//		cout<<"col: "<<j<<endl;
		for (int k = 0; k < length; k++){
//			cout<<output[k]<<" <> "<<desired[j][k]<<endl;
			if (output[k] == desired[j][k])
				matches++;
		}
		if (matches >= length){
			column = j;
			break;
		}
	}
	return column;
}
/****************************************
 * returns a vector of errors for each output bit
 * 0 - if a proper correspondence was found
 * >0 and 1< if not.
 * **************************************/
double compareColumns(char **output, char **desired, int length, int width){
	int matches, alloc_error_col;
	int *column = new int[width];
	double min_error = -1.0;
	double total_error, col_error;
	for (int i = 0; i < width; i++)
		column[i] = -1;
	total_error = 0.0;
	//check for all columns of the output circuit
	for (int i = 0; i < width; i++){
		alloc_error_col = 0;
		//against all columns of the target circuit
		for (int j = 0; j < width; j++){
			col_error = 0.0;
			matches = 0;
			for (int k = 0; k < length; k++){
				if (output[i][k] == desired[j][k])
					matches++;
				else col_error++;
			}
			if (min_error == -1.0){
				min_error = col_error;
				alloc_error_col = j;
			}else if(min_error > col_error){
				min_error = col_error;
				alloc_error_col = j;
			} else if (matches >= length && !findInArray(j, column, i)){
				column[i] = j;
				break;
			}
		}
		if (column[i] == -1.0)
			total_error += min_error/length;	
	}
//	cout<<total_error<<endl;
	return total_error;
}
/****************************************
 * returns a vector of errors for each output bit
 * 0 - if a proper correspondence was found
 * >0 and 1< if not.
 * this function is analogous to the one above but we assume
 * that it is computed using the complemented variables thus 
 * there are more possibilities to generate the correct output
 * **************************************/
double compareColumns(char **output, char **desired, int length, int width, int goal_width){
	int matches, alloc_error_col;
	int *column = new int[width];
	double min_error = -1.0;
	double total_error, col_error;
	for (int i = 0; i < width; i++)
		column[i] = -1;
	total_error = 0.0;
	//check for all columns of the output circuit
	for (int i = 0; i < width; i++){
		alloc_error_col = 0;
		//against all columns of the target circuit
		for (int j = 0; j < goal_width; j++){
			col_error = 0.0;
			matches = 0;
			for (int k = 0; k < length; k++){
				if (output[i][k] == desired[j][k])
					matches++;
				else col_error++;
			}
			if (min_error == -1.0){
				min_error = col_error;
				alloc_error_col = j;
			}else if(min_error > col_error){
				min_error = col_error;
				alloc_error_col = j;
			} else if (matches >= length && !findInArray(j, column, i)){
				column[i] = j;
				break;
			}
		}
		if (column[i] == -1.0)
			total_error += min_error/length;	
	}
//	cout<<total_error<<endl;
	return total_error;
}
/****************************************
 * returns a connection id from the connection input index
 * **************************************/
bool findInArray(int needle, int *haystack, int range){
	for (int i = 0; i < range; i++)
		if(needle == haystack[i])
			return true;
	return false;
}
/****************************************
 * returns a connection id from the connection input index
 * **************************************/
int getConnByInput(aConnection **cons, int count, int in_id){
	for (int n = 0; n < count; n++)
		if (cons[n]->input_index == in_id) return n;
	return -1;
}
/****************************************
 * returns a connection id from the connection output index
 * **************************************/
int getConnByOutput(aConnection **cons, int count, int out_id){
	for (int n = 0; n < count; n++)
		if (cons[n]->output_index == out_id) return n;
	return -1;
}
/****************************************
 * returns a connection id from the connection output index
 * **************************************/
char getConnValByOutput(aConnection **cons, int count, int out_id){
	for (int n = 0; n < count; n++)
		if (cons[n]->output_index == out_id) 
			if(cons[n]->current_value == 0) 
				return '0';
			else if (cons[n]->current_value == 1) return '1';
	return -1;
}/****************************************
 * Dot Multiply two vectors
 * **************************************/
void vectorScalarMult(int w, cuComplex scalar, cuComplex *a_Vec, cuComplex *b_Vec){
	cuComplex result = make_cuFloatComplex(0, 0);
        for (int m = 0; m < w; m++)
                b_Vec[m] =  cuCmulf(a_Vec[m], scalar);
}
/****************************************
 * Creates the gate and set the appropriate parameters
 * **************************************/
void initGate(aGate *A, int I, int O, int cost, int val) {
	A->valuedness = val;
	A -> Cost = cost;
	A -> numI = I;
	A -> numO = O;
	A->inLU = new int[val^A->numI*val^A->numI];
	if (I > O)
		A->outLU = new int[val^A->numI*val^A->numI];
	else
		A->outLU = new int[val^A->numO*val^A->numO];
	A -> representation = new char[REPSIZE];
	A -> coordI = new int[I];
	A -> coordO = new int[O];
}
/****************************************
 * Initiate connection to default values
 * **************************************/
void initConnection(aConnection *con){
	con->input_index = -1;
	con->output_index = -1;
	con->current_value = -1;
	con->circuit_output = false;
	con->circuit_input = false;
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
	binary[n-1] = '\0';

}
int bin2dec(char *binary, int radix)
{
	int k = radix-1;
	int decimal = 0;
	// take care of negative input

//	cout<<endl;
	for (int p = 0; p < radix; p++){
//		cout<<binary[p]<<",  ";
		if (binary[p] == '1'){
			decimal  += pow(2.0, double(k));
		}
		k--;
	}
//	cout<<": "<<decimal<<"  ";
//	cout<<endl;
	return decimal;
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
void getGateRepfromStr(int pos, string str, char *representation){
	representation[0] = str.at(pos);
	representation[1] = str.at(pos+1);
	representation[2] = str.at(pos+2);
	representation[REPSIZE] = '\0';
	
}
/**********************************
* copy array of integers to another array of same length
 **********************************/
void copyIntArray(int *one, int *two, int length){
	for (int h = 0; h < length; h++)
		two[h] = one[h];
}
/**********************************
* copy array of connections to another array of conenctions of the same size
 **********************************/
void copyConnections(aConnection **one, aConnection **two, int length){
	for (int h = 0; h < length; h++){
		two[h]->input_index = one[h]->input_index;
		two[h]->output_index = one[h]->output_index;
		two[h]->current_value = one[h]->current_value;
		two[h]->circuit_input = one[h]->circuit_input;
		two[h]->circuit_output = one[h]->circuit_output;
	}
}
