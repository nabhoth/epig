#define __LINUX__     // must be __LINUX__ or __WINDOWS__ depending on platform

extern "C" {

#ifdef __LINUX__

#define VERBOSE 0
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdint.h>
#include <fpu_control.h>
#include "qcost.c"

#include "QMDDinclude.h"


static inline int memReadStat(int field)
{
   char    name[256];
   pid_t pid = getpid();
   sprintf(name, "/proc/%d/statm", pid);
   FILE*   in = fopen(name, "rb");
   if (in == NULL) return 0;
   int     value;
   for (; field >= 0; field--)
     fscanf(in, "%d", &value);
   fclose(in);
   return value;}

   uint64_t memUsed() { return (uint64_t)memReadStat(0) * (uint64_t)getpagesize(); }


#endif

}

#include "EpiG.h"

/**********************************
* decodes the name of the gate into a RevLib compatible format
**********************************/
void decodeGateToRevLib(qGate *A){
	string name = "";
	int c_count = 0;
	int m, u_count_name;
	string u;
	char buff[10];
	int wire_counter = 0;
	int controls[A->my_string.length()];
	int u_placement = -1;
	u_count_name = 0;
	string u_name = "";
	for (int m = 0; m < A->my_string.length(); m++){
		name += A->my_string.at(m);
		if (A->my_string.at(m) == 'C'){
			controls[c_count] = wire_counter;
			c_count++;
			wire_counter++;
			u_count_name = 0;
		} else if (A->my_string.at(m) == 'I'){
			wire_counter++;
			u_count_name = 0;
		} else {
			if (u_placement == -1){
				u_placement = m;
				wire_counter++;
			}
			u_name += A->my_string.at(m);
		}
	}

	if (u_name == "NOT"){
		if (c_count == 0)
			A->rev_lib_char = "not";
		if (c_count == 1)
			A->rev_lib_char = "c";
		if (c_count == 2)
			A->rev_lib_char = "t3";
		if (c_count == 3)
			A->rev_lib_char = "t4";
	} else if (u_name == "V"){
			A->rev_lib_char = "v";
	}  else if (u_name == "V+"){
			A->rev_lib_char = "v+";
	}

//	cout << "Gate: "<< A->my_string<<endl;
	A->rev_lib_string = "";
	for (int r = 0; r < c_count; r++){
		sprintf(buff, "%d", controls[r]);
		m = 0;
		while(buff[m] > 0)
		A->rev_lib_string += buff[m++];
		A->rev_lib_string += ' ';
	}
	sprintf(buff, "%d", u_placement);
	m = 0;
	while(buff[m] > 0)
	A->rev_lib_string += buff[m++];
//	cout<<"revlibs representation: "<<A->rev_lib_char<<"  "<<A->rev_lib_string<<endl;
}


/****************************************
 * Linking to the DMDD package 
 * **************************************/
void doQMDDFitness(Individual *ind, Individual *best, bool display) {
	//two dummy gates and one pointer for operations
	QMDDrevlibDescription circ, bestcirc, result;
	QMDDedge trans, mult;
	QMDDinit(0);
	//complex val;
	cuComplex d;

	Individual *indi = ind;
	qGate* myfinalG;
	int rowcount, rowcounter,maxcount;
	int errorcounter = 0;
	float re,im,equal;

	int num_entries = (int)(pow(pow((float)indi->valuedness, (float)indi->ioNumber),2));
	indi->Error = 0;
	indi->fitness = 0;


//cout<<"Building: "<<indi->rev_lib_string<<endl;
	const char *p = indi->rev_lib_string.c_str();
    	circ=QMDDcircuitRevlib(&p,circ,0);
if (display)
	QMDDmatrixPrint2(circ.e);
	cuComplex x = make_cuFloatComplex(QMDDsumDiagonalDoubleR(circ.e)/num_entries, QMDDsumDiagonalDoubleI(circ.e)/num_entries);
//	cout<<cuCrealf(x)<<", "<<cuCimagf(x)<<endl;
//cout<<"END: "<<endl;
	const char *b = best->rev_lib_string.c_str();
    	bestcirc=QMDDspecRevlib(&b);
if (display)
	QMDDmatrixPrint2(bestcirc.e);
//cout<<"END: "<<endl;
//exit(0);
	trans=QMDDconjugateTranspose(bestcirc.e);
	//trans=QMDDconjugateTranspose(circ.e);
        mult=QMDDmultiply(circ.e,trans);
if (display)
	QMDDmatrixPrint2(mult);
	if(mult.p->ident) equal = 1; else equal = 0;
//	QMDDprint(mult,100);
//	printf("sum on diagonal ");
	cuComplex a = make_cuFloatComplex(QMDDsumDiagonalDoubleR(mult), QMDDsumDiagonalDoubleI(mult));
//	cout<<cuCrealf(a)<<", "<<cuCimagf(a)<<endl;
	re = QMDDsumDiagonalDoubleR(mult)/num_entries;
	im = QMDDsumDiagonalDoubleI(mult)/num_entries;
	d = cuCmulf(make_cuFloatComplex(re,im),make_cuFloatComplex(re,-im));
//	cout<<cuCrealf(d)<<", "<<cuCimagf(d)<<":: "<<re<<", "<<im<<endl;
	if(mult.p->ident){ 
		//if ((re == abs(1)))
			indi->fitness = 1;
		cout<<"Building: "<<indi->rev_lib_string<<endl;
		cout<<"Building: "<<best->rev_lib_string<<endl;
	}else {
		indi->fitness = cuCrealf(d);
	//	if (re > 1)
	//		indi->fitness = 1/re;
	//	else if (re == 1) indi->fitness = re/2;
	//		else indi->fitness = re;
	}
	indi->Error = 1 - indi->fitness;
//	cout<<indi->fitness<<endl;

}
