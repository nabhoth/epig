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
#include "EpiG.h"
#include "string.h"
//some methods for linked list
//inserting before the first element of the list
group* QE::insertB(group* A, group* m_List)
{
	group *n_head;
	A->next = m_List;
	m_List = A;
	*n_head = *A;
	return n_head;
}
//inserting after the current pointer
void QE::insertA(group* A, group* m_List)
{
	A->next = m_List->next;
	m_List->next = A;
}
//remove first element
group* QE::removeFirst(group* m_List)
{
	group *removed;
	*removed = *m_List;
	*m_List = *m_List->next;
	removed->next = NULL;
	return removed;
}
//shift all elements to the left(dump the leftmost)
group* QE::moveL(group* m_List)
{
	group *n_ext  = m_List;

	//changing order
	m_List = n_ext->next;
	//assigning next
	n_ext->next = m_List->next;
	m_List->next= n_ext;
	return m_List;
}

//return the size of the list
int QE::sizeL(group* m_List)
{
	int c = 0;
	while(m_List->next != 0)
	{
		m_List = m_List->next;
		c++;
	}
	return c;
}
//no more of List methods ... new can be added
//generate combinations for groups of gates
//return position of the gate
int QE::getGate(char *b)
{
//	for (int a = 0;a<numofgates;a++)
//		if ((*gateArray[a]).representation == b) return a;

}

//calculates whreas the search has to be finnished
bool QE::setTerminatingCondition(int type)
{
	if (type != 2){

	for (int a=0;a<=Segments;a++)
		if (gateCounter[a] != numofgates-1)
			return false;
	out_stream<<" generation max reached "<<endl<<endl<<endl;
	return true;
	} else {
		return end;
	}
}

//when reading constatns from file -> transform into real shit
cuComplex QE::assignVal(float x)
{
	if (x == 2)
		return make_cuFloatComplex(0, 1);
//		return complex<float>(0, 1);
	else if (x == -2)
		return make_cuFloatComplex(0, -1);
//		return complex<float>(0, -1);
	else if (x == 4)
		return make_cuFloatComplex(real(exp(complex<float>(0, (3.1416/4)))), imag(exp(complex<float>(0, (3.1416/4)))));
//		return exp(complex<float>(0, (3.1416/4)));
	else if (x == 5)
		return make_cuFloatComplex(0.5, 0.5);
//		return complex<float>(0.5, 0.5);
	else if (x == -5)
		return make_cuFloatComplex(0.5, -0.5);
//		return complex<float>(0.5, -0.5);
	else 	return make_cuFloatComplex(x,0);
}

void QE::addGateSet(qGate *array[]){
	for (int a = 0; a < MAXNUMOFGATES; a++){
		gateArray[a] = new qGate;
		gateArray[a]->numIO = array[a]->numIO;
		gateArray[a]->Cost = array[a]->Cost;
		//gateArray[a]->representation = array[a]->representation;
	int rows = (int(pow(float(2),gateArray[a]->numIO)));
        for (int j = 0; j < rows; j++)
         {
             for (int k = 0; k < rows; k++)
             {
				gateArray[a]->gateMatrix1[j*rows+k] = array[a]->gateMatrix1[j*rows+k];
             }
         }
	}
}
/*******************************************************
*the array of gates IS normalized to the number of wires of the circuit
*******************************************************/
string** QE::setgateNormArray(qGate** array, qGate ideal, int segments, float err){
	cout<<"starting exhaustive search"<<endl;
	int solutionscount = 0, c = 0;
	error = err;
	clock_t time = 0;
	bool solution = false;
	char a;
	int ga, gb;
	string t = string("");

	//initiate bestlocal initial values
	bestlocal.error = 1;
	bestlocal.fitness = 0;
	bestlocal.my_string = string("");

	//in case the first n is not filled jump
	while (array[c] != NULL)
		c++;

	//get the array of pointers to the gates
	for (int g = 0; g<c;g++){
		gateArray[g] = *array++;
		//cout<< gateArray[g]->representation<<endl;
	}

	//initiate the solution array empty
	for (int m = 0; m <MAXNUMOFGATES; m++)
		solutionArray[m] = NULL;

	numofgates = c;
	bestOne = ideal;
	Segments= segments;
	if (INDV != NULL)
		delete (INDV);
	INDV = new qGate;
	INDV->numIO = bestOne.numIO;
	generationCondition = 0;
	//end of initiation

	if (numofgates <= Segments)
		Segments = numofgates - 1;
	//cout<<"starting exhaustive search with the following settings:"<<endl;
	cout<<numofgates<<" gates, "<<Segments<<" segs"<<endl;
	//get the next modified string -- the 0th element --long and painful
	//init the 0th element short and fun
	//synthPermutations(true);
	synthAlterations(true);
	for (int u = 0; u < INDV->my_string.size()/2; u++){
		a = INDV->my_string.at(u);
		INDV->my_string.at(u) = INDV->my_string.at(INDV->my_string.size()-u-1);
		INDV->my_string.at(INDV->my_string.size()-u-1) = a;
	}
	//evaluate it` -- the 0th element
	Matrix(false);
	// life cycle
	solutionscount = 0;
	time = clock();
	//start the cycle of search. yet exhaustive on gates or segments but other improvements to come
	cout<<"starting exhaustive search"<<endl;
	while (!setTerminatingCondition(0))
	{
		Matrix(false);
		if (compareIndividual(INDV))
			++solutionscount;
		if (solutionscount >=  SolutionCount)
			break;
		//checking if the searc is taking too longjmp -- more than 2 minutes
		if (abs(((double)clock()-time))/CLOCKS_PER_SEC >10){
			break;
		} else if (!synthAlterations(false)){
			//altering the INDV string -- inverse it
		cout<<"inverting: "<<endl;
			ga = rand()%INDV->my_string.size()/2;
			do {
				gb = rand()%INDV->my_string.size();
			}while(gb <= ga);


			for (int u = ga; u < gb; u++){
				a = INDV->my_string.at(u);
				INDV->my_string.at(u) = INDV->my_string.at(INDV->my_string.size()-u-1);
				INDV->my_string.at(INDV->my_string.size()-u-1) = a;
			}

			cout<<"inverted: "<<INDV->my_string<<endl;
//			cout<<numofgates<<" gates, "<<Segments<<" segsments, "<<abs(clock()-time)/CLOCKS_PER_SEC<<" secs." <<endl;
		} else {
		}

	}
	//closeStream();
	for (int g = 0; g<solutionscount;g++){
		cout<<*solutionArray[g]<<endl;
		//cout<< gateArray[g]->representation<<endl;
	}
	if (solutionscount <= 0)
		return NULL;
	return solutionArray;
}
/*******************************************************
*compare the currently evaluated individual with the best one
*returns 0 if the bestindividual is better
*returns 1 if the currrent individual is better
*******************************************************/
bool QE::compareIndividual(qGate *A)
{
	float x,y,err= 0;
	int errorcounter = 0;
	int rows = (int(pow((float)2,(float)A->numIO)));
	for (int a = 0; a<rows;a++)
    {
		for (int b = 0; b<rows;b++)
		{
			//simple comparison -- allowing don't cares as value "10"
			if ((cuCrealf(bestOne.gateMatrix1[a*rows+b]) != 10 && cuCimagf(bestOne.gateMatrix1[a*rows+b]) !=  0)
				&& (cuCrealf(bestOne.gateMatrix1[a*rows+b]) != 0 && cuCimagf(bestOne.gateMatrix1[a*rows+b]) != 0))
			{
				errorcounter++;
				err += cuCabsf(cuCsubf(cuCmulf(bestOne.gateMatrix1[a*rows+b],cuConjf(bestOne.gateMatrix1[b*rows+a])), cuCmulf(A->gateMatrix1[a*rows+b], cuConjf(A->gateMatrix1[b*rows+a]))));
				}
		}
	}

    //normalizing error
	if (errorcounter != 0)
		err =  err/errorcounter;

	cout<<"error: "<<err<<" :: ideal: "<<error<<endl;
	if (err < error)
	{
		if (solutionArray[0] != NULL)
		{
			if (err <= bestlocal.error){
				for (int t = SolutionCount-1; t>0 ; t--)
					if (solutionArray[t-1] != NULL){
						if (solutionArray[t] == NULL)
						solutionArray[t] = new string("");
						else{
							delete (solutionArray[t]);
							solutionArray[t] = new string("");
						}

						*solutionArray[t] = *solutionArray[t-1];
					}
				if (solutionArray[0] != NULL)
					delete (solutionArray[0]);
				solutionArray[0]  = new string("");
				for (int h = 0; h < A->my_string.size(); h++)
						*solutionArray[0] += A->my_string.at(h);
					//cout<<"error: "<<err<<" :: ideal: "<<error<<endl;
					//cout<<"string: "<<*solutionArray[0]<<endl;
				bestlocal.my_string = string(A->my_string);
				bestlocal.error = err;
				return true;
			}
			return false;
		}else{
			if (solutionArray[0] != NULL)
					delete (solutionArray[0]);
			solutionArray[0] = new string("");
			for (int h = 0; h < A->my_string.size(); h++)
				*solutionArray[0] += A->my_string.at(h);
					//cout<<"error: "<<err<<" :: ideal: "<<error<<endl;
					//cout<<"new string: "<<*solutionArray[0]<<endl;
				bestlocal.my_string = string(A->my_string);
				bestlocal.error = err;
				return true;
		}
	}
	return false;
}


//close streams
void QE::closeStream()
{
     in_stream.close();
     out_stream.close();
}

/*******************************************************
*This function generates a small binary counter ordered
*quantum circuits.
*
*
********************************************************/
//synthetize segment for the exhaustive search with all gates
void QE::synthSegment(bool begin)
{
	int carry, a;

	if (!begin)
	{
		if (gateCounter[0] < numofgates-1)
		{
			gateCounter[0]++;
		}
		else
		{
			gateCounter[0] = 0;
			carry = 1;
			a = 1;
			while(true)
			{
				if 	(gateCounter[a] == numofgates-1)
				{
					carry = 1;
					gateCounter[a] = 0;
					if (segPointer == a)
					{
						segPointer++;
						break;
					}

				}
				else
				{
					gateCounter[a]++;
					carry = 0;
					break;
				}

				if (a < numofgates-1) a++;
				else break;

			}
		}
	}

	INDV->my_string = "";
	char *temp = new char[2];
	carry = 0;
	for (a = segPointer;a>=0;a--)
	{
			//temp[0] = gateArray[gateCounter[a]]->representation;
			temp[1] = '\0';
			INDV->my_string.insert(0, temp);
	}
	delete temp;

	cout<<" String is "<<INDV->my_string<<endl;
}

/*********************************************************
*Function generating permutationsin a given set of gates.
*All is controlled by the set of gates this just generates
*all permutations without repetition. Ex gates:0,1,2,3;
* Segments: 3; Thus we will have 012, 013, 023, 123.
********************************************************/
//the set of input gate MUST be biger than the size of circuit
bool QE::synthPermutations(bool begin){

	int carry, a, hash1 = 0, hash2 = 0;
//cout<<" String is "<<endl;
	if (!begin)
	{
		if (gateCounter[0] < numofgates-1)
		{
			gateCounter[0]++;
		}
		else
		{
			//gateCounter[0] = 0;
			a = 1;
			while(true)
			{
				//increment if all indexes have to be updated
				if 	(gateCounter[a] < gateCounter[a-1]-1)
				{
					gateCounter[a]++;
					for (int t = a-1; t >= 0; t--)
						gateCounter[t] = gateCounter[t+1]+1;
					break;
				}
				//shifts the whole sequence to the left/right
				else if (a == Segments-1)
				{
					if (gateCounter[0] == numofgates -1){
						end = true;
						break;
					}
				}

				if (a < Segments-1) a++;
				else break;
			}
		}
	} else {
		cout<<" String is "<<Segments<<endl;

		for (a = Segments-1; a>=0; a--)
			gateCounter[a] = Segments-1-a;
	}
		cout<<" String is: "<<Segments<<endl;

	for (int m = 0; m< INDV->my_string.size(); m++)
		hash1 += 2^INDV->my_string.at(m);


	//setting the result in the memeory
	INDV->my_string = string("");
	char *temp = new char[2];
	carry = 0;
	//cout<<" String is "<<Segments<<endl;
	for (a = Segments-1;a>=0;a--)
	{
		//cout<<" String is "<<Segments-1<<" "<<gateCounter[a]<<" "<<numofgates<<endl;
//			temp[0] = gateArray[gateCounter[a]]->representation;
			temp[1] = '\0';
			INDV->my_string.insert(0, temp);
	}
	//cout<<" String is "<<Segments<<endl;
	delete temp;

	for (int m = 0; m< INDV->my_string.size(); m++)
		hash2 += 2^INDV->my_string.at(m);

	//if we are at the end of the sequence we must give up anyway
		cout<<" String is INDV: "<<INDV->my_string<<endl;
	if (hash2 - hash1 != 0){
		cout<<hash2 - hash1<<" hash "<<endl;
		return true;
	}else return false;

}


/*********************************************************
*Function generating alterations in a given set of gates.
*All is controlled by the set of gates this just generates
*all permutations without repetition. Ex gates:0,1,2,3;
* Segments: 3; Thus we will have 0123, 1023, 1203, 1230.
********************************************************/
//the set of input gate MUST be biger than the size of circuit
bool QE::synthAlterations(bool begin){
	int carry, a, hash1 = 0, hash2 = 0;

	if (!begin)
	{
		a = segmentstate->currentpos;
		//bring the circuit segment one step towards the end
		if (a+1 < Segments-1){
			carry = gateCounter[a];
			gateCounter[a] = gateCounter[a+1];
			gateCounter[a+1] = carry;
			segmentstate->currentpos += 1;
		}else{
			//check rules to move the strings or blocks
			if (segmentstate->initpos == Segments-1){
				carry = gateCounter[a];
				gateCounter[a] = gateCounter[a+1];
				gateCounter[a+1] = carry;
				segmentstate->currentpos = 0;
				segmentstate->initpos += 1;
				////break;
			} else {
//				carry = gateCounter[a];
//				gateCounter[a] = gateCounter[a+1];
//				gateCounter[a+1] = carry;
				segmentstate->currentpos = 0;
				segmentstate->initpos += 1;
			}
		}
	} else {
		if (segmentstate == NULL){
			segmentstate = new qState();
		}
		segmentstate->initpos = 0;
		segmentstate->currentpos = 0;

//		cout<<" String is "<<Segments<<endl;
		segmentstate->initstate = string(INDV->my_string);

		for (a = Segments-1; a>=0; a--)
			gateCounter[a] = Segments-1-a;
	}
	//new segment generated
//	cout<<" String is: "<<Segments<<endl;
	for (int m = 0; m< INDV->my_string.size(); m++)
		hash1 += 2^INDV->my_string.at(m);


	//setting the result in the memeory
	INDV->my_string = string("");
	char *temp = new char[2];
	carry = 0;
	//cout<<" String is "<<Segments<<endl;
	for (a = Segments-1;a>=0;a--)
	{
		//cout<<" String is "<<Segments-1<<" "<<gateCounter[a]<<" "<<numofgates<<endl;
			//temp[0] = gateArray[gateCounter[a]]->representation;
			temp[1] = '\0';
			INDV->my_string.insert(0, temp);
	}
	//cout<<" String is "<<Segments<<endl;
	delete temp;

	for (int m = 0; m< INDV->my_string.size(); m++)
		hash2 += 2^INDV->my_string.at(m);

	//if we are at the end of the sequence we must give up anyway
	cout<<" String is INDV: "<<INDV->my_string<<endl;
	if (hash2 - hash1 != 0){
		cout<<hash2 - hash1<<" hash "<<endl;
		return true;
	}else return false;
}

void QE::synthRestric(bool begin)
{
	group *curr, *head;
	head = NULL;

	//creating linked list
	for (int a = 0;a<=Segments;a++)
	{
		curr = (group *)malloc(sizeof(group));
		curr->V = a;
		curr->next = head;
		head = curr;
	}
	curr = head;


}


group* QE::genPerms(group *head, bool perm)
{
	group *h_ead;

	if (sizeL(head) >1)
	{
		h_ead = removeFirst(head);
		genPerms(h_ead, true);
	}else
	{
		if(perm)
			return moveL(head);
		else return head;
	}

}

//set the segments counter to the specified value
void QE::setSegments(int begin)
{

	for (int a = 0;a<MAXSEGNUM;a++)
		 gateCounter[a] = 0;

	if (begin > 0)
		for (int a = 0;a<begin-2;a++)
			gateCounter[a] = 0;

		segPointer = begin-1;

}

void QE::restricSegment()
{
	char segs[2], grps[100], *pt;
	int a, b, counter;
	string iN;

			 cout<<" input the number of segments you want to begin the syntheis for: ";
			 cin>>segs[0];
			 if ( segs[0] >= '0' && segs[0] <='9' )
			 {
				 cout<<" setting segments to "<<*segs<<endl;
				 setSegments(atoi(segs));
				 Segments = segPointer;

				 cout<<" gates for selection are: "<<endl;
				 for (a = 0;a<numofgates;a++)
				 {
					 cout<<a<<" - "<<nameArray[a]<<endl;
				 }
				 cout<<" type group or single gates separated by space "<<endl;
				 cout<<" and then the number of segments they should allocate "<<endl<<endl<<endl;

				 a = 0;
				 b = 0;
				 counter = 0;
				 while(a != (Segments+1))
				 {
					 pt = grps;
					 cout<<" type the number of segments they should allocate "<<endl;
					 cout<<" should be smaller or equal than "<<(Segments - a+1)<<endl;
					 cout<<" and in the same line type please a group of gates --only numbers-- "<<endl;
					 cin>>counter;
					 cin.getline(pt,100, '\n');
					 if (counter <= (Segments - a+1))
					 {
						 a += counter;
						 int c = 0;
						 for (c = 0;c<counter;c++)
						 	restriGate[b+c] = new int[MAXNUMOFGATES];
						 pt = &grps[0];
						 c = 0;
						 while(pt < (pt+strlen(pt)))
						 {
							 int w =  atoi(pt);
							 for (int t = 0;t<counter;t++)
							 	restriGate[b+t][c] = w;
							 pt++;
							 while (pt[0] != ' ' && pt < (pt+strlen(pt)))
							 {
								 pt++;
							 }
							 c++;
							  if (*pt== '\n') break;

						 }
					 for (int t = 0;t<counter;t++)
						 restriGate[b+t][c] = -1;

					 }
					 else cout<<" WRONG !!  "<<endl;

					 b+= counter;
				 }

			 }else
			 {
				 setSegments(0);
			 }

			 for (int g = 0;g<=Segments;g++)
			 {
				 a = 0;
				 while(true)
					 {
						 if (restriGate[g][a] == -1) break;
						 cout<<restriGate[g][a]<<" ";
						a++;
					}
				 cout<<endl;
		 }

		 synthRestric(true);
}


//compare two gates
bool QE::compareGate(qGate *A)
{
	bool result;

	int rows = (int(pow(float(2),A->numIO)));
	for (int a = 0;a<population;a++)
	{
		result = true;
		for (int  j = 0; j<rows;j++)
			{
				for (int k = 0; k<rows;k++)
					{
						if ((cuCrealf(resultCircuits[a]->gateMatrix1[j*rows+k]) != cuCrealf(A->gateMatrix1[j*rows+k])) &&(cuCimagf(resultCircuits[a]->gateMatrix1[j*rows+k]) != cuCimagf(A->gateMatrix1[j*rows+k])))
							result = false;
					}
			}
			if (result) break;
	}
	return result;
}

void QE::Matrix(bool display)
{    //temporary storage values
	 float Cost = 0;
	 //three free gates for operations
	 //finalG is the mtrix represenation for fitness evaluation
     qGate *temp0, *temp1, *temp2;
	int rows;
	 string S = INDV->my_string;
	 int j = 0;
	temp2 = new qGate;
	for (int a = 0;a<= numofgates;a++)
/*        if (gateArray[a]->representation == S.at(j))
		{
			temp0 = gateArray[a];
			Cost+= (*gateArray[a]).Cost;
			break;
		}
*/	 j++;

	while (j<=S.length()-1)
	 {
		 for (int  a = 0;a<= numofgates;a++)
/*			 if (gateArray[a]->representation == S.at(j))
			 {
				 temp1 = gateArray[a];
				 Cost+= gateArray[a]->Cost;
				 break;
			}
*/			 initMatrix(temp2, temp0->numIO);
			 cblasMatrixProduct(temp0, temp1, temp2);
		temp0 = ccGate(temp2, temp0);
		j++;
	 }
	 //assignement
		rows = (int(pow(float(2),temp0->numIO)));
		for (int  a = 0; a<rows;a++)
				for (int b = 0; b<rows;b++)
					INDV->gateMatrix1[a*rows+b] = temp2->gateMatrix1[a*rows+b];
	delete temp2;
	//cout<< "work done on "<<INDV->my_string<<endl;
}

//initialize all from another class with preset field of data or no
void QE::subInit()
{
	int command = 0;
	 INDV = new qGate;
	 generationCondition = 0;
	 cout<<" complete search, restricted search or permutation search?(0/1/2) ";
	 cin>>command;


			 cout<<" enter the number of maximum allowed segments to search ";
			 cin>>Segments;
			 synthPermutations(true);


		Matrix(false);

	// life cycle
	while (!setTerminatingCondition(command))
	{
		if (command == 2)
			synthPermutations(false);
		else
			synthSegment(false);
		Matrix(false);

	}
	closeStream();

}
