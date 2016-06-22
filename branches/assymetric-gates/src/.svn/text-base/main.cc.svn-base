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
#include <pthread.h>
#include "string.h"
/********************************************
* The one and only application object
********************************************/

int main(int argc, char* argv[])
{
	char str[128];
	char out_path[120];
	char command[7]  = "mkdir ";;
	char s[120];
	if ( argc != 3 ){ // argc should be 2 for correct execution
		// We print argv[0] assuming it is the program name
		cout<<argc<<" usage: "<< argv[0] <<" <filename> <output_path>\n";
		exit(0);
	}else {

		FILE *fp = fopen(argv[1],"r");
		if( fp ) {
			// exists
			 fclose(fp);
			strcpy (str, argv[1] );
		} else {
			 // doesnt exist
			cout<<"Could not open input file\n";
			exit(0);
		}
		fp = fopen(argv[2],"r");
		if( fp ) {
			// exists
			 fclose(fp);
			strcpy ( out_path, argv[2] );
		} else {
			 // doesnt exist
			strcpy ( out_path, argv[2] );
			strcpy ( s, command );
			strcat(s, out_path);;
			cout<<"Could not open output directory! Creating...\n";
			cout<<"File location: "<<s<<endl;
			system(s);
		}
	}
	cout<<"\nOutput path: "<<out_path<<endl;
	//ignition
	GA g;
	srand( time(NULL) );
	int generationCondition = 0;
	pthread_t *childs = new pthread_t[MAXNUMBER];
	int status, died;
	int datadisplaycounter;
	int pipefd[2];
	float results[MAXNUMBER];
	int bestcounter = 0;
	//the circuits passed form the GA to various optimization and search modules
	string *transmitted[MAXNUMOFGATES];
	//start the GA
#ifdef __TIMING__
	long time = clock();
#endif
	g.initializePop(0, str, out_path);
#ifdef __TIMING__
	time = clock() -time ;
	cout<<"Initialization Time: "<<time<<",  "<<(double)(time)/(double)(CLOCKS_PER_SEC)<<"s"<< (double)((double)(time)/(double)(CLOCKS_PER_SEC))/(60)<<"min"<<endl;
#endif
	g.setTerminatingCondition(generationCondition);
#ifdef __TIMING__
	time = clock();
#endif
	g.evaluatePopulation();
#ifdef __TIMING__
	time = clock() -time ;
	cout<<"Initial Evaluation Time: "<<time<<",  "<<(double)(time)/(double)(CLOCKS_PER_SEC)<<"s"<< (double)((double)(time)/(double)(CLOCKS_PER_SEC))/(60)<<"min"<<endl;
#endif
	cout<< "INIT Stage  done: GA Ready"<<endl;
	//evaluate the population of the circuits for the first time
	//start to iterate the GA with modules
	int counter = 0;
	g.finalResult.counter = counter;
	//life cycle
	datadisplaycounter = 0;
	while (!g.terminatingCondition())
	{
#ifdef __TIMING__
	time = clock() ;
#endif
		datadisplaycounter++;
		if (generationCondition%100 == 0)
			cout<<" current generation is "<<generationCondition<<endl;
//		cout<<" replication start"<<endl;
		g.applyReplication();
//		cout<<" replication done"<<endl;
		for (int a = 0;a<g.populationNumber;a++)
		{
			g.doMutation(g.population[a]);
		}
//		cout<<" mutation done"<<endl;
		bestcounter = 0;
		for (int i = 0; i < g.populationNumber; i++)
		{
	            //pthread_create(&childs[i], NULL, &GA::startSub, &g);
	            g.makeFitness(g.population[i]);
		}	
//		cout<<" fitness done"<<endl;
/*		for (int n = 0; n < g.populationNumber; n++){
			void* thr_retval;
			pthread_join(childs[n], &thr_retval);
		}
*/		g.evaluatePopulation();
//		cout<<" evaluation done"<<endl;
//		if (datadisplaycounter >= 10){
//			g.calcAverage();
//			datadisplaycounter = 0;
//		}


		generationCondition++;
		g.setTerminatingCondition(generationCondition);
//		if (g.terminatingCondition()){
//			g.outputBest(true);
//		}

		counter++;
#ifdef __TIMING__
	time = clock() -time ;
	cout<<"Generation Time: "<<time<<",  "<<(double)(time)/(double)(CLOCKS_PER_SEC)<<"s"<< (double)((double)(time)/(double)(CLOCKS_PER_SEC))/(60)<<"min"<<endl;
#endif
	}
	
	//---------closing the output stream---------------
//	g.closeStream();
	return 0;
}
