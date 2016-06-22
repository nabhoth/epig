//Size of the population: 
#define ga_populationNumber 100
//Number of segments in each circuit (approx):             
#define ga_segments 25 
//Minimal Number of segments in each circuit (approx):     
#define ga_mincircuitsize 18
//Estimated minimum cost of the goal circuit:              
#define ga_divider 15 
//Number of GA generations before EX search starts:        
#define ga_alterations 10000000
//Total number of GA generations:                          
#define ga_generations 500000
//Mutation probability is:                                 
#define ga_proba_Mutation 0.1
//Crossover probability is:                                
#define ga_proba_Crossover 0.8
//Factor alpha is (for complex fitness):                   
#define ga_alpha 0.98
//Factor beta is (for complex fitness):                    
#define ga_beta 0.02
//Factor alpha1 is (for pareto replication):               
#define ga_alpha1 1
//Factor beta1 is  (for pareto replication):               
#define ga_beta1 1
//Phase switch (0 - off, 1 - on):                          
#define ga_phase 0
//Display level:                                           
#define ga_displaylevel 1
//Type of GA (0 - normal, 1 - Darwinian):                  
#define ga_Ga 1
//Use Elitism (0 - no, 1 - yes):                           
#define ga_elitism 1
//Type of mutation (0 - normal, 1 - bitwise):              
#define ga_mutation 1
//Type of crossover (0 - 1point, 1 - 2point):              
#define ga_crossover 1
//Type of replication (0 - RW, 1 - SUS):                   
#define ga_replicator 1
//Type of fitness (0 - simplest, 3 - complex):             
#define ga_fitness 3
//Fitness calculation (0 - individual, 1 - grouped):       
#define ga_grouped 0
//Pareto optimization (0 - no, 1 - yes):                   
#define ga_pareto 0
//Max-Fitness Threshold for Replication (0 - no, 1< >0):   
#define ga_threshold 0.05
//The number of wires the final circuit has:               
#define ga_resultnum 5
//Valuedness of designed logic:                            
#define ga_valuedness 2
//Error adjustment for Entanglement design:                         
#define ga_error_adjust 0.6
//Design of Sequential Logic/ Sequence Detector:           
#define ga_seq_detect_enabled 1
//Measurement used:                                        
#define ga_measurement 0
//Measurement used:                                        
#define ga_poly_search 0

