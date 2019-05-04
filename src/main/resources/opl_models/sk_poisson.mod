/*********************************************
 * OPL 12.8.0.0 Model
 * Author: Roberto Rossi
 * Creation Date: 11 Sep 2018 at 09:48:21
 *********************************************/

int N = ...; // Number of objects
range objects = 1..N;

float expectedValues[objects] = ...;      
float expectedWeights[objects] = ...;

int C = ...;   // Knapsack capacity
float c = ...; // Penalty cost for exceeding capacity


/**
 * Loss function piecewise linearization - Poisson distribution
 */
int nbpartitions = ...;
range partitions = 1..nbpartitions;
int maxWeight = ...;
range weight = 0..maxWeight;
float prob[partitions] = ...;
float means[weight,partitions] = ...;
float error[weight] = ...;

dvar boolean X[objects]; // Object selector
dvar float+ M;			 // Expected knapsack weight
dvar float+ P;			 // Expected capacity shortage

maximize sum(i in objects) X[i]*expectedValues[i] - c*P; // Maximize profit minus expected capacity shortage cost

subject to{
	M == sum(i in objects) X[i]*expectedWeights[i];      // Expected knapsack weight
	
	/**
     * Piecewise linearization of expected capacity shortage (complementary loss)
    *
    forall(w in weight){    
		forall(p in partitions)
		    (M >= w) => P >= -(C-w) + sum(k in 1..p)prob[k]*C - sum(k in 1..p)prob[k]*means[w][k] + error[w];
  		(M >= w) => P >= - (C-w) + error[w];
 	}  	
 	(M >= maxWeight) => P >= - (C-M) + error[maxWeight];
 	*/
 	
 	/**
     * Piecewise linearization of expected capacity shortage (loss)
    */
    forall(w in weight){    
		forall(p in partitions)
		    (M >= w) => P >= sum(k in p..nbpartitions)prob[k]*means[w][k] - C*sum(k in p..nbpartitions)prob[k]  + error[w];
  		(M >= w) => P >= error[w];
 	}  	
 	(M >= maxWeight) => P >= - (C-M) + error[maxWeight];
}

