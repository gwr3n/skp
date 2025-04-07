/*********************************************
 * OPL 12.8.0.0 Model
 * Author: Roberto Rossi
 * Creation Date: 11 Sep 2018 at 09:48:21
 *********************************************/

int N = ...; // Number of objects
range objects = 1..N;

float expectedValues[objects] = ...;
float expectedWeights[objects] = ...;
float varianceCovarianceWeights[objects][objects] = ...;

float C = ...;  // Knapsack capacity
float c = ...; // Penalty cost for exceeding capacity


/**
 * Loss function piecewise linearization - standard normal distribution
 */

int nbpartitions = ...;
range partitions = 1..nbpartitions;
float prob[partitions] = ...;
float means[partitions] = ...;
float error = ...;
float s = ...;           // linearisation step          (possible value 1)
float x0 = ...;          // value of sqrt function at 0	(possible value 1/4)

dvar boolean X[objects]; // Object selector
dvar float+ M;			 // Expected knapsack weight
dvar float+ V;           // Knapsack weight variance
dvar float+ S;           // Knapsack weight standard deviation
dvar float+ P;			 // Expected capacity shortage

/**
 * Piecewise linearisation of sqrt()
 */
range breakpoints = 1..ftoi(ceil((sum(i in objects, j in objects)varianceCovarianceWeights[i][j])/s));
float slopes[i in breakpoints] = (sqrt(i*s) - sqrt((i - 1)*s))/s; 

maximize sum(i in objects) X[i]*expectedValues[i] - c*P; // Maximize profit minus expected capacity shortage cost

subject to{
	M == sum(i in objects) X[i]*expectedWeights[i];                                       // Expected knapsack weight
	V >= sum(j in objects) (sum(i in objects) X[i]*varianceCovarianceWeights[i][j])*X[j]; // Knapsack weight variance
	S == piecewise(b in breakpoints){slopes[b] -> b*s; 0}(0, x0)(V); // sqrt(V)
	
	/**
     * Piecewise linearization of expected capacity shortage (complementary loss)
    *
	forall(p in partitions)
  		P >= -(C-M) + sum(k in 1..p)prob[k]*(C-M) - (sum(k in 1..p)prob[k]*means[k])*S + error*S;
  	P >= - (C-M) + error*S;
  	*/
  	
  	/**
     * Piecewise linearization of expected capacity shortage (loss)
     */
	forall(p in partitions)
  		P >= (sum(k in p..nbpartitions)prob[k]*means[k])*S - sum(k in p..nbpartitions)prob[k]*(C-M) + error*S;
  	P >= error*S;
}

