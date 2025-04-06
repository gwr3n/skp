/*********************************************
 * OPL 12.8.0.0 Model
 * Author: gwren
 * Creation Date: 4 May 2019 at 04:10:58
 *********************************************/

int N = ...; // Number of objects
range objects = 1..N;
int S = ...; // Number of scenarios
range scenarios = 1..S;

float expectedValues[objects] = ...;
float weights[scenarios][objects] = ...;

float C = ...;  // Knapsack capacity
float c = ...;  // Penalty cost for exceeding capacity

dvar boolean X[objects];  // Object selector
dvar float+ P[scenarios]; // Expected capacity shortage

// Maximize profit minus expected capacity shortage cost
maximize sum(i in objects) X[i]*expectedValues[i] - c/S* sum(s in scenarios) P[s]; 

subject to{
  	forall(s in scenarios)
		P[s] >= sum(i in objects) X[i]*weights[s][i] - C;
}