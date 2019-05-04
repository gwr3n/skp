/*********************************************
 * OPL 12.8.0.0 Model
 * Author: gwren
 * Creation Date: 4 May 2019 at 04:10:58
 *********************************************/

int N = ...; // Number of objects
range objects = 1..N;

float values[objects] = ...;
float weights[objects] = ...;

float C = ...;  // Knapsack capacity
float c = ...; // Penalty cost for exceeding capacity

dvar boolean X[objects]; // Object selector
dvar float+ P;			 // Expected capacity shortage

maximize sum(i in objects) X[i]*values[i] - c*P; // Maximize profit minus expected capacity shortage cost

subject to{
	P == maxl( sum(i in objects) X[i]*weights[i] - C, 0);
}