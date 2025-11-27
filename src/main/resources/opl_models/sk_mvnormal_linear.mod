
int N = ...; // Number of objects
range objects = 1..N;

float expectedValues[objects] = ...;
float expectedWeights[objects] = ...;
float varianceCovarianceWeights[objects][objects] = ...;

float C = ...;  // Knapsack capacity
float c = ...;  // Penalty cost for exceeding capacity

int nbpartitions = ...;
range partitions = 1..nbpartitions;
float prob[partitions] = ...;
float means[partitions] = ...;
float error = ...;
float s = ...;           // linearisation step
float x0 = ...;          // sqrt function at 0

dvar boolean X[objects]; // Object selector
dvar float+ M;           // Expected knapsack weight
dvar float+ V;           // Knapsack weight variance
dvar float+ S;           // Knapsack weight standard deviation
dvar float+ P;           // Expected capacity shortage

// Auxiliary variables for linearization
dvar boolean Y[i in objects, j in objects: i <= j]; // Represents X[i]*X[j]

range breakpoints = 1..ftoi(ceil((sum(i in objects, j in objects)varianceCovarianceWeights[i][j])/s));
float slopes[b in breakpoints] = (sqrt(b*s) - sqrt((b - 1)*s))/s;

maximize sum(i in objects) X[i]*expectedValues[i] - c*P;

subject to {
    // Expected knapsack weight
    M == sum(i in objects) X[i]*expectedWeights[i];

    // Linearized variance constraint using Y[i,j]
    V >= sum(i in objects, j in objects: i <= j) varianceCovarianceWeights[i][j] * (i == j ? Y[i,j] : 2*Y[i,j]);

    // Linking constraints for Y[i,j]
    forall(i in objects, j in objects: i <= j) {
        Y[i,j] <= X[i];
        Y[i,j] <= X[j];
        Y[i,j] >= X[i] + X[j] - 1;
    }

    // Piecewise linearization of sqrt(V)
    S == piecewise(b in breakpoints){slopes[b] -> b*s; 0}(0, x0)(V);

    // Piecewise linearization of expected capacity shortage (loss)
    forall(p in partitions)
        P >= (sum(k in p..nbpartitions)prob[k]*means[k])*S
             - sum(k in p..nbpartitions)prob[k]*(C-M) + error*S;
    P >= error*S;
}
