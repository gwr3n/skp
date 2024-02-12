package skp.milp;

import com.gurobi.gurobi.GRB;
import com.gurobi.gurobi.GRBEnv;
import com.gurobi.gurobi.GRBException;
import com.gurobi.gurobi.GRBLinExpr;
import com.gurobi.gurobi.GRBModel;
import com.gurobi.gurobi.GRBVar;

import skp.instance.SKPGenericDistribution;
import skp.milp.instance.SKPGenericDistributionMILPSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.sim.SimulateNormal;
import skp.utililities.gson.GSONUtility;

public class SKPGenericDistributionMILP{   
   int linearizationSamples;
   
   int[] optimalKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   SKPGenericDistribution instance;
   
   public SKPGenericDistributionMILP(SKPGenericDistribution instance, int linearizationSamples){
      this.instance = instance;
      this.linearizationSamples = linearizationSamples;
   }
   
   double gurobiSolutionTimeMs;
   int simplexIterations;
   int exploredNodes;
   int cuts = 0;

   public int[] getOptimalKnapsack() {
      return this.optimalKnapsack;
   }
   
   public double getMILPSolutionValue() {
      return this.milpSolutionValue;
   }
   
   public double getMILPOptimalityGap() {
      return this.milpOptimalityGap;
   }
   
   public double getMILPMaxLinearizationError() {
      return this.milpMaxLinearizationError;
   }

   void computeMILPMaxLinearizationError() {
      // TODO Auto-generated method stub
      // use PiecewiseFirstOrderLossFunction.getMaxApproximationError 
      // by feeding the optimal knapsack weight distributions to the constructor.
      System.err.print("Unimplemented method");
   }
   
   public SKPGenericDistributionMILPSolvedInstance solve(int simulationRuns) {
      try {
         GRBEnv env = new GRBEnv();
         GRBModel model = new GRBModel(env);
         
         // Must set LazyConstraints parameter when using lazy constraints
         //model.set(GRB.IntParam.LazyConstraints, 1);
         
         model.set(GRB.StringAttr.ModelName, "SKPGenericDistribution");
         
         // Create decision variables
         
         // Object selectors
         GRBVar[] X = new GRBVar[this.instance.getItems()];
         for (int i = 0; i < this.instance.getItems(); ++i) {
            X[i] = model.addVar(0.0, 1.0, this.instance.getExpectedValues()[i], GRB.BINARY, "X"+i);
         }
         // Expected knapsack weight
         GRBVar M = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "M");
         // Expected capacity shortage
         GRBVar P = model.addVar(0, GRB.INFINITY, -instance.getShortageCost(), GRB.CONTINUOUS, "P");
         
         // Create constraints
         
         // Expected knapsack weight
         {
            GRBLinExpr expr = new GRBLinExpr();
            for (int i = 0; i < this.instance.getItems(); i++)
               expr.addTerm(this.instance.getWeights()[i].getMean(), X[i]);
            model.addConstr(expr, GRB.EQUAL, M, "Expected knapsack weight");
         }
         
         // Expected capacity shortage (loss): P >= - (C-M)
         //                                    P-M >= -C
         {
            GRBLinExpr expr = new GRBLinExpr();
            expr.addTerm(1.0, P);
            expr.addTerm(-1.0, M);
            model.addConstr(expr, GRB.EQUAL, -this.instance.getCapacity(), "Expected knapsack weight");
         }
         
         // The objective is to maximize the profit minus the expected capacity shortage
         //
         // Note: The objective coefficients are set during the creation of
         //       the decision variables above.
         model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
         
         // Solve
         model.optimize();
         
         if (model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {
         
         this.milpSolutionValue = model.get(GRB.DoubleAttr.ObjVal);
         this.milpOptimalityGap = model.get(GRB.DoubleAttr.MIPGap);
         this.gurobiSolutionTimeMs = model.get(GRB.DoubleAttr.Runtime);
         this.simplexIterations = (int) Math.round(model.get(GRB.DoubleAttr.IterCount));
         this.exploredNodes = (int) Math.round(model.get(GRB.DoubleAttr.NodeCount));
         
         this.optimalKnapsack = new int[instance.getItems()];
         for(int i = 0; i < instance.getItems(); i++){
            this.optimalKnapsack[i] = (int) Math.round(X[i].get(GRB.DoubleAttr.X));
         }
         
         this.milpMaxLinearizationError = 0;
         
         } else {
            System.err.println("No solution");
         }
         
         System.out.println("JSON solution:" + model.getJSONSolution());
         
         // Dispose of model and environment
         model.dispose();
         env.dispose();
         
         SimulateGenericDistribution sim = new SimulateGenericDistribution(instance);
         double simulatedSolutionValue = sim.simulate(optimalKnapsack, simulationRuns);
         
         double milpMaxLinearizationError = 0; // cuts are tight
         double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;
         
         SKPGenericDistributionMILPSolvedInstance solvedInstance = new SKPGenericDistributionMILPSolvedInstance(
               instance,
               this.getOptimalKnapsack(),
               simulatedSolutionValue,
               simulationRuns,
               this.getMILPSolutionValue(),
               this.getMILPOptimalityGap(),
               cuts,
               this.linearizationSamples,
               milpMaxLinearizationError,
               simulatedLinearizationError,
               gurobiSolutionTimeMs,
               simplexIterations,
               exploredNodes
               );
         
         return solvedInstance;
         
      } catch (GRBException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      
      return null;
   }
   
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstance();
      
      try {
         SKPGenericDistributionMILP sskp = new SKPGenericDistributionMILP(instance, 100000);
         
         System.out.println(GSONUtility.<SKPGenericDistributionMILPSolvedInstance>printInstanceAsJSON(sskp.solve(100000)));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
