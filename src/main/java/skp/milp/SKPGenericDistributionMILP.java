package skp.milp;

import java.util.Arrays;

import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.opl.*;

import skp.folf.PiecewiseFirstOrderLossFunction;
import skp.instance.SKPGenericDistribution;
import skp.milp.instance.SKPGenericDistributionMILPSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utililities.gson.GSONUtility;

import umontreal.ssj.probdist.Distribution;

public class SKPGenericDistributionMILP { 
   private static long[] seed = {1,2,3,4,5,6};
   
   int linearizationSamples;
   
   int[] optimalKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   SKPGenericDistribution instance;
   
   double cplexSolutionTimeMs;
   int simplexIterations;
   int exploredNodes;
   int cuts = 0; // cuts counter
   
   public SKPGenericDistributionMILP(SKPGenericDistribution instance, int linearizationSamples){
      this.instance = instance;
      this.linearizationSamples = linearizationSamples;
   }
   
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
   
   // Piecewise linearization callback.  Whenever a feasible solution is found,
   // find the cut that corresponds to the current partial assignment.
   /*protected void callback() {
      try {
         if (where == GRB.CB_MIPSOL) {
            this.cuts++;
            
            double[] x = this.getNodeRel(X);

            Distribution[] weights = this.instance.getWeights();
            Distribution[] reducedWeights = IntStream.iterate(0, i -> i + 1)
                  .limit(weights.length)
                  .filter(i -> x[i] == 1.0)
                  .mapToObj(i -> weights[i])
                  .toArray(Distribution[]::new);

            PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(reducedWeights, seed);

            double a = pwfolf.getFirstOrderLossFunctionValue(this.instance.getCapacity(), linearizationSamples)-
                  (pwfolf.getEmpiricalDistribution(linearizationSamples).cdf(this.instance.getCapacity())-1)*this.instance.getCapacity();
            double b = pwfolf.getEmpiricalDistribution(linearizationSamples).cdf(this.instance.getCapacity())-1; //slope

            // Expected capacity shortage (loss): P >= bC+a 
            //                                    
            GRBLinExpr expr = new GRBLinExpr();
            expr.addTerm(1.0, this.P);
            addLazy(expr, GRB.GREATER_EQUAL, b*this.instance.getCapacity()+a);
         }
      } catch (GRBException e) {
         System.out.println("Error code: " + e.getErrorCode() + ". " +
               e.getMessage());
         e.printStackTrace();
      }
   }*/
   
   public SKPGenericDistributionMILPSolvedInstance solve(int simulationRuns) {
      try {
         IloCplex cplex = new IloCplex();
         
         // Create decision variables
         
         // Object selectors
         IloNumVar[] X = cplex.numVarArray(this.instance.getItems(), 0, 1, IloNumVarType.Int);
         
         // Expected knapsack weight
         IloNumVar M = cplex.numVar(0, Double.POSITIVE_INFINITY);
         
         // Expected capacity shortage
         IloNumVar P = cplex.numVar(0, Double.POSITIVE_INFINITY);
         
         // Create constraints
         
         // Expected knapsack weight
         cplex.addEq(cplex.scalProd(Arrays.stream(this.instance.getWeights()).mapToDouble(w -> w.getMean()).toArray(), X), M);
         
         // Expected capacity shortage (loss): P >= - (C-M)
         //                                    P-M >= -C
         cplex.addGe(cplex.sum(P,cplex.prod(M, -1)), -this.instance.getCapacity());
         
         
         // The objective is to maximize the profit minus the expected capacity shortage
         cplex.addMaximize(cplex.sum(
               cplex.scalProd(this.instance.getExpectedValues(), X),
               cplex.prod(-instance.getShortageCost(), P))
               );
         
         
         cplex.setParam(IloCplex.Param.MIP.Strategy.Search,
                        IloCplex.MIPSearch.Traditional);
         
         double start = cplex.getCplexImpl().getCplexTime();
         boolean status =  cplex.solve();
         double end = cplex.getCplexImpl().getCplexTime();
         if ( status ) {   
            this.milpSolutionValue = cplex.getObjValue();
            this.milpOptimalityGap = cplex.getMIPRelativeGap();
            this.cplexSolutionTimeMs = (end - start)*1000;
            this.simplexIterations = cplex.getNiterations();
            this.exploredNodes = cplex.getNnodes();
            
            this.optimalKnapsack = new int[instance.getItems()];
            for(int i = 0; i < instance.getItems(); i++){
               this.optimalKnapsack[i] = (int) Math.round(cplex.getValue(X[i]));
            }
            System.out.println("M: "+cplex.getValue(M));
            System.out.println("P: "+cplex.getValue(P));
            
            this.milpMaxLinearizationError = 0; // cuts are tight
            
         } else {
            System.out.println("No solution!");
         } 
         cplex.end();
         System.gc();
         
         SimulateGenericDistribution sim = new SimulateGenericDistribution(instance);
         double simulatedSolutionValue = sim.simulate(optimalKnapsack, simulationRuns);
         double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;
         
         SKPGenericDistributionMILPSolvedInstance solvedInstance = new SKPGenericDistributionMILPSolvedInstance(
               instance,
               this.getOptimalKnapsack(),
               simulatedSolutionValue,
               simulationRuns,
               this.getMILPSolutionValue(),
               this.getMILPOptimalityGap(),
               this.cuts,
               this.linearizationSamples,
               milpMaxLinearizationError,
               simulatedLinearizationError,
               cplexSolutionTimeMs,
               simplexIterations,
               exploredNodes
               );
         
         return solvedInstance;
         
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      
      return null;
   }
   
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstance();
      
      try {
         SKPGenericDistributionMILP sskp = new SKPGenericDistributionMILP(instance, 100000);
         
         SKPGenericDistributionMILPSolvedInstance solvedInstance = sskp.solve(100000);
         System.out.println(GSONUtility.<SKPGenericDistributionMILPSolvedInstance>printInstanceAsJSON(solvedInstance));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
