package skp.milp;

import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.opl.IloCplex;
import ilog.opl.IloOplFactory;
import skp.folf.FirstOrderLossFunctionScalarProduct;
import skp.instance.SKPGenericDistribution;
import skp.milp.instance.SKPGenericDistributionCutsSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.Distribution;

public class SKPGenericDistributionLazyCuts {
   
   int linearizationSamples;
   int simulationRuns;
   int lazyCuts;
   
   int[] optimalKnapsack;
   int[] lastKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   SKPGenericDistribution instance;
   
   double cplexSolutionTimeMs = 0;
   int simplexIterations = 0;
   int exploredNodes = 0;
   
   IloNumVar[] knapsack;
   IloNumVar P;
   IloCplex cplex;
   
   public SKPGenericDistributionLazyCuts(SKPGenericDistribution instance, int linearizationSamples, int simulationRuns){
      this.instance = instance;
      this.linearizationSamples = linearizationSamples;
      this.simulationRuns = simulationRuns;
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
   
   private static final double step = 0.01;
   
   private static double[] computeDirectionalDerivative(SKPGenericDistribution instance, double[] knapsack, int linearizationSamples) {
      Distribution[] weights = instance.getWeights();
      FirstOrderLossFunctionScalarProduct folfsp = new FirstOrderLossFunctionScalarProduct(weights);
      double[] dd = new double[weights.length];
      for(int i = 0; i < dd.length; i++) {
         double[] kp = Arrays.copyOf(knapsack, knapsack.length);
         kp[i] = kp[i] + step;
         dd[i] = (folfsp.getFirstOrderLossFunctionValue(instance.getCapacity(), kp, linearizationSamples) -
                  folfsp.getFirstOrderLossFunctionValue(instance.getCapacity(), knapsack, linearizationSamples))/step;
      }
      return dd;
   }
   
   private static double computeLX(SKPGenericDistribution instance, double[] knapsack, int linearizationSamples) {
      Distribution[] weights = instance.getWeights();
      FirstOrderLossFunctionScalarProduct folfsp = new FirstOrderLossFunctionScalarProduct(weights);
      return folfsp.getFirstOrderLossFunctionValue(instance.getCapacity(), knapsack, linearizationSamples);
   }
   
   private static long time_limit = 60*10; //10 minutes
   private static double tolerance = 1e-3; // // Equivalent to CPLEX https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-relative-mip-gap-tolerance

   public SKPGenericDistributionCutsSolvedInstance solve() throws IloException {
      this.milpSolutionValue = Double.MAX_VALUE;
      this.lastKnapsack = null;
      this.lazyCuts = 0;
      
      SKPGenericDistributionCutsSolvedInstance solvedInstance = null;
      IloOplFactory.setDebugMode(false);
      
      cplex = new IloCplex();
      cplex.setParam(IloCplex.Param.TimeLimit, time_limit);
      cplex.setParam(IloCplex.Param.Threads, 1); // Lazy cuts do not allow multithreading
      cplex.setParam(IloCplex.Param.MIP.Display, 2);
      cplex.setParam(IloCplex.Param.MIP.Strategy.Search, IloCplex.MIPSearch.Traditional);
      cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, tolerance);
      cplex.setParam(IloCplex.Param.Preprocessing.Presolve, false); // we must deactivate the dual presolve
      cplex.setOut(null);

      // Create decision variables

      // Object selectors
      IloNumVar[] X = cplex.numVarArray(this.instance.getItems(), 0, 1, IloNumVarType.Int);
      knapsack = X;

      // Expected knapsack weight
      IloNumVar M = cplex.numVar(0, Double.POSITIVE_INFINITY);

      // Expected capacity shortage
      P = cplex.numVar(0, Double.POSITIVE_INFINITY);

      // Create constraints

      // Expected knapsack weight
      cplex.addEq(cplex.scalProd(Arrays.stream(this.instance.getWeights()).mapToDouble(w -> w.getMean()).toArray(), X), M);

      // Expected capacity shortage (loss): P >= - (C-M)
      //                                    P-M >= -C
      cplex.addGe(cplex.sum(P,cplex.prod(M, -1)), -this.instance.getCapacity());

      //Lazy Cuts
      cplex.use(new LPNLPLazyConstraintCallback());

      // Bound objective function value
      cplex.add(cplex.le(cplex.sum(
            cplex.scalProd(this.instance.getExpectedValues(), X),
            cplex.prod(-instance.getShortageCost(), P)), this.milpSolutionValue));

      // The objective is to maximize the profit minus the expected capacity shortage
      cplex.addMaximize(cplex.sum(
            cplex.scalProd(this.instance.getExpectedValues(), X),
            cplex.prod(-instance.getShortageCost(), P))
            );

      
      double start = cplex.getCplexImpl().getCplexTime();
      boolean status =  cplex.solve();
      double end = cplex.getCplexImpl().getCplexTime();
      if ( status ) {   
         this.milpSolutionValue = cplex.getObjValue();
         this.milpOptimalityGap = cplex.getMIPRelativeGap();
         this.cplexSolutionTimeMs += (end - start)*1000;
         this.simplexIterations += cplex.getNiterations();
         this.exploredNodes += cplex.getNnodes();

         this.optimalKnapsack = new int[instance.getItems()];
         for(int i = 0; i < instance.getItems(); i++){
            this.optimalKnapsack[i] = (int) Math.round(cplex.getValue(X[i]));
         }
         //DEBUG
         //System.out.println("X: "+Arrays.toString(this.optimalKnapsack));
         //System.out.println("M: "+cplex.getValue(M));
         //System.out.println("P: "+cplex.getValue(P));

         this.milpMaxLinearizationError = 0; // cuts are tight

         SimulateGenericDistribution sim = new SimulateGenericDistribution(instance);
         double simulatedSolutionValue = sim.simulate(optimalKnapsack, this.simulationRuns);
         double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;

         solvedInstance = new SKPGenericDistributionCutsSolvedInstance(
               instance,
               this.optimalKnapsack,
               simulatedSolutionValue,
               this.simulationRuns,
               this.getMILPSolutionValue(),
               this.getMILPOptimalityGap(),
               this.lazyCuts,
               this.linearizationSamples,
               milpMaxLinearizationError,
               simulatedLinearizationError,
               cplexSolutionTimeMs,
               simplexIterations,
               exploredNodes
               );     
      } else {
         System.out.println("No solution!");
      } 
      cplex.end();
      System.gc();
      
      return solvedInstance;
   }
   
   /**
    * https://rma350.github.io/2012/06/16/from-separation-to-optimization-in-cplex.html
    */
   private class LPNLPLazyConstraintCallback extends IloCplex.LazyConstraintCallback{
      @Override
      protected void main() throws IloException {
         double[] kp = this.getValues(knapsack);
         kp = Arrays.stream(kp).map(v -> Math.abs(v)).toArray();
         LPNLPLazyCut cut = new LPNLPLazyCut(kp, 
               computeLX(instance, kp, linearizationSamples), 
               computeDirectionalDerivative(instance, kp, linearizationSamples));
         IloNumExpr[] sum = new IloNumExpr[kp.length];
         for(int i = 0; i < cut.X.length; i++)
            sum[i] = cplex.sum(-cut.getDirectionalDerivative()[i]*cut.X[i], cplex.prod(cut.getDirectionalDerivative()[i], knapsack[i]));
         IloNumExpr rhs = cplex.sum(cut.getKnapsackLoss(), cplex.sum(sum));
         add(cplex.ge(P, rhs), IloCplex.CutManagement.UseCutForce);
         lazyCuts++;
      }
   }
  
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstanceLarge();
      
      int linearizationSamples = 1000;
      int simulationRuns = 100000;
      
      try {
         SKPGenericDistributionLazyCuts sskp = new SKPGenericDistributionLazyCuts(instance, linearizationSamples, simulationRuns);
         
         SKPGenericDistributionCutsSolvedInstance solvedInstance = sskp.solve();
         System.out.println(GSONUtility.<SKPGenericDistributionCutsSolvedInstance>printInstanceAsJSON(solvedInstance));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}

class LPNLPLazyCut {
   double[] X;  // Knapsack
   double LX;   // Value of loss function at X
   double[] dd; // Directional derivative
   
   public LPNLPLazyCut(double[] X, double LX, double[] dd) {
      this.X = X;
      this.LX = LX;
      this.dd = dd;
   }
   
   public double[] getKnapsack() {
      return X;
   }
   
   public double getKnapsackLoss() {
      return this.LX;
   }
   
   public double[] getDirectionalDerivative() {
      return dd;
   }
}
