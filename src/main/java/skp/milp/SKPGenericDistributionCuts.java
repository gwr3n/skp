package skp.milp;

import java.util.ArrayList;
import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.opl.*;

import skp.folf.FirstOrderLossFunctionScalarProduct;
import skp.instance.SKPGenericDistribution;
import skp.milp.instance.SKPGenericDistributionCutsSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.Distribution;

public class SKPGenericDistributionCuts {    
   
   int linearizationSamples;
   int simulationRuns;
   int maxCuts;
   
   int[] optimalKnapsack;
   int[] lastKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   SKPGenericDistribution instance;
   
   double cplexSolutionTimeMs = 0;
   int simplexIterations = 0;
   int exploredNodes = 0;
   
   ArrayList<LPNLPCut> cutList = new ArrayList<LPNLPCut>();
   
   public SKPGenericDistributionCuts(SKPGenericDistribution instance, int linearizationSamples, int maxCuts, int simulationRuns){
      this.instance = instance;
      this.linearizationSamples = linearizationSamples;
      this.maxCuts = maxCuts;
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
   
   private static final double step = 1e-4;
   
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
   
   private static long time_limitMs = 60*10*1000; //10 minutes
   private static double tolerance = 1e-4; // Equivalent to CPLEX https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-relative-mip-gap-tolerance
   
   public SKPGenericDistributionCutsSolvedInstance solve() throws IloException {
      long startGlobal = System.currentTimeMillis();
      this.milpSolutionValue = Double.MAX_VALUE;
      this.lastKnapsack = null;
      
      //long start = System.currentTimeMillis();
      //generateInitialCut();
      //System.out.println(System.currentTimeMillis() - start);
      
      SKPGenericDistributionCutsSolvedInstance solvedInstance = null;
      IloOplFactory.setDebugMode(false);
      boolean stop = false;
      while(!stop) {
         IloCplex cplex = new IloCplex();
         cplex.setParam(IloCplex.Param.TimeLimit, 60*10); //10 minutes
         cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, tolerance); 
         cplex.setParam(IloCplex.Param.Threads, 8);
         cplex.setParam(IloCplex.Param.MIP.Display, 2);
         
         /*
          * Disable all cuts and preprocessing
          * 
         cplex.setParam(IloCplex.Param.MIP.PolishAfter.AbsMIPGap, 0);
         cplex.setParam(IloCplex.Param.MIP.PolishAfter.MIPGap, 0);
         
         cplex.setParam(IloCplex.Param.Preprocessing.Presolve, false);
         cplex.setParam(IloCplex.Param.Preprocessing.Aggregator, 1);
         cplex.setParam(IloCplex.Param.Preprocessing.Symmetry, 0);
         
         cplex.setParam(    IloCplex.Param.MIP.Cuts.BQP, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.Cliques, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.Covers, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.Disjunctive, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.FlowCovers, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.Gomory, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.GUBCovers, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.Implied, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.LiftProj, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.LocalImplied, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.MCFCut, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.MIRCut, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.Nodecuts, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.PathCut, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.RLT, -1);
         cplex.setParam(    IloCplex.Param.MIP.Cuts.ZeroHalfCut, -1);
         */
         
         //cplex.setParam(IloCplex.Param.MIP.Strategy.Search, IloCplex.MIPSearch.Traditional);
         
         cplex.setOut(null);

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
         
         //Cuts
         for(LPNLPCut cut: this.cutList) {
            IloNumExpr[] sum = new IloNumExpr[this.instance.getItems()];
            for(int i = 0; i < cut.X.length; i++)
               sum[i] = cplex.sum(-cut.getDirectionalDerivative()[i]*cut.X[i], cplex.prod(cut.getDirectionalDerivative()[i], X[i]));
            IloNumExpr rhs = cplex.sum(cut.getKnapsackLoss(), cplex.sum(sum));
            cplex.add(cplex.ge(P, rhs));
         }

         // Bound objective function value
         cplex.add(cplex.le(cplex.sum(
               cplex.scalProd(this.instance.getExpectedValues(), X),
               cplex.prod(-instance.getShortageCost(), P)), this.milpSolutionValue));
         
         // The objective is to maximize the profit minus the expected capacity shortage
         cplex.addMaximize(cplex.sum(
               cplex.scalProd(this.instance.getExpectedValues(), X),
               cplex.prod(-instance.getShortageCost(), P))
               );

         boolean status =  cplex.solve();
         if ( status ) {   
            this.milpSolutionValue = cplex.getObjValue();
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
            
            //Update optimality gap
            double objValue = 0;
            for(int i = 0; i < instance.getItems(); i++){
               objValue += this.optimalKnapsack[i]*instance.getExpectedValues()[i];
            }
            objValue -= instance.getShortageCost()*computeLX(instance, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(), linearizationSamples);
            this.milpOptimalityGap = (cplex.getObjValue() - objValue)/(1e-10 + objValue);

            if(Arrays.equals(this.optimalKnapsack, this.lastKnapsack) || 
                  this.cutList.size() > this.maxCuts ||
                  System.currentTimeMillis() - startGlobal >= time_limitMs ||
                  this.milpOptimalityGap < tolerance) {
               stop = true;
               double endGlobal = System.currentTimeMillis();
               this.cplexSolutionTimeMs = endGlobal - startGlobal;
               
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
                     this.cutList.size(),
                     this.linearizationSamples,
                     milpMaxLinearizationError,
                     simulatedLinearizationError,
                     this.cplexSolutionTimeMs,
                     simplexIterations,
                     exploredNodes
                     );
            } else {               
               this.lastKnapsack = this.optimalKnapsack;
               //New cut
               LPNLPCut cut = new LPNLPCut(Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(), 
                                           computeLX(instance, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(), linearizationSamples), 
                                           computeDirectionalDerivative(instance, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(), linearizationSamples));
               this.cutList.add(cut);
               //if(this.cutList.size() % 10 == 0) System.out.println("Cuts: "+this.cutList.size()+"\tObj:"+this.milpSolutionValue);
            }
         } else {
            System.out.println("No solution!");
            stop = true;
         } 
         cplex.end();
         System.gc();
      }
      return solvedInstance;
   }

   @SuppressWarnings("unused")
   private void generateInitialCut() {
      this.optimalKnapsack = new int[instance.getItems()];
      int low = 0;
      int high = instance.getItems() - 1;

      while (low <= high) {
          int mid = (low + high) / 2;
          this.optimalKnapsack = new int[instance.getItems()];
          Arrays.fill(this.optimalKnapsack, 0, mid + 1, 1);
          SimulateGenericDistribution sim = new SimulateGenericDistribution(instance);
          double simulatedSolutionValueR = sim.simulate(optimalKnapsack, this.simulationRuns);
          this.optimalKnapsack[mid] = 0;
          double simulatedSolutionValueL = sim.simulate(optimalKnapsack, this.simulationRuns);
          
          if (simulatedSolutionValueR > simulatedSolutionValueL) {
              low = mid + 1; // Move to the right half
          } else {
              high = mid - 1; // Move to the left half
          }
      }
      this.lastKnapsack = this.optimalKnapsack;
      // New cut
      LPNLPCut cut = new LPNLPCut(
          Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(),
          computeLX(instance, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(), linearizationSamples),
          computeDirectionalDerivative(instance, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(), linearizationSamples)
      );
      this.cutList.add(cut);
   }
   
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstanceLarge();
      
      int linearizationSamples = 1000;
      int maxCuts = 1000;
      int simulationRuns = 100000;
      
      try {
         SKPGenericDistributionCuts sskp = new SKPGenericDistributionCuts(instance, linearizationSamples, maxCuts, simulationRuns);
         
         SKPGenericDistributionCutsSolvedInstance solvedInstance = sskp.solve();
         System.out.println(GSONUtility.<SKPGenericDistributionCutsSolvedInstance>printInstanceAsJSON(solvedInstance));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}

class LPNLPCut {
   double[] X;  // Knapsack
   double LX;   // Value of loss function at X
   double[] dd; // Directional derivative
   
   public LPNLPCut(double[] X, double LX, double[] dd) {
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
