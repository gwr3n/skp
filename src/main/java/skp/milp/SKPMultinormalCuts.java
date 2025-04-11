package skp.milp;

import java.util.ArrayList;
import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.opl.*;

import skp.folf.FirstOrderLossFunctionScalarProductMVN;
import skp.instance.SKPMultinormal;
import skp.milp.instance.SKPGenericDistributionCutsMVNSolvedInstance;
import skp.sim.SimulateMultinormal;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdistmulti.MultiNormalDist;

public class SKPMultinormalCuts {    
   
   int simulationRuns;
   int maxCuts;
   
   int[] optimalKnapsack;
   int[] lastKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   SKPMultinormal instance;
   boolean independentDemand;
   
   double cplexSolutionTimeMs = 0;
   int simplexIterations = 0;
   int exploredNodes = 0;
   
   ArrayList<LPNLPCutMVN> cutList = new ArrayList<LPNLPCutMVN>();
   
   public SKPMultinormalCuts(SKPMultinormal instance, int maxCuts, int simulationRuns){
      this.instance = instance;
      this.independentDemand = independentDemand(instance);
      this.maxCuts = maxCuts;
      this.simulationRuns = simulationRuns;
   }
   
   static boolean independentDemand(SKPMultinormal instance) {
      double[][] cov = instance.getWeights().getCovariance();
      for(int i = 0; i < cov.length; i++) {
         for(int j = 0; j < cov[i].length; j++) {
            if(i != j && cov[i][j] != 0)
               return false;
         }
      }
      return true;
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
   
   private static double[] computeDirectionalDerivative(SKPMultinormal instance, boolean independentDemand, double[] knapsack) {
      MultiNormalDist weights = instance.getWeights();
      FirstOrderLossFunctionScalarProductMVN folfsp = new FirstOrderLossFunctionScalarProductMVN(weights, independentDemand);
      double[] dd = new double[weights.getMean().length];
      for(int i = 0; i < dd.length; i++) {
         double[] kp = Arrays.copyOf(knapsack, knapsack.length);
         kp[i] = kp[i] + step;
         dd[i] = (folfsp.getFirstOrderLossFunctionValue(instance.getCapacity(), kp) -
                  folfsp.getFirstOrderLossFunctionValue(instance.getCapacity(), knapsack))/step;
      }
      return dd;
   }
   
   private static double computeLX(SKPMultinormal instance, boolean independentDemand, double[] knapsack) {
      MultiNormalDist weights = instance.getWeights();
      FirstOrderLossFunctionScalarProductMVN folfsp = new FirstOrderLossFunctionScalarProductMVN(weights, independentDemand);
      return folfsp.getFirstOrderLossFunctionValue(instance.getCapacity(), knapsack);
   }
   
   private static long time_limitMs = 60*10*1000; //10 minutes
   private static double tolerance = 1e-2;
   
   public SKPGenericDistributionCutsMVNSolvedInstance solve() throws IloException {
      long startGlobal = System.currentTimeMillis();
      this.milpSolutionValue = Double.MAX_VALUE;
      this.lastKnapsack = null;
      
      SKPGenericDistributionCutsMVNSolvedInstance solvedInstance = null;
      IloOplFactory.setDebugMode(false);
      boolean stop = false;
      while(!stop) {
         IloCplex cplex = new IloCplex();
         cplex.setParam(IloCplex.Param.TimeLimit, 60*10); //10 minutes
         cplex.setParam(IloCplex.Param.Threads, 8);
         cplex.setParam(IloCplex.Param.MIP.Display, 2);
         //cplex.setParam(IloCplex.Param.Simplex.DGradient, 1);
         //cplex.setParam(IloCplex.Param.Preprocessing.Presolve, false);
         //cplex.setParam(IloCplex.Param.Preprocessing.Aggregator, 1);
         //cplex.setParam(IloCplex.Param.Preprocessing.Symmetry, 0);
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
         cplex.addEq(cplex.scalProd(this.instance.getWeights().getMean(), X), M);

         // Expected capacity shortage (loss): P >= - (C-M)
         //                                    P-M >= -C
         cplex.addGe(cplex.sum(P,cplex.prod(M, -1)), -this.instance.getCapacity());
         
         //Cuts
         for(LPNLPCutMVN cut: this.cutList) {
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
            objValue -= instance.getShortageCost()*computeLX(instance, independentDemand, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray());
            this.milpOptimalityGap = (cplex.getObjValue() - objValue)/objValue;

            if(Arrays.equals(this.optimalKnapsack, this.lastKnapsack) || 
                  this.cutList.size() > this.maxCuts ||
                  System.currentTimeMillis() - startGlobal >= time_limitMs ||
                  this.milpOptimalityGap < tolerance) {
               stop = true;
               double endGlobal = System.currentTimeMillis();
               this.cplexSolutionTimeMs = endGlobal - startGlobal;
               
               SimulateMultinormal sim = new SimulateMultinormal(instance);
               double simulatedSolutionValue = sim.simulate(optimalKnapsack, this.simulationRuns);
               double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;

               solvedInstance = new SKPGenericDistributionCutsMVNSolvedInstance(
                     instance,
                     this.optimalKnapsack,
                     simulatedSolutionValue,
                     this.simulationRuns,
                     this.getMILPSolutionValue(),
                     this.getMILPOptimalityGap(),
                     this.cutList.size(),
                     milpMaxLinearizationError,
                     simulatedLinearizationError,
                     this.cplexSolutionTimeMs,
                     simplexIterations,
                     exploredNodes
                     );
            } else {               
               this.lastKnapsack = this.optimalKnapsack;
               //New cut
               LPNLPCutMVN cut = new LPNLPCutMVN(Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray(), 
                                           computeLX(instance, independentDemand, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray()), 
                                           computeDirectionalDerivative(instance, independentDemand, Arrays.stream(this.optimalKnapsack).asDoubleStream().toArray()));
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
   
   public static void main(String args[]) {

      SKPMultinormal instance = SKPMultinormal.getTestInstance();
      
      int maxCuts = 1000;
      int simulationRuns = 10000;
      
      try {
         SKPMultinormalCuts sskp = new SKPMultinormalCuts(instance, maxCuts, simulationRuns);
         
         SKPGenericDistributionCutsMVNSolvedInstance solvedInstance = sskp.solve();
         System.out.println(GSONUtility.<SKPGenericDistributionCutsMVNSolvedInstance>printInstanceAsJSON(solvedInstance));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}

class LPNLPCutMVN {
   double[] X;  // Knapsack
   double LX;   // Value of loss function at X
   double[] dd; // Directional derivative
   
   public LPNLPCutMVN(double[] X, double LX, double[] dd) {
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
