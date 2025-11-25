package skp.milp;

import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloRange;
import ilog.opl.IloCplex;
import ilog.opl.IloOplFactory;
import skp.folf.FirstOrderLossFunctionScalarProductMVN;
import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.SKPMultinormal;
import skp.milp.instance.SKPMultinormalCutsSolvedInstance;
import skp.sim.SimulateMultinormal;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdistmulti.MultiNormalDist;

public class SKPMultinormalLazyCuts {
   
   int simulationRuns;
   int lazyCuts;
   
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
   
   IloNumVar[] knapsack;
   IloNumVar P;
   IloCplex cplex;
   
   public SKPMultinormalLazyCuts(SKPMultinormal instance, int simulationRuns){
      this.instance = instance;
      this.independentDemand = independentDemand(instance);
      this.simulationRuns = simulationRuns;
      this.warmStart = false;
   }
   
   public SKPMultinormalLazyCuts(SKPMultinormal instance, int simulationRuns, int warmStartPartitions){
      this.instance = instance;
      this.independentDemand = independentDemand(instance);
      this.simulationRuns = simulationRuns;
      this.warmStart = true;
      this.warmStartPartitions = warmStartPartitions;
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
   
   private static final double step = 1e-4;
   
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
   
   private static long time_limit = 60*10; //10 minutes
   private static double tolerance = 1e-4; // // Equivalent to CPLEX https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-relative-mip-gap-tolerance
   private boolean warmStart;
   private int warmStartPartitions;
   
   public SKPMultinormalCutsSolvedInstance solve() throws IloException {
      long startGlobal = System.currentTimeMillis();
      this.milpSolutionValue = Double.MAX_VALUE;
      this.lastKnapsack = null;
      this.lazyCuts = 0;
      
      SKPMultinormalCutsSolvedInstance solvedInstance = null;
      IloOplFactory.setDebugMode(false);
      
      cplex = new IloCplex();
      cplex.setParam(IloCplex.Param.TimeLimit, time_limit);
      cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, tolerance);
      cplex.setParam(IloCplex.Param.Threads, 1); // Lazy cuts do not allow multithreading      
      cplex.setParam(IloCplex.Param.MIP.Strategy.Search, IloCplex.MIPSearch.Traditional); // we must deactivate dynamic search (done automatically)
      cplex.setParam(IloCplex.Param.Preprocessing.Presolve, false); // we must deactivate the dual presolve
      cplex.setParam(IloCplex.Param.MIP.Strategy.PresolveNode, 0); // we must deactivate the dual presolve
      cplex.setParam(IloCplex.Param.MIP.Display, 2);
      cplex.setOut(null);

      // Create decision variables

      // Object selectors
      IloNumVar[] X = cplex.numVarArray(this.instance.getItems(), 0, 1, IloNumVarType.Bool);
      knapsack = X;

      // Expected knapsack weight
      IloNumVar M = cplex.numVar(0, Double.POSITIVE_INFINITY);
      
      // Knapsack weight variance
      IloNumVar V = cplex.numVar(0, Double.POSITIVE_INFINITY);
      
      // Expected capacity shortage
      P = cplex.numVar(0, Double.POSITIVE_INFINITY);

      // Create constraints

      // Expected knapsack weight
      cplex.addEq(cplex.scalProd(this.instance.getWeights().getMean(), X), M);

      // Expected capacity shortage (loss): P >= - (C-M)
      //                                    P-M >= -C
      cplex.addGe(cplex.sum(P,cplex.prod(M, -1)), -this.instance.getCapacity());
      
      
      // ************** Warm start (add Jensens cut) START ********************
      
      if(warmStart) {
         // Knapsack weight variance
         double[][] cov = this.instance.getWeights().getCovariance();
         IloNumExpr varExpr;
         int n = cov.length;

         if (independentDemand) {
             // Variance simplifies to sum_i Var_i * X_i (diagonal only)
             double[] diag = new double[n];
             for (int i = 0; i < n; i++) {
                 diag[i] = cov[i][i];
             }
             varExpr = cplex.scalProd(diag, X);
         } else {
             // Full quadratic form: X^T * Cov * X = sum_{i,j} cov[i][j] * X_i * X_j
             varExpr = cplex.numExpr();
             for (int i = 0; i < n; i++) {
                 for (int j = i + 1; j < n; j++) {
                     double cij = cov[i][j];
                     if (cij != 0.0) {
                         varExpr = cplex.sum(varExpr, cplex.prod(2 * cij, cplex.prod(X[i], X[j])));
                     }
                 }
             }
         }

         // V = variance of selected weights
         cplex.addGe(V, varExpr);
         
         // Knapsack weight standard deviation via piecewise linearization of sqrt(V)
         double s = 10e-2;                     // linearisation step
         double x0 = Math.sqrt(s)/4;           // value at 0 (same role as x0 in .mod)

         // Upper bound on V
         double totalVarianceUpper = 0.0;
         for (int i = 0; i < cov.length; i++) {
             for (int j = 0; j < cov[i].length; j++) {
                 totalVarianceUpper += cov[i][j];
             }
         }
         int nbPartitions = (int) Math.ceil(totalVarianceUpper / s);

         // Breakpoints and function values
         double[] xBreaks = new double[nbPartitions + 1];
         double[] yVals   = new double[nbPartitions + 1];
         xBreaks[0] = 0.0;
         yVals[0]   = x0;
         for (int b = 1; b <= nbPartitions; b++) {
             xBreaks[b] = b * s;
             yVals[b]   = x0 + Math.sqrt(xBreaks[b]);
         }

         // Standard deviation variable
         IloNumVar S = cplex.numVar(0.0, Double.POSITIVE_INFINITY, "S");

         // Extreme slopes (left, right)
         double slopeLeft = (yVals[1] - yVals[0]) / (xBreaks[1] - xBreaks[0]); // = 1 / sqrt(s)
         double slopeRight = 0.0;

         // Piecewise linear approximation: S == pwl(V)
         IloNumExpr pwS = cplex.piecewiseLinear(V, xBreaks, yVals, slopeLeft, slopeRight);
         cplex.addEq(S, pwS);
         
         double[] prob = PiecewiseStandardNormalFirstOrderLossFunction.getProbabilities(warmStartPartitions);    // size = nbPartitions
         double[] means = PiecewiseStandardNormalFirstOrderLossFunction.getMeans(warmStartPartitions);           // size = nbPartitions (or item means per partition)
         double error = 0;// scalar error term
         
         // Auxiliary scaling variable S (define bounds as needed)
         // IloNumVar S = cplex.numVar(0.0, Double.POSITIVE_INFINITY, "S");
         
         // For each partition p build: P >= (sum_{k=p..nbPartitions-1} prob[k]*means[k]) * S
         //    - sum_{k=p..nbPartitions-1} prob[k]*(C - M) + error * S
         for (int p = 0; p < warmStartPartitions; p++) {
            double sumProbMeans = 0.0;
            double sumProb = 0.0;
            for (int k = p; k < warmStartPartitions; k++) {
               sumProbMeans += prob[k] * means[k];
               sumProb += prob[k];
            }
            IloLinearNumExpr rhs = cplex.linearNumExpr();
            // (sum prob[k]*means[k]) * S
            rhs.addTerm(sumProbMeans, S);
            // - sum prob[k]*(C - M) = -sumProb*C + sumProb*M
            rhs.addTerm(sumProb, M); // + sumProb * M
            double constantPart = -sumProb * instance.getCapacity(); // - sumProb * C
            // + error * S
            rhs.addTerm(error, S);
            
            // Final: P >= rhs + constantPart
            cplex.addGe(P, cplex.sum(rhs, constantPart));
         }
         
         //Global lower bound: P >= error * S
         cplex.addGe(P, cplex.prod(error, S));
      }
      
      // ************** Warm start (add Jensens cut) END ********************
      

      // The objective is to maximize the profit minus the expected capacity shortage
      
      IloLinearNumExpr objExpr = cplex.linearNumExpr();
      objExpr.add(cplex.scalProd(instance.getExpectedValues(), X));
      objExpr.addTerm(-instance.getShortageCost(), P);
      
      // incumbent bounding constraint :  objExpr  <=  bestKnown
      milpSolutionValue = Double.POSITIVE_INFINITY;   // nothing known yet
      IloRange incumbentBound = cplex.addLe(objExpr, milpSolutionValue);
      
      cplex.addMaximize(objExpr);
      
      //Lazy Cuts
      cplex.use(new LPNLPLazyConstraintCallback());
      //cplex.use(new LPNLPUserCutCallback());
      cplex.use(new IncumbentUpdater(incumbentBound));

      boolean status =  cplex.solve();      
      if ( status ) {   
         this.milpSolutionValue = cplex.getObjValue();
         this.milpOptimalityGap = cplex.getMIPRelativeGap();
         double endGlobal = System.currentTimeMillis();
         this.cplexSolutionTimeMs = endGlobal - startGlobal;
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

         SimulateMultinormal sim = new SimulateMultinormal(instance);
         double simulatedSolutionValue = sim.simulate(optimalKnapsack, this.simulationRuns);
         double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;

         solvedInstance = new SKPMultinormalCutsSolvedInstance(
               instance,
               this.optimalKnapsack,
               simulatedSolutionValue,
               this.simulationRuns,
               this.getMILPSolutionValue(),
               this.getMILPOptimalityGap(),
               this.lazyCuts,
               milpMaxLinearizationError,
               simulatedLinearizationError,
               this.cplexSolutionTimeMs,
               simplexIterations,
               exploredNodes
               );     
      } else {
         System.out.println("**** No solution! ****");
         this.milpSolutionValue = Double.NaN;
         this.milpOptimalityGap = Double.NaN;
         double endGlobal = System.currentTimeMillis();
         this.cplexSolutionTimeMs = endGlobal - startGlobal;
         this.simplexIterations += cplex.getNiterations();
         this.exploredNodes += cplex.getNnodes();

         this.optimalKnapsack = new int[instance.getItems()];
         //DEBUG
         //System.out.println("X: "+Arrays.toString(this.optimalKnapsack));
         //System.out.println("M: "+cplex.getValue(M));
         //System.out.println("P: "+cplex.getValue(P));

         this.milpMaxLinearizationError = 0; // cuts are tight

         double simulatedSolutionValue = Double.NaN;
         double simulatedLinearizationError = Double.NaN;

         solvedInstance = new SKPMultinormalCutsSolvedInstance(
               instance,
               this.optimalKnapsack,
               simulatedSolutionValue,
               this.simulationRuns,
               this.getMILPSolutionValue(),
               this.getMILPOptimalityGap(),
               this.lazyCuts,
               milpMaxLinearizationError,
               simulatedLinearizationError,
               this.cplexSolutionTimeMs,
               simplexIterations,
               exploredNodes
               );
      } 
      cplex.end();
      System.gc();
      
      return solvedInstance;
   }
   
   /**
    * Lazy-constraint callback that closes the linearisation gap at integer solutions.
    * https://rma350.github.io/2012/06/16/from-separation-to-optimization-in-cplex.html
    */
   private class LPNLPLazyConstraintCallback extends IloCplex.LazyConstraintCallback{
      @Override
      protected void main() throws IloException {
         double[] kp = this.getValues(knapsack);
         kp = Arrays.stream(kp).map(v -> Math.abs(v)).toArray();
         LPNLPLazyCut cut = new LPNLPLazyCut(kp, 
               computeLX(instance, independentDemand, kp), 
               computeDirectionalDerivative(instance, independentDemand, kp));
         
         double pVal = getValue(P);
         if(cut.getKnapsackLoss() - pVal <= tolerance) return;
         
         IloNumExpr[] sum = new IloNumExpr[kp.length];
         for(int i = 0; i < cut.X.length; i++)
            sum[i] = cplex.sum(-cut.getDirectionalDerivative()[i]*cut.X[i], cplex.prod(cut.getDirectionalDerivative()[i], knapsack[i]));
         IloNumExpr rhs = cplex.sum(cut.getKnapsackLoss(), cplex.sum(sum));
         add(cplex.ge(P, rhs), IloCplex.CutManagement.UseCutForce);
         lazyCuts++;
      }
   }
   
   /********************************************************************
    *  Outer-approximation cuts at every LP solution (fractional or not)
    ********************************************************************/
   @SuppressWarnings("unused")
   private class LPNLPUserCutCallback extends IloCplex.UserCutCallback {
       @Override
       protected void main() throws IloException {
          double[] kp = this.getValues(knapsack);
          kp = Arrays.stream(kp).map(v -> Math.abs(v)).toArray();
          LPNLPLazyCut cut = new LPNLPLazyCut(kp, 
                computeLX(instance, independentDemand, kp), 
                computeDirectionalDerivative(instance, independentDemand, kp));
          
          double pVal = getValue(P);
          if(cut.getKnapsackLoss() - pVal <= tolerance) return;
          
          IloNumExpr[] sum = new IloNumExpr[kp.length];
          for(int i = 0; i < cut.X.length; i++)
             sum[i] = cplex.sum(-cut.getDirectionalDerivative()[i]*cut.X[i], cplex.prod(cut.getDirectionalDerivative()[i], knapsack[i]));
          IloNumExpr rhs = cplex.sum(cut.getKnapsackLoss(), cplex.sum(sum));
          add(cplex.ge(P, rhs), IloCplex.CutManagement.UseCutForce);
          lazyCuts++;
       }
   }
   
   /**************************************************************************
    *  Updates the global best value and tightens the bounding constraint
    *  whenever CPLEX finds a new incumbent.
    **************************************************************************/
   private class IncumbentUpdater extends IloCplex.IncumbentCallback {

       private final IloRange bound;

       IncumbentUpdater(IloRange incumbentBound) {
           this.bound = incumbentBound;
       }

       @Override
       protected void main() throws IloException {

           double newObj = getObjValue();

           /* strict improvement ?  (use a small tolerance to avoid loops) */
           if (newObj > milpSolutionValue + 1e-6) { // Same as CPLEX absolute MIP tolerance
               milpSolutionValue = newObj;

               /* tighten the right–hand side of  objExpr <= bestInc        */
               bound.setUB(newObj);

               /* optional: also tell CPLEX to cut off nodes whose best
                  bound cannot beat the new incumbent                       */
               cplex.setParam(IloCplex.Param.MIP.Tolerances.LowerCutoff, newObj - 1e-6);
           }
       }
   }
  
   public static void main(String args[]) {

      SKPMultinormal instance = SKPMultinormal.getTestInstance();
      
      int simulationRuns = 100000;
      
      try {
         SKPMultinormalLazyCuts sskp = new SKPMultinormalLazyCuts(instance, simulationRuns);
         
         SKPMultinormalCutsSolvedInstance solvedInstance = sskp.solve();
         System.out.println(GSONUtility.<SKPMultinormalCutsSolvedInstance>printInstanceAsJSON(solvedInstance));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}

class LPNLPLazyCutMVN {
   double[] X;  // Knapsack
   double LX;   // Value of loss function at X
   double[] dd; // Directional derivative
   
   public LPNLPLazyCutMVN(double[] X, double LX, double[] dd) {
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
