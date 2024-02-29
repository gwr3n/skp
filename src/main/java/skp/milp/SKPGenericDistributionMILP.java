package skp.milp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.opl.*;

import skp.folf.PiecewiseFirstOrderLossFunction;
import skp.instance.SKPGenericDistribution;
import skp.milp.instance.SKPGenericDistributionMILPSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.Distribution;

public class SKPGenericDistributionMILP { 
   private static final long[] seed = {12345,54321,21435,53412,54321,14235};
   
   int linearizationSamples;
   
   int[] optimalKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   SKPGenericDistribution instance;
   
   double cplexSolutionTimeMs = 0;
   int simplexIterations = 0;
   int exploredNodes = 0;
   int[] lastKnapsack = null;
   
   ArrayList<Cut> cutList = new ArrayList<Cut>();
   
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
   
   private double computeCutRHS(int[] knapsack) {
      Distribution[] weights = this.instance.getWeights();
      Distribution[] reducedWeights = IntStream.iterate(0, i -> i + 1)
            .limit(weights.length)
            .filter(i -> knapsack[i] == 1)
            .mapToObj(i -> weights[i])
            .toArray(Distribution[]::new);

      PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(reducedWeights, seed);

      double a = pwfolf.getFirstOrderLossFunctionValue(this.instance.getCapacity(), linearizationSamples)-
            (pwfolf.getEmpiricalDistribution(linearizationSamples).cdf(this.instance.getCapacity())-1)*this.instance.getCapacity();
      double b = pwfolf.getEmpiricalDistribution(linearizationSamples).cdf(this.instance.getCapacity())-1; //slope
      
      return b*this.instance.getCapacity()+a;
   }
   
   public void solve() throws IloException {
      IloOplFactory.setDebugMode(false);
      
      boolean stop = false;
      while(lastKnapsack == null || !stop) {
         IloCplex cplex = new IloCplex();
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
         for(int i = 0; i < this.cutList.size(); i++) {
            Cut cut = this.cutList.get(i);
            // if all objects in the knapsack are selected, then P should be greater or equal than cut RHS
            cplex.add(cplex.ifThen(cplex.ge(cplex.scalProd(X, cut.getKnapsack()), Arrays.stream(cut.getKnapsack()).sum()), cplex.ge(P, cut.getRHS())));
         }

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

            if(Arrays.equals(this.optimalKnapsack, this.lastKnapsack)) {
               stop = true;
            } else {
               this.lastKnapsack = this.optimalKnapsack;
               this.cutList.add(new Cut(this.optimalKnapsack, this.computeCutRHS(this.optimalKnapsack)));
            }

         } else {
            System.out.println("No solution!");
         } 
         cplex.end();
         System.gc();
      }
   }
   
   public SKPGenericDistributionMILPSolvedInstance solve(int simulationRuns) throws IloException {
      SKPGenericDistributionMILPSolvedInstance solvedInstance = null;
      IloOplFactory.setDebugMode(false);
      boolean stop = false;
      while(lastKnapsack == null || !stop) {
         IloCplex cplex = new IloCplex();
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
         for(int i = 0; i < this.cutList.size(); i++) {
            Cut cut = this.cutList.get(i);
            // if all objects in the knapsack are selected, then P should be greater or equal than cut RHS
            cplex.add(cplex.ifThen(cplex.ge(cplex.scalProd(X, cut.getKnapsack()), Arrays.stream(cut.getKnapsack()).sum()), cplex.ge(P, cut.getRHS())));
         }

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

            if(Arrays.equals(this.optimalKnapsack, this.lastKnapsack)) {
               stop = true;
               
               SimulateGenericDistribution sim = new SimulateGenericDistribution(instance);
               double simulatedSolutionValue = sim.simulate(optimalKnapsack, simulationRuns);
               double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;

               solvedInstance = new SKPGenericDistributionMILPSolvedInstance(
                     instance,
                     this.getOptimalKnapsack(),
                     simulatedSolutionValue,
                     simulationRuns,
                     this.getMILPSolutionValue(),
                     this.getMILPOptimalityGap(),
                     this.cutList.size(),
                     this.linearizationSamples,
                     milpMaxLinearizationError,
                     simulatedLinearizationError,
                     cplexSolutionTimeMs,
                     simplexIterations,
                     exploredNodes
                     );
            } else {
               this.lastKnapsack = this.optimalKnapsack;
               this.cutList.add(new Cut(this.optimalKnapsack, this.computeCutRHS(this.optimalKnapsack)));
            }

         } else {
            System.out.println("No solution!");
         } 
         cplex.end();
         System.gc();
      }
      return solvedInstance;
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

class Cut {
   int[] X;
   double rhs;
   
   public Cut(int[] X, double rhs) {
      this.X = X;
      this.rhs = rhs;
   }
   
   public int[] getKnapsack() {
      return X;
   }
   
   public double getRHS() {
      return rhs;
   }
}
