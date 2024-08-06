package skp.milp;

import java.util.stream.IntStream;

import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.opl.IloCplex;

import skp.folf.PiecewiseFirstOrderLossFunction;
import skp.instance.SKPGenericDistribution;
import skp.milp.instance.SKPGenericDistributionBandBSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.probdist.Distribution;

import java.util.Arrays;
import java.util.Stack;

public class SKPGenericDistributionBandB {

   int linearizationSamples;
   int simulationRuns;
   
   int[] optimalKnapsack;
   double maxProfit = 0;
   
   double BandBSolutionTimeMs = 0;
   int exploredNodes = 0;
   
   SKPGenericDistribution instance;
   
   public SKPGenericDistributionBandB(SKPGenericDistribution instance, int linearizationSamples, int simulationRuns){
      this.instance = instance;
      this.linearizationSamples = linearizationSamples;
      this.simulationRuns = simulationRuns;
      this.optimalKnapsack = new int[instance.getItems()];
   }
   
   public int[] getOptimalKnapsack() {
      return this.optimalKnapsack;
   }
   
   public double getMILPSolutionValue() {
      return this.maxProfit;
   }
   
   private static double computeCutRHS(SKPGenericDistribution instance, int[] knapsack, int linearizationSamples) {
      Distribution[] weights = instance.getWeights();
      Distribution[] reducedWeights = IntStream.iterate(0, i -> i + 1)
            .limit(weights.length)
            .filter(i -> knapsack[i] == 1)
            .mapToObj(i -> weights[i])
            .toArray(Distribution[]::new);

      PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(reducedWeights);

      double a = pwfolf.getFirstOrderLossFunctionValue(instance.getCapacity(), linearizationSamples)-
            (pwfolf.getEmpiricalDistribution(linearizationSamples).cdf(instance.getCapacity())-1)*instance.getCapacity();
      double b = pwfolf.getEmpiricalDistribution(linearizationSamples).cdf(instance.getCapacity())-1; //slope
      
      return b*instance.getCapacity()+a;
   }
   
   
   // Calculate the bound value for a given node
   private static void bound(Node u, SKPGenericDistribution instance, int linearizationSamples, int simulationRuns) throws IloException {
      IloCplex cplex = new IloCplex();
      cplex.setParam(IloCplex.Param.TimeLimit, 60);
      cplex.setParam(IloCplex.Param.Threads, 8);
      cplex.setParam(IloCplex.Param.MIP.Display, 2);
      cplex.setOut(null);

      // Create decision variables

      // Object selectors
      IloNumVar[] X = new IloNumVar[instance.getItems()];
      for(int i = 0; i <= u.level; i++)
         X[i] = cplex.numVar(u.knapsack[i], u.knapsack[i], IloNumVarType.Int);
      for(int i = u.level + 1; i < instance.getItems(); i++)
         X[i] = cplex.numVar(0, 1, IloNumVarType.Int);

      // Expected knapsack weight
      IloNumVar M = cplex.numVar(0, Double.POSITIVE_INFINITY);

      // Expected capacity shortage
      IloNumVar P = cplex.numVar(0, Double.POSITIVE_INFINITY);

      // Create constraints

      // Expected knapsack weight
      cplex.addEq(cplex.scalProd(Arrays.stream(instance.getWeights()).mapToDouble(w -> w.getMean()).toArray(), X), M);

      // Expected capacity shortage (loss): P >= - (C-M)
      //                                    P-M >= -C
      cplex.addGe(cplex.sum(P,cplex.prod(M, -1)), -instance.getCapacity());

      //Cut
      Cut cut = new Cut(u.knapsack, computeCutRHS(instance, u.knapsack, linearizationSamples));
      cplex.add(cplex.ge(P, cut.getRHS()));
      
      // The objective is to maximize the profit minus the expected capacity shortage
      cplex.addMaximize(cplex.sum(
            cplex.scalProd(instance.getExpectedValues(), X),
            cplex.prod(-instance.getShortageCost(), P))
            );

      cplex.setParam(IloCplex.Param.MIP.Strategy.Search, IloCplex.MIPSearch.Traditional);

      boolean status =  cplex.solve(); 
      if(status) {
         u.bound = cplex.getObjValue();
         int[] optimalKnapsack = new int[instance.getItems()];
         for(int i = 0; i < instance.getItems(); i++){
            optimalKnapsack[i] = (int) Math.round(cplex.getValue(X[i]));
         }
         u.bestKnapsack = optimalKnapsack;
         SimulateGenericDistribution sim = new SimulateGenericDistribution(instance);
         u.profit = sim.simulate(optimalKnapsack, simulationRuns);
      } else {
         System.out.println("No solution!");
      } 
      // Release all resources
      cplex.end();
      System.gc();
   }
   
   private static long time_limitMs = 60000; 
   
   // Branch and bound algorithm
   public SKPGenericDistributionBandBSolvedInstance solve() throws IloException {
       long start = System.currentTimeMillis();
       Stack<Node> stack = new Stack<>();
       Node u, v;

       int[] knapsack = new int[instance.getItems()];
       Node rootNode = new Node(-1, knapsack);
       bound(rootNode, this.instance, this.linearizationSamples, this.simulationRuns);
       stack.push(rootNode);
       exploredNodes++;

       while (!stack.isEmpty()) {
           u = stack.pop();

           if (u.level == instance.getItems() - 1) {
               continue;
           }
           
           int[] branch_1 = Arrays.copyOf(knapsack, instance.getItems());
           branch_1[u.level + 1] = 1;
           v = new Node(u.level + 1, branch_1);
           bound(v, this.instance, this.linearizationSamples, this.simulationRuns);

           if (v.profit > maxProfit) {
               maxProfit = v.profit;
               optimalKnapsack = v.bestKnapsack;
           }

           if (v.bound > maxProfit) {
               stack.push(v);
               exploredNodes++;
           }

           int[] branch_0 = Arrays.copyOf(knapsack, instance.getItems());
           branch_0[u.level + 1] = 0;
           v = new Node(u.level + 1, branch_0);
           bound(v, this.instance, this.linearizationSamples, this.simulationRuns);

           if (v.bound > maxProfit) {
               stack.push(v);
               exploredNodes++;
           }
           
           if(exploredNodes % 10 == 0) {
              System.out.println("Explored nodes: "+exploredNodes);
              System.out.println("Solution time: "+(System.currentTimeMillis() - start));
           }
           if(System.currentTimeMillis() - start >= time_limitMs)
              break;
       }
       long end = System.currentTimeMillis();
       return new SKPGenericDistributionBandBSolvedInstance(this.instance, this.optimalKnapsack, this.maxProfit, this.simulationRuns, end-start, this.exploredNodes);
   }
   
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstance();
      
      int linearizationSamples = 10000;
      int simulationRuns = 100000;
      
      try {
         SKPGenericDistributionBandB sskp = new SKPGenericDistributionBandB(instance, linearizationSamples, simulationRuns);
         
         SKPGenericDistributionBandBSolvedInstance solvedInstance = sskp.solve();
         System.out.println(GSONUtility.<SKPGenericDistributionBandBSolvedInstance>printInstanceAsJSON(solvedInstance));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}

class Node {
   int level;
   int[] knapsack;
   
   int[] bestKnapsack;
   double bound;
   double profit;

   public Node(int level, int[] knapsack) {
       this.level = level;
       this.knapsack = knapsack;
   }
}

