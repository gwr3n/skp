package skp.milp;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;

import ilog.opl.IloCplex;

import skp.folf.FirstOrderLossFunctionScalarProduct;
import skp.instance.SKPGenericDistribution;
import skp.milp.instance.SKPGenericDistributionBandBSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.probdist.Distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Stack;

public class SKPGenericDistributionBandB {

   int linearizationSamples;
   int simulationRuns;
   
   int[] optimalKnapsack;
   double maxProfit;
   double bestUB;
   double optGap;
   
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
   
   private static final double step = 0.1;
   
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
   
   // Calculate the bound value for a given node
   private static void bound(Node u, SKPGenericDistribution instance, int linearizationSamples, int simulationRuns) throws IloException {
      IloCplex cplex = new IloCplex();
      cplex.setParam(IloCplex.Param.TimeLimit, 60*10);
      cplex.setParam(IloCplex.Param.Threads, 8);
      cplex.setParam(IloCplex.Param.MIP.Display, 2);
      cplex.setParam(IloCplex.Param.MIP.Strategy.Search, IloCplex.MIPSearch.Traditional);
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
      
      //Existing cuts
      for(LPNLPCut cut: u.cutList) {
         IloNumExpr[] sum = new IloNumExpr[u.knapsack.length];
         for(int i = 0; i < cut.X.length; i++)
            sum[i] = cplex.sum(-cut.getDirectionalDerivative()[i]*cut.X[i], cplex.prod(cut.getDirectionalDerivative()[i], X[i]));
         IloNumExpr rhs = cplex.sum(cut.getKnapsackLoss(), cplex.sum(sum));
         cplex.add(cplex.ge(P, rhs));
      }
         
      //New cut
      LPNLPCut cut = new LPNLPCut(Arrays.stream(u.knapsack).asDoubleStream().toArray(), 
                                  computeLX(instance, Arrays.stream(u.knapsack).asDoubleStream().toArray(), linearizationSamples), 
                                  computeDirectionalDerivative(instance, Arrays.stream(u.knapsack).asDoubleStream().toArray(), linearizationSamples));
      u.cutList.add(cut);
      IloNumExpr[] sum = new IloNumExpr[u.knapsack.length];
      for(int i = 0; i < u.knapsack.length; i++)
         sum[i] = cplex.sum(-cut.getDirectionalDerivative()[i]*u.knapsack[i], cplex.prod(cut.getDirectionalDerivative()[i], X[i]));
      IloNumExpr rhs = cplex.sum(cut.getKnapsackLoss(), cplex.sum(sum));
      cplex.add(cplex.ge(P, rhs));
      
      // The objective is to maximize the profit minus the expected capacity shortage
      cplex.addMaximize(cplex.sum(
            cplex.scalProd(instance.getExpectedValues(), X),
            cplex.prod(-instance.getShortageCost(), P))
            );

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
   
   private static long time_limitMs = 60*10*1000; 
   boolean logging = false;
   
   // Branch and bound algorithm
   public SKPGenericDistributionBandBSolvedInstance solve() throws IloException {
       long start = System.currentTimeMillis();
       Stack<Node> stack = new Stack<>();
       Node u, v;
       this.maxProfit = 0; 

       Node rootNode = new Node(-1, new int[instance.getItems()]);
       bound(rootNode, this.instance, this.linearizationSamples, this.simulationRuns);
       this.bestUB = rootNode.bound;
       this.optGap = (this.bestUB - this.maxProfit)/this.bestUB;
       stack.push(rootNode);
       exploredNodes++;

       while (!stack.isEmpty()) {
           u = stack.pop();

           if (u.level == instance.getItems() - 1) {
               continue;
           }
           
           int[] branch_1 = Arrays.copyOf(u.knapsack, instance.getItems());
           branch_1[u.level + 1] = 1;
           v = new Node(u.level + 1, branch_1);
           bound(v, this.instance, this.linearizationSamples, this.simulationRuns);

           if (v.profit > maxProfit) {
               this.maxProfit = v.profit;
               this.optimalKnapsack = v.bestKnapsack;
               this.optGap = (this.bestUB - this.maxProfit)/this.bestUB;
           }

           if (v.bound > maxProfit) {
               stack.push(v);
               this.exploredNodes++;
           }

           int[] branch_0 = Arrays.copyOf(u.knapsack, instance.getItems());
           branch_0[u.level + 1] = 0;
           v = new Node(u.level + 1, branch_0);
           bound(v, this.instance, this.linearizationSamples, this.simulationRuns);

           if (v.bound > maxProfit) {
               stack.push(v);
               this.exploredNodes++;
           }
           
           if(exploredNodes % 10 == 0 && logging) {
              System.out.println(instance.getInstanceID());
              System.out.println("Explored nodes: "+exploredNodes);
              System.out.println("Solution time: "+(System.currentTimeMillis() - start));
              System.out.println("UB: "+this.bestUB);
              System.out.println("Max profit: "+this.maxProfit);
              System.out.printf("Opt gap: %.2f%%\n", 100*this.optGap);
              System.out.println();
           }
           
           if(System.currentTimeMillis() - start >= time_limitMs)
              break;
       }
       long end = System.currentTimeMillis();
       if(stack.isEmpty()) {
          this.bestUB = this.maxProfit;
          this.optGap = (this.bestUB - this.maxProfit)/this.bestUB;
       }
       return new SKPGenericDistributionBandBSolvedInstance(this.instance, 
                                                            this.optimalKnapsack, 
                                                            this.maxProfit, 
                                                            this.simulationRuns, 
                                                            end-start, 
                                                            this.exploredNodes, 
                                                            this.optGap);
   }
   
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstanceLarge();
      
      int linearizationSamples = 1000;
      int simulationRuns = 1000;
      
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
   
   ArrayList<LPNLPCut> cutList;

   public Node(int level, int[] knapsack) {
       this.level = level;
       this.knapsack = knapsack;
       this.cutList = new ArrayList<LPNLPCut>();
   }
   
   @SuppressWarnings("unchecked")
   public Node(int level, int[] knapsack, ArrayList<Cut> cutList) {
      this.level = level;
      this.knapsack = knapsack;
      this.cutList = (ArrayList<LPNLPCut>) cutList.clone();
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

