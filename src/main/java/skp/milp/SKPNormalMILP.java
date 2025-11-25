package skp.milp;

import java.util.Arrays;

import ilog.concert.IloException;
import ilog.opl.IloCplex;
import ilog.opl.IloCustomOplDataSource;
import ilog.opl.IloOplDataHandler;
import ilog.opl.IloOplDataSource;
import ilog.opl.IloOplFactory;
import ilog.opl.IloOplModel;
import skp.folf.LinearisationFactory;
import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.SKPNormal;
import skp.milp.instance.SKPNormalMILPSolvedInstance;
import skp.sim.SimulateNormal;
import skp.utilities.gson.GSONUtility;

/**
 * To run from Mac OS
 * 
 * -Djava.library.path=/Applications/CPLEX_Studio1210/opl/bin/x86-64_osx/
 * 
 * Environment variable
 * 
 * DYLD_LIBRARY_PATH /Applications/CPLEX_Studio1210/opl/bin/x86-64_osx/
 * 
 * @author gwren
 *
 */

public class SKPNormalMILP extends SKPMILP{
   SKPNormal instance;
   
   PWAPPROXIMATION pwa;        // sqrt approximation       
   double s = 1e-2;            // sqrt approximation step    
   double x0 = Math.sqrt(s)/4; // sqrt approximation x0   
   
   public SKPNormalMILP(SKPNormal instance, int partitions, PWAPPROXIMATION pwa){
      this.instance = instance;
      this.partitions = partitions;
      linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model =  "sk_normal";
      this.pwa = pwa;
   }
   
   public SKPNormalMILP(SKPNormal instance, int partitions, PWAPPROXIMATION pwa, double s, double x0){
      this.instance = instance;
      this.partitions = partitions;
      linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model =  "sk_normal";
      this.pwa = pwa;   
      this.s = s;       
      this.x0 = x0;     
   }
   
   public SKPNormalMILP(SKPNormal instance, PWAPPROXIMATION pwa, double epsilon){
      this.instance = instance;
      linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model =  "sk_normal";
      this.pwa = pwa;
      
      int[] param = chooseLinearisationParameters(epsilon);
      
      this.partitions = param[0] - 1;
      double Vmax = Arrays.stream(this.instance.getWeights()).mapToDouble(w -> w.getVariance()).sum();
      this.s = Vmax/param[1];       
      this.x0 = Math.sqrt(this.s)/4; 
      
      System.out.println("Linearisation parameters: W = " + this.partitions + "; s = " + this.s);
   }
   
   IloOplDataSource getDataSource(IloOplFactory oplF) {
      return new SKPNormalMILP.MyData(oplF);
   }
   
   void computeMILPMaxLinearizationError(IloOplModel opl, IloCplex cplex) throws IloException{
      this.milpMaxLinearizationError = this.instance.getShortageCost()*
                                       cplex.getValue(opl.getElement("S").asNumVar())*
                                       PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions);
   }
   
   /*********************************************************************
    * Choose the number of linear segments so that the absolute gap of
    * both MILP relaxations is ≤ epsilon.
    *
    * return int[2]
    *   [0] : W + 1   (number of FOLF linear segments)
    *   [1] : Q       (number of √–segments)
    *
    *  Pre-requisites
    *  ---------------
    *  •  instance.getWeights() returns objects with method getVariance().
    *  •  instance.getShortageCost() returns c > 0.
    *  •  JensenMinimaxPartitioner.compute(m) returns a Partition object
    *     for a loss approximation with  m  linear segments (= m–1
    *     partitions).  The object exposes
    *        double   error   ––   the worst–case  ẑe_{m-1};
    *        double[] prob    ––   probability masses  p_k;
    *        double[] mean    ––   conditional means   E[ω|Ω_k].
    ********************************************************************/
   public int[] chooseLinearisationParameters(final double epsilon)
   {
       /* ---------- constants that depend only on the instance ---------- */
       double Vmax = Arrays.stream(instance.getWeights())
                           .mapToDouble(w -> w.getVariance())
                           .sum();
       double c    = instance.getShortageCost();
       
       /* ================================================================
        * Delegate to helper method
        * ================================================================ */
       return LinearisationFactory.chooseLinearisationParameters(epsilon, Vmax, c);
   }
   

   public static SKPNormalMILPSolvedInstance solve(SKPNormal instance, int partitions, int simulationRuns) throws IloException {      
      SKPNormalMILPSolvedInstance solvedJensens = new SKPNormalMILP(instance, partitions, PWAPPROXIMATION.JENSENS).solve(simulationRuns);
      SKPNormalMILPSolvedInstance solvedEM = new SKPNormalMILP(instance, partitions, PWAPPROXIMATION.EDMUNDSON_MADANSKI).solve(simulationRuns);
      
      double optimality_gap = 0.0;
      SKPNormalMILPSolvedInstance opt = solvedJensens;
      if(!Arrays.equals(solvedJensens.optimalKnapsack, solvedEM.optimalKnapsack)) {
         optimality_gap = (solvedJensens.milpSolutionValue - solvedEM.milpSolutionValue)/(1e-10 + solvedEM.milpSolutionValue);
         if(solvedJensens.simulatedSolutionValue < solvedEM.simulatedSolutionValue) {
            opt = solvedEM;
         } 
      }
      
      return new SKPNormalMILPSolvedInstance(
            instance,
            opt.optimalKnapsack,
            opt.simulatedSolutionValue,
            simulationRuns,
            opt.milpSolutionValue,
            optimality_gap,
            partitions,
            PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples(),
            opt.milpMaxLinearizationError,
            opt.simulatedLinearizationError,
            solvedJensens.cplexSolutionTimeMs+solvedEM.cplexSolutionTimeMs,
            solvedJensens.simplexIterations+solvedEM.simplexIterations,
            solvedJensens.exploredNodes+solvedEM.exploredNodes
            );
   }
   
   public void solve() throws IloException {
      this.solveMILP(model, instance);
   }
   
   public SKPNormalMILPSolvedInstance solve(int simulationRuns) throws IloException {
      this.solveMILP(model, instance);
      
      SimulateNormal sim = new SimulateNormal(instance);
      double simulatedSolutionValue = sim.simulate(optimalKnapsack, simulationRuns);
      
      double milpMaxLinearizationError = 100*this.getMILPMaxLinearizationError()/simulatedSolutionValue;
      double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;
      
      SKPNormalMILPSolvedInstance solvedInstance = new SKPNormalMILPSolvedInstance(
            instance,
            this.getOptimalKnapsack(),
            simulatedSolutionValue,
            simulationRuns,
            this.getMILPSolutionValue(),
            this.getMILPOptimalityGap(),
            this.partitions,
            this.linearizationSamples,
            milpMaxLinearizationError,
            simulatedLinearizationError,
            cplexSolutionTimeMs,
            simplexIterations,
            exploredNodes
            );
      
      return solvedInstance;
   }

   class MyData extends IloCustomOplDataSource
   {
      MyData(IloOplFactory oplF){
         super(oplF);
      }

      public void customRead(){
         IloOplDataHandler handler = getDataHandler();

         handler.startElement("N");
         handler.addIntItem(instance.getItems());
         handler.endElement();

         handler.startElement("expectedValues");
         handler.startArray();
         for (int j = 0 ; j<instance.getExpectedValues().length ; j++)
            handler.addNumItem(instance.getExpectedValues()[j]);
         handler.endArray();
         handler.endElement();

         handler.startElement("expectedWeights");
         handler.startArray();
         for (int j = 0 ; j<instance.getWeights().length ; j++)
            handler.addNumItem(instance.getWeights()[j].getMean());
         handler.endArray();
         handler.endElement();

         handler.startElement("varianceWeights");
         handler.startArray();
         for (int j = 0 ; j<instance.getWeights().length ; j++)
            handler.addNumItem(instance.getWeights()[j].getVariance());
         handler.endArray();
         handler.endElement();

         handler.startElement("C");
         handler.addNumItem(instance.getCapacity());
         handler.endElement();

         handler.startElement("c");
         handler.addNumItem(instance.getShortageCost());
         handler.endElement();

         handler.startElement("nbpartitions");
         handler.addIntItem(partitions);
         handler.endElement();
         
         double[] probabilities = PiecewiseStandardNormalFirstOrderLossFunction.getProbabilities(partitions);
         handler.startElement("prob");
         handler.startArray();
         for (int j = 0 ; j<probabilities.length ; j++)
            handler.addNumItem(probabilities[j]);
         handler.endArray();
         handler.endElement();         

         double[] means = PiecewiseStandardNormalFirstOrderLossFunction.getMeans(partitions);
         handler.startElement("means");
         handler.startArray();
         for (int j = 0 ; j<means.length ; j++)
            handler.addNumItem(means[j]);
         handler.endArray();
         handler.endElement();

         double error = (pwa == PWAPPROXIMATION.EDMUNDSON_MADANSKI ? PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions) : 0);
         handler.startElement("error");
         handler.addNumItem(error);
         handler.endElement();
         
         double sqrt_step = s;
         handler.startElement("s");
         handler.addNumItem(sqrt_step);
         handler.endElement();
         
         double sqrt_x0 = (pwa == PWAPPROXIMATION.EDMUNDSON_MADANSKI ? x0 : 0);
         handler.startElement("x0");
         handler.addNumItem(sqrt_x0);
         handler.endElement();
      }
   };
   
   public static void main(String args[]) {

      SKPNormal instance = SKPNormal.getTestInstance();
      double epsilon = 0.1;
      int simulationRuns = 100000;
      
      try {
         SKPNormalMILP sskp = new SKPNormalMILP(instance, PWAPPROXIMATION.JENSENS, epsilon);
         
         System.out.println(GSONUtility.<SKPNormalMILPSolvedInstance>printInstanceAsJSON(sskp.solve(simulationRuns)));
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
