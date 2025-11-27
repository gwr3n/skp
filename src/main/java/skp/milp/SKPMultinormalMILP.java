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
import skp.instance.SKPMultinormal;
import skp.milp.instance.SKPMultinormalMILPSolvedInstance;
import skp.sim.SimulateMultinormal;
import skp.utilities.gson.GSONUtility;

public class SKPMultinormalMILP extends SKPMILP{
   SKPMultinormal instance;
   boolean ignoreCorrelation;
   
   PWAPPROXIMATION pwa;        // sqrt approximation       
   double s = 1e-2;            // sqrt approximation step    
   double x0 = Math.sqrt(s)/4; // sqrt approximation x0  

   /*public SKPMultinormalMILP(SKPMultinormal instance, int partitions, PWAPPROXIMATION pwa)  throws IloException{
      this.instance = instance;
      this.partitions = partitions;
      this.linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model = "sk_mvnormal_linear";
      this.pwa = pwa;
   }*/
   
   public SKPMultinormalMILP(SKPMultinormal instance, int partitions, double s, PWAPPROXIMATION pwa)  throws IloException{
      this.instance = instance;
      this.partitions = partitions;
      this.linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model = "sk_mvnormal_linear";
      this.pwa = pwa;
      this.s = s;       
      this.x0 = Math.sqrt(s)/4;
   }
   
   public SKPMultinormalMILP(SKPMultinormal instance, PWAPPROXIMATION pwa, double epsilon)  throws IloException{
      this.instance = instance;
      this.linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model = "sk_mvnormal_linear";
      this.pwa = pwa;
      
      int[] param = chooseLinearisationParameters(epsilon);
      
      this.partitions = param[0] - 1;
      double[][] covMatrix = instance.getWeights().getCovariance();
      double Vmax = 0.0;
      for (int i = 0; i < covMatrix.length; i++) {
          for (int j = 0; j < covMatrix[i].length; j++) {
              Vmax += covMatrix[i][j];
          }
      }
      this.s = Vmax/param[1];       
      this.x0 = Math.sqrt(this.s)/4; 
      
      System.out.println("Linearisation parameters: W = " + this.partitions + "; s = " + this.s);
   }
   
   IloOplDataSource getDataSource(IloOplFactory oplF) {
      return new SKPMultinormalMILP.MyData(oplF);
   }
   
   void computeMILPMaxLinearizationError(IloOplModel opl, IloCplex cplex) throws IloException{
      this.milpMaxLinearizationError = this.instance.getShortageCost()*
                                       cplex.getValue(opl.getElement("S").asNumVar())*
                                       PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions);
   }
   
   public int[] chooseLinearisationParameters(final double epsilon)
   {
       /* ---------- constants that depend only on the instance ---------- */
       double[][] covMatrix = instance.getWeights().getCovariance();
       double Vmax = 0.0;
       for (int i = 0; i < covMatrix.length; i++) {
           for (int j = 0; j < covMatrix[i].length; j++) {
               Vmax += covMatrix[i][j];
           }
       }
       double c    = instance.getShortageCost();
       
       /* ================================================================
        * Delegate to helper method
        * ================================================================ */
       return LinearisationFactory.chooseLinearisationParameters(epsilon, Vmax, c);
   }
   
   public static SKPMultinormalMILPSolvedInstance solve(SKPMultinormal instance, int partitions, double s, int simulationRuns) throws IloException {
      SKPMultinormalMILPSolvedInstance solvedJensens = new SKPMultinormalMILP(instance, partitions, s, PWAPPROXIMATION.JENSENS).solve(simulationRuns);
      SKPMultinormalMILPSolvedInstance solvedEM = new SKPMultinormalMILP(instance, partitions, s, PWAPPROXIMATION.EDMUNDSON_MADANSKI).solve(simulationRuns);
      
      double optimality_gap = 0.0;
      SKPMultinormalMILPSolvedInstance opt = solvedJensens;
      if(!Arrays.equals(solvedJensens.optimalKnapsack, solvedEM.optimalKnapsack)) {
         optimality_gap = (solvedJensens.milpSolutionValue - solvedEM.milpSolutionValue)/(1e-10 + solvedEM.milpSolutionValue);
         if(solvedJensens.simulatedSolutionValue < solvedEM.simulatedSolutionValue) {
            opt = solvedEM;
         } 
      }
      
      return new SKPMultinormalMILPSolvedInstance(
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
   
   public SKPMultinormalMILPSolvedInstance solve(int simulationRuns) throws IloException {
      this.solve();
      
      SimulateMultinormal sim = new SimulateMultinormal(instance);
      double simulatedSolutionValue = sim.simulate(optimalKnapsack, simulationRuns);
      
      double milpMaxLinearizationError = 100*this.getMILPMaxLinearizationError()/simulatedSolutionValue;
      double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;
      
      SKPMultinormalMILPSolvedInstance solvedInstance = new SKPMultinormalMILPSolvedInstance(instance,
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
         for (int j = 0 ; j<instance.getWeights().getMean().length ; j++)
            handler.addNumItem(instance.getWeights().getMean()[j]);
         handler.endArray();
         handler.endElement();

         handler.startElement("varianceCovarianceWeights");
         handler.startArray();
         for (int i = 0 ; i<instance.getItems() ; i++) {
            handler.startArray();
            for (int j = 0 ; j<instance.getItems() ; j++) {
               handler.addNumItem(instance.getWeights().getCovariance()[i][j]);
            }
            handler.endArray();
         }
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
      
      SKPMultinormal instance = SKPMultinormal.getTestInstance();
      double epsilon = 0.1;
      int simulationRuns = 1000000;
      
      try {
         SKPMultinormalMILP sskp_J = new SKPMultinormalMILP(instance, PWAPPROXIMATION.JENSENS, epsilon);
         
         SKPMultinormalMILPSolvedInstance solved_J = sskp_J.solve(simulationRuns);
         System.out.println(GSONUtility.<SKPMultinormalMILPSolvedInstance>printInstanceAsJSON(solved_J));
         
         boolean test = true;
         if(test) {
            SKPMultinormalMILP sskp_EM = new SKPMultinormalMILP(instance, PWAPPROXIMATION.EDMUNDSON_MADANSKI, epsilon);
            
            SKPMultinormalMILPSolvedInstance solved_EM = sskp_EM.solve(simulationRuns);
            System.out.println(GSONUtility.<SKPMultinormalMILPSolvedInstance>printInstanceAsJSON(solved_EM));
            
            System.out.println("Jensen's MILP solution value: " + solved_J.milpSolutionValue);
            System.out.println("Edmundson-Madanski MILP solution value: " + solved_EM.milpSolutionValue);
            System.out.println("Difference: " + (solved_J.milpSolutionValue - solved_EM.milpSolutionValue));
            System.out.println("Allowed epsilon: " + epsilon);
            if(solved_J.milpSolutionValue - solved_EM.milpSolutionValue <= epsilon + 1e-12)
               System.out.println("Difference is within allowed epsilon.");
            else
               System.out.println("Difference exceeds allowed epsilon!");
         }
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
   }
}
