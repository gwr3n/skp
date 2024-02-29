package skp.milp;

import java.util.Arrays;

import ilog.concert.IloException;
import ilog.opl.IloCplex;
import ilog.opl.IloCustomOplDataSource;
import ilog.opl.IloOplDataHandler;
import ilog.opl.IloOplDataSource;
import ilog.opl.IloOplFactory;
import ilog.opl.IloOplModel;

import skp.folf.PiecewiseFirstOrderLossFunction;
import skp.instance.SKPPoisson;
import skp.milp.instance.SKPPoissonMILPSolvedInstance;
import skp.sim.SimulatePoisson;
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

public class SKPPoissonMILP extends SKPMILP{
   SKPPoisson instance;
   
   double[] probabilityMasses;
   
   public SKPPoissonMILP(SKPPoisson instance, int partitions, int linearizationSamples) {
      this.instance = instance;
      this.partitions = partitions;
      this.linearizationSamples = linearizationSamples;
      this.probabilityMasses = new double[partitions];
      Arrays.fill(this.probabilityMasses, 1.0/partitions);
      this.model = "sk_poisson";
   }
   
   IloOplDataSource getDataSource(IloOplFactory oplF) {
      return new SKPPoissonMILP.MyData(oplF);
   }
   
   void computeMILPMaxLinearizationError(IloOplModel opl, IloCplex cplex) throws IloException {
      double[] errors = PiecewiseFirstOrderLossFunction.poissonKnapsackPiecewiseFOLFApproximationErrors(instance.getCapacity(), probabilityMasses, linearizationSamples);
      this.milpMaxLinearizationError = instance.getShortageCost()*
            errors[(int)Math.round(cplex.getValue(opl.getElement("M").asNumVar()))];
   }
   
   public void solve() throws IloException {
      this.solveMILP(model, instance);
   }
   
   public SKPPoissonMILPSolvedInstance solve(int simulationRuns) throws IloException {
      this.solveMILP(model, instance);
      
      SimulatePoisson sim = new SimulatePoisson(instance);
      double simulatedSolutionValue = sim.simulate(optimalKnapsack, simulationRuns);
      
      double milpMaxLinearizationError = 100*this.getMILPMaxLinearizationError()/simulatedSolutionValue;
      double simulatedLinearizationError = 100*(simulatedSolutionValue-milpSolutionValue)/simulatedSolutionValue;
      
      SKPPoissonMILPSolvedInstance solvedInstance = new SKPPoissonMILPSolvedInstance(
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

         handler.startElement("C");
         handler.addIntItem(instance.getCapacity());
         handler.endElement();

         handler.startElement("c");
         handler.addNumItem(instance.getShortageCost());
         handler.endElement();

         handler.startElement("nbpartitions");
         handler.addIntItem(partitions);
         handler.endElement();
         
         int maxWeight = PiecewiseFirstOrderLossFunction.poissonKnapsackPiecewiseFOLFMaxWeight(instance.getCapacity(), probabilityMasses, linearizationSamples);
         handler.startElement("maxWeight");
         handler.addIntItem(maxWeight);
         handler.endElement();
         
         handler.startElement("prob");
         handler.startArray();
         for (int j = 0 ; j<probabilityMasses.length ; j++)
            handler.addNumItem(probabilityMasses[j]);
         handler.endArray();
         handler.endElement();

         double[][] means = PiecewiseFirstOrderLossFunction.poissonKnapsackPiecewiseFOLFConditionalExpectations(instance.getCapacity(), probabilityMasses, linearizationSamples);
         handler.startElement("means");
         handler.startArray();
         for (int i = 0 ; i < maxWeight + 1; i++){
            handler.startArray();
            for (int j = 0 ; j < partitions ; j++){
               handler.addNumItem(means[i][j]);
            }
            handler.endArray();
         }
         handler.endArray();
         handler.endElement();

         double[] errors = PiecewiseFirstOrderLossFunction.poissonKnapsackPiecewiseFOLFApproximationErrors(instance.getCapacity(), probabilityMasses, linearizationSamples);
         handler.startElement("error");
         handler.startArray();
         for (int j = 0 ; j<errors.length ; j++)
            handler.addNumItem(errors[j]);
         handler.endArray();
         handler.endElement();
      }
   };
   
   public static void main(String args[]) {
      SKPPoisson instance = SKPPoisson.getTestInstance();
      
      try {
         SKPPoissonMILP sskp = new SKPPoissonMILP(instance, 10, 100000);
         
         System.out.println(GSONUtility.<SKPPoissonMILPSolvedInstance>printInstanceAsJSON(sskp.solve(100000)));
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
