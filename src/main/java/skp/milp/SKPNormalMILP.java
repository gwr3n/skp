package skp.milp;

import ilog.concert.IloException;
import ilog.opl.IloCplex;
import ilog.opl.IloCustomOplDataSource;
import ilog.opl.IloOplDataHandler;
import ilog.opl.IloOplDataSource;
import ilog.opl.IloOplFactory;
import ilog.opl.IloOplModel;
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
   
   public SKPNormalMILP(SKPNormal instance, int partitions) throws IloException{
      this.instance = instance;
      this.partitions = partitions;
      linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model =  "sk_normal";
   }
   
   IloOplDataSource getDataSource(IloOplFactory oplF) {
      return new SKPNormalMILP.MyData(oplF);
   }
   
   void computeMILPMaxLinearizationError(IloOplModel opl, IloCplex cplex) throws IloException{
      this.milpMaxLinearizationError = this.instance.getShortageCost()*
                                       cplex.getValue(opl.getElement("S").asNumVar())*
                                       PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions);
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
         handler.addIntItem(instance.getExpectedValuesPerUnit().length);
         handler.endElement();

         handler.startElement("expectedValues");
         handler.startArray();
         for (int j = 0 ; j<instance.getExpectedValuesPerUnit().length ; j++)
            handler.addNumItem(instance.getExpectedValuesPerUnit()[j]*instance.getWeights()[j].getMean());
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

         double error = PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions);
         handler.startElement("error");
         handler.addNumItem(error);
         handler.endElement();
      }
   };
   
   public static void main(String args[]) {

      SKPNormal instance = SKPNormal.getTestInstance();
      
      try {
         SKPNormalMILP sskp = new SKPNormalMILP(instance, 10);
         
         System.out.println(GSONUtility.<SKPNormalMILPSolvedInstance>printInstanceAsJSON(sskp.solve(100000)));
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
