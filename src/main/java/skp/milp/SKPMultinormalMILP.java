package skp.milp;

import ilog.concert.IloException;
import ilog.opl.IloCplex;
import ilog.opl.IloCustomOplDataSource;
import ilog.opl.IloOplDataHandler;
import ilog.opl.IloOplDataSource;
import ilog.opl.IloOplFactory;
import ilog.opl.IloOplModel;
import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.SKPMultinormal;
import skp.milp.instance.SKPMultinormalMILPSolvedInstance;
import skp.sim.SimulateMultinormal;

public class SKPMultinormalMILP extends SKPMILP{
   SKPMultinormal instance;

   public SKPMultinormalMILP(SKPMultinormal instance, int partitions)  throws IloException{
      this.instance = instance;
      this.partitions = partitions;
      this.linearizationSamples = PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples();
      this.model = "sk_mvnormal";
   }
   
   IloOplDataSource getDataSource(IloOplFactory oplF) {
      return new SKPMultinormalMILP.MyData(oplF);
   }
   
   void computeMILPMaxLinearizationError(IloOplModel opl, IloCplex cplex) throws IloException{
      this.milpMaxLinearizationError = this.instance.getShortageCost()*
                                       cplex.getValue(opl.getElement("S").asNumVar())*
                                       PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions);
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
         handler.addIntItem(instance.getExpectedValuesPerUnit().length);
         handler.endElement();

         handler.startElement("expectedValues");
         handler.startArray();
         for (int j = 0 ; j<instance.getExpectedValuesPerUnit().length ; j++)
            handler.addNumItem(instance.getExpectedValuesPerUnit()[j]*instance.getWeights().getMean()[j]);
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

         double error = PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions);
         handler.startElement("error");
         handler.addNumItem(error);
         handler.endElement();
      }
   };
}
