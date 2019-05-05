package skp.milp;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import ilog.concert.IloException;
import ilog.opl.IloCplex;
import ilog.opl.IloCustomOplDataSource;
import ilog.opl.IloOplDataHandler;
import ilog.opl.IloOplDataSource;
import ilog.opl.IloOplErrorHandler;
import ilog.opl.IloOplFactory;
import ilog.opl.IloOplModel;
import ilog.opl.IloOplModelDefinition;
import ilog.opl.IloOplModelSource;
import ilog.opl.IloOplSettings;

import skp.folf.PiecewiseFirstOrderLossFunction;
import skp.instance.SKPPoisson;
import skp.sim.SimulatePoisson;

public class SKPPoissonMILP {
   SKPPoisson instance;
   
   int partitions;
   double[] probabilityMasses;
   int linearizationSamples;
   
   int[] optimalKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   private final String model = "sk_poisson";
   double milpMaxLinearizationError;
   
   double cplexSolutionTimeMs;
   int simplexIterations;
   int exploredNodes;
   
   public SKPPoissonMILP(SKPPoisson instance, int partitions, int linearizationSamples) {
      this.instance = instance;
      this.partitions = partitions;
      this.linearizationSamples = linearizationSamples;
      this.probabilityMasses = new double[partitions];
      Arrays.fill(this.probabilityMasses, 1.0/partitions);
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
   
   private InputStream getMILPModelStream(File file){
      FileInputStream is = null;
      try{
         is = new FileInputStream(file);
      }catch(IOException e){
         e.printStackTrace();
      }
      return is;
   }
   
   public SKPPoissonMILPSolvedInstance solve(int simulationRuns) throws IloException {
      this.solveMILP(model);
      
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
   
   private void solveMILP(String model_name) throws IloException{
      IloOplFactory.setDebugMode(false);
      IloOplFactory oplF = new IloOplFactory();
      IloOplErrorHandler errHandler = oplF.createOplErrorHandler(System.out);
      IloCplex cplex = oplF.createCplex();
      ClassLoader classLoader = getClass().getClassLoader();
      IloOplModelSource modelSource=oplF.createOplModelSourceFromStream(getMILPModelStream(new File(classLoader.getResource("./opl_models/"+model_name+".mod").getFile())),model_name);
      IloOplSettings settings = oplF.createOplSettings(errHandler);
      IloOplModelDefinition def=oplF.createOplModelDefinition(modelSource,settings);
      IloOplModel opl=oplF.createOplModel(def,cplex);
      cplex.setParam(IloCplex.IntParam.Threads, 8);
      cplex.setParam(IloCplex.IntParam.MIPDisplay, 2);
      /*cplex.setParam(IloCplex.IntParam.VarSel, 1);
      cplex.setParam(IloCplex.IntParam.ZeroHalfCuts, 2);
      cplex.setParam(IloCplex.IntParam.ImplBd, 2);
      cplex.setParam(IloCplex.IntParam.FracCuts, 2);
      cplex.setParam(IloCplex.IntParam.GUBCovers, 2);
      cplex.setParam(IloCplex.IntParam.DisjCuts, 2);
      cplex.setParam(IloCplex.IntParam.Covers, 2);
      cplex.setParam(IloCplex.IntParam.Cliques, 2);
      cplex.setParam(IloCplex.IntParam.FlowCovers, 2);
      cplex.setParam(IloCplex.IntParam.FlowPaths, 2);
      cplex.setParam(IloCplex.IntParam.MIRCuts, 2);
      cplex.setParam(IloCplex.IntParam.MIPEmphasis, 3);
       */

      IloOplDataSource dataSource = new SKPPoissonMILP.MyData(oplF);
      opl.addDataSource(dataSource);
      opl.generate();

      cplex.setOut(null);

      double start = cplex.getCplexImpl().getCplexTime();
      boolean status =  cplex.solve();
      double end = cplex.getCplexImpl().getCplexTime();
      if ( status ) {   
         this.milpSolutionValue = cplex.getObjValue();
         this.milpOptimalityGap = cplex.getMIPRelativeGap();
         cplexSolutionTimeMs = (end - start)*1000;
         simplexIterations = cplex.getNiterations();
         exploredNodes = cplex.getNnodes();
         
         this.optimalKnapsack = new int[instance.getItems()];
         for(int i = 0; i < instance.getItems(); i++){
            this.optimalKnapsack[i] = (int) Math.round(cplex.getValue(opl.getElement("X").asIntVarMap().get(i+1)));
         }
         
         double[] errors = PiecewiseFirstOrderLossFunction.poissonKnapsackPiecewiseFOLFApproximationErrors(instance.getCapacity(), probabilityMasses, linearizationSamples);
         this.milpMaxLinearizationError = this.instance.getShortageCost()*
                                      errors[(int)Math.round(cplex.getValue(opl.getElement("M").asNumVar()))];
      } else {
         System.out.println("No solution!");
         opl.end();
         oplF.end();//
         errHandler.end();
         cplex.end();
         System.gc(); 
      } 
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
}
