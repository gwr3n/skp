package skp.milp;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

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

import skp.instance.SKP;
import skp.instance.SKPMultinormal;
import skp.milp.instance.SKPMultinormalScenarioBasedSolvedInstance;
import skp.sim.SimulateMultinormal;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

public class SKPMultinormalScenarioBased {
   
   static final long[] seed = {12345, 24513, 24531, 42531, 35124, 32451};
   RandomStream randGenerator;
   
   SKPMultinormal instance;
   int numberOfScenarios;
   public double[][] scenarios;
   
   final String model =  "sk_scenario_based";
   
   int[] optimalKnapsack;
   double milpSolutionValue;
   
   double milpOptimalityGap;
   double cplexSolutionTimeMs;
   int simplexIterations;
   int exploredNodes;
   
   public SKPMultinormalScenarioBased(SKPMultinormal instance, int numberOfScenarios, RandomStream randGenerator) {
      this.randGenerator = randGenerator;
      this.randGenerator.resetNextSubstream();
      this.instance = instance;
      this.numberOfScenarios = numberOfScenarios;
      SimulateMultinormal sm = new SimulateMultinormal(instance);
      scenarios = sm.sampleWeights(numberOfScenarios, this.randGenerator);
   }
   
   public SKPMultinormalScenarioBased(SKPMultinormal instance, int numberOfScenarios) {
      MRG32k3aL rnd = new MRG32k3aL();
      rnd.setSeed(seed);
      this.randGenerator = rnd;
      this.randGenerator.resetStartStream();
      
      this.instance = instance;
      this.numberOfScenarios = numberOfScenarios;
      SimulateMultinormal sm = new SimulateMultinormal(instance);
      scenarios = sm.sampleWeights(numberOfScenarios, this.randGenerator);
   }
   
   InputStream getMILPModelStream(File file){
      FileInputStream is = null;
      try{
         is = new FileInputStream(file);
      }catch(IOException e){
         e.printStackTrace();
      }
      return is;
   }
   
   IloOplDataSource getDataSource(IloOplFactory oplF) {
      return new MyData(oplF);
   }
   
   public SKPMultinormalScenarioBasedSolvedInstance solve(int simulationRuns) throws IloException {
      this.solveMILP(model, instance);
      
      SimulateMultinormal sim = new SimulateMultinormal(instance);
      double[] simulatedSolutionValue1 = sim.simulateMeanVariance(optimalKnapsack, simulationRuns, this.randGenerator);
      double[] simulatedSolutionValue2 = sim.simulateMeanVariance(optimalKnapsack, simulationRuns, this.randGenerator);
      
      SKPMultinormalScenarioBasedSolvedInstance solvedInstance = new SKPMultinormalScenarioBasedSolvedInstance(
            instance,
            optimalKnapsack,
            simulatedSolutionValue1[0],
            simulatedSolutionValue2[0],
            simulatedSolutionValue1[1],
            simulatedSolutionValue2[1],
            simulationRuns,
            milpSolutionValue,
            milpOptimalityGap,
            cplexSolutionTimeMs,
            simplexIterations,
            exploredNodes);
      
      return solvedInstance;
   }
   
   void solveMILP(String model_name, SKP instance) throws IloException{
      double startGlobal = System.currentTimeMillis();
      IloOplFactory.setDebugMode(false);
      IloOplFactory oplF = new IloOplFactory();
      IloOplErrorHandler errHandler = oplF.createOplErrorHandler(System.out);
      IloCplex cplex = oplF.createCplex();
      ClassLoader classLoader = getClass().getClassLoader();
      InputStream isModel = getMILPModelStream(new File(classLoader.getResource("./opl_models/"+model_name+".mod").getFile()));
      IloOplModelSource modelSource=oplF.createOplModelSourceFromStream(isModel,model_name);
      IloOplSettings settings = oplF.createOplSettings(errHandler);
      IloOplModelDefinition def=oplF.createOplModelDefinition(modelSource,settings);
      IloOplModel opl=oplF.createOplModel(def,cplex);
      cplex.setParam(IloCplex.Param.Threads, 8);
      cplex.setParam(IloCplex.Param.MIP.Display, 2);
      cplex.setParam(IloCplex.Param.TimeLimit, 60*10);
      //cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-4); //this is the default value
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

      IloOplDataSource dataSource = this.getDataSource(oplF);
      opl.addDataSource(dataSource);
      opl.generate();

      cplex.setOut(null);

      //double start = cplex.getCplexImpl().getCplexTime();
      boolean status =  cplex.solve();
      //double end = cplex.getCplexImpl().getCplexTime();
      if ( status ) {   
         double endGlobal = System.currentTimeMillis();
         this.milpSolutionValue = cplex.getObjValue();
         this.milpOptimalityGap = cplex.getMIPRelativeGap();
         this.cplexSolutionTimeMs = endGlobal - startGlobal; // (end - start)*1000;
         this.simplexIterations = cplex.getNiterations();
         this.exploredNodes = cplex.getNnodes();
         
         this.optimalKnapsack = new int[instance.getItems()];
         for(int i = 0; i < instance.getItems(); i++){
            this.optimalKnapsack[i] = (int) Math.round(cplex.getValue(opl.getElement("X").asIntVarMap().get(i+1)));
         }
      } else {
         System.out.println("No solution!");
      } 
      // Release all resources
      opl.end();
      def.end();
      settings.end();
      modelSource.end();
      try {
         isModel.close();
      } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      cplex.end();
      errHandler.end();
      oplF.end();
      System.gc();
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
         
         handler.startElement("S");
         handler.addIntItem(numberOfScenarios);
         handler.endElement();

         handler.startElement("expectedValues");
         handler.startArray();
         for (int j = 0 ; j<instance.getExpectedValues().length ; j++)
            handler.addNumItem(instance.getExpectedValues()[j]);
         handler.endArray();
         handler.endElement();
         
         handler.startElement("weights");
         handler.startArray();
         for (int i = 0 ; i < numberOfScenarios ; i++) {
            handler.startArray();
            for (int j = 0 ; j < instance.getItems(); j++) {
               handler.addNumItem(scenarios[i][j]);
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
      }
   };
   
   public static void main(String args[]) {

      SKPMultinormal instance = SKPMultinormal.getTestInstance();
      
      try {
         int numberOfScenarios = 1000;
         SKPMultinormalScenarioBased sskp = new SKPMultinormalScenarioBased(instance, numberOfScenarios);
         
         int simulationRuns = 10000;
         System.out.println(GSONUtility.<SKPMultinormalScenarioBasedSolvedInstance>printInstanceAsJSON(sskp.solve(simulationRuns)));
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
