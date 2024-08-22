package skp.milp;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import ilog.concert.IloException;
import ilog.opl.IloCplex;
import ilog.opl.IloOplDataSource;
import ilog.opl.IloOplErrorHandler;
import ilog.opl.IloOplFactory;
import ilog.opl.IloOplModel;
import ilog.opl.IloOplModelDefinition;
import ilog.opl.IloOplModelSource;
import ilog.opl.IloOplSettings;

import skp.instance.SKP;

public abstract class SKPMILP {
   
   int partitions;
   int linearizationSamples;
   
   int[] optimalKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   String model;
   
   double cplexSolutionTimeMs;
   int simplexIterations;
   int exploredNodes;

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
   
   InputStream getMILPModelStream(File file){
      FileInputStream is = null;
      try{
         is = new FileInputStream(file);
      }catch(IOException e){
         e.printStackTrace();
      }
      return is;
   }
   
   abstract IloOplDataSource getDataSource(IloOplFactory oplF);
   
   abstract void computeMILPMaxLinearizationError(IloOplModel opl, IloCplex cplex) throws IloException;
   
   void solveMILP(String model_name, SKP instance) throws IloException{
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

      double start = cplex.getCplexImpl().getCplexTime();
      boolean status =  cplex.solve();
      double end = cplex.getCplexImpl().getCplexTime();
      if ( status ) {   
         this.milpSolutionValue = cplex.getObjValue();
         this.milpOptimalityGap = cplex.getMIPRelativeGap();
         this.cplexSolutionTimeMs = (end - start)*1000;
         this.simplexIterations = cplex.getNiterations();
         this.exploredNodes = cplex.getNnodes();
         
         this.optimalKnapsack = new int[instance.getItems()];
         for(int i = 0; i < instance.getItems(); i++){
            this.optimalKnapsack[i] = (int) Math.round(cplex.getValue(opl.getElement("X").asIntVarMap().get(i+1)));
         }
         
         computeMILPMaxLinearizationError(opl, cplex);
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
}
