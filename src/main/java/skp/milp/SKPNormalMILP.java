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
import skp.SKPNormal;
import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;

public class SKPNormalMILP {

   SKPNormal instance;
   int partitions;
   int[] knapsack;
   double solutionValue;
   private final String model = "sk_normal";
   double maxLinearizationError;

   public SKPNormalMILP(SKPNormal instance, int partitions) throws IloException{
      this.instance = instance;
      this.partitions = partitions;
      this.solve(model);
   }
   
   public int[] getKnapsack() {
      return this.knapsack;
   }
   
   public double getSolutionValue() {
      return this.solutionValue;
   }
   
   public double getMaxLinearizationError() {
      return this.maxLinearizationError;
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

   public void solve(String model_name) throws IloException{
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

      IloOplDataSource dataSource = new SKPNormalMILP.MyData(oplF);
      opl.addDataSource(dataSource);
      opl.generate();

      cplex.setOut(null);

      double start = cplex.getCplexImpl().getCplexTime();
      boolean status =  cplex.solve();
      double end = cplex.getCplexImpl().getCplexTime();
      if ( status ) {   
         this.solutionValue = cplex.getObjValue();
         @SuppressWarnings("unused")
         double time = end - start;
         
         this.knapsack = new int[instance.getItems()];
         for(int i = 0; i < instance.getItems(); i++){
            this.knapsack[i] = (int) Math.round(cplex.getValue(opl.getElement("X").asIntVarMap().get(i+1)));
         }
         
         this.maxLinearizationError = this.instance.getShortageCost()*
                                      cplex.getValue(opl.getElement("S").asNumVar())*
                                      PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions);
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
}
