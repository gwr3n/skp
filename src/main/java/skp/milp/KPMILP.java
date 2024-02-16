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
import skp.instance.KP;

public class KPMILP {
   KP instance;
   int[] knapsack;
   double solutionValue;
   private final String model;
   
   boolean selector;
   int selectFirstObject;
   
   public KPMILP(KP instance) throws IloException{
      this.instance = instance;
      this.selector = false;
      this.model = "kp";
      this.solve(model);
   }
   
   public KPMILP(KP instance, int selectFirstObject) throws IloException{
      this.instance = instance;
      this.selector = true;
      this.selectFirstObject = selectFirstObject;
      this.model = "kp_selector_x1";
      this.solve(model);
   }
   
   public int[] getKnapsack() {
      return this.knapsack;
   }
   
   public double getSolutionValue() {
      return this.solutionValue;
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
      InputStream isModel = getMILPModelStream(new File(classLoader.getResource("./opl_models/"+model_name+".mod").getFile()));
      IloOplModelSource modelSource=oplF.createOplModelSourceFromStream(isModel,model_name);
      IloOplSettings settings = oplF.createOplSettings(errHandler);
      IloOplModelDefinition def=oplF.createOplModelDefinition(modelSource,settings);
      IloOplModel opl=oplF.createOplModel(def,cplex);
      cplex.setParam(IloCplex.Param.Threads, 8);
      cplex.setParam(IloCplex.Param.MIP.Display, 2);
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

      IloOplDataSource dataSource = selector ? new KPMILP.MyDataSelector(oplF) : new KPMILP.MyData(oplF);
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
      } else {
         System.out.println("No solution!");
      } 
      try {
         isModel.close();
      } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      opl.end();
      errHandler.end();
      cplex.end();
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

         handler.startElement("values");
         handler.startArray();
         for (int j = 0 ; j<instance.getValues().length ; j++)
            handler.addNumItem(instance.getValues()[j]*instance.getWeights()[j]);
         handler.endArray();
         handler.endElement();

         handler.startElement("weights");
         handler.startArray();
         for (int j = 0 ; j<instance.getWeights().length ; j++)
            handler.addNumItem(instance.getWeights()[j]);
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
   
   class MyDataSelector extends IloCustomOplDataSource
   {
      MyDataSelector(IloOplFactory oplF){
         super(oplF);
      }

      public void customRead(){
         IloOplDataHandler handler = getDataHandler();

         handler.startElement("N");
         handler.addIntItem(instance.getItems());
         handler.endElement();

         handler.startElement("values");
         handler.startArray();
         for (int j = 0 ; j<instance.getValues().length ; j++)
            handler.addNumItem(instance.getValues()[j]*instance.getWeights()[j]);
         handler.endArray();
         handler.endElement();

         handler.startElement("weights");
         handler.startArray();
         for (int j = 0 ; j<instance.getWeights().length ; j++)
            handler.addNumItem(instance.getWeights()[j]);
         handler.endArray();
         handler.endElement();

         handler.startElement("C");
         handler.addNumItem(instance.getCapacity());
         handler.endElement();

         handler.startElement("c");
         handler.addNumItem(instance.getShortageCost());
         handler.endElement();
         
         handler.startElement("selectFirstObject");
         handler.addIntItem(selectFirstObject);
         handler.endElement();
      }
   };
}  
