package skp.folf;

import umontreal.ssj.probdist.EmpiricalDist;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.probdist.Distribution;

public class PiecewiseFirstOrderLossFunction extends FirstOrderLossFunction{
   
   public PiecewiseFirstOrderLossFunction(Distribution[] distributions){
      super(distributions);
   }
   
   /**
    * Piecewise linearization parameters
    */
   
   public double[] getConditionalExpectations(double[] probabilityMasses, int nbSamples){
      double[] conditionalExpectations = new double[probabilityMasses.length];
      EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples);
      double probabilityMass = 0;
      int conditionalExpectationIndex = 0;
      for(int i = 0; i < empDistribution.getN(); i++){
         if(probabilityMass < 1 && probabilityMass < probabilityMasses[conditionalExpectationIndex]){
            conditionalExpectations[conditionalExpectationIndex] += empDistribution.getObs(i)/empDistribution.getN();
            probabilityMass += 1.0/empDistribution.getN();
         }else{
            conditionalExpectations[conditionalExpectationIndex] /= probabilityMasses[conditionalExpectationIndex];
            probabilityMass = 0;
            conditionalExpectationIndex++;
         }
      }
      conditionalExpectations[conditionalExpectationIndex] /= probabilityMasses[conditionalExpectationIndex];
      return conditionalExpectations;
   }
   
   private double[] getApproximationErrors(double[] probabilityMasses, int nbSamples){
      double[] conditionalExpectations = this.getConditionalExpectations(probabilityMasses, nbSamples);
      double[] approximationErrors = new double[conditionalExpectations.length];

      for(int i = 0; i < probabilityMasses.length; i++){
         approximationErrors[i] = getComplementaryFirstOrderLossFunctionValue(conditionalExpectations[i], nbSamples)-
               getPiecewiseComplementaryFirstOrderLossFunctionValue(i, conditionalExpectations[i], probabilityMasses, conditionalExpectations);
      }

      return approximationErrors;
   }

   double getMaxApproximationError(double[] probabilityMasses, int nbSamples){
      double[] approximationErrors = this.getApproximationErrors(probabilityMasses, nbSamples);
      double maxApproximationError = 0;

      for(int i = 0; i < probabilityMasses.length; i++){
         maxApproximationError = Math.max(maxApproximationError, approximationErrors[i]);
      }

      return maxApproximationError;
   }
   
   private double getPiecewiseComplementaryFirstOrderLossFunctionValue(int segmentIndex, double x, double[] probabilityMasses, double[] conditionalExpectations){
      double value = 0;
      for(int i = 0; i < segmentIndex; i++){
         value += (x-conditionalExpectations[i])*probabilityMasses[i];
      }
      return value;
   }
   
   private double getPiecewiseFirstOrderLossFunctionValue(int segmentIndex, double x, double[] probabilityMasses, double[] conditionalExpectations){
      double value = 0;
      for(int i = segmentIndex; i < probabilityMasses.length; i++){
         value += (conditionalExpectations[i]-x)*probabilityMasses[i];
      }
      return value;
   }
   
   private double getPiecewiseComplementaryFirstOrderLossFunctionErrorValue(double x, int nbSamples, double[] probabilityMasses, double[] conditionalExpectations){
      double lossFunctionValue = this.getComplementaryFirstOrderLossFunctionValue(x, nbSamples);

      double maxValue = 0;
      for(int j = 0; j <= probabilityMasses.length; j++){
         double value = 0;
         for(int i = 0; i < j; i++){
            value += (x-conditionalExpectations[i])*probabilityMasses[i];
         }
         maxValue = Math.max(maxValue, value);
      }
      
      return lossFunctionValue-maxValue;
   }
   
   private double getPiecewiseFirstOrderLossFunctionErrorValue(double x, int nbSamples, double[] probabilityMasses, double[] conditionalExpectations){
      double lossFunctionValue = this.getFirstOrderLossFunctionValue(x, nbSamples);

      double maxValue = 0;
      for(int j = 0; j < probabilityMasses.length; j++){
         double value = 0;
         for(int i = j; i < probabilityMasses.length; i++){
            value += (conditionalExpectations[i]-x)*probabilityMasses[i];
         }
         maxValue = Math.max(maxValue, value);
      }
      
      return lossFunctionValue-maxValue;
   }
   
   /**
    * Charting
    */
   
   private XYSeries getPiecewiseComplementaryFirstOrderLossFunctionXYSeriesForSegment(int segmentIndex, 
         double[] probabilityMasses, 
         double[] conditionalExpectations, 
         double min, 
         double max, 
         double minYValue,
         double precision){
      XYSeries series = new XYSeries("Piecewise Complementary First Order Loss Function");
      for(double x = min; x <= max; x+= precision){
         double value = getPiecewiseComplementaryFirstOrderLossFunctionValue(segmentIndex, x, probabilityMasses, conditionalExpectations);
         if(value >= minYValue) series.add(x, value);
      }
      return series;
   }
   
   private XYSeries getPiecewiseFirstOrderLossFunctionXYSeriesForSegment(int segmentIndex, 
         double[] probabilityMasses, 
         double[] conditionalExpectations, 
         double min, 
         double max, 
         double minYValue,
         double precision){
      XYSeries series = new XYSeries("Piecewise First Order Loss Function");
      for(double x = min; x <= max; x+= precision){
         double value = getPiecewiseFirstOrderLossFunctionValue(segmentIndex, x, probabilityMasses, conditionalExpectations);
         if(value >= minYValue) series.add(x, value);
      }
      return series;
   }
   
   private XYSeries getPiecewiseComplementaryFirstOrderLossFunctionErrorXYSeries(double min, double max, int nbSamples, double[] probabilityMasses, double[] conditionalExpectations, double precision){
      XYSeries series = new XYSeries("Piecewise Complementary First Order Loss Function Error");
      for(double x = min; x <= max; x+= precision){
         series.add(x, getPiecewiseComplementaryFirstOrderLossFunctionErrorValue(x, nbSamples, probabilityMasses, conditionalExpectations));
      }
      return series;
   }
   
   private XYSeries getPiecewiseFirstOrderLossFunctionErrorXYSeries(double min, double max, int nbSamples, double[] probabilityMasses, double[] conditionalExpectations, double precision){
      XYSeries series = new XYSeries("Piecewise First Order Loss Function Error");
      for(double x = min; x <= max; x+= precision){
         series.add(x, getPiecewiseFirstOrderLossFunctionErrorValue(x, nbSamples, probabilityMasses, conditionalExpectations));
      }
      return series;
   }
   
   private void plotPiecewiseComplementaryFirstOrderLossFunction(double min, double max, double minYValue, double[] probabilityMasses, int nbSamples, double precision, boolean saveToDisk){
      int segments = probabilityMasses.length + 1;
      double[] conditionalExpectations = this.getConditionalExpectations(probabilityMasses, nbSamples);

      XYSeriesCollection xyDataset = new XYSeriesCollection();

      xyDataset.addSeries(this.getComplementaryFirstOrderLossFunctionXYSeries(min, max, nbSamples, precision));

      for(int i = 0; i < segments; i++)
         xyDataset.addSeries(this.getPiecewiseComplementaryFirstOrderLossFunctionXYSeriesForSegment(i, probabilityMasses, conditionalExpectations, min, max, minYValue, precision));

      xyDataset.addSeries(this.getPiecewiseComplementaryFirstOrderLossFunctionErrorXYSeries(min, max, nbSamples, probabilityMasses, conditionalExpectations, precision));
      
      JFreeChart chart = ChartFactory.createXYLineChart("Complementary First Order Loss Function", "x", "CL(x)",
            xyDataset, PlotOrientation.VERTICAL, false, true, false);
      ChartFrame frame = new ChartFrame("Empirical complementary loss function",chart);
      frame.setVisible(true);
      frame.setSize(500,400);

      if(saveToDisk){

         XYLineChart lc = new XYLineChart("Piecewise linearization", "x", "CL(x)", xyDataset);

         try {
            File latexFolder = new File("./latex");
            if(!latexFolder.exists()){
               latexFolder.mkdir();
            }
            Writer file = new FileWriter("./latex/pw_cfolf_graph.tex");
            file.write(lc.toLatex(8, 5));
            file.close();
         } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }
   
   private void plotPiecewiseFirstOrderLossFunction(double min, double max, double minYValue, double[] probabilityMasses, int nbSamples, double precision, boolean saveToDisk){
      int segments = probabilityMasses.length + 1;
      double[] conditionalExpectations = this.getConditionalExpectations(probabilityMasses, nbSamples);

      XYSeriesCollection xyDataset = new XYSeriesCollection();

      xyDataset.addSeries(this.getFirstOrderLossFunctionXYSeries(min, max, nbSamples, precision));

      for(int i = 0; i < segments; i++)
         xyDataset.addSeries(this.getPiecewiseFirstOrderLossFunctionXYSeriesForSegment(i, probabilityMasses, conditionalExpectations, min, max, minYValue, precision));

      xyDataset.addSeries(this.getPiecewiseFirstOrderLossFunctionErrorXYSeries(min, max, nbSamples, probabilityMasses, conditionalExpectations, precision));
      
      JFreeChart chart = ChartFactory.createXYLineChart("First Order Loss Function", "x", "CL(x)",
            xyDataset, PlotOrientation.VERTICAL, false, true, false);
      ChartFrame frame = new ChartFrame("Empirical complementary loss function",chart);
      frame.setVisible(true);
      frame.setSize(500,400);

      if(saveToDisk){

         XYLineChart lc = new XYLineChart("Piecewise linearization", "x", "CL(x)", xyDataset);

         try {
            File latexFolder = new File("./latex");
            if(!latexFolder.exists()){
               latexFolder.mkdir();
            }
            Writer file = new FileWriter("./latex/pw_folf_graph.tex");
            file.write(lc.toLatex(8, 5));
            file.close();
         } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }
   
   /**
    * Poisson knapsack
    */
   
   /**
    *  Memoization
    */
   static Hashtable<String,double[]> poisson_ce = new Hashtable<String,double[]>();
   static Hashtable<String,Double> poisson_ae = new Hashtable<String,Double>();
   
   public static double[][] poissonKnapsackPiecewiseFOLFConditionalExpectations(int capacity, double[] probabilityMasses, int nbSamples) {
      ArrayList<double[]> conditionalExpectations = new ArrayList<double[]>();
      ArrayList<Double> approximationErrors = new ArrayList<Double>();
      
      int demand = 0;
      while(true) {
         String hash = demand + "_" + Arrays.toString(probabilityMasses) + "_" + nbSamples;
         Distribution[] distributions = new Distribution[1];
         distributions[0] = new PoissonDist(demand);
         PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
         double[] ce = null;
         if(poisson_ce.containsKey(hash)) {
            ce = poisson_ce.get(hash);
         } else {
            ce = pwfolf.getConditionalExpectations(probabilityMasses, nbSamples);
            poisson_ce.put(hash, ce);
         }
         conditionalExpectations.add(ce);
         double ae;
         if(poisson_ae.containsKey(hash)) {
            ae = poisson_ae.get(hash);
         } else {
            ae = pwfolf.getMaxApproximationError(probabilityMasses, nbSamples);
            poisson_ae.put(hash, ae);
         }
         approximationErrors.add(ae);
         if(ce[0] > capacity) break;
         else demand++;
      }
      return conditionalExpectations.toArray(new double[conditionalExpectations.size()][]);
   }
   
   public static double[] poissonKnapsackPiecewiseFOLFApproximationErrors(int capacity, double[] probabilityMasses, int nbSamples) {
      ArrayList<double[]> conditionalExpectations = new ArrayList<double[]>();
      ArrayList<Double> approximationErrors = new ArrayList<Double>();
      
      int demand = 0;
      while(true) {
         String hash = demand + "_" + Arrays.toString(probabilityMasses) + "_" + nbSamples;
         Distribution[] distributions = new Distribution[1];
         distributions[0] = new PoissonDist(demand);
         PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
         double[] ce = null;
         if(poisson_ce.containsKey(hash)) {
            ce = poisson_ce.get(hash);
         } else {
            ce = pwfolf.getConditionalExpectations(probabilityMasses, nbSamples);
            poisson_ce.put(hash, ce);
         }
         conditionalExpectations.add(ce);
         double ae;
         if(poisson_ae.containsKey(hash)) {
            ae = poisson_ae.get(hash);
         } else {
            ae = pwfolf.getMaxApproximationError(probabilityMasses, nbSamples);
            poisson_ae.put(hash, ae);
         }
         approximationErrors.add(ae);
         if(ce[0] > capacity) break;
         else demand++;
      }
      return toPrimitive(approximationErrors.toArray(new Double[conditionalExpectations.size()]));
   }
   
   public static int poissonKnapsackPiecewiseFOLFMaxWeight(int capacity, double[] probabilityMasses, int nbSamples) {
      ArrayList<double[]> conditionalExpectations = new ArrayList<double[]>();
      ArrayList<Double> approximationErrors = new ArrayList<Double>();
      
      int demand = 0;
      while(true) {
         String hash = demand + "_" + Arrays.toString(probabilityMasses) + "_" + nbSamples;
         Distribution[] distributions = new Distribution[1];
         distributions[0] = new PoissonDist(demand);
         PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
         double[] ce = null;
         if(poisson_ce.containsKey(hash)) {
            ce = poisson_ce.get(hash);
         } else {
            ce = pwfolf.getConditionalExpectations(probabilityMasses, nbSamples);
            poisson_ce.put(hash, ce);
         }
         conditionalExpectations.add(ce);
         double ae;
         if(poisson_ae.containsKey(hash)) {
            ae = poisson_ae.get(hash);
         } else {            
            ae = pwfolf.getMaxApproximationError(probabilityMasses, nbSamples);
            poisson_ae.put(hash, ae);
         }
         approximationErrors.add(ae);
         if(ce[0] > capacity) break;
         else demand++;
      }
      return demand;
   }
   
   /**
    * Main
    */
   
   public static void main(String[] args){
      //testPiecewiseFirstOrderLossFunction();
      //printLinearizationParameters();
      //generatePoissonPiecewiseTables();
   }
   
   @SuppressWarnings("unused")
   private static void testPiecewiseFirstOrderLossFunction(){
      Distribution[] distributions = new Distribution[1];
      distributions[0] = new NormalDist(0,1);
      PiecewiseFirstOrderLossFunction pwcfolf = new PiecewiseFirstOrderLossFunction(distributions);
      double[] probabilityMasses = {0.5,0.5};
      
      /* Plot parameters */
      int min = -2;
      int max = 2;
      int minYvalue = -1;
      int samples = 1000;
      double precision = 0.1;
      boolean generateLatex = true;
      
      pwcfolf.plotPiecewiseComplementaryFirstOrderLossFunction(min, max, minYvalue, probabilityMasses, samples, precision, generateLatex);
      pwcfolf.plotPiecewiseFirstOrderLossFunction(min, max, minYvalue, probabilityMasses, samples, precision, generateLatex);
   }
   
   @SuppressWarnings("unused")
   private static void printLinearizationParameters() {
      Distribution[] distributions = new Distribution[3];
      double lambda[] = {20,5,50};
      distributions[0] = new PoissonDist(lambda[0]);
      distributions[1] = new PoissonDist(lambda[1]);
      distributions[2] = new PoissonDist(lambda[2]);
      PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
      
      int partitions = 5;
      int nbSamples = 1000;
      double[] probabilityMasses = new double[partitions]; 
      Arrays.fill(probabilityMasses, 1.0/partitions);
      String out = "float prob[partitions] = ";
      out += Arrays.toString(probabilityMasses) + ";\n";
      out += "float means[partitions] = ";
      out += Arrays.toString(pwfolf.getConditionalExpectations(probabilityMasses, nbSamples)) + ";\n";
      out += "float error = "+pwfolf.getMaxApproximationError(probabilityMasses, nbSamples) + ";\n";
      System.out.print(out);
   }
   
   @SuppressWarnings("unused")
   private static void generatePoissonPiecewiseTables() {
      int capacity = 100;
      int partitions = 5;
      int nbSamples = 10000;
      double[] probabilityMasses = new double[partitions];
      Arrays.fill(probabilityMasses, 1.0/partitions);
      System.out.println("prob = " + Arrays.toString(probabilityMasses) + ";");
      System.out.println("means = " + Arrays.deepToString(poissonKnapsackPiecewiseFOLFConditionalExpectations(capacity, probabilityMasses, nbSamples)) + ";");
      System.out.println("error = " + Arrays.toString(poissonKnapsackPiecewiseFOLFApproximationErrors(capacity, probabilityMasses, nbSamples)) + ";");
      System.out.println("maxWeight = " + poissonKnapsackPiecewiseFOLFMaxWeight(capacity, probabilityMasses, nbSamples) + ";");
   }
   
   /**
    * Auxiliary methods
    */
   
   private static double[] toPrimitive(Double[] array) {
      if (array == null) {
        return null;
      } else if (array.length == 0) {
        return new double[0];
      }
      final double[] result = new double[array.length];
      for (int i = 0; i < array.length; i++) {
        result[i] = array[i].doubleValue();
      }
      return result;
    }
}
