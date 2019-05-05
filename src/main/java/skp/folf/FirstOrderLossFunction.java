package skp.folf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import umontreal.ssj.probdist.EmpiricalDist;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.rng.MRG32k3aL;

public class FirstOrderLossFunction {
   Distribution[] distributions;
   long[] seed;
   MRG32k3aL randGenerator;
   
   FirstOrderLossFunction(Distribution[] distributions, long[] seed){
      this.distributions = distributions;
      this.randGenerator = new MRG32k3aL();
      this.randGenerator.setSeed(seed);
   }
   
   private double[][] sample(int nbSamples){
      this.randGenerator.resetStartStream();
      double[][] sampleMatrix = new double[nbSamples][this.distributions.length];
      for(int i = 0; i < sampleMatrix.length; i++){
         for(int j = 0; j < sampleMatrix[i].length; j++){
            sampleMatrix[i][j] = distributions[j].inverseF(UniformGen.nextDouble(this.randGenerator, 0, 1));
         }
      }
      return sampleMatrix;
   }
   
   protected EmpiricalDist getEmpiricalDistribution(int nbSamples){
      double[][] sampleMatrix = this.sample(nbSamples);
      double[] observations = new double[nbSamples];
      for(int i = 0; i < sampleMatrix.length; i++){
         for(int j = 0; j < sampleMatrix[i].length; j++){
            observations[i] += sampleMatrix[i][j];
         }
      }
      Arrays.sort(observations);
      EmpiricalDist empDistribution = new EmpiricalDist(observations);
      return empDistribution;
   }
   
   public double getComplementaryFirstOrderLossFunctionValue(double x, int nbSamples){
      EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples);
      double value = 0;
      for(int i = 0; i < empDistribution.getN(); i++){
         value += Math.max(x-empDistribution.getObs(i),0)/empDistribution.getN();
      }
      return value;
   }
   
   public double getFirstOrderLossFunctionValue(double x, int nbSamples){
      EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples);
      double value = 0;
      for(int i = 0; i < empDistribution.getN(); i++){
         value += Math.max(empDistribution.getObs(i)-x,0)/empDistribution.getN();
      }
      return value;
   }
   
   /**
    * Charting
    */
   
   private XYSeries getDistributionXYSeries(int nbSamples, double precision){
      XYSeries series = new XYSeries("Empirical distribution");
      EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples);
      for(int i = 0; i < empDistribution.getN(); i++){
         while(i>0 && i<empDistribution.getN() && empDistribution.getObs(i)==empDistribution.getObs(i-1))i++;
         series.add(empDistribution.getObs(i),empDistribution.cdf(empDistribution.getObs(i)));
      }
      return series;
   }
   
   private void plotEmpiricalDistribution(int nbSamples, double precision, boolean saveToDisk){
      XYSeriesCollection xyDataset = new XYSeriesCollection(this.getDistributionXYSeries(nbSamples, precision));
      JFreeChart chart = ChartFactory.createXYLineChart("Empirical distribution", "Support", "Frequency",
             xyDataset, PlotOrientation.VERTICAL, false, true, false);
      ChartFrame frame = new ChartFrame("Empirical distribution",chart);
      frame.setVisible(true);
      frame.setSize(500,400);
      
      if(saveToDisk){

         XYLineChart lc = new XYLineChart("Empirical Distribution", "x", "f(x)", xyDataset);

         try {
            File latexFolder = new File("./latex");
            if(!latexFolder.exists()){
               latexFolder.mkdir();
            }
            Writer file = new FileWriter("./latex/emp_graph.tex");
            file.write(lc.toLatex(8, 5));
            file.close();
         } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }
   
   protected XYSeries getComplementaryFirstOrderLossFunctionXYSeries(double min, double max, int nbSamples, double precision){
      XYSeries series = new XYSeries("Empirical complementary loss function");
      for(double x = min; x <= max; x+= precision){
         series.add(x, getComplementaryFirstOrderLossFunctionValue(x, nbSamples));
      }
      return series;
   }
   
   protected XYSeries getFirstOrderLossFunctionXYSeries(double min, double max, int nbSamples, double precision){
      XYSeries series = new XYSeries("Empirical loss function");
      for(double x = min; x <= max; x+= precision){
         series.add(x, getFirstOrderLossFunctionValue(x, nbSamples));
      }
      return series;
   }
   
   private void plotEmpiricalComplementaryFirstOrderLossFunction(double min, double max, int nbSamples, double precision, boolean saveToDisk){
      XYSeriesCollection xyDataset = new XYSeriesCollection(this.getComplementaryFirstOrderLossFunctionXYSeries(min, max, nbSamples, precision));
      JFreeChart chart = ChartFactory.createXYLineChart("Empirical Complementary First Order Loss Function", "x", "CL(x)",
             xyDataset, PlotOrientation.VERTICAL, false, true, false);
      ChartFrame frame = new ChartFrame("Empirical complementary first order loss function",chart);
      frame.setVisible(true);
      frame.setSize(500,400);
      
      if(saveToDisk){

         XYLineChart lc = new XYLineChart("Empirical Complementary First Order Loss Function", "x", "CL(x)", xyDataset);

         try {
            File latexFolder = new File("./latex");
            if(!latexFolder.exists()){
               latexFolder.mkdir();
            }
            Writer file = new FileWriter("./latex/cfolf_graph.tex");
            file.write(lc.toLatex(8, 5));
            file.close();
         } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }
   
   private void plotEmpiricalFirstOrderLossFunction(double min, double max, int nbSamples, double precision, boolean saveToDisk){
      XYSeriesCollection xyDataset = new XYSeriesCollection(this.getFirstOrderLossFunctionXYSeries(min, max, nbSamples, precision));
      JFreeChart chart = ChartFactory.createXYLineChart("Empirical First Order Loss Function", "x", "L(x)",
             xyDataset, PlotOrientation.VERTICAL, false, true, false);
      ChartFrame frame = new ChartFrame("Empirical complementary loss function",chart);
      frame.setVisible(true);
      frame.setSize(500,400);
      
      if(saveToDisk){

         XYLineChart lc = new XYLineChart("Empirical First Order Loss Function", "x", "L(x)", xyDataset);

         try {
            File latexFolder = new File("./latex");
            if(!latexFolder.exists()){
               latexFolder.mkdir();
            }
            Writer file = new FileWriter("./latex/folf_graph.tex");
            file.write(lc.toLatex(8, 5));
            file.close();
         } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }
   
   /**
    * Main
    */
   
   public static void main(String[] args){
      testDistributionPlot();
      testLossFunctionPlot();
   }
   
   private static void testDistributionPlot(){
      long[] seed = {1,2,3,4,5,6};
      Distribution[] distributions = new Distribution[3];
      double lambda[] = {20,5,50};
      distributions[0] = new PoissonDist(lambda[0]);
      distributions[1] = new PoissonDist(lambda[1]);
      distributions[2] = new PoissonDist(lambda[2]);
      FirstOrderLossFunction lf = new FirstOrderLossFunction(distributions, seed);
      
      /* Plot parameters */
      int samples = 1000;
      double precision = 1;
      boolean generateLatex = true;
      
      lf.plotEmpiricalDistribution(samples, precision, generateLatex);
   }
   
   private static void testLossFunctionPlot(){
      long[] seed = {1,2,3,4,5,6};
      Distribution[] distributions = new Distribution[3];
      double lambda[] = {20,5,50};
      distributions[0] = new PoissonDist(lambda[0]);
      distributions[1] = new PoissonDist(lambda[1]);
      distributions[2] = new PoissonDist(lambda[2]);
      FirstOrderLossFunction lf = new FirstOrderLossFunction(distributions, seed);
      
      /* Plot parameters */
      int min = 50;
      int max = 100;
      int samples = 1000;
      double precision = 1;
      boolean generateLatex = true;
      
      lf.plotEmpiricalComplementaryFirstOrderLossFunction(min, max, samples, precision, generateLatex);
      lf.plotEmpiricalFirstOrderLossFunction(min, max, samples, precision, generateLatex);
   }
}
