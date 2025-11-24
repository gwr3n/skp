package skp.folf;

import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public final class PiecewiseStandardNormalFirstOrderLossFunction {
   
   private static final int linearizationSamples = 0;
   
   private static final Map<Integer, Double> errorCache = new ConcurrentHashMap<>();
   private static final Map<Integer, double[]> probabilitiesCache = new ConcurrentHashMap<>(); 
   private static final Map<Integer, double[]> meansCache = new ConcurrentHashMap<>();
   
   public static int getLinearizationSamples() {
      return linearizationSamples;
   }
   
   public static double getError(int partitions){ 
      if(partitions <= 0) { 
         throw new IllegalArgumentException("partitions must be >= 1"); 
      } 
      return errorCache.computeIfAbsent(partitions, p -> JensenMinimaxPartitioner.compute(p + 1).error); 
   }
   
   /*public static double getError(int partitions){
      if(partitions <= 10) {
         double[] errors = {
            0.3989422804014327,
            0.1206560496714961,
            0.05784405029198253,
            0.033905164962384104,
            0.022270929512393414,
            0.01574607463566398,
            0.011721769576577057,
            0.00906528789647753,
            0.007219916411227892,
            0.005885974956458359
            };
         return errors[partitions-1];
      }else {
         Distribution[] distributions = new Distribution[1];
         distributions[0] = new NormalDist(0,1);
         PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
         double[] probabilityMasses = new double[partitions];
         Arrays.fill(probabilityMasses, 1.0/partitions);
         double error = pwfolf.getMaxApproximationError(probabilityMasses, linearizationSamples);
         return error;
      }
   }*/
   
   public static double[] getProbabilities(int partitions){ 
      if(partitions <= 0) { 
         throw new IllegalArgumentException("partitions must be >= 1"); 
      } 
      double[] cached = probabilitiesCache.computeIfAbsent(partitions, p -> { 
         double[] src = JensenMinimaxPartitioner.compute(p + 1).p; 
         return Arrays.copyOf(src, src.length); // store defensive copy 
         }); 
      return Arrays.copyOf(cached, cached.length); // return defensive copy 
   }
   
   /*public static double[] getProbabilities(int partitions){
      switch(partitions){
      case 1:
         return new double[] {1};
      case 2:
         return new double[] {0.5, 0.5};  
      case 3:
         return new double[] {0.28783338731597996, 0.4243332253680401, 0.28783338731597996};
      case 4:
         return new double[] {0.18755516774758485, 0.31244483225241515, 0.31244483225241515, 0.18755516774758485};
      case 5:
         return new double[] {0.1324110437406592, 0.23491250409192982, 0.26535290433482195, 0.23491250409192987, 0.13241104374065915};
      case 6:
         return new double[] {0.09877694599482933, 0.18223645973091096, 0.2189865942742597, 0.2189865942742597, 0.18223645973091096, 0.09877694599482933};
      case 7:
         return new double[] {0.07669891602586965, 0.14538182479573014, 0.18144834397745296, 0.19294183040189444, 0.18144834397745302, 0.14538182479573014, 0.07669891602586965};
      case 8:
         return new double[] {0.061394553470121016, 0.11872108750901467, 0.15205073490726895, 0.16783362411359537, 0.16783362411359537, 0.15205073490726895, 0.11872108750901467, 0.061394553470121016};
      case 9:
         return new double[] {0.05033061540430428, 0.09884442068481658, 0.12900389832341652, 0.1460368860056377, 0.15156835916364986, 0.1460368860056377, 0.12900389832341652, 0.09884442068481658, 0.05033061540430428};
      case 10:
         return new double[] {0.04206108420763477, 0.0836356495308449, 0.11074334596058821, 0.1276821455299152, 0.13587777477101692, 0.13587777477101692, 0.1276821455299152, 0.11074334596058821, 0.0836356495308449, 0.04206108420763477};
      default:
         {
            double[] probabilities = new double[partitions];
            Arrays.fill(probabilities, 1.0/partitions);
            return probabilities;
         }
      }
   }*/
   
   public static double[] getMeans(int partitions){ 
      if(partitions <= 0) { 
         throw new IllegalArgumentException("partitions must be >= 1"); 
      } 
      double[] cached = meansCache.computeIfAbsent(partitions, p -> { 
         double[] src = JensenMinimaxPartitioner.compute(p + 1).expect; 
         return Arrays.copyOf(src, src.length); // store defensive copy 
      }); 
      return Arrays.copyOf(cached, cached.length); // return defensive copy 
   }
   
   /*public static double[] getMeans(int partitions){
      switch(partitions){
      case 1:
         return new double[] {0};
      case 2:
         return new double[] {-0.7978845608028654, 0.7978845608028654};  
      case 3:
         return new double[] {-1.1850544278068644, 0, 1.1850544278068644};
      case 4:
         return new double[] {-1.4353532729205845, -0.41522324304905966, 0.41522324304905966, 1.4353532729205845};
      case 5:
         return new double[] {-1.6180463502161044, -0.6914240068499904, 0, 0.6914240068499903, 1.6180463502161053};
      case 6:
         return new double[] {-1.7608020666235031, -0.8960107374480083, -0.28188851144130117, 0.28188851144130117, 0.8960107374480083, 1.7608020666235031};
      case 7:
         return new double[] {-1.8773528492652836, -1.0572304450884658, -0.4934048390251067, 0, 0.4934048390251065, 1.0572304450884658, 1.8773528492652836};
      case 8:
         return new double[] {-1.9754694729585056, -1.1895340795157716, -0.6615516528578579, -0.213586638906901, 0.213586638906901, 0.6615516528578579, 1.1895340795157716, 1.9754694729585056};
      case 9:
         return new double[] {-2.059957433491476, -1.30127090280595, -0.8004000560466271, -0.3845969617811554, 0, 0.3845969617811554, 0.8004000560466271, 1.30127090280595, 2.059957433491476};
      case 10:
         return new double[] {-2.133986195498256, -1.3976822972668839, -0.918199946431143, -0.5265753462727588, -0.17199013069262026, 0.17199013069262026, 0.5265753462727588, 0.918199946431143, 1.3976822972668839, 2.133986195498256};
      default:
         {
            Distribution[] distributions = new Distribution[1];
            distributions[0] = new NormalDist(0,1);
            PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
            double[] probabilityMasses = new double[partitions];
            Arrays.fill(probabilityMasses, 1.0/partitions);
            double[] means = pwfolf.getConditionalExpectations(probabilityMasses, linearizationSamples);
            return means;
         }
      }
   }*/
}
