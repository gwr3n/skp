package skp.instance;

import java.util.Arrays;
import java.util.stream.IntStream;

import umontreal.ssj.probdist.NormalDist;

public class SKPNormal extends SKP {

   double[] stdWeights;
   double capacity;
   
   public SKPNormal(double[] expectedValues, double[] expectedWeights, double[] stdWeights, double capacity, double shortageCost) {
      super(expectedValues, expectedWeights, shortageCost);
      this.stdWeights = stdWeights;
      this.capacity = capacity;
      
      generateInstanceID();
   }
   
   public SKPNormal(double[] expectedValues, double[] expectedWeights, double coefficientOfVariation, double capacity, double shortageCost) {
      this(expectedValues, 
            expectedWeights,
            IntStream.iterate(0, i -> i + 1)
                     .limit(expectedWeights.length)
                     .mapToDouble(i -> expectedWeights[i]*coefficientOfVariation)
                     .toArray(),
            capacity, 
            shortageCost);
   }
   
   public SKPNormal(double[] expectedValues, NormalDist[] weights, double capacity, double shortageCost) {
      this(expectedValues,
            IntStream.iterate(0, i -> i + 1)
                     .limit(weights.length)
                     .mapToDouble(i -> weights[i].getMean())
                     .toArray(),
            IntStream.iterate(0, i -> i + 1)
                     .limit(weights.length)
                     .mapToDouble(i -> weights[i].getSigma())
                     .toArray(),
            capacity,
            shortageCost);
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public NormalDist[] getWeights() {
      return IntStream.iterate(0, i -> i + 1).limit(expectedWeights.length)
                      .mapToObj(i -> new NormalDist(expectedWeights[i], stdWeights[i]))
                      .toArray(NormalDist[]::new);
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.expectedValues),
                                        Arrays.hashCode(this.expectedWeights),
                                        Arrays.hashCode(this.stdWeights),
                                        Double.hashCode(this.capacity),
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPNormal) {
         SKPNormal o = (SKPNormal) obj;
         return Arrays.equals(this.expectedValues, o.expectedValues) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                Arrays.equals(this.stdWeights, o.stdWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPNormal getTestInstance() {
      double[] expectedValues = {111, 111, 21, 117, 123, 34, 3, 121, 112, 12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      int capacity = 100;
      int shortageCost = 10;
      return new SKPNormal(expectedValues, expectedWeights, cv, capacity, shortageCost);
   }
   
   public static SKPNormal getTestOptimalityGapInstance() {
      double[] expectedValues = {30.0, 81.0, 66.0, 24.0, 69.0, 36.0, 93.0, 39.0, 9.0,  45.0, 66.0, 39.0, 18.0, 39.0, 42.0, 78.0, 78.0, 84.0, 18.0, 42.0, 51.0, 87.0, 30.0, 18.0, 96.0};
      double[] expectedWeights = {27.374614896234103, 78.4130955596743, 64.65813191442086, 23.266107454511886, 67.4093146834377, 33.59693259866023, 92.16485054308757, 37.79821214150361, 8.236141621162524, 42.6355367987863, 65.54323036758973, 36.508548229182615, 15.968166079227474, 38.83708767781832, 39.062391091831344, 75.03159125698056, 77.3931165367291, 82.83151230238252, 17.34237700822168, 40.044541913379156, 49.92889765445393, 85.56454225080667, 27.426648865161223, 17.446202319304014, 94.46271171962006};
      double cv = 0.2;
      double capacity = 112.6313185;
      int shortageCost = 10;
      return new SKPNormal(expectedValues, expectedWeights, cv, capacity, shortageCost);
   }
}
