package skp.sim.instance;

import skp.instance.SKPNormal;

public class SKPNormalRecedingSolvedInstance {
public SKPNormal instance;
   
   public double simulatedSolutionMean;
   public double simulatedSolutionStd;
   public int simulationRuns;
   
   public double evp;
   public double evwpi;
   public double evwpi_obj_2_n;
   
   public int piecewisePartitions;
   public int piecewiseSamples;
   
   public SKPNormalRecedingSolvedInstance(
         SKPNormal instance,
         double simulatedSolutionMean,
         double simulatedSolutionStd,
         int simulationRuns,
         double evp,
         double evwpi,
         double evwpi_obj_2_n,
         int piecewisePartitions,
         int piecewiseSamples) {
      this.instance = instance;
      this.simulatedSolutionMean = simulatedSolutionMean;
      this.simulatedSolutionStd = simulatedSolutionStd;
      this.simulationRuns = simulationRuns;
      this.evp = evp;
      this.evwpi = evwpi;
      this.evwpi_obj_2_n = evwpi_obj_2_n;
      this.piecewisePartitions = piecewisePartitions;
      this.piecewiseSamples = piecewiseSamples;
   }
}
