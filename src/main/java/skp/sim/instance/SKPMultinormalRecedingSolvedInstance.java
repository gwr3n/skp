package skp.sim.instance;

import skp.instance.SKPMultinormal;

public class SKPMultinormalRecedingSolvedInstance {
   public SKPMultinormal instance;
   
   public double simulatedSolutionValue;
   public int simulationRuns;
   
   public double evp;
   public double evwpi;
   public double evwpi_obj_2_n;
   
   public int piecewisePartitions;
   public int piecewiseSamples;
   
   public SKPMultinormalRecedingSolvedInstance(
         SKPMultinormal instance,
         double simulatedSolutionValue,
         int simulationRuns,
         double evp,
         double evwpi,
         double evwpi_obj_2_n,
         int piecewisePartitions,
         int piecewiseSamples) {
      this.instance = instance;
      this.simulatedSolutionValue = simulatedSolutionValue;
      this.simulationRuns = simulationRuns;
      this.evp = evp;
      this.evwpi = evwpi;
      this.evwpi_obj_2_n = evwpi_obj_2_n;
      this.piecewisePartitions = piecewisePartitions;
      this.piecewiseSamples = piecewiseSamples;
   }
}
