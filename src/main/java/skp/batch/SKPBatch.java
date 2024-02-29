package skp.batch;

import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.rng.MRG32k3aL;

public abstract class SKPBatch {

   protected static final long[] seed = {1,2,3,4,5,6};
   protected static final MRG32k3aL randGenerator = new MRG32k3aL();

   protected Distribution expectedValue;
   protected Distribution expectedWeight;
   protected DiscreteDistributionInt capacity;
   protected Distribution shortageCost;
   
   public SKPBatch(
         Distribution expectedValue,
         Distribution expectedWeight,
         DiscreteDistributionInt capacity,
         Distribution shortageCost
         ) {
      this.expectedValue = expectedValue;
      this.expectedWeight = expectedWeight;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
   }
}