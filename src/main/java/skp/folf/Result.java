package skp.folf;

/**
 * Result container for the minimax partitioning.
 * Holds the segment probabilities, conditional means, and the maximum error.
 */
public class Result {
    /** Segment probabilities \( p_i \) */
    public final double[] p;
    /** Conditional means \( E[Z|\Omega_i] \) */
    public final double[] expect;
    /** Maximum error of the Jensen lower bound */
    public final double error;

    public Result(double[] p, double[] e, double err) {
        this.p = p;
        this.expect = e;
        this.error = err;
    }
}
