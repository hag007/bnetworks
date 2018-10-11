package org.mskcc.netbox.test.stats;

import junit.framework.TestCase;
import org.apache.commons.math.distribution.HypergeometricDistribution;
import org.apache.commons.math.distribution.HypergeometricDistributionImpl;
import org.mskcc.netbox.stats.BenjaminiHochbergFDR;

/**
 * Tests the Benjamin Hochberg FDR Correction.
 *
 * @author Ethan Cerami.
 */
public final class TestBenjaminiHochbergFDR extends TestCase {

    /**
     * Test FDR Correction.
     */
    public void testFDRCorrection() {
        double[] pValues = new double[5];
        pValues[0] = 0.001;
        pValues[1] = 0.002;
        pValues[2] = 0.004;
        pValues[3] = 0.04;
        pValues[4] = 0.1;
        BenjaminiHochbergFDR fdrCorrector = new BenjaminiHochbergFDR(pValues);
        fdrCorrector.calculate();

        double[] adjustedPValues = fdrCorrector.getAdjustedPvalues();
        double[] unAdjustedPValues = fdrCorrector.getOrdenedPvalues();

        // Output from R.  Verify that Java output matches the R Output
        //[1,] 0.001 0.005000000
        //[2,] 0.002 0.005000000
        //[3,] 0.004 0.006666667
        //[4,] 0.040 0.050000000
        //[5,] 0.100 0.100000000
        assertEquals(0.001, unAdjustedPValues[0]);
        assertEquals(0.002, unAdjustedPValues[1]);
        assertEquals(0.004, unAdjustedPValues[2]);
        assertEquals(0.040, unAdjustedPValues[3]);
        assertEquals(0.100, unAdjustedPValues[4]);

        assertEquals(0.005, adjustedPValues[0]);
        assertEquals(0.005, adjustedPValues[1]);
        assertEquals(0.006666667, adjustedPValues[2], 0.00001);
        assertEquals(0.050, adjustedPValues[3]);
        assertEquals(0.100, adjustedPValues[4]);

        pValues = new double[7];
        pValues[0] = 0.01;
        pValues[1] = 0.002;
        pValues[2] = 0.004;
        pValues[3] = 0.04;
        pValues[4] = 0.001;
        pValues[5] = 0.9;
        pValues[6] = 0.1;
        fdrCorrector = new BenjaminiHochbergFDR(pValues);
        fdrCorrector.calculate();

        adjustedPValues = fdrCorrector.getAdjustedPvalues();
        unAdjustedPValues = fdrCorrector.getOrdenedPvalues();

        assertEquals(0.116666667, adjustedPValues[5], 0.0001);
        //        for (int i=0; i<adjustedPValues.length; i++) {
        //            System.out.println (unAdjustedPValues[i] + " " + adjustedPValues[i]);
        //        }

        HypergeometricDistribution hyper = new HypergeometricDistributionImpl(9264, 14, 222);
        double pValue = hyper.probability(2);
        // System.out.println (pValue);

    }
}
