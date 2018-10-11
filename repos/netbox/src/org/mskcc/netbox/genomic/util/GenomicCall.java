package org.mskcc.netbox.genomic.util;

import org.mskcc.netbox.util.GlobalConfig;

/**
 * Utility Class for Reading in Mutation Calls.
 */
public final class GenomicCall {

    private GenomicCall() {
    }

    /**
     * Translates a Mutation String into a Mutation Call.
     *
     * @param value a value, such as:  R787*
     * @return 0 or 1.
     */
    public static byte getMutationCall(String value) {
        byte byteValue;
        if (value.equals("0") || value.equals("NA")) {
            byteValue = 0;
        } else {
            byteValue = 1;
        }
        return byteValue;
    }

    /**
     * Translates an Expression Value into a Byte Call.
     *
     * @param s a value, such as:  0.0001
     * @return 0 or 1.
     */
    public static byte getExpressionCall(String s) {
        double value = Math.abs(Double.parseDouble(s));
        GlobalConfig config = GlobalConfig.getInstance();
        if (value > config.getMRNAZScoreThreshold()) {
            return 1;
        } else {
            return 0;
        }
    }
}
