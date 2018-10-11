package org.mskcc.netbox.util;

import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Global Formatter Classes.
 *
 * @author Ethan Cerami.
 */
public final class Formatter {

    /**
     * Private Constructor to prevent instantiation.
     */
    private Formatter() {
    }

    /**
     * Gets the Decimal Formatter.
     *
     * @return Number Format Object.
     */
    public static NumberFormat getDecimalFormat() {
        NumberFormat formatter = new DecimalFormat("#0.00000");
        return formatter;
    }

    /**
     * Gets the Short Decimal Formatter.
     *
     * @return Number Format Object.
     */
    public static NumberFormat getShortDecimalFormat() {
        NumberFormat formatter = new DecimalFormat("#0.00");
        return formatter;
    }

    /**
     * Gets the P-Value Formatter.
     *
     * @return Number Format Object.
     */
    public static NumberFormat getPValueFormat() {
        NumberFormat formatter = new DecimalFormat("0.#####E0");
        return formatter;
    }

}
