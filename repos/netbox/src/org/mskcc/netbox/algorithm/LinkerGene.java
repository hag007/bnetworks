package org.mskcc.netbox.algorithm;

/**
 * Encapsulates Information Regrarding a Linker Gene.
 *
 * @author Ethan Cerami.
 */
public final class LinkerGene {
    private String gene;
    private int localDegree;
    private int globalDegree;
    private double unadjustedPValue;
    private double fdrAdjustedPValue;
    private boolean inNetwork;

    /**
     * Constructor.
     *
     * @param g       Gene Name.
     * @param lDegree Local Degree.
     * @param gDegree Global Degree.
     * @param p       Unadjusted P-Value.
     */
    public LinkerGene(String g, int lDegree, int gDegree, double p) {
        this.gene = g;
        this.localDegree = lDegree;
        this.globalDegree = gDegree;
        this.unadjustedPValue = p;
    }

    /**
     * Gets the Gene Symbol.
     *
     * @return Gene Symbol.
     */
    public String getGene() {
        return gene;
    }

    /**
     * Gets the Local Degree.
     *
     * @return local degree.
     */
    public int getLocalDegree() {
        return localDegree;
    }

    /**
     * Gets the Global Degree.
     *
     * @return global degree.
     */
    public int getGlobalDegree() {
        return globalDegree;
    }

    /**
     * Is the linker currently is the network or pruned from the network.
     *
     * @return true or false.
     */
    public boolean isInNetwork() {
        return inNetwork;
    }

    /**
     * Gets the Unadjusted P-Value.
     *
     * @return unadjusted p-value.
     */
    public double getUnadjustedPValue() {
        return unadjustedPValue;
    }

    /**
     * Gets the FDR Adjusted P-Value.
     *
     * @return FDR Adjusted P-Value.
     */
    public double getFdrAdjustedPValue() {
        return fdrAdjustedPValue;
    }

    /**
     * Sets the FDR Adjusted P-Value.
     *
     * @param p FDR Adjusted P-Value.
     */
    public void setFdrAdjustedPValue(double p) {
        this.fdrAdjustedPValue = p;
    }

    /**
     * Sets the status of whether the linker is in the network or pruned from the network.
     *
     * @param in true or false.
     */
    public void setInNetwork(boolean in) {
        this.inNetwork = in;
    }
}
