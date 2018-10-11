package org.mskcc.netbox.graph;

/**
 * Encapsulates a Specific Network Parition State.
 *
 * @author Ethan Cerami.
 */
public final class NetworkPartitionState {
    private int numEdgesRemoved;
    private double networkModularity;
    private int numModules;

    /**
     * Gets Number of Edges Removed.
     *
     * @return Number of Edges Removed.
     */
    public int getNumEdgesRemoved() {
        return numEdgesRemoved;
    }

    /**
     * Sets Number of Edges Removed.
     *
     * @param n Number of Edges Removed.
     */
    public void setNumEdgesRemoved(int n) {
        this.numEdgesRemoved = n;
    }

    /**
     * Gets the Network Modularity.
     *
     * @return network modularity.
     */
    public double getNetworkModularity() {
        return networkModularity;
    }

    /**
     * Sets the Network Modularity.
     *
     * @param q Network modularity.
     */
    public void setNetworkModularity(double q) {
        this.networkModularity = q;
    }

    /**
     * Gets the number of modules.
     *
     * @return number of modules.
     */
    public int getNumModules() {
        return numModules;
    }

    /**
     * Sets the number of modules.
     *
     * @param n number of modules.
     */
    public void setNumModules(int n) {
        this.numModules = n;
    }
}
