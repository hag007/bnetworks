package org.mskcc.netbox.algorithm;

/**
 * Encapsulates Information Regarding Number of Nodes / Number of Edges in a Network
 * or a Network Component.
 *
 * @author Ethan Cerami.
 */
public final class VertexEdgePair {
    private int numVertices;
    private int numEdges;

    /**
     * Gets the number of vertices.
     *
     * @return number of vertices.
     */
    public int getNumVertices() {
        return numVertices;
    }

    /**
     * Sets the number of vertices.
     *
     * @param n number of vertices.
     */
    public void setNumVertices(int n) {
        this.numVertices = n;
    }

    /**
     * Gets the number of Edges.
     *
     * @return number of edges.
     */
    public int getNumEdges() {
        return numEdges;
    }

    /**
     * Sets the number of edges.
     *
     * @param n number of edges.
     */
    public void setNumEdges(int n) {
        this.numEdges = n;
    }
}
