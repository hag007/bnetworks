package org.mskcc.netbox.graph;

/**
 * Encapsulates Minimal Information Regarding an Edge.
 *
 * @author Ethan Cerami.
 */
public final class MiniEdge {
    private String nodeA;
    private String nodeB;
    private String hash;

    /**
     * Constructor.
     *
     * @param a Vertex A.
     * @param b Vertex B.
     */
    public MiniEdge(String a, String b) {
        nodeA = a;
        nodeB = b;

        if (nodeA.compareTo(nodeB) > 0) {
            hash = nodeA + ":" + nodeB;
        } else {
            hash = nodeB + ":" + nodeA;
        }
    }

    /**
     * Gets Node A.
     *
     * @return Node A.
     */
    public String getNodeA() {
        return nodeA;
    }

    /**
     * Gets Node B.
     *
     * @return Node B.
     */
    public String getNodeB() {
        return nodeB;
    }

    /**
     * Gets Opposite of Specified Node.
     *
     * @param node Node.
     * @return opposite node.
     */
    public String getOpposite(String node) {
        if (node.equals(nodeA)) {
            return nodeB;
        } else {
            return nodeA;
        }
    }

    /**
     * Hash Case.
     *
     * @return hase code.
     */
    public int hashCode() {
        return hash.hashCode();
    }

    /**
     * Override equals methods.
     *
     * @param o Object.
     * @return true or false.
     */
    public boolean equals(Object o) {
        MiniEdge e = (MiniEdge) o;
        if (e.hashCode() == this.hashCode()) {
            return true;
        } else {
            return false;
        }
    }
}
