package org.mskcc.netbox.algorithm;

import edu.uci.ics.jung.algorithms.cluster.ClusterSet;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;

import java.util.Set;

/**
 * Utility to Identify the Size of the Largest Component within a Network.
 *
 * @author Ethan Cerami.
 */
public final class LargestComponentUtil {

    private LargestComponentUtil() {
    }

    /**
     * Determines the Size of the Largest Component within a Network.
     * @param g Graph.
     * @return VertexEdgePair Object.
     */
    public static VertexEdgePair determineSizeOfLargestComponent(Graph g) {
        WeakComponentClusterer wcSearch = new WeakComponentClusterer();
        ClusterSet clusterSet = wcSearch.extract(g);

        int numVertices = 0;
        int numEdges = 0;
        for (int i = 0; i < clusterSet.size(); i++) {
            Set<Vertex> set = clusterSet.getCluster(i);
            Graph subGraph = clusterSet.getClusterAsNewSubGraph(i);
            if (set.size() > numVertices) {
                numVertices = set.size();
                numEdges = subGraph.numEdges();
            }
        }
        VertexEdgePair p = new VertexEdgePair();
        p.setNumVertices(numVertices);
        p.setNumEdges(numEdges);
        return p;
    }
}
