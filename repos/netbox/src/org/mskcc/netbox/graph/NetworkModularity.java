package org.mskcc.netbox.graph;

import edu.uci.ics.jung.algorithms.cluster.ClusterSet;
import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Calculates Network Modularity.
 * <p/>
 * Based on M. E. J. Newman, M. Girvan, Finding and evaluating community structure in networks
 * Phys. Rev. E 69, 026113 (2004)
 * http://arxiv.org/abs/cond-mat/0308217
 */
public final class NetworkModularity {

    /**
     * Private Constructor to prevent instantiation.
     */
    private NetworkModularity() {
    }

    /**
     * Calculates Network Modularity of Given Graph and Cluster Set.
     *
     * @param graph     Graph graph.
     * @param moduleSet Module Set Object.
     * @param debug     Debug Flag.
     * @return network modularity score.
     */
    public static double calculateNetworkModularity(Graph graph, ClusterSet moduleSet,
                                                    boolean debug) {
        double q = 0.0;
        int numEdges = graph.numEdges();

        if (debug) {
            System.out.println("L = " + numEdges);
        }

        //  Iterate through all modules
        for (int i = 0; i < moduleSet.size(); i++) {
            Graph g = moduleSet.getClusterAsNewSubGraph(i);
            int l = g.numEdges();
            int d = 0;
            Set<Vertex> vertexSet = g.getVertices();
            for (Vertex vertex : vertexSet) {
                //  Get Vertex Degree in the Original Graph
                Vertex vertexInOriginalGraph = (Vertex) vertex.getEqualVertex(graph);
                d += vertexInOriginalGraph.degree();
            }
            double currentQ = ((l / (double) numEdges))
                    - (Math.pow(d / ((double) 2 * numEdges), 2));
            if (debug) {
                System.out.println("l = " + l);
                System.out.println("d = " + d);
                System.out.println("currentQ = " + currentQ);
            }

            q += currentQ;
        }
        return q;
    }

    /**
     * Calculates Network Modularity of Given Graph and Cluster Set.
     * This version returns identical results, but is more efficient than
     * calculateNetworkModularity(Graph graph, ClusterSet moduleSet) because it does not
     * require the use or creation of Cluster Set Objects.
     *
     * @param graph Graph graph.
     * @param moduleList module list.
     * @param moduleMap module map.
     * @param debug Debug flag.
     * @return network modularity.
     */
    public static double calculateNetworkModularity(Graph graph, ArrayList<String> moduleList,
                    HashMap<String, String> moduleMap, boolean debug) {
        int numEdges = graph.numEdges();
        StringLabeller labeller = StringLabeller.getLabeller(graph);
        if (debug) {
            System.out.println("L = " + numEdges);
        }

        //  Create an inverse look-up map
        HashMap<String, ArrayList<String>> nodesInModule =
                new HashMap<String, ArrayList<String>>();
        for (String gene : moduleMap.keySet()) {
            String moduleId = moduleMap.get(gene);
            ArrayList<String> nodeList = nodesInModule.get(moduleId);
            if (nodeList == null) {
                nodeList = new ArrayList<String>();
                nodeList.add(gene);
                nodesInModule.put(moduleId, nodeList);
            } else {
                nodeList.add(gene);
            }
        }

        double q = 0.0;

        //  Iterate through all modules
        for (String moduleId : moduleList) {
            if (debug) {
                System.out.println("Module ID:  " + moduleId);
            }
            HashSet<Edge> edgesInModule = new HashSet<Edge>();
            int l = 0;
            int d = 0;

            ArrayList<String> nodeList = nodesInModule.get(moduleId);
            if (nodeList != null) {
                //  Iterate through all nodes in module
                for (String node : nodeList) {
                    Vertex v = labeller.getVertex(node);

                    //  Now iterate through all adjacent vertices.
                    d += v.degree();
                    Set<Edge> edgeSet = v.getIncidentEdges();
                    for (Edge edge : edgeSet) {
                        Vertex other = edge.getOpposite(v);
                        String neighborLabel = GraphUtil.getVertexLabel(labeller, other);
                        if (nodeList.contains(neighborLabel)) {
                            edgesInModule.add(edge);
                        }
                    }
                }
            }

            l = edgesInModule.size();
            double currentQ = ((l / (double) numEdges))
                    - (Math.pow(d / ((double) 2 * numEdges), 2));
            if (debug) {
                System.out.println("l = " + l);
                System.out.println("d = " + d);
                System.out.println("q = " + currentQ);
            }
            q += currentQ;
        }
        return q;
    }
}
