package org.mskcc.netbox.graph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.utils.Pair;

import java.util.Iterator;

/**
 * Utility Class for Converting a JUNG Graph to Cytoscape SIF Format.
 *
 * @author Ethan Cerami.
 */
public final class JungToSif {

    /**
     * Private Constructor to prevent instantiation.
     */
    private JungToSif() {
    }

    /**
     * Converts the Specified JUNG Graph to a Cytoscape SIF Format.
     *
     * @param g JUNG Graph.
     * @return SIF Format.
     */
    public static String convertToSif(Graph g) {
        StringBuffer buf = new StringBuffer();
        StringLabeller labeller = StringLabeller.getLabeller(g);
        Iterator edgeIterator = g.getEdges().iterator();
        while (edgeIterator.hasNext()) {
            Edge edge = (Edge) edgeIterator.next();
            Pair pair = edge.getEndpoints();
            Vertex vertexA = (Vertex) pair.getFirst();
            Vertex vertexB = (Vertex) pair.getSecond();
            String geneA = GraphUtil.getVertexLabel(labeller, vertexA);
            String geneB = GraphUtil.getVertexLabel(labeller, vertexB);
            buf.append(geneA + " INTERACTS " + geneB + "\n");
        }
        return buf.toString();
    }
}
