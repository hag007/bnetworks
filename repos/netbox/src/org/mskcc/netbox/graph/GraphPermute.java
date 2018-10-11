package org.mskcc.netbox.graph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.impl.UndirectedSparseEdge;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Generates Locally Permuted Graphs.
 */
public final class GraphPermute {
    private Graph graph;
    private Random randomGenerator;

    /**
     * Constructor.
     *
     * @param g Graph Object.
     */
    public GraphPermute(Graph g) {
        randomGenerator = new Random(10);
        this.graph = (Graph) g.copy();
    }

    /**
     * Gets the Next Random Graph.
     *
     * @return Random Graph.
     */
    public Graph next() {
        //  Copy all Vertices
        Set<Edge> edgeSet = graph.getEdges();
        ArrayList<Edge> edgeList = new ArrayList<Edge>();
        edgeList.addAll(edgeSet);

        int i = 0;
        while (i < edgeSet.size()) {

            //  Randomly Pick two Edges
            int randomIndex0 = randomGenerator.nextInt(edgeList.size());
            int randomIndex1 = randomGenerator.nextInt(edgeList.size());

            // Edge 0:  A <--> B
            Edge edge0 = edgeList.get(randomIndex0);

            // Edge 1:  C <--> D
            Edge edge1 = edgeList.get(randomIndex1);

            if (swapWillSucceed(edge0, edge1)) {
                Vertex vertexA = (Vertex) edge0.getEndpoints().getFirst();
                Vertex vertexB = (Vertex) edge0.getEndpoints().getSecond();
                Vertex vertexC = (Vertex) edge1.getEndpoints().getFirst();
                Vertex vertexD = (Vertex) edge1.getEndpoints().getSecond();

                //  Rewire:  A <--> D
                Edge newEdge0 = graph.addEdge(new UndirectedSparseEdge(vertexA, vertexD));
                edgeList.add(newEdge0);

                //  Rewire:  C <--> B
                Edge newEdge1 = graph.addEdge(new UndirectedSparseEdge(vertexC, vertexB));
                edgeList.add(newEdge1);

                //  Remove the old edges
                graph.removeEdge(edge0);
                edgeList.remove(edge0);

                graph.removeEdge(edge1);
                edgeList.remove(edge1);

                i++;
            }
        }

        return graph;
    }

    private boolean swapWillSucceed(Edge edge0, Edge edge1) {
        Vertex vertexA = (Vertex) edge0.getEndpoints().getFirst();
        Vertex vertexB = (Vertex) edge0.getEndpoints().getSecond();
        Vertex vertexC = (Vertex) edge1.getEndpoints().getFirst();
        Vertex vertexD = (Vertex) edge1.getEndpoints().getSecond();

        HashSet<Vertex> vertexSet = new HashSet<Vertex>();
        vertexSet.add(vertexA);
        vertexSet.add(vertexB);
        vertexSet.add(vertexC);
        vertexSet.add(vertexD);

        if (vertexSet.size() == 4) {
            //  If the edges that we want to create already exist, the swap will not succeed.
            if (vertexA.isNeighborOf(vertexD) || vertexC.isNeighborOf(vertexB)) {
                return false;
            } else {
                return true;
            }
        } else {
            return false;
        }
    }
}
