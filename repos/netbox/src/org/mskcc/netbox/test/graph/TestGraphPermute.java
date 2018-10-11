package org.mskcc.netbox.test.graph;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.graph.impl.SparseVertex;
import edu.uci.ics.jung.graph.impl.UndirectedSparseEdge;
import edu.uci.ics.jung.graph.impl.UndirectedSparseGraph;
import junit.framework.TestCase;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.GraphPermute;

/**
 * JUnit Tests for the Graph Permute Class.
 *
 * @author Ethan Cerami.
 */
public class TestGraphPermute extends TestCase {

    /**
     * Tests the Graph Permute Functionality.
     *
     * @throws GraphCreationException Graph Creation Error.
     */
    public final void testGraphPermute() throws GraphCreationException {

        //  Create a sample graph with simple community structure
        Graph g = new UndirectedSparseGraph();

        Vertex vertexA = createVertex(g, "A");  // V0
        Vertex vertexB = createVertex(g, "B");  // V1
        Vertex vertexX = createVertex(g, "X");  // V2
        Vertex vertexY = createVertex(g, "Y");  // V3
        Vertex vertexZ = createVertex(g, "Z");  // V4

        g.addEdge(new UndirectedSparseEdge(vertexA, vertexX));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexY));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexZ));

        g.addEdge(new UndirectedSparseEdge(vertexX, vertexB));
        g.addEdge(new UndirectedSparseEdge(vertexY, vertexB));
        g.addEdge(new UndirectedSparseEdge(vertexZ, vertexB));

        GraphPermute graphPermute = new GraphPermute(g);
        Graph randomGraph = graphPermute.next();
        assertEquals(6, g.getEdges().size());
        assertEquals(6, randomGraph.getEdges().size());
    }

    private Vertex createVertex(Graph g, String label) throws GraphCreationException {
        try {
            StringLabeller labeller = StringLabeller.getLabeller(g);
            Vertex vertex = new SparseVertex();
            g.addVertex(vertex);
            labeller.setLabel(vertex, label);
            return vertex;
        } catch (StringLabeller.UniqueLabelException e) {
            throw new GraphCreationException(e);
        }
    }
}
