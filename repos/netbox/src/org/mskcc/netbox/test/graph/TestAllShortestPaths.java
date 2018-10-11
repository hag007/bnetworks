package org.mskcc.netbox.test.graph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.graph.impl.SparseVertex;
import edu.uci.ics.jung.graph.impl.UndirectedSparseEdge;
import edu.uci.ics.jung.graph.impl.UndirectedSparseGraph;
import junit.framework.TestCase;
import org.mskcc.netbox.graph.AllShortestPaths;
import org.mskcc.netbox.graph.GraphCreationException;

import java.util.ArrayList;
import java.util.List;

/**
 * JUnit Tests for AllShortestPaths.
 *
 * @author Ethan Cerami.
 */
public class TestAllShortestPaths extends TestCase {
    private static final int THREE = 3;

    /**
     * Tests All Shortest Paths on a Sample Graph.
     *
     * @throws GraphCreationException Graph Creation Error.
     */
    public final void testAllShortestPaths() throws GraphCreationException {

        //  Create a sample graph with simple community structure
        Graph g = new UndirectedSparseGraph();

        Vertex vertexA = createVertex(g, "A");
        Vertex vertexB = createVertex(g, "B");
        Vertex vertexX = createVertex(g, "X");
        Vertex vertexY = createVertex(g, "Y");
        Vertex vertexZ = createVertex(g, "Z");

        g.addEdge(new UndirectedSparseEdge(vertexA, vertexX));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexY));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexZ));

        g.addEdge(new UndirectedSparseEdge(vertexX, vertexB));
        g.addEdge(new UndirectedSparseEdge(vertexY, vertexB));
        g.addEdge(new UndirectedSparseEdge(vertexZ, vertexB));

        AllShortestPaths util = new AllShortestPaths(g, "A", "X", 1);
        ArrayList<List<Edge>> paths = util.getAllShortestPaths();
        assertEquals(1, paths.size());

        util = new AllShortestPaths(g, "A", "B", 2);
        paths = util.getAllShortestPaths();
        assertEquals(THREE, paths.size());

        //for (int i=0; i<paths.size(); i++) {
        //    System.out.println (util.getPathDescription(i));
        //    System.out.println (util.getLinkerLabel(i));
        //}
        assertEquals("A --> X --> B", util.getPathDescription(0));
        assertEquals("X", util.getLinkerLabel(0));
        assertEquals("A --> Y --> B", util.getPathDescription(1));
        assertEquals("Y", util.getLinkerLabel(1));
        assertEquals("A --> Z --> B", util.getPathDescription(2));
        assertEquals("Z", util.getLinkerLabel(2));
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
