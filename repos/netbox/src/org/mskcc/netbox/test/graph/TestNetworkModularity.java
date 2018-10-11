package org.mskcc.netbox.test.graph;

import edu.uci.ics.jung.algorithms.cluster.ClusterSet;
import edu.uci.ics.jung.algorithms.cluster.EdgeBetweennessClusterer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.graph.impl.SparseVertex;
import edu.uci.ics.jung.graph.impl.UndirectedSparseEdge;
import edu.uci.ics.jung.graph.impl.UndirectedSparseGraph;
import junit.framework.TestCase;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.MiniEdge;
import org.mskcc.netbox.graph.NetworkModularity;
import org.mskcc.netbox.graph.NetworkModularityFast;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Tests the NetworkModularity Utility Class.
 *
 * @author Ethan Cerami
 */
public final class TestNetworkModularity extends TestCase {
    private static final double PRECISION = 0.00001;
    private static final double EXPECTED_1 = 0.42307;
    private static final double EXPECTED_2 = 0.5;
    private static final double EXPECTED_3 = 0.25;
    private static final double EXPECTED_4 = 0.1171875;
    private static final double ZERO = 0.0;
    private static final int ONE = 1;
    private static final int TWO = 2;
    private static final int FIVE = 5;
    private static final int SIX = 6;
    private static final int EIGHT = 8;

    /**
     * Test Network Modularity, 1.
     */
    public void testNetworkModularity1() {

        //  Create a sample graph with simple community structure
        Graph g = new UndirectedSparseGraph();

        //  Community 1
        Vertex vertexA = new SparseVertex();
        Vertex vertexB = new SparseVertex();
        Vertex vertexC = new SparseVertex();
        Vertex vertexD = new SparseVertex();

        //  Community 2
        Vertex vertexE = new SparseVertex();
        Vertex vertexF = new SparseVertex();
        Vertex vertexG = new SparseVertex();
        Vertex vertexH = new SparseVertex();

        g.addVertex(vertexA);
        g.addVertex(vertexB);
        g.addVertex(vertexC);
        g.addVertex(vertexD);
        g.addVertex(vertexE);
        g.addVertex(vertexF);
        g.addVertex(vertexG);
        g.addVertex(vertexH);

        //  Community 1 Edges
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexB));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexC));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexD));
        g.addEdge(new UndirectedSparseEdge(vertexB, vertexC));
        g.addEdge(new UndirectedSparseEdge(vertexC, vertexD));
        g.addEdge(new UndirectedSparseEdge(vertexD, vertexB));

        //  Community 2 Edges
        g.addEdge(new UndirectedSparseEdge(vertexE, vertexF));
        g.addEdge(new UndirectedSparseEdge(vertexE, vertexG));
        g.addEdge(new UndirectedSparseEdge(vertexE, vertexH));
        g.addEdge(new UndirectedSparseEdge(vertexG, vertexF));
        g.addEdge(new UndirectedSparseEdge(vertexG, vertexH));
        g.addEdge(new UndirectedSparseEdge(vertexH, vertexF));

        //  Single edge connecting the two communities
        g.addEdge(new UndirectedSparseEdge(vertexD, vertexE));

        //  Cluster with Newman EdgeBetweenness Clusterer:  remove exactly 1 edge
        EdgeBetweennessClusterer clusterer = new EdgeBetweennessClusterer(1);
        ClusterSet clusterSet = clusterer.extract(g);

        double q = NetworkModularity.calculateNetworkModularity(g, clusterSet, false);
        assertEquals(EXPECTED_1, q, PRECISION);
    }

    /**
     * Test Network Modularity, 2.
     *
     * @throws GraphCreationException Graph Creation Error.
     */
    public void testNetworkModularity2() throws GraphCreationException {

        //  Create a sample graph with simple community structure
        Graph g = new UndirectedSparseGraph();

        //  Community 1
        Vertex vertexA = createVertex(g, "A");
        Vertex vertexB = createVertex(g, "B");
        Vertex vertexC = createVertex(g, "C");

        //  Community 2
        Vertex vertexD = createVertex(g, "D");
        Vertex vertexE = createVertex(g, "E");
        Vertex vertexF = createVertex(g, "F");

        //  Community 1 Edges
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexB));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexC));
        g.addEdge(new UndirectedSparseEdge(vertexB, vertexC));

        //  Community 2 Edges
        g.addEdge(new UndirectedSparseEdge(vertexD, vertexE));
        g.addEdge(new UndirectedSparseEdge(vertexD, vertexF));
        g.addEdge(new UndirectedSparseEdge(vertexE, vertexF));

        //  Cluster with Newman EdgeBetweenness Clusterer:  remove 0 edges
        EdgeBetweennessClusterer clusterer = new EdgeBetweennessClusterer(0);
        ClusterSet clusterSet = clusterer.extract(g);

        double q = NetworkModularity.calculateNetworkModularity(g, clusterSet, false);
        assertEquals(EXPECTED_2, q, PRECISION);

        //  Test second Network Modularity Method.  Should return identical results
        ArrayList<String> moduleList = new ArrayList<String>();
        moduleList.add("MODULE_1");
        moduleList.add("MODULE_2");
        HashMap<String, String> moduleMap = new HashMap<String, String>();
        moduleMap.put("A", "MODULE_1");
        moduleMap.put("B", "MODULE_1");
        moduleMap.put("C", "MODULE_1");
        moduleMap.put("D", "MODULE_2");
        moduleMap.put("E", "MODULE_2");
        moduleMap.put("F", "MODULE_2");

        //  This should return the same exact q as above.
        q = NetworkModularity.calculateNetworkModularity(g, moduleList, moduleMap, false);
        assertEquals(EXPECTED_2, q, PRECISION);

        //  Test third Network Modularity Method.  Should return identical results
        NetworkModularityFast fastUtil = new NetworkModularityFast(g);
        q = fastUtil.calculateNetworkModularity(moduleList, moduleMap, false, false);
        assertEquals(EXPECTED_2, q, PRECISION);

        //  Test Fourth Network Modularity Method.
        moduleMap.put("B", "MODULE_2");
        q = fastUtil.calculateNetworkModularity(moduleList, moduleMap, false, false);
        assertEquals(0.11111111, q, PRECISION);

        moduleMap.put("C", "MODULE_2");
        q = fastUtil.calculateNetworkModularity(moduleList, moduleMap, false, false);
        assertEquals(-0.0555555, q, PRECISION);

        moduleMap.put("B", "MODULE_1");
        moduleMap.put("C", "MODULE_1");
        fastUtil.setBaseLine(moduleList, moduleMap);
        q = fastUtil.updateNetworkModularityOneNodeMove(moduleList, "MODULE_1", "MODULE_2", "B");
        fastUtil.accept();
        assertEquals(0.11111111, q, PRECISION);
        q = fastUtil.updateNetworkModularityOneNodeMove(moduleList, "MODULE_1", "MODULE_2", "C");
        assertEquals(-0.0555555, q, PRECISION);
        fastUtil.accept();
        q = fastUtil.updateNetworkModularityOneNodeMove(moduleList, "MODULE_2", "MODULE_1", "C");
        fastUtil.accept();
        q = fastUtil.updateNetworkModularityOneNodeMove(moduleList, "MODULE_2", "MODULE_1", "B");
        assertEquals(0.5, q, PRECISION);
    }

    /**
     * Test Network Modularity, 3.
     */
    public void testNetworkModularity3() {

        //  Create a sample graph with simple community structure
        Graph g = new UndirectedSparseGraph();

        //  Community 1
        Vertex vertexA = new SparseVertex();
        Vertex vertexB = new SparseVertex();
        Vertex vertexC = new SparseVertex();

        //  Community 2
        Vertex vertexD = new SparseVertex();
        Vertex vertexE = new SparseVertex();
        Vertex vertexF = new SparseVertex();

        g.addVertex(vertexA);
        g.addVertex(vertexB);
        g.addVertex(vertexC);
        g.addVertex(vertexD);
        g.addVertex(vertexE);
        g.addVertex(vertexF);

        //  Community 1 Edges
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexB));
        g.addEdge(new UndirectedSparseEdge(vertexA, vertexC));
        g.addEdge(new UndirectedSparseEdge(vertexB, vertexC));

        //  Community 2 Edges
        g.addEdge(new UndirectedSparseEdge(vertexD, vertexE));
        g.addEdge(new UndirectedSparseEdge(vertexD, vertexF));
        g.addEdge(new UndirectedSparseEdge(vertexE, vertexF));

        //  Two edges connecting the two communities
        g.addEdge(new UndirectedSparseEdge(vertexB, vertexE));
        g.addEdge(new UndirectedSparseEdge(vertexC, vertexD));

        //  Cluster with Newman EdgeBetweenness Clusterer:  remove 2 edges
        EdgeBetweennessClusterer clusterer = new EdgeBetweennessClusterer(0);
        ClusterSet clusterSet = clusterer.extract(g);
        assertEquals(1, clusterSet.size());

        double q = NetworkModularity.calculateNetworkModularity(g, clusterSet, false);
        assertEquals(ZERO, q, PRECISION);

        clusterer = new EdgeBetweennessClusterer(1);
        clusterSet = clusterer.extract(g);
        assertEquals(ONE, clusterSet.size());

        q = NetworkModularity.calculateNetworkModularity(g, clusterSet, false);
        assertEquals(ZERO, q, PRECISION);

        clusterer = new EdgeBetweennessClusterer(2);
        clusterSet = clusterer.extract(g);
        assertEquals(TWO, clusterSet.size());

        q = NetworkModularity.calculateNetworkModularity(g, clusterSet, false);
        assertEquals(EXPECTED_3, q, PRECISION);

        clusterer = new EdgeBetweennessClusterer(FIVE);
        clusterSet = clusterer.extract(g);
        assertEquals(3, clusterSet.size());

        q = NetworkModularity.calculateNetworkModularity(g, clusterSet, false);
        assertEquals(EXPECTED_4, q, PRECISION);

        clusterer = new EdgeBetweennessClusterer(EIGHT);
        clusterSet = clusterer.extract(g);
        assertEquals(SIX, clusterSet.size());

        q = NetworkModularity.calculateNetworkModularity(g, clusterSet, false);
        assertEquals(-0.171875, q, PRECISION);
    }

    /**
     * More modularity tests.
     *
     * @throws GraphCreationException Graph Creation Error.
     */
    public void testNetworkModularity4() throws GraphCreationException {
        //  Create a sample graph with simple community structure
        Graph g = new UndirectedSparseGraph();

        Vertex ppara = createVertex(g, "PPARA");
        Vertex frfr1op = createVertex(g, "FGFR1OP");
        Vertex ptpre = createVertex(g, "PTPRE");
        Vertex cct2 = createVertex(g, "CCT2");
        Vertex a2m = createVertex(g, "A2M");
        Vertex dctn2 = createVertex(g, "DCTN2");
        Vertex cct6a = createVertex(g, "CCT6A");
        Vertex ncam1 = createVertex(g, "NCAM1");
        Vertex kcna5 = createVertex(g, "KCNA5");
        Vertex tubgcp2 = createVertex(g, "TUBGCP2");
        Vertex stac3 = createVertex(g, "STAC3");
        Vertex tubgcp6 = createVertex(g, "TUBGCP6");
        Vertex cntn2 = createVertex(g, "CNTN2");
        Vertex ptprd = createVertex(g, "PTPRD");
        Vertex lyz = createVertex(g, "LYZ");


        g.addEdge(new UndirectedSparseEdge(lyz, a2m));
        g.addEdge(new UndirectedSparseEdge(stac3, ppara));
        g.addEdge(new UndirectedSparseEdge(cct6a, cct2));
        g.addEdge(new UndirectedSparseEdge(cntn2, ncam1));
        g.addEdge(new UndirectedSparseEdge(kcna5, ptpre));
        g.addEdge(new UndirectedSparseEdge(ptpre, ptprd));

        g.addEdge(new UndirectedSparseEdge(dctn2, tubgcp2));
        g.addEdge(new UndirectedSparseEdge(dctn2, frfr1op));
        g.addEdge(new UndirectedSparseEdge(dctn2, tubgcp6));
        g.addEdge(new UndirectedSparseEdge(tubgcp6, tubgcp2));
        g.addEdge(new UndirectedSparseEdge(tubgcp6, frfr1op));
        g.addEdge(new UndirectedSparseEdge(frfr1op, tubgcp2));


        ArrayList<String> moduleList = new ArrayList<String>();
        moduleList.add("MODULE_0");

        HashMap<String, String> moduleMap = new HashMap<String, String>();
        moduleMap.put("FGFR1OP", "MODULE_0");
        moduleMap.put("PTPRE", "MODULE_0");
        moduleMap.put("A2M", "MODULE_0");
        moduleMap.put("DCTN2", "MODULE_0");
        moduleMap.put("CCT6A", "MODULE_0");
        moduleMap.put("NCAM1", "MODULE_0");
        moduleMap.put("KCNA5", "MODULE_0");
        moduleMap.put("TUBGCP2", "MODULE_0");
        moduleMap.put("STAC3", "MODULE_0");
        moduleMap.put("TUBGCP6", "MODULE_0");
        moduleMap.put("CNTN2", "MODULE_0");
        moduleMap.put("PTPRD", "MODULE_0");
        moduleMap.put("LYZ", "MODULE_0");

        double q = NetworkModularity.calculateNetworkModularity(g, moduleList, moduleMap, false);
        assertEquals(-0.006944, q, PRECISION);

        moduleList = new ArrayList<String>();
        moduleList.add("MODULE_0");
        moduleList.add("MODULE_1");
        moduleMap = new HashMap<String, String>();
        moduleMap.put("FGFR1OP", "MODULE_0");
        moduleMap.put("PTPRE", "MODULE_0");
        moduleMap.put("DCTN2", "MODULE_0");
        moduleMap.put("CCT6A", "MODULE_0");
        moduleMap.put("NCAM1", "MODULE_0");
        moduleMap.put("KCNA5", "MODULE_0");
        moduleMap.put("TUBGCP2", "MODULE_0");
        moduleMap.put("STAC3", "MODULE_0");
        moduleMap.put("TUBGCP6", "MODULE_0");
        moduleMap.put("CNTN2", "MODULE_0");
        moduleMap.put("PTPRD", "MODULE_0");
        moduleMap.put("A2M", "MODULE_1");
        moduleMap.put("LYZ", "MODULE_1");
        q = NetworkModularity.calculateNetworkModularity(g, moduleList, moduleMap, false);
        assertEquals(0.13194, q, PRECISION);
    }

    /**
     * Tests the Mini Edge Equals / Hashcode Functionality.
     */
    public void testMiniEdge() {
        HashSet<MiniEdge> edgeSet = new HashSet<MiniEdge>();
        MiniEdge edge0 = new MiniEdge("A", "B");
        MiniEdge edge1 = new MiniEdge("B", "A");
        edgeSet.add(edge0);
        edgeSet.add(edge1);
        assertTrue(edge0.equals(edge1));
        assertEquals(1, edgeSet.size());
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
