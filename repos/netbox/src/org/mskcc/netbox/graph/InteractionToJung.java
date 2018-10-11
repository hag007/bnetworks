package org.mskcc.netbox.graph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.graph.impl.SparseVertex;
import edu.uci.ics.jung.graph.impl.UndirectedSparseEdge;
import edu.uci.ics.jung.graph.impl.UndirectedSparseGraph;
import edu.uci.ics.jung.utils.Pair;
import edu.uci.ics.jung.utils.UserDataContainer;
import org.mskcc.netbox.model.Interaction;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * Utility Class for converting Interaction Objects to JUNG Graphs.
 *
 * @author Ethan Cerami.
 */
public final class InteractionToJung {
    /**
     * Gene Key.
     */
    public static final String GENE_KEY = "GENE_KEY";

    /**
     * Private constructor to prevent instantiation.
     */
    private InteractionToJung() {
    }

    /**
     * Creates a JUNG Graph from the Specified List of Interactions.
     *
     * @param interactionList ArrayList of Interaction Objects.
     * @return JUNG Graph.
     * @throws GraphCreationException Graph Creation Error.
     */
    public static Graph createGraph(ArrayList<Interaction> interactionList)
            throws GraphCreationException {
        Graph g = new UndirectedSparseGraph();
        StringLabeller labeller = StringLabeller.getLabeller(g);

        HashSet<String> edgeSet = new HashSet<String>();
        for (Interaction interaction : interactionList) {
            String geneA = interaction.getGeneA();
            String geneB = interaction.getGeneB();

            //  Do not add duplicate edges
            String key = createKey(geneA, geneB);
            if (!edgeSet.contains(key)) {
                Vertex vertexA = getVertex(g, labeller, geneA);
                Vertex vertexB = getVertex(g, labeller, geneB);
                g.addEdge(new UndirectedSparseEdge(vertexA, vertexB));
                edgeSet.add(key);
            }
        }
        return g;
    }

    /**
     * Adds the Specified List of Interactions to the JUNG Graph.
     *
     * @param g Graph.
     * @param interactionList ArrayList of Interaction Objects.
     * @return JUNG Graph.
     * @throws GraphCreationException Graph Creation Error.
     */
    public static Graph appendGraph(Graph g, ArrayList<Interaction> interactionList)
            throws GraphCreationException {
        StringLabeller labeller = StringLabeller.getLabeller(g);

        HashSet<String> edgeSet = new HashSet<String>();
        Set<Edge> originalEdgeSet = g.getEdges();
        for (Edge edge : originalEdgeSet) {
            Pair p = edge.getEndpoints();
            Vertex v1 = (Vertex) p.getFirst();
            Vertex v2 = (Vertex) p.getSecond();
            String node1 = labeller.getLabel(v1);
            String node2 = labeller.getLabel(v2);
            String key = createKey(node1, node2);
            edgeSet.add(key);
        }

        for (Interaction interaction : interactionList) {
            String geneA = interaction.getGeneA();
            String geneB = interaction.getGeneB();

            //  Do not add duplicate edges
            String key = createKey(geneA, geneB);
            if (!edgeSet.contains(key)) {
                Vertex vertexA = getVertex(g, labeller, geneA);
                Vertex vertexB = getVertex(g, labeller, geneB);
                g.addEdge(new UndirectedSparseEdge(vertexA, vertexB));
                edgeSet.add(key);
            }
        }
        return g;
    }

    private static String createKey(String geneA, String geneB) {
        String first = null;
        String second = null;
        if (geneA.compareTo(geneB) < 0) {
            first = geneA;
            second = geneB;
        } else {
            first = geneB;
            second = geneA;
        }
        return first + ":" + second;
    }

    private static Vertex getVertex(Graph g, StringLabeller labeller, String gene)
            throws GraphCreationException {
        try {
            Vertex vertex = labeller.getVertex(gene);
            if (vertex == null) {
                vertex = new SparseVertex();
                vertex = g.addVertex(vertex);
                labeller.setLabel(vertex, gene);
                vertex.setUserDatum(GENE_KEY, gene, new UserDataContainer.CopyAction.Shared());
            }
            return vertex;
        } catch (StringLabeller.UniqueLabelException e) {
            throw new GraphCreationException(e);
        }
    }
}
