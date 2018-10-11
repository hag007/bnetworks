package org.mskcc.netbox.graph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * Gets All Shortest Paths Between A and B.
 *
 * @author Ethan Cerami.
 */
public final class AllShortestPaths {
    private ArrayList<List<Edge>> paths;
    private StringLabeller labeller;
    private Vertex startVertex;
    private Vertex endVertex;
    private String start;
    private String end;

    /**
     * Constructor.
     *
     * @param g                     Graph
     * @param s                     Start
     * @param e                     End
     * @param shortestPathThreshold Shortest Path Threshold. 1 or 2.
     */
    public AllShortestPaths(Graph g, String s, String e,
                            int shortestPathThreshold) {
        this.start = s;
        this.end = e;
        paths = new ArrayList<List<Edge>>();
        labeller = StringLabeller.getLabeller(g);
        if (shortestPathThreshold < 1 || shortestPathThreshold > 2) {
            throw new IllegalArgumentException("shortestPathThreshold must be set to 1 or 2.");
        }
        startVertex = labeller.getVertex(start);
        if (startVertex == null) {
            return;
        }
        endVertex = labeller.getVertex(end);
        if (endVertex == null) {
            return;
        }
        Set incidentEdgeSet = startVertex.getIncidentEdges();
        if (incidentEdgeSet.size() > 0) {
            Iterator iterator = incidentEdgeSet.iterator();
            while (iterator.hasNext()) {
                Edge edge = (Edge) iterator.next();
                Vertex neighborVertex = edge.getOpposite(startVertex);
                String neighborLabel = labeller.getLabel(neighborVertex);
                if (shortestPathThreshold == 1) {
                    if (neighborLabel.equalsIgnoreCase(end)) {
                        List<Edge> edgeList = new ArrayList<Edge>();
                        edgeList.add(edge);
                        paths.add(edgeList);
                    }
                } else {
                    Set twoHopsEdgeSet = neighborVertex.getIncidentEdges();
                    if (twoHopsEdgeSet.size() > 0) {
                        Iterator twoHopsEdgeIterator = twoHopsEdgeSet.iterator();
                        while (twoHopsEdgeIterator.hasNext()) {
                            Edge twoHopEdge = (Edge) twoHopsEdgeIterator.next();
                            Vertex twoHopNeighborVertex = twoHopEdge.getOpposite(neighborVertex);
                            String twoHopNeighborLabel = labeller.getLabel(twoHopNeighborVertex);
                            if (twoHopNeighborLabel.equalsIgnoreCase(end)
                                    && !twoHopNeighborLabel.equalsIgnoreCase(start)) {
                                List<Edge> edgeList = new ArrayList<Edge>();
                                edgeList.add(edge);
                                edgeList.add(twoHopEdge);
                                paths.add(edgeList);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Gets a list of all the shortest paths.
     *
     * @return ArrayList of List of Edges.
     */
    public ArrayList<List<Edge>> getAllShortestPaths() {
        return paths;
    }

    /**
     * Gets a Description of Specified Path.
     *
     * @param i Index value.
     * @return String description, e.g. A--> B --> C.
     */
    public String getPathDescription(int i) {
        List<Edge> path = paths.get(i);
        if (path.size() == 1) {
            return (start + " --> " + end);
        } else {
            Edge edge = path.get(0);
            Vertex linker = edge.getOpposite(startVertex);
            return (start + " --> " + labeller.getLabel(linker) + " --> " + end);
        }
    }

    /**
     * Gets the Linker Vertex for the Specified Path.
     *
     * @param i Index Value.
     * @return Linker Vertex.
     */
    public Vertex getLinkerVertex(int i) {
        List<Edge> path = paths.get(i);
        if (path.size() == 1) {
            return null;
        } else {
            Edge edge = path.get(0);
            Vertex linker = edge.getOpposite(startVertex);
            return linker;
        }
    }

    /**
     * Gets the Linker Label for the Specified Path.
     *
     * @param i Index Value.
     * @return Linker Label.
     */
    public String getLinkerLabel(int i) {
        Vertex linkerVertex = getLinkerVertex(i);
        return labeller.getLabel(linkerVertex);
    }
}
