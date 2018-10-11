package org.mskcc.netbox.algorithm;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.utils.UserDataContainer;
import org.apache.commons.math.distribution.HypergeometricDistribution;
import org.apache.commons.math.distribution.HypergeometricDistributionImpl;
import org.apache.log4j.Logger;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.InteractionToJung;
import org.mskcc.netbox.graph.InteractionUtil;
import org.mskcc.netbox.graph.NetworkStatsUtil;
import org.mskcc.netbox.model.Interaction;
import org.mskcc.netbox.stats.BenjaminiHochbergFDR;
import org.mskcc.netbox.util.Formatter;
import org.mskcc.netbox.util.ProgressMonitor;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;

/**
 * Algorithm for Connecting a list of Genes into a Network.
 *
 * @author Ethan Cerami.
 */
public final class GeneConnector {
    /**
     * Gene Status.
     */
    public static final String STATUS = "STATUS";

    /**
     * Gene is Altered.
     */
    public static final String ALTERED = "ALTERED";

    /**
     * Gene is a Linker.
     */
    public static final String LINKER = "LINKER";

    private Graph g;
    private ArrayList<String> alteredGeneList;
    private HashSet<String> alteredGeneSet = new HashSet<String>();
    private ArrayList<LinkerGene> linkerList;
    private HashSet<String> alteredGenesInNetwork = new HashSet<String>();
    private NetworkStatsUtil networkStatsUtil;
    private static NumberFormat numberFormat = Formatter.getPValueFormat();
    private double pValueCutOff;
    private int shortestPathThreshold;
    private static Logger log = Logger.getLogger(GeneConnector.class);
    private ProgressMonitor pMonitor = ProgressMonitor.getInstance();
    private int numAlteredGenes;
    private int numLinkerGenes;

    /**
     * Constructor.
     *
     * @param list Gene List.
     * @param sp   Shortest Path Threshold
     * @param p    P-Value Threshold.
     * @throws GraphCreationException Graph Creation Error.
     */
    public GeneConnector(ArrayList<String> list, int sp, double p)
            throws GraphCreationException {
        this.alteredGeneList = new ArrayList<String>();
        this.linkerList = new ArrayList<LinkerGene>();
        this.pValueCutOff = p;
        this.shortestPathThreshold = sp;
        if (shortestPathThreshold < 1 || shortestPathThreshold > 2) {
            throw new IllegalArgumentException("Invalid value for shortestPathThreshold:  "
                    + shortestPathThreshold);
        }
        this.networkStatsUtil = NetworkStatsUtil.getInstance();

        pMonitor.setCurrentMessage("Input gene list consists of " + list.size() + " genes.");
        for (String gene : list) {
            int degree = networkStatsUtil.getGeneDegree(gene);
            if (degree > 0) {
                alteredGeneList.add(gene);
            }
        }
        this.alteredGeneSet.addAll(alteredGeneList);
        pMonitor.setCurrentMessage("Of these genes, " + alteredGeneList.size() + " are in "
                + "the reference network.");

        //  Get all interactions associated with all altered genes and create a JUNG Graph
        ArrayList<Interaction> interactionList =
                InteractionUtil.getAllInteractions(alteredGeneList);
        log.info("With an input list of " + alteredGeneList.size() + ", I am starting "
                + "out with " + interactionList.size() + " interactions.");

        g = InteractionToJung.createGraph(interactionList);

        //  Prune the graph of any un-needed linkers
        StringLabeller labeller = StringLabeller.getLabeller(g);
        pruneGraph(g, shortestPathThreshold, labeller);

        //  Interconnect the remaining linkers
        ArrayList<String> linkerStr = new ArrayList<String>();
        for (LinkerGene linker : linkerList) {
            if (linker.isInNetwork()) {
                linkerStr.add(linker.getGene());
            }
        }
        if (linkerStr.size() > 0) {
            ArrayList<Interaction> linkerInteractionList =
                    InteractionUtil.connectGenesNoLinkers(linkerStr);
            InteractionToJung.appendGraph(g, linkerInteractionList);
        }
    }

    /**
     * Get the Graph.
     *
     * @return Graph Object.
     */
    public Graph getGraph() {
        return g;
    }

    /**
     * Get the gene list.
     *
     * @return ArrayList of Gene Symbols.
     */
    public ArrayList<String> getAlteredGeneList() {
        return alteredGeneList;
    }

    /**
     * Get the linker gene list.
     *
     * @return ArrayList of Linker Gene Symbols.
     */
    public ArrayList<LinkerGene> getLinkerList() {
        return linkerList;
    }

    private void pruneGraph(Graph graph, int sp, StringLabeller labeller) {
        int totalNumGenes = 0;
        if (sp == 2) {
            totalNumGenes = (int) networkStatsUtil.getNetworkStats().getNumGenes();
        }
        Set<Vertex> vertexSet = graph.getVertices();

        //  Mark all vertices as ALTERED or LINKER
        for (Vertex vertex : vertexSet) {
            String label = labeller.getLabel(vertex);
            if (alteredGeneSet.contains(label)) {
                alteredGenesInNetwork.add(label);
                vertex.setUserDatum(STATUS, ALTERED, new UserDataContainer.CopyAction.Shared());
            } else {
                vertex.setUserDatum(STATUS, LINKER, new UserDataContainer.CopyAction.Shared());
            }
        }
        int numAlteredGenesInNetwork = alteredGenesInNetwork.size();

        ArrayList<Vertex> removalList = new ArrayList<Vertex>();
        if (sp == 2) {

            //  Prune Linker Genes that do not pass threshold
            for (Vertex vertex : vertexSet) {
                String status = (String) vertex.getUserDatum(STATUS);
                String label = labeller.getLabel(vertex);
                if (status.equals(LINKER)) {
                    Set<Vertex> neighborSet = vertex.getNeighbors();
                    int counter = 0;

                    //  Count number of neighbors for this linker in observed network
                    for (Vertex neighborVertex : neighborSet) {
                        String neighborLabel = labeller.getLabel(neighborVertex);
                        if (alteredGeneList.contains(neighborLabel)) {
                            counter++;
                        }
                    }

                    if (counter > 1) {
                        //  Get the global degree for this linker
                        int globalDegree = networkStatsUtil.getGeneDegree(label);

                        //  Calculate the Hypergeometric distribution for the linker
                        //  (population size, number of successes, sample size)
                        HypergeometricDistribution hyper = new HypergeometricDistributionImpl(
                                totalNumGenes, globalDegree, numAlteredGenesInNetwork);
                        double pValue = hyper.probability(counter);
                        log.info("Gene: " + label + " Local:  " + counter
                                + " global:  " + globalDegree + ", pvalue:  " + pValue);
                        log.info("new HypergeometricDistributionImpl("
                                + totalNumGenes + "," + globalDegree + ","
                                + numAlteredGenesInNetwork + ");");
                        LinkerGene linker = new LinkerGene(label, counter, globalDegree, pValue);
                        linkerList.add(linker);
                    } else {
                        //  If the vertex is only connected to 1 altered node, then it's not a
                        //  linker.  Remove it.
                        removeVertex(removalList, vertex);
                    }
                }
            }

            //  FDR Adjustment for Linkers
            Collections.sort(linkerList, new LinkerGeneComparator());
            double[] pValues = new double[linkerList.size()];
            for (int i = 0; i < linkerList.size(); i++) {
                pValues[i] = linkerList.get(i).getUnadjustedPValue();
            }
            BenjaminiHochbergFDR fdr = new BenjaminiHochbergFDR(pValues);
            fdr.calculate();
            double[] adjustedPValues = fdr.getAdjustedPvalues();
            for (int i = 0; i < linkerList.size(); i++) {
                linkerList.get(i).setFdrAdjustedPValue(adjustedPValues[i]);
            }

            for (LinkerGene linker : linkerList) {
                if (linker.getFdrAdjustedPValue() > pValueCutOff) {
                    Vertex vertex = labeller.getVertex(linker.getGene());
                    removeVertex(removalList, vertex);
                    linker.setInNetwork(false);
                } else {
                    linker.setInNetwork(true);
                }
            }
        } else {
            //  At shortest path = 1, remove all linkers.
            for (Vertex vertex : vertexSet) {
                String status = (String) vertex.getUserDatum(STATUS);
                if (status.equals(LINKER)) {
                    removeVertex(removalList, vertex);
                }
            }
        }

        //  Now perform the actual removal of the linkers that don't make the cut
        for (Vertex vertex : removalList) {
            graph.removeVertex(vertex);
        }

        //  Remove Any Remaining Isolates
        removalList = new ArrayList<Vertex>();
        numAlteredGenes = 0;
        numLinkerGenes = 0;
        for (Vertex vertex : vertexSet) {
            if (vertex.numNeighbors() == 0) {
                removalList.add(vertex);
            } else {
                String status = (String) vertex.getUserDatum(STATUS);
                if (status.equals(ALTERED)) {
                    numAlteredGenes++;
                } else if (status.equals(LINKER)) {
                    numLinkerGenes++;
                }
            }
        }
        for (Vertex vertex : removalList) {
            graph.removeVertex(vertex);
        }

        if (sp == 1) {
            pMonitor.setCurrentMessage("At shortest path threshold of:  "
                    + sp
                    + ", I can "
                    + "connect " + numAlteredGenes + " genes.");
        } else {
            pMonitor.setCurrentMessage("At shortest path threshold of:  "
                    + sp + " and p-value cut-off of:  "
                    + numberFormat.format(pValueCutOff)
                    + ", I can "
                    + "connect " + numAlteredGenes + " genes with "
                    + numLinkerGenes + " linker genes.");
        }
    }

    /**
     * Gets Number of Altered Genes within the Network.
     * @return number of altered genes.
     */
    public int getNumAlteredGenes() {
        return numAlteredGenes;
    }

    /**
     * Gets Number of Linker Genes within the Network.
     * @return number of linker genes.
     */
    public int getNumLinkerGenes() {
        return numLinkerGenes;
    }

    private void removeVertex(ArrayList<Vertex> removalList, Vertex vertex) {
        //  Remove Incident Edges
        Set<Edge> edgeSet = vertex.getIncidentEdges();
        for (Edge edge : edgeSet) {
            g.removeEdge(edge);
        }
        //  Mark vertex for removal;  can't remove them now;  otherwise, you
        //  get ConcurrentModificationException
        removalList.add(vertex);
    }
}

/**
 * Compares Linker Genes, based on P-Values.
 */
class LinkerGeneComparator implements Comparator {

    public int compare(Object o1, Object o2) {
        LinkerGene linker1 = (LinkerGene) o1;
        LinkerGene linker2 = (LinkerGene) o2;
        return Double.compare(linker1.getUnadjustedPValue(), linker2.getUnadjustedPValue());
    }
}
