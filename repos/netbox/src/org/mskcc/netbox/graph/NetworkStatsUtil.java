package org.mskcc.netbox.graph;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import org.mskcc.netbox.model.NetworkStats;
import org.mskcc.netbox.query.InteractionQuery;
import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.ProgressMonitor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

/**
 * Network Stats Utility Class.
 */
public final class NetworkStatsUtil {
    private HashMap<String, Integer> degreeMap = new HashMap<String, Integer>();
    private ArrayList<String> genesInNetwork = new ArrayList<String>();
    private NetworkStats networkStats;
    private static NetworkStatsUtil util;

    private NetworkStatsUtil() throws GraphCreationException {
        ProgressMonitor pMonitor = ProgressMonitor.getInstance();
        pMonitor.setCurrentMessage("Getting Network Stats");
        Graph g = InteractionQuery.getGlobalGraph();
        int numNodes = g.numVertices();
        int numEdges = g.numEdges();
        networkStats = new NetworkStats();
        networkStats.setNumGenes(numNodes);
        networkStats.setNumEdges(numEdges);

        Set<Vertex> vertexSet = g.getVertices();
        StringLabeller labeller = StringLabeller.getLabeller(g);

        pMonitor.setMaxValue(vertexSet.size());
        for (Vertex vertex : vertexSet) {
            String geneSymbol = labeller.getLabel(vertex);
            degreeMap.put(geneSymbol, vertex.degree());
            genesInNetwork.add(geneSymbol);
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }
    }

    /**
     * Gets Instance of Singleton.
     *
     * @return NetworkStatsUtil Object.
     * @throws GraphCreationException Graph Creation Error.
     */
    public static NetworkStatsUtil getInstance() throws GraphCreationException {
        if (util == null) {
            util = new NetworkStatsUtil();
        }
        return util;
    }

    /**
     * Gets the Global Degree of the Specified Gene.
     *
     * @param geneSymbol Gene Symbol.
     * @return global degree.
     */
    public int getGeneDegree(String geneSymbol) {
        if (degreeMap.containsKey(geneSymbol)) {
            return degreeMap.get(geneSymbol);
        } else {
            return 0;
        }
    }

    /**
     * Gets all genes in the reference network.
     *
     * @return ArrayList of Genes.
     */
    public ArrayList<String> getGenesInNetwork() {
        return genesInNetwork;
    }

    /**
     * Gets the Global Network Stats.
     *
     * @return Network Stats Object.
     */
    public NetworkStats getNetworkStats() {
        return networkStats;
    }
}
