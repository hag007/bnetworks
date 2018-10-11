package org.mskcc.netbox.graph;

import edu.uci.ics.jung.algorithms.cluster.ClusterSet;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.importance.BetweennessCentrality;
import edu.uci.ics.jung.algorithms.importance.EdgeRanking;
import edu.uci.ics.jung.graph.ArchetypeGraph;
import edu.uci.ics.jung.graph.Graph;
import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.Formatter;
import org.mskcc.netbox.util.ProgressMonitor;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Newman-Girvan Module Detector.
 *
 * @author Ethan Cerami, and JUNG Community Authors.
 */
public final class NewmanGirvanModuleDetector {
    private List mEdgesRemoved;
    private ArrayList<NetworkPartitionState> networkPartionStateList =
            new ArrayList<NetworkPartitionState>();
    private double maxModularityScore = 0.0;
    private NetworkPartitionState optimalPartitionState = null;
    private ClusterSet optimalClusterSet = null;
    private Graph optimalGraph;
    private boolean initialStateRepresentsSingleConnectedComponent = false;

    /**
     * Constructs a new clusterer for the specified graph.
     *
     * @param g Graph Object.
     */
    public NewmanGirvanModuleDetector(ArchetypeGraph g) {
        mEdgesRemoved = new ArrayList();
        extract(g);
    }

    /**
     * Finds the set of clusters which have the strongest "community structure".
     * The more edges removed the smaller and more cohesive the clusters.
     *
     * @param g Graph Object.
     */
    private void extract(ArchetypeGraph g) {
        if (!(g instanceof Graph)) {
            throw new IllegalArgumentException("Argument must be of type Graph.");
        }

        Graph originalGraph = (Graph) g;

        //  First, init set to remove all edges
        int mNumEdgesToRemove = originalGraph.getEdges().size();
        mEdgesRemoved.clear();

        //  Calculate Base Line Stats Before Any Partitioning
        WeakComponentClusterer wcSearch = new WeakComponentClusterer();
        ClusterSet clusterSet = wcSearch.extract(originalGraph);
        NetworkPartitionState state = calculateNetworkModularity(originalGraph, clusterSet, 0);
        if (clusterSet.size() == 1) {
            initialStateRepresentsSingleConnectedComponent = true;
        }
        // outputState(state);
        networkPartionStateList.add(state);

        //  Set baseline
        maxModularityScore = state.getNetworkModularity();
        optimalPartitionState = state;
        optimalClusterSet = clusterSet;
        optimalGraph = (Graph) originalGraph.copy();

        //  Execute algorithm, remove all edges
        executeAlgorithm(originalGraph, mNumEdgesToRemove, state.getNumModules());
        mEdgesRemoved.clear();
    }

    private void executeAlgorithm(Graph originalGraph, int numEdgesToRemove,
                                  int currentNumModules) {
        Graph graphCopy = (Graph) originalGraph.copy();

        ProgressMonitor pMonitor = ProgressMonitor.getInstance();
        pMonitor.setMaxValue(numEdgesToRemove);
        for (int k = 0; k < numEdgesToRemove; k++) {
            BetweennessCentrality bc = new BetweennessCentrality(graphCopy, false);
            bc.setRemoveRankScoresOnFinalize(false);
            bc.evaluate();

            //  Determine the edge with highest betweenness score and remove it
            EdgeRanking highestBetweenness = (EdgeRanking) bc.getRankings().get(0);

            //  Output the edge info.
            mEdgesRemoved.add(highestBetweenness.edge.getEqualEdge(graphCopy));
            graphCopy.removeEdge(highestBetweenness.edge);

            //  After removal of so many edges, we have multiple weak components
            WeakComponentClusterer wcSearch = new WeakComponentClusterer();
            ClusterSet clusterSet = wcSearch.extract(graphCopy);
            if (currentNumModules != clusterSet.size()) {
                currentNumModules = clusterSet.size();

                // Calculate Network Modularity using *Original Graph*,
                // not the Newly Partitioned Graph!
                NetworkPartitionState state = calculateNetworkModularity(originalGraph,
                        clusterSet, k + 1);
                //  outputState(state);
                networkPartionStateList.add(state);
                if (state.getNetworkModularity() > maxModularityScore) {
                    maxModularityScore = state.getNetworkModularity();
                    optimalPartitionState = state;
                    optimalClusterSet = clusterSet;
                    optimalGraph = (Graph) graphCopy.copy();
                }
            }
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }
    }

    private NetworkPartitionState calculateNetworkModularity(Graph graph, ClusterSet clusterSet,
                                                             int numEdgesRemoved) {
        double networkModularityScore =
                NetworkModularity.calculateNetworkModularity(graph, clusterSet, false);
        NetworkPartitionState state = new NetworkPartitionState();
        state.setNumEdgesRemoved(numEdgesRemoved);
        state.setNumModules(clusterSet.size());
        state.setNetworkModularity(networkModularityScore);
        return state;
    }

    private void outputState(NetworkPartitionState state) {
        NumberFormat formatter = Formatter.getDecimalFormat();
        ProgressMonitor.getInstance().setCurrentMessage(state.getNumModules() + "\t"
                + state.getNumEdgesRemoved() + "\t"
                + formatter.format(state.getNetworkModularity()));
    }

    /**
     * Gets the Max Network Modularity Score.
     *
     * @return max network modularity score.
     */
    public NetworkPartitionState getOptimalPartitionState() {
        return optimalPartitionState;
    }

    /**
     * Gets the Cluster Set corresponding to the Max Modularity Score.
     *
     * @return Cluster Set.
     */
    public ClusterSet getOptimalClusterSet() {
        return optimalClusterSet;
    }

    /**
     * Gets the Graph corresponding to the Max Modularity Score.
     *
     * @return Graph Object.
     */
    public Graph getOptimalGraph() {
        return optimalGraph;
    }

    /**
     * Indicates if the initival graph is represented by a single connected component,
     * prior to partitioning.
     *
     * @return true or false.
     */
    public boolean isInitialStateRepresentedBySingleConnectedComponent() {
        return initialStateRepresentsSingleConnectedComponent;
    }
}

