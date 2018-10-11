package org.mskcc.netbox.netcarto;

import edu.uci.ics.jung.algorithms.cluster.ClusterSet;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import org.apache.log4j.Logger;
import org.mskcc.netbox.graph.GraphUtil;
import org.mskcc.netbox.graph.NetworkModularity;
import org.mskcc.netbox.graph.NetworkModularityFast;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Implementation of the NetCarto Simulated Annealing Algorithm.
 * <p/>
 * Based on algorithm description in: Guimerˆ R, Nunes Amaral LA.
 * Functional cartography of complex metabolic networks.  Nature. 2005 Feb 24;433(7028):895-900.
 * http://www.ncbi.nlm.nih.gov/pubmed/15729348
 *
 * @author Ethan Cerami
 */
public final class NetCartoAnnealing {
    private static Logger logger = Logger.getLogger(NetCartoAnnealing.class);
    private double c = 0.995;
    private double f = 1;
    private Random randomGenerator = new Random(500);
    private Graph graph;
    private int currentId = 0;
    private Set<String> ignoreSet = new HashSet<String>();
    private NetworkModularityFast fastModularityUtil;
    private NetCartoState currentState;
    private double finalModularity;
    private ArrayList<String> finalModuleList;
    private HashMap<String, String> finalModuleMap;
    private HashMap<String, ArrayList<String>> globalModuleMap;

    /**
     * Constructor.
     *
     * @param originalGraph Original Graph.
     */
    public NetCartoAnnealing(Graph originalGraph) {
        this.graph = originalGraph;
        this.fastModularityUtil = new NetworkModularityFast(graph);
    }

    /**
     * Execute Simulated Annealing.
     */
    public void execute() {
        // Do not include disconnected component of size <= 4 in the simulated annealing
        WeakComponentClusterer wcSearch = new WeakComponentClusterer();
        ClusterSet clusterSet = wcSearch.extract(graph);
        StringLabeller labeller = StringLabeller.getLabeller(graph);

        finalModuleList = new ArrayList<String>();
        finalModuleMap = new HashMap<String, String>();
        int counter = 0;
        for (int i = 0; i < clusterSet.size(); i++) {
            Set<Vertex> set = clusterSet.getCluster(i);
            if (set.size() <= 4) {
                String moduleId = "" + counter;
                finalModuleList.add(moduleId);
                counter++;
                for (Vertex v : set) {
                    String nodeLabel = GraphUtil.getVertexLabel(labeller, v);
                    ignoreSet.add(nodeLabel);
                    finalModuleMap.put(nodeLabel, moduleId);
                }
            }
        }
        executeAnnealing();

        //  Assemble the final modules
        ArrayList<String> moduleList = currentState.getModuleList();
        HashMap<String, String> moduleMap = currentState.getModuleMap();
        for (String moduleId : moduleList) {
            String newModuleId = "" + counter;
            counter++;
            finalModuleList.add(newModuleId);

            for (String node : moduleMap.keySet()) {
                String moduleAssignmentId = moduleMap.get(node);
                if (moduleAssignmentId.equals(moduleId)) {
                    finalModuleMap.put(node, newModuleId);
                }
            }
        }

        finalModularity = NetworkModularity.calculateNetworkModularity(graph,
                finalModuleList, finalModuleMap, false);

        globalModuleMap = new HashMap<String, ArrayList<String>>();
        for (String moduleId : finalModuleList) {
            ArrayList<String> list = new ArrayList<String>();
            for (String node : finalModuleMap.keySet()) {
                String moduleAssignmentId = finalModuleMap.get(node);
                if (moduleAssignmentId.equals(moduleId)) {
                    list.add(node);
                }
            }
            globalModuleMap.put(moduleId, list);
        }
    }

    /**
     * Gets the Final Module List Identified By SA.
     * @return List of Modules.
     */
    public ArrayList<String> getFinalModuleList() {
        return finalModuleList;
    }

    /**
     * Gets the Global Module Map.
     * @return Module Map.
     */
    public HashMap<String, ArrayList<String>> getGlobalModuleMap() {
        return globalModuleMap;
    }

    /**
     * Gets the Final Network Modularity.
     * @return Network Modularity.
     */
    public double getFinalModularity() {
        return finalModularity;
    }

    private void executeAnnealing() {
        int counter = 0;
        boolean continueAnnealing = true;

        currentState = initState();

        //  Set start temperature
        double temperature = 2.0 / graph.numVertices();

        //  Execute Simulated Annealing
        while (continueAnnealing) {
            outputState(temperature, currentState);

            //  Update the State
            NetCartoState updatedState = updateState(temperature, currentState);

            //  According to RGraph / NetCarto documentation, "program will stop when the
            //  modularity remains unchanged during 25 different temperatures."
            double currentModularity = currentState.getNetworkModularity();
            double updatedModularity = updatedState.getNetworkModularity();
            if (Math.abs(currentModularity - updatedModularity) < 1E-6) {
                counter++;
            } else {
                counter = 0;
            }

            if (counter == 100) {
                continueAnnealing = false;
            }

            //  Cool down
            currentState = updatedState;
            temperature = temperature * c;
        }
    }

    /**
     * Initial State is Defined by N Modules.  All nodes in separate modules.
     */
    private NetCartoState initState() {
        NetCartoState state = new NetCartoState();
        ArrayList<String> moduleList = new ArrayList<String>();
        HashMap<String, String> moduleMap = new HashMap<String, String>();
        StringLabeller labeller = StringLabeller.getLabeller(graph);
        Set<Vertex> vertexSet = graph.getVertices();
        for (Vertex vertex : vertexSet) {
            String label = labeller.getLabel(vertex);
            if (!ignoreSet.contains(label)) {
                String moduleId = getNextId();
                moduleList.add(moduleId);
                moduleMap.put(label, moduleId);
            }
        }
        double networkModularity = fastModularityUtil.calculateNetworkModularity(moduleList,
                moduleMap, false, false);
        state.setModuleList(moduleList);
        state.setModuleMap(moduleMap);
        state.setNetworkModularity(networkModularity);
        return state;
    }

    private void outputState(double temperature, NetCartoState state) {
        logger.warn("1/T:  " + (1 / temperature)
                + "\tTemp:  " + temperature + "\tModularity:  "
                + (state.getNetworkModularity() * -1) + "\t"
                + "Num_Modules:  " + state.getModuleList().size());
    }

    private NetCartoState updateState(double temperature, NetCartoState state) {
        ArrayList<String> moduleList = (ArrayList<String>) state.getModuleList().clone();

        //  From original algorithm description:
        //  "we propose n1 = fS^2 individual node movements from one module to another, where S is
        //  the number of nodes in the network."
        int s = graph.numVertices();
        int n1 = (int) (f * Math.pow(s, 2));

        fastModularityUtil.setBaseLine(moduleList, state.getModuleMap());
        HashMap<String, String> moduleMap = state.getModuleMap();

        ArrayList<String> nodeList = new ArrayList<String>();
        nodeList.addAll(moduleMap.keySet());

        if (moduleList.size() > 1) {
            for (int i = 0; i < n1; i++) {
                //  At each iteration, randomly pick a vertex from G
                int randomIndex = randomGenerator.nextInt(nodeList.size());
                String nodeLabel = nodeList.get(randomIndex);
                //  and, assign it to a random module
                String originalModuleId = moduleMap.get(nodeLabel);
                int randomModuleIndex = randomGenerator.nextInt(moduleList.size());
                String randomModuleId = moduleList.get(randomModuleIndex);
                if (!originalModuleId.equals(randomModuleId)) {
                    logger.info("Randomly re-assigning " + nodeLabel + " from module: "
                            + originalModuleId + " to module:  " + randomModuleId + ".");

                    //  Create updated State, and Decide to Accept / Reject it
                    double updatedQ = fastModularityUtil.updateNetworkModularityOneNodeMove(
                            state.getModuleList(), originalModuleId, randomModuleId,
                            nodeLabel);
                    double dE = updatedQ - state.getNetworkModularity();
                    double p = randomGenerator.nextDouble();
                    if (dE > 0.0 && Math.exp((-1 * dE) / temperature) < p) {
                        logger.info("Rejecting updated state");
                    } else {
                        fastModularityUtil.accept();
                        state.setNetworkModularity(updatedQ);
                        moduleMap.put(nodeLabel, randomModuleId);
                    }
                }
            }
        }

        //  Before proceeding, prune the module list of empty modules.
        pruneModuleList(state);

        //  From original algorithm description:
        //  "We also propose nc = fS collective movements, which involve either merging two modules
        //  or splitting a module.
        int nc = (int) (f * s);
        for (int i = 0; i < nc; i++) {
            moduleList = (ArrayList<String>) state.getModuleList().clone();
            moduleMap = (HashMap<String, String>)
                    state.getModuleMap().clone();

            int action = randomGenerator.nextInt(2);
            if (action == 0) {
                //  Merge Two Modules
                mergeTwoModules(moduleList, moduleMap);
            } else {
                //  Split a Module
                splitModule(moduleList, moduleMap);
            }
            //  Create updated State, and Decide to Accept / Reject it
            NetCartoState updatedState = createUpdatedState(moduleList, moduleMap);
            if (acceptUpdatedState(state, updatedState, temperature)) {
                state = updatedState;
                //  Again, before proceeding, prune the module list of empty modules.
                pruneModuleList(state);
            }
        }

        return state;
    }

    private void pruneModuleList(NetCartoState state) {
        ArrayList<String> moduleList;
        moduleList = state.getModuleList();
        HashMap<String, String> moduleMap = state.getModuleMap();

        //  Before proceeding, prune the module list of any empty modules.
        //  To do so, first, create an inverse look-up map
        HashMap<String, ArrayList<String>> nodesInModule =
                new HashMap<String, ArrayList<String>>();
        for (String gene : moduleMap.keySet()) {
            String moduleId = moduleMap.get(gene);
            ArrayList<String> nodeList = nodesInModule.get(moduleId);
            if (nodeList == null) {
                nodeList = new ArrayList<String>();
                nodeList.add(gene);
                nodesInModule.put(moduleId, nodeList);
            } else {
                nodeList.add(gene);
            }
        }
        ArrayList<String> markedForRemoval = new ArrayList<String>();
        for (String moduleId : moduleList) {
            ArrayList nodeList = nodesInModule.get(moduleId);
            if (nodeList == null || nodeList.size() == 0) {
                markedForRemoval.add(moduleId);
            }
        }
        for (String moduleId : markedForRemoval) {
            moduleList.remove(moduleId);
        }
    }

    private NetCartoState createUpdatedState(ArrayList<String> moduleList,
                                             HashMap<String, String> moduleMap) {
        NetCartoState updatedState = new NetCartoState(moduleList, moduleMap);
        updatedState.setNetworkModularity(NetworkModularity.calculateNetworkModularity(graph,
                moduleList, moduleMap, false));
        return updatedState;
    }

    private boolean acceptUpdatedState(NetCartoState state, NetCartoState updatedState,
                                       double temperature) {
        double dE = updatedState.getNetworkModularity()
                - state.getNetworkModularity();
        double p = randomGenerator.nextDouble();
        logger.info("dE = " + dE + ", random:  " + p);
        if (dE > 0.0 && Math.exp((-1 * dE) / temperature) < p) {
            //  Reject
            logger.info("Rejecting updated state");
            return false;
        } else {
            //  Accept
            logger.info("Accepting updated state");
            return true;
        }
    }

    private void mergeTwoModules(ArrayList<String> moduleList,
                                 HashMap<String, String> moduleMap) {
        if (moduleList.size() <= 1) {
            return;
        }
        //  Randomly pick two different modules
        int moduleAIndex = randomGenerator.nextInt(moduleList.size());
        int moduleBIndex = randomGenerator.nextInt(moduleList.size());
        while (moduleAIndex == moduleBIndex) {
            moduleBIndex = randomGenerator.nextInt(moduleList.size());
        }

        String moduleAId = moduleList.get(moduleAIndex);
        String moduleBId = moduleList.get(moduleBIndex);
        logger.info("Merging modules:  " + moduleAId + " and " + moduleBId + " together.");

        //  Remove module B from list
        moduleList.remove(moduleBIndex);

        // Update all nodes that previously pointed to B, and now point them to A
        for (String node : moduleMap.keySet()) {
            String moduleId = moduleMap.get(node);
            if (moduleId.equals(moduleBId)) {
                moduleMap.put(node, moduleAId);
            }
        }
    }

    private void splitModule(ArrayList<String> moduleList,
                             HashMap<String, String> moduleMap) {

        //  Return if we have reached the upper limit of modules
        if (graph.numVertices() == moduleList.size()) {
            return;
        }

        //  Randomly pick a module to split
        int moduleIdIndex = randomGenerator.nextInt(moduleList.size());
        String moduleId = moduleList.get(moduleIdIndex);

        //  Get all nodes in that module
        ArrayList<String> nodeList = new ArrayList<String>();
        for (String node : moduleMap.keySet()) {
            String currentModuleId = moduleMap.get(node);
            if (moduleId.equals(currentModuleId)) {
                nodeList.add(node);
            }
        }

        //  Create a new module
        String newModuleId = getNextId();
        moduleList.add(newModuleId);
        logger.info("Splitting module " + moduleId + " into two.");

        //  Shuffle all the nodes
        Collections.shuffle(nodeList, randomGenerator);

        //  Pick a random split index
        int splitIndex = randomGenerator.nextInt(nodeList.size());

        //  Randomly assign genes between the modules
        for (int i = 0; i < splitIndex; i++) {
            String nodeLabel = nodeList.get(i);
            moduleMap.put(nodeLabel, moduleId);
        }
        for (int i = splitIndex; i < nodeList.size(); i++) {
            String nodeLabel = nodeList.get(i);
            moduleMap.put(nodeLabel, newModuleId);
        }
    }

    private String getNextId() {
        String id = "MODULE_" + currentId;
        currentId++;
        return id;
    }
}
