package org.mskcc.netbox.graph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.log4j.Logger;

/**
 * Utility to Calculate Network Modularity (Faster).
 *
 * @author Ethan Cerami.
 */
public final class NetworkModularityFast {
    private int totalNumEdges;
    private HashMap<String, Integer> degreeMap = new HashMap<String, Integer>();
    private HashMap<String, ArrayList<MiniEdge>> edgeMap
            = new HashMap<String, ArrayList<MiniEdge>>();
    private HashMap<String, ModuleStats> qMap = new HashMap<String, ModuleStats>();
    private HashMap<String, ArrayList<String>> storedNodesInModule;
    private String sourceModuleId;
    private String targetModuleId;
    private String nodeId;
    private ModuleStats sourceModuleStats;
    private ModuleStats targetModuleStats;
    private static Logger logger = Logger.getLogger(NetworkModularityFast.class);


    /**
     * Constructor.
     * @param g Graph.
     */
    public NetworkModularityFast(Graph g) {
        this.totalNumEdges = g.numEdges();
        StringLabeller labeller = StringLabeller.getLabeller(g);

        //  Store Degree of all vertices in hash map
        //  Store Mini Edges for all vertices in another hash map
        Set<Vertex> vertexSet = g.getVertices();
        for (Vertex vertexA : vertexSet) {
            String a = GraphUtil.getVertexLabel(labeller, vertexA);
            int degree = vertexA.degree();
            degreeMap.put(a, degree);

            Set<Edge> edgeSet = vertexA.getIncidentEdges();
            ArrayList<MiniEdge> edgeList = new ArrayList<MiniEdge>();
            for (Edge edge : edgeSet) {
                Vertex vertexB = edge.getOpposite(vertexA);
                String b = GraphUtil.getVertexLabel(labeller, vertexB);
                MiniEdge miniEdge = new MiniEdge(a, b);
                edgeList.add(miniEdge);
            }
            edgeMap.put(a, edgeList);
        }
    }

    /**
     * Sets the Base Line.
     * @param moduleList    Module List.
     * @param moduleMap     Module Map.
     */
    public void setBaseLine(ArrayList<String> moduleList,
                            HashMap<String, String> moduleMap) {
        calculateNetworkModularity(moduleList, moduleMap, false, true);
    }

    /**
     * Calculates Network Modularity.
     * @param moduleList    Module List.
     * @param moduleMap     Module Map.
     * @param debug         Debug Flag.
     * @param store         Store Flag.
     * @return network modularity.
     */
    public double calculateNetworkModularity(ArrayList<String> moduleList,
                                             HashMap<String, String> moduleMap,
                                             boolean debug, boolean store) {

        if (store) {
            qMap = new HashMap<String, ModuleStats>();
        }

        //  Create an inverse look-up map
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

        if (store) {
            storedNodesInModule = nodesInModule;
        }

        double q = 0.0;

        //  Iterate through all modules
        for (String moduleId : moduleList) {
            HashSet<MiniEdge> edgesInModule = new HashSet<MiniEdge>();
            int d = 0;

            ArrayList<String> nodeList = nodesInModule.get(moduleId);
            if (nodeList != null) {
                //  Iterate through all genes in module
                for (String nodeA : nodeList) {

                    //  Now iterate through all adjacent vertices.
                    d += degreeMap.get(nodeA);
                    ArrayList<MiniEdge> edgeList = edgeMap.get(nodeA);
                    for (MiniEdge edge : edgeList) {
                        String nodeB = edge.getOpposite(nodeA);
                        if (nodeList.contains(nodeB)) {
                            edgesInModule.add(edge);
                        }
                    }
                }
            }

            ModuleStats moduleStats = new ModuleStats(totalNumEdges);
            moduleStats.setD(d);
            moduleStats.setEdgesInModule(edgesInModule);

            if (store) {
                qMap.put(moduleId, moduleStats);
            }
            q += moduleStats.getModularity();
        }
        return q;
    }

    /**
     * Updates the Network Modularity, Based on a Single Movement of a Node from an old module
     * to a new module.
     *
     * @param moduleList   List of Modules.
     * @param sourceModule Source module where the node is coming from.
     * @param targetModule Target module where the node is going to.
     * @param nodeA        Identifier of the node.
     * @return network modularity.
     */
    public double updateNetworkModularityOneNodeMove(ArrayList<String> moduleList,
         String sourceModule, String targetModule, String nodeA) {

        //  If the modules are the same, throw an exception.
        if (targetModule.equals(sourceModule)) {
            throw new IllegalArgumentException("modules must be different");
        }

        //  Store these arguments for the accept() method.
        this.sourceModuleId = sourceModule;
        this.targetModuleId = targetModule;
        this.nodeId = nodeA;

        //  Get All the Nodes in the Target Module
        ArrayList<String> nodeListTo = storedNodesInModule.get(targetModule);

        double q = 0.0;

        //  Get Topological info about the node that we just moved
        int degree = degreeMap.get(nodeA);
        ArrayList<MiniEdge> edgeList = edgeMap.get(nodeA);

        try {
            //  Update a Clone of the Source Module
            //  We work with Clones, just in case we decide not to accept this change,
            //  and need to roll back.
            sourceModuleStats = (ModuleStats) qMap.get(sourceModule).clone();
            sourceModuleStats.setD(sourceModuleStats.getD() - degree);
            HashSet<MiniEdge> miniEdgesSourceList = sourceModuleStats.getEdgesInModule();

            //  Update a Clone of the Target Module
            targetModuleStats = (ModuleStats) qMap.get(targetModule).clone();
            targetModuleStats.setD(targetModuleStats.getD() + degree);
            HashSet<MiniEdge> miniEdgesTargetList = targetModuleStats.getEdgesInModule();

            //  Update Both Clones of the Source and Target Modules
            for (MiniEdge edge : edgeList) {
                // If node A has any existing within module edges in the source module, remove them
                miniEdgesSourceList.remove(edge);

                //  If node A results in any new within modules in the target module, add them
                String nodeB = edge.getOpposite(nodeA);
                if (nodeListTo.contains(nodeB)) {
                    miniEdgesTargetList.add(edge);
                }
            }

            q += sourceModuleStats.getModularity();
            q += targetModuleStats.getModularity();

            //  Iterate through all other modules via quick look up
            for (String moduleId : moduleList) {
                if (!moduleId.equals(sourceModule) && !moduleId.equals(targetModule)) {
                    ModuleStats moduleStats = qMap.get(moduleId);
                    q += moduleStats.getModularity();
                }
            }
        } catch (CloneNotSupportedException e) {
            logger.warn(e);
        }
        return q;
    }

    /**
     * Accepts the Last Single Node Movement.
     */
    public void accept() {
        qMap.put(sourceModuleId, sourceModuleStats);
        qMap.put(targetModuleId, targetModuleStats);
        ArrayList<String> nodeListSource = storedNodesInModule.get(sourceModuleId);
        nodeListSource.remove(nodeId);
        ArrayList<String> nodeListTarget = storedNodesInModule.get(targetModuleId);
        nodeListTarget.add(nodeId);
    }
}

/**
 * Encapsulates Module Stats.
 *
 * @author Ethan Cerami.
 */
class ModuleStats {
    private int totalNumEdges;
    private HashSet<MiniEdge> edgesInModule;
    private int d;

    /**
     * Constructor.
     * @param n number of edges.
     */
    public ModuleStats(int n) {
        totalNumEdges = n;
    }

    /**
     * Gets All Edges within the Module.
     * @return HashSet of all MiniEdges.
     */
    public HashSet<MiniEdge> getEdgesInModule() {
        return edgesInModule;
    }

    /**
     * Sets All Edges within the Module.
     * @param edgesInM HashSet of all MiniEdges.
     */
    public void setEdgesInModule(HashSet<MiniEdge> edgesInM) {
        this.edgesInModule = edgesInM;
    }

    /**
     * Gets D.
     * @return d.
     */
    public int getD() {
        return d;
    }

    /**
     * Sets D.
     * @param dValue d.
     */
    public void setD(int dValue) {
        this.d = dValue;
    }

    /**
     * Gets Network Modularity.
     * @return Network Modularity.
     */
    public double getModularity() {
        int l = edgesInModule.size();
        return ((l / (double) totalNumEdges))
                - (Math.pow(d / ((double) 2 * totalNumEdges), 2));
    }

    /**
     * Overrides clone().
     * @return Cloned object.
     * @throws CloneNotSupportedException Clone Error.
     */
    protected Object clone() throws CloneNotSupportedException {
        ModuleStats clone = new ModuleStats(totalNumEdges);
        clone.setD(d);
        clone.setEdgesInModule((HashSet<MiniEdge>) edgesInModule.clone());
        return clone;
    }
}
