package org.mskcc.netbox.netcarto;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Encapsulates information regarding a Graph Partitioning State.
 *
 * @author Ethan Cerami.
 */
public final class NetCartoState {
    private ArrayList<String> moduleList;
    private HashMap<String, String> moduleMap;
    private double networkModularity;

    /**
     * Constructor.
     */
    public NetCartoState() {
    }

    /**
     * Constructor.
     * @param mList Module List.
     * @param mMap  Module Map.
     */
    public NetCartoState(ArrayList<String> mList, HashMap<String, String> mMap) {
        this.moduleList = mList;
        this.moduleMap = mMap;
    }

    /**
     * Gets the Module List.
     * @return module list.
     */
    public ArrayList<String> getModuleList() {
        return moduleList;
    }

    /**
     * Sets the Module List.
     * @param ml module list.
     */
    public void setModuleList(ArrayList<String> ml) {
        this.moduleList = ml;
    }

    /**
     * Gets the Module Map.
     * @return Gets the Module Map.
     */
    public HashMap<String, String> getModuleMap() {
        return moduleMap;
    }

    /**
     * Gets the Module Map.
     * @param mm Module Map.
     */
    public void setModuleMap(HashMap<String, String> mm) {
        this.moduleMap = mm;
    }

    /**
     * Gets the Network Modularity.
     * @return network modularity.
     */
    public double getNetworkModularity() {
        return networkModularity * -1;
    }

    /**
     * Sets the Network Modularity.
     * @param q Network Modularity.
     */
    public void setNetworkModularity(double q) {
        this.networkModularity = q;
    }
}
