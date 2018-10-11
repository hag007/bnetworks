package org.mskcc.netbox.algorithm;

import org.mskcc.netbox.graph.Module;

import java.util.ArrayList;

/**
 * Interface for Module Detector.
 *
 * @author Ethan Cerami.
 */
public interface ModuleDetector {

    /**
     * Gets all Modules Detected.
     * @return Module List.
     */
    ArrayList<Module> getModules();
}
