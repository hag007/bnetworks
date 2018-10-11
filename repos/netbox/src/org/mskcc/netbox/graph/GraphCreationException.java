package org.mskcc.netbox.graph;

/**
 * Graph Creation Error.
 *
 * @author Ethan Cerami.
 */
public class GraphCreationException extends Exception {

    /**
     * Constructor.
     *
     * @param e Exception Object.
     */
    public GraphCreationException(Exception e) {
        super(e);
    }
}
