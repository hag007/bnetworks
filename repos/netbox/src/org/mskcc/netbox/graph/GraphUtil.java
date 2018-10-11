package org.mskcc.netbox.graph;

import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;

/**
 * Misc. Graph Utils.
 *
 * @author Ethan Cerami.
 */
public final class GraphUtil {

    /**
     * Private constructor, prevents instantiation.
     */
    private GraphUtil() {
    }

    /**
     * Gets the Label of the Specified Vertex.
     *
     * @param labeller StringLabeller Object.
     * @param vertex   Vertex Object.
     * @return Vertex Label.
     */
    public static String getVertexLabel(StringLabeller labeller, Vertex vertex) {
        String gene = labeller.getLabel(vertex);
        if (gene == null) {
            gene = (String) vertex.getUserDatum(InteractionToJung.GENE_KEY);
        }
        return gene;
    }
}
