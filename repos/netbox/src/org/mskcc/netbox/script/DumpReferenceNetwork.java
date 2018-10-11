package org.mskcc.netbox.script;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;

import java.io.IOException;
import java.util.Set;
import java.util.HashMap;

import org.mskcc.netbox.query.InteractionQuery;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.model.Gene;

/**
 * Dumps the Human Interaction Network.
 *
 * @author Ethan Cerami.
 */
public class DumpReferenceNetwork {

    /**
     * Private constructor to prevent instantiation.
     */
    private DumpReferenceNetwork() {
    }

    /**
     * Command Line Tool to Dump the Reference Network.
     *
     * @param args Command Line Arguments.
     * @throws java.io.IOException IO Error.
     */
    public static void main(String[] args) throws IOException, GraphCreationException {
        HashMap<String, Gene> geneMap = GeneQuery.getGeneMapBySymbol();
        Graph g = InteractionQuery.getGlobalGraph();
        Set<Vertex> vertexSet = g.getVertices();
        StringLabeller labeller = StringLabeller.getLabeller(g);
        for (Vertex vertex: vertexSet) {
            String geneSymbol = labeller.getLabel(vertex);
            Gene gene = geneMap.get(geneSymbol);
            System.out.println(gene.getEntrezGeneId());
        }
    }
}
