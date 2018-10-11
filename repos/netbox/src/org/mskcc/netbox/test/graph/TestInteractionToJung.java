package org.mskcc.netbox.test.graph;

import edu.uci.ics.jung.graph.Graph;
import junit.framework.TestCase;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.InteractionToJung;
import org.mskcc.netbox.model.Interaction;

import java.util.ArrayList;

/**
 * Tests the InteractionToJung Class.
 *
 * @author Ethan Cerami.
 */
public class TestInteractionToJung extends TestCase {

    /**
     * Test the InteractionToJung Class.
     *
     * @throws GraphCreationException Graph Creation Error.
     */
    public final void testJungInteractionUtil() throws GraphCreationException {
        ArrayList<Interaction> interactionList = new ArrayList<Interaction>();

        Interaction interaction = new Interaction();
        interaction.setGeneA("A");
        interaction.setGeneB("B");
        interactionList.add(interaction);

        interaction = new Interaction();
        interaction.setGeneA("B");
        interaction.setGeneB("A");
        interactionList.add(interaction);

        Graph g = InteractionToJung.createGraph(interactionList);
        assertEquals(2, g.getVertices().size());
        assertEquals(1, g.getEdges().size());
    }
}
