package org.mskcc.netbox.test.query;

import junit.framework.TestCase;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.util.GlobalSession;

/**
 * Tests the Gene Query Class.
 *
 * @author Ethan Cerami.
 */
public class TestGeneQuery extends TestCase {
    private static final String AGAP2 = "AGAP2";
    private static final long AGAP2_ID = 116986;

    /**
     * Tests the Basic Gene Query.
     *
     * @throws Exception All Errors.
     */
    public final void testGeneQuery() throws Exception {

        SessionFactory sessionFactory = GlobalSession.getInstance().getSessionFactory();

        // Ask for a session using the JDBC information we've configured
        storeGene(sessionFactory);

        Session session = sessionFactory.openSession();
        try {
            //  Verify different mechanisms of gene lookup.
            Gene gene = GeneQuery.getGeneBySymbol(AGAP2);
            assertEquals(AGAP2, gene.getGeneSymbol());
            assertEquals(AGAP2_ID, gene.getEntrezGeneId());

            gene = GeneQuery.getGeneByEntrezGeneId(AGAP2_ID);
            assertEquals(AGAP2, gene.getGeneSymbol());
            assertEquals(AGAP2_ID, gene.getEntrezGeneId());

        } finally {
            // No matter what, close the session
            session.close();
        }

        // Clean up after ourselves
        sessionFactory.close();
    }

    /**
     * Stores a Sample Gene.
     *
     * @param sessionFactory Session Factor Object.
     */
    private void storeGene(SessionFactory sessionFactory) {
        Session session = sessionFactory.openSession();
        Transaction tx = null;
        //  Start with clean slate
        GeneQuery.deleteAllGenes();
        try {
            // Create some model and persist it
            tx = session.beginTransaction();

            //  Store AGAP2 Gene
            Gene gene = new Gene(AGAP2, AGAP2_ID);
            session.save(gene);

            // We're done; make our changes permanent
            tx.commit();
        } catch (Exception e) {
            if (tx != null) {
                // Something went wrong; discard all partial changes
                tx.rollback();
            }
            fail("Could not save gene");
        } finally {
            // No matter what, close the session
            session.close();
        }
    }
}
