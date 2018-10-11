package org.mskcc.netbox.script;

import org.hibernate.Session;
import org.hibernate.Transaction;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.query.InteractionQuery;
import org.mskcc.netbox.query.NetworkStatsQuery;
import org.mskcc.netbox.util.GlobalSession;

import java.io.IOException;

/**
 * Command Line Tool to Clear the Database of All Content.
 *
 * @author Ethan Cerami.
 */
public final class ClearDatabase {

    /**
     * Private Constructor to prevent instantiation.
     */
    private ClearDatabase() {
    }

    /**
     * Command Line Tool to Clear the Database of All Content.
     *
     * @param args Command Line Arguments.  None Expected.
     * @throws IOException IO Error.
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Initializing Database.");
        Session session = GlobalSession.getInstance().getSession();
        Transaction tx = session.beginTransaction();
        System.out.println("Deleting all Genes.");
        GeneQuery.deleteAllGenes();
        System.out.println("Deleting all Interactions.");
        InteractionQuery.deleteAllInteractions();
        System.out.println("Deleting all Network Stats.");
        NetworkStatsQuery.deleteAllNetworkStats();
        tx.commit();
        GlobalSession.getInstance().closeAll();
    }
}
