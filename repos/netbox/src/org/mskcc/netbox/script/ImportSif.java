package org.mskcc.netbox.script;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.model.Interaction;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.util.GlobalSession;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Command Line Tool to Import Interactions in SIF Format.
 *
 * @author Ethan Cerami.
 */
public final class ImportSif {
    private static final int NUM_INTERACTIONS_PER_BATCH = 1000;

    /**
     * Private constructor to prevent instantiation.
     */
    private ImportSif() {
    }

    /**
     * Command Line Tool to Import Interactions in SIF Format.
     *
     * @param args Command Line Arguments.
     * @throws IOException IO Error.
     */
    public static void main(String[] args) throws IOException {
        File file = new File(args[0]);
        System.out.println("Reading in File:  " + file.getAbsolutePath());
        System.out.println("Initializing Database...");
        Session session = GlobalSession.getInstance().getSession();
        SessionFactory sessionFactory = GlobalSession.getInstance().getSessionFactory();
        Transaction tx = session.beginTransaction();
        HashMap<String, Gene> geneMap = GeneQuery.getGeneMapBySymbol();

        String dataSource = args[1];

        int pass = 0;
        int fail = 0;
        int counter = 0;
        HashSet<String> invalidGeneSet = new HashSet<String>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line = reader.readLine();
            while (line != null) {
                String[] parts = line.split("\\s");
                String geneSymbolA = parts[0];
                String interactionType = parts[1];
                String geneSymbolB = parts[2];
                if (geneSymbolA != null && geneSymbolB != null) {
                    if (geneSymbolA.length() > 0 && geneSymbolB.length() > 0) {
                        if (!geneSymbolA.equals(geneSymbolB)) {
                            Gene geneA = geneMap.get(geneSymbolA);
                            Gene geneB = geneMap.get(geneSymbolB);
                            if (geneA != null && geneB != null) {
                                Interaction interaction = createInteraction(dataSource,
                                        geneSymbolA, geneSymbolB,
                                        interactionType);
                                session.save(interaction);
                                pass++;
                            } else {
                                if (geneA == null) {
                                    invalidGeneSet.add(geneSymbolA);
                                }
                                if (geneB == null) {
                                    invalidGeneSet.add(geneSymbolB);
                                }
                                fail++;
                            }
                        }
                    }
                }
                line = reader.readLine();
                counter++;
                if ((counter % NUM_INTERACTIONS_PER_BATCH) == 0) {
                    tx.commit();
                    session.close();

                    System.out.println("Interactions processed:  " + counter);
                    session = sessionFactory.openSession();
                    tx = session.beginTransaction();
                }
            }
        } finally {
            if (tx.isActive()) {
                tx.commit();
            }
            if (session.isOpen()) {
                session.close();
            }
            sessionFactory.close();
        }
        System.out.println("Number of Interactions Stored:  " + pass);
        System.out.println("Number of Interactions that could not be stored, due to "
                + "invalid gene symbols:  " + fail);
        System.out.println("Total number of invalid gene symbols:  " + invalidGeneSet.size());
    }

    /**
     * Creates a New Interaction Object.
     *
     * @param dataSource      Data Source.
     * @param geneSymbolA     Symbol for Gene A.
     * @param geneSymbolB     Symbol for Gene B.
     * @param interactionType Interaction Type.
     * @return Interaction Object.
     */
    private static Interaction createInteraction(String dataSource, String geneSymbolA,
                                                 String geneSymbolB, String interactionType) {
        Interaction interaction = new Interaction();
        interaction.setGeneA(geneSymbolA);
        interaction.setGeneB(geneSymbolB);
        interaction.setSource(dataSource);
        interaction.setInteractionType(interactionType);
        return interaction;
    }
}
