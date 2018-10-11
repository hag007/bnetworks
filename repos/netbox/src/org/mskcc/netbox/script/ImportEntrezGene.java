package org.mskcc.netbox.script;

import org.hibernate.Session;
import org.hibernate.Transaction;
import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.util.GlobalSession;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * Command Line Tool to Import Entrez Gene Information.
 *
 * @author Ethan Cerami.
 */
public final class ImportEntrezGene {

    /**
     * Private Constructor to prevent instantiation.
     */
    private ImportEntrezGene() {
    }

    /**
     * Command Line Tool to Import Entrez Gene Information.
     *
     * @param args Command Line Arguments.  None Expected.
     * @throws IOException IO Error.
     */
    public static void main(String[] args) throws IOException {
        File file = new File(args[0]);
        System.out.println("Reading in File:  " + file.getAbsolutePath());
        System.out.println("Initializing Database...");

        Session session = GlobalSession.getInstance().getSession();
        Transaction tx = session.beginTransaction();
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = reader.readLine();
        int counter = 0;
        while (line != null) {
            String[] parts = line.split("\t");
            String taxId = parts[0];
            String entrezGeneId = parts[1];
            String geneSymbol = parts[2];
            if (taxId.equals("9606")) {
                Gene gene = new Gene(geneSymbol, Long.parseLong(entrezGeneId));
                session.save(gene);
                counter++;
            }
            line = reader.readLine();
        }
        tx.commit();
        System.out.println("Total number of genes saved:  " + counter);
        System.out.println("Verifying that genes were loaded...");
        HashMap<String, Gene> geneMap = GeneQuery.getGeneMapBySymbol();
        System.out.println("Total number of genes retrieved:  " + geneMap.size());
        GlobalSession.getInstance().closeAll();
    }
}
