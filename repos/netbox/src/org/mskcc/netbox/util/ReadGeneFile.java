package org.mskcc.netbox.util;

import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.query.GeneQuery;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Utility Class for Reading in Text Files that Contain Gene Symbols.
 * These files can contain comments, beginning with #.
 *
 * @author Ethan Cerami.
 */
public final class ReadGeneFile {

    /**
     * Private Constructor to prevent instantiation.
     */
    private ReadGeneFile() {
    }

    /**
     * Reads in a Text File Containing Gene Symbols.
     *
     * @param geneSymbolMap A HashMap of Valid Human Gene Symbols.
     * @param file          File to Read.
     * @return ArrayList of Gene Symbols.
     * @throws IOException IO Error.
     */
    public static ArrayList<String> readGeneFile(HashMap<String, Gene> geneSymbolMap, File file)
            throws IOException {
        HashSet<String> invalidGeneSet = new HashSet<String>();
        ArrayList<String> geneList = new ArrayList<String>();
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = reader.readLine();
        while (line != null) {
            line = line.trim();
            if (!line.startsWith("#")) {

                //  First, assume it's an entrez gene Id.
                //  If this fails, assume it's a gene symbol.
                String ensemblGeneId = null;
                Gene currentGene = null;
                try {
                	
                    ensemblGeneId = line;
                    currentGene = GeneQuery.getGeneByEnsemblGeneId(ensemblGeneId);
                } catch (NumberFormatException e) {
                    currentGene = geneSymbolMap.get(line);
                }
                if (currentGene == null) {
                    invalidGeneSet.add(line);
                } else {
                    geneList.add(currentGene.getEnsemblGeneId());
                }
            }
            line = reader.readLine();
        }
        if (invalidGeneSet.size() > 0) {
        	System.err.println("File contains "+ invalidGeneSet.size() +" invalid gene symbols or "
                    + "invalid Entrez Gene Ids");
            StringBuffer buf = new StringBuffer("File contains"+ invalidGeneSet.size() +"invalid gene symbols or "
                    + "invalid Entrez Gene Ids:  ");
            for (String gene : invalidGeneSet) {
                buf.append(gene + " ");
            }
            buf.append("\nPlease check these genes or comment them out before trying again.");
            // throw new IllegalArgumentException(buf.toString());
        } else {
            if (geneList.size() == 0) {
                throw new IllegalArgumentException("Could not find any genes in file:  "
                        + file.getAbsolutePath());
            }
        }
        return geneList;
    }
}
