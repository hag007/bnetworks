package org.mskcc.netbox.util;

import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.model.GeneSet;
import org.mskcc.netbox.query.GeneQuery;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Reads In Gene Sets.
 *
 * @author Ethan Cerami
 */
public final class ReadGeneSets {
    private ArrayList<GeneSet> geneSetList = new ArrayList<GeneSet>();

    /**
     * Constructor.
     *
     * @param file  File Containing Gene Sets.
     *
     * @throws IOException  IO Error.
     */
    public ReadGeneSets(File file) throws IOException {
        FileReader reader = new FileReader(file);
        BufferedReader buf = new BufferedReader(reader);
        String line = buf.readLine();
        HashMap<String, Gene> geneMap = GeneQuery.getGeneMapBySymbol();
        while (line != null) {
            String [] parts = line.split("\t");
            String name = parts[0];
            ArrayList<Gene> geneList = new ArrayList<Gene>();
            for (int i = 2; i < parts.length; i++) {
                String geneSymbol = parts[i];
                Gene gene = geneMap.get(geneSymbol);
                if (gene != null) {
                    geneList.add(gene);
                }
            }
            GeneSet geneSet = new GeneSet();
            geneSet.setName(name);
            geneSet.setGeneList(geneList);
            geneSetList.add(geneSet);
            line = buf.readLine();
        }
    }

    /**
     * Gets ArrayList of All Gene Sets.
     * @return List of All Gene Sets.
     */
    public ArrayList<GeneSet> getGeneSetList() {
        return geneSetList;
    }
}
