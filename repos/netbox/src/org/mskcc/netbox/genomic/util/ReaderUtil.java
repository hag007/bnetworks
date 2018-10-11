package org.mskcc.netbox.genomic.util;

import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.util.ProgressMonitor;

/**
 * Read Utility Class.
 *
 * @author Ethan Cerami.
 */
public final class ReaderUtil {

    /**
     * Private Constructor to prevent instantiation.
     */
    private ReaderUtil() {
    }

    /**
     * Given an Entrez Gene ID, Returns the Corresponding Gene Symbol.
     * However, when validation of genes is turned off, this method simply returns the
     * entrez gene ID.
     *
     * @param entrezGeneId  Entrez Gene ID.
     * @param validateGenes Validate Gene ID.
     * @param pMonitor      Progress Monitor to Record Warning Messages.
     * @return gene symbol.
     */
    public static String getGeneId(String entrezGeneId, boolean validateGenes,
                                   ProgressMonitor pMonitor) {
        Gene gene;
        String id = null;
        if (validateGenes) {
            try {
                gene = GeneQuery.getGeneByEntrezGeneId(Long.parseLong(entrezGeneId));
            } catch (NumberFormatException e) {
                pMonitor.logWarning("Cannot parse Entrez Gene ID:  " + id);
                return null;
            }
            if (gene == null) {
                pMonitor.logWarning("Do not know gene with Entrez Gene ID:  " + entrezGeneId);
                return null;
            } else {
                return gene.getGeneSymbol();
            }
        } else {
            return entrezGeneId;
        }
    }

}
