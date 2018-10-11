package org.mskcc.netbox.model;

import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;

/**
 * Encapsulate Basic Network Stats.
 *
 * @author Ethan Cerami.
 */
@Entity
@NamedQueries({
        @NamedQuery(name = "org.mskcc.netbox.getNetworkStats",
                query = "from NetworkStats as networkstats"),
        @NamedQuery(name = "org.mskcc.netbox.deleteAllNetworkStats",
                query = "delete from NetworkStats")
})
public final class NetworkStats {
    /**
     * The database key associated with this gene once it is
     * persisted to a database.
     */
    @Id
    private long networkStatsId;

    private long numGenes;
    private long numEdges;

    /**
     * Gets Number of Genes in Database.
     *
     * @return Number of Genes.
     */
    public long getNumGenes() {
        return numGenes;
    }

    /**
     * Sets Number of Genes in Database.
     *
     * @param n Number of Genes.
     */
    public void setNumGenes(long n) {
        this.numGenes = n;
    }

    /**
     * Gets Number of Edges in Database.
     *
     * @return number of edges.
     */
    public long getNumEdges() {
        return numEdges;
    }

    /**
     * Sets Number of Edges in Database.
     *
     * @param n number of edges.
     */
    public void setNumEdges(long n) {
        this.numEdges = n;
    }
}
