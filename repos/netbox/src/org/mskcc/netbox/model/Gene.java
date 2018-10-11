package org.mskcc.netbox.model;

import org.hibernate.annotations.Cache;
import org.hibernate.annotations.CacheConcurrencyStrategy;
import org.hibernate.annotations.Index;

import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import java.io.Serializable;

/**
 * Gene Object.
 */
@Entity
@Cache(usage = CacheConcurrencyStrategy.READ_WRITE)
@NamedQueries({
	@NamedQuery(name = "org.mskcc.netbox.getGeneBySymbol",
			query = "from Gene as gene where gene.geneSymbol = :symbol"),
	@NamedQuery(name = "org.mskcc.netbox.getGeneMapBySymbol",
	query = "from Gene as gene"),
	@NamedQuery(name = "org.mskcc.netbox.deleteAllGenes",
	query = "delete from Gene"),
	@NamedQuery(name = "org.mskcc.netbox.getGeneByEntrezGeneId",
	query = "from Gene as gene where gene.entrezGeneId = :geneId")
})
public final class Gene implements Serializable {

	/**
	 * The database key associated with this gene once it is
	 * persisted to a database.
	 */
	@Id
	private long entrezGeneId;

	@Index(name = "geneSymbolIndex")
	private String geneSymbol;
	private String refSeqId;

	private String ensemblGeneId;

	/**
	 * Constructor.
	 *
	 * @param s  Gene Symbol.
	 * @param id Entrez Gene ID.
	 */
	public Gene(String s, long id) {
		this.geneSymbol = s;
		this.entrezGeneId = id;
	}

	/**
	 * Constructor.
	 *
	 * @param s  Gene Symbol.
	 * @param id Entrez Gene ID.
	 */
	public Gene(String s, String id) {
		if (id ==null)
		{
			id="";
		}

		this.geneSymbol = s;
		this.ensemblGeneId = id;
	}

	/**
	 * No-Arg Constructor.
	 */
	public Gene() {
	}

	/**
	 * Gets the Official Gene Symbol.
	 *
	 * @return Official Gene Symbol.
	 */
	public String getGeneSymbol() {
		return geneSymbol;
	}

	/**
	 * Sets the Official Gene Symbol.
	 *
	 * @param s Gene Symbol.
	 */
	public void setGeneSymbol(String s) {
		this.geneSymbol = s;
	}

	/**
	 * Gets the Entrez Gene ID.
	 *
	 * @return Entrez Gene ID.
	 */
	public long getEntrezGeneId() {
		return entrezGeneId;
	}

	/**
	 * Sets the Entrez Gene ID.
	 *
	 * @param id Entrez Gene ID.
	 */
	public void setEntrezGeneId(long id) {
		this.entrezGeneId = id;
	}

	/**
	 * Gets the Entrez Gene ID.
	 *
	 * @return Entrez Gene ID.
	 */
	public String getEnsemblGeneId() {
		return ensemblGeneId;
	}

	/**
	 * Sets the Entrez Gene ID.
	 *
	 * @param id Entrez Gene ID.
	 */
	public void setEnsemblGeneId(String id) {
		this.ensemblGeneId= id;
	}

	
	/**
	 * Gets the RefSeq Accession ID.
	 *
	 * @return RefSeq Accession ID.
	 */
	public String getRefSeqId() {
		return refSeqId;
	}

	/**
	 * Sets the RefSeq Accession ID.
	 *
	 * @param id RefSeq Accession ID.
	 */
	public void setRefSeqId(String id) {
		this.refSeqId = id;
	}

	@Override
	/**
	 * toString() override.
	 */
	public String toString() {
		return ensemblGeneId;
	}
}
