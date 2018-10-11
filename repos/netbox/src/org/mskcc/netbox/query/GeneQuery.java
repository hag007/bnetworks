package org.mskcc.netbox.query;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

import org.apache.log4j.Logger;
import org.hibernate.ObjectNotFoundException;
import org.hibernate.Query;
import org.hibernate.Session;
import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.util.GlobalSession;

/**
 * Gene Queries.
 *
 * @author Ethan Cerami.
 */
public final class GeneQuery {
	private static Logger logger = Logger.getLogger(GeneQuery.class);
	private static HashMap<String, Gene> geneMapByGeneSymbol = new HashMap<String, Gene>();
	private static HashMap<String, String> ensembl2gs = new HashMap<String, String>();
	private static HashMap<String, String> gs2ensembl = new HashMap<String, String>();
	private static HashMap<String, String> entrez2ensembl = new HashMap<String, String>();
	private static HashMap<String, String> ensembl2entrez = new HashMap<String, String>();

	static {

		try {
			String fname = "/home/hag007/bnet/dictionaries/ensembl2gene_symbol.txt";
			FileReader fr = new FileReader(fname);
			BufferedReader br = new BufferedReader(fr);
			String line=br.readLine();
			while((line = br.readLine()) != null) {
				String[] ls = line.trim().split("\t");
				if (ls.length == 2) {
					ensembl2gs.put(ls[0].split("\\.")[0], ls[1]);
					gs2ensembl.put(ls[1], ls[0].split("\\.")[0]);
				}
			}

			fname = "/home/hag007/bnet/dictionaries/ensembl2entrez.txt";
			fr = new FileReader(fname);
			br = new BufferedReader(fr);

			while((line = br.readLine()) != null) {
				String[] ls = line.trim().split("\t");
				if (ls.length == 2) {
					ensembl2entrez.put(ls[0], ls[1]);
					entrez2ensembl.put(ls[1], ls[0]);
				}
			}	

		}
		catch(IOException e){
			System.out.println("error while loading dictionaries: "+ e.getMessage());
		}



	}


	/**
	 * Private Constructor to Prevent Instantiation.
	 */
	private GeneQuery() {
	}

	/**
	 * Look up Gene by Gene Symbol.
	 *
	 * @param geneSymbol Gene Symbol.
	 * @return Gene Object or Null if Not Found.
	 */
	public static Gene getGeneBySymbol(String geneSymbol) {
		Session session = GlobalSession.getInstance().getSession();
		Query query = session.getNamedQuery("org.mskcc.netbox.getGeneBySymbol");
		query.setString("symbol", geneSymbol);
		Gene found = (Gene) query.uniqueResult();
		if (found != null) {
			return found;
		}
		return found;
	}

	/**
	 * Look up Gene by Entrez Gene ID.
	 *
	 * @param entrezGeneId Entrez Gene ID.
	 * @return Gene Object or Null if Not Found.
	 */
	public static Gene getGeneByEntrezGeneId(long entrezGeneId) {
		Session session = GlobalSession.getInstance().getSession();
		try {
			return (Gene) session.load(Gene.class, entrezGeneId);
		} catch (ObjectNotFoundException e) {
			return null;
		}
	}

	/**
	 * Look up Gene by Entrez Gene ID.
	 *
	 * @param entrezGeneId Entrez Gene ID.
	 * @return Gene Object or Null if Not Found.
	 */
	public static Gene getGeneByEnsemblGeneId(String ensemblGeneId) {
		Session session = GlobalSession.getInstance().getSession();
		try {
			return (Gene) new Gene(ensembl2gs.get(ensemblGeneId), ensemblGeneId);
		} catch (ObjectNotFoundException e) {
			return null;
		}
	}

	/**
	 * Gets a HashMap of all Genes Indexed by Gene Symbol.
	 *
	 * @return HashMap of Gene Objects indexed by Gene Symbol.
	 */
	public static HashMap<String, Gene> getGeneMapBySymbol() {
		lookUp();
		return geneMapByGeneSymbol;
	}

	/**
	 * Delete all Genes in the Database.
	 */
	public static void deleteAllGenes() {
		Session session = GlobalSession.getInstance().getSession();
		Query query = session.getNamedQuery("org.mskcc.netbox.deleteAllGenes");
		query.executeUpdate();
	}

	/**
	 * Internal Lookup Used to Populate the HashMaps.
	 */
	private static void lookUp() {
		if (geneMapByGeneSymbol.size() == 0) {
			logger.info("Initialize Gene HashMaps");
			Session session = GlobalSession.getInstance().getSession();
			Query query = session.getNamedQuery("org.mskcc.netbox.getGeneMapBySymbol");
			Iterator iterator = query.iterate();
			while (iterator.hasNext()) {
				Gene gene = (Gene) iterator.next();
				geneMapByGeneSymbol.put(gene.getGeneSymbol(), gene);
			}
			logger.info("Gene HashMaps have a size of:  " + geneMapByGeneSymbol.size());
		}
	}
}
