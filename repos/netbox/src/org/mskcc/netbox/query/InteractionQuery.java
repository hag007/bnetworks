package org.mskcc.netbox.query;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.hibernate.Query;
import org.hibernate.Session;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.InteractionToJung;
import org.mskcc.netbox.model.Interaction;
import org.mskcc.netbox.util.GlobalSession;
import org.mskcc.netbox.util.ProgressMonitor;

import edu.uci.ics.jung.graph.Graph;

/**
 * Interaction Queries.
 *
 * @author Ethan Cerami.
 */
public final class InteractionQuery {

	private static HashMap<String, ArrayList<Interaction>> interactionMap;
	private static ArrayList<Interaction> globalInteractionList;
	private static Graph globalGraph;



	static {
		globalInteractionList= new ArrayList<>(); 
				
		try {
		String fname = "/home/hag007/bnet/networks/dip.sif";
		FileReader fr = new FileReader(fname);
		BufferedReader br = new BufferedReader(fr);
		String line;
		while((line = br.readLine()) != null) {
			String[] ls = line.trim().split("\t");
			if (ls.length == 3) {
				Interaction inter = new Interaction();
				inter.setGeneA(ls[0]);
				inter.setGeneB(ls[2]);
				inter.setInteractionType(ls[1]);
				inter.setSource("dip");
				globalInteractionList.add(inter);
			}
		}
		}
		catch(IOException e)
		{
			System.out.println("error reading edges file");
		}


	}



	/**
	 * Private Constructor to prevent instatiation.
	 */
	private InteractionQuery() {
	}

	/**
	 * Deletes all Interactions in the Database.
	 */
	public static void deleteAllInteractions() {
		Session session = GlobalSession.getInstance().getSession();
		Query query = session.getNamedQuery("org.mskcc.netbox.deleteAllInteractions");
		query.executeUpdate();
	}

	/**
	 * Gets a HashMap of All Interactions in the Database.
	 *
	 * @return HashMap of All Interactions in the Database, indexed by Interaction Key.
	 * @throws GraphCreationException Graph Creation Error.
	 */
	public static HashMap<String, ArrayList<Interaction>> getAllInteractions()
			throws GraphCreationException {
		if (interactionMap == null) {
			init();
		}
		return interactionMap;
	}

	/**
	 * Gets the Global Graph of All Interactions in the Database.
	 *
	 * @return Global Graph of all Interactions in the Database.
	 * @throws GraphCreationException Graph Creation Error.
	 */
	public static Graph getGlobalGraph() throws GraphCreationException {
		if (globalGraph == null) {
			initGraph();
		}
		return globalGraph;
	}

	/**
	 * Gets all Interactions Assocatied with the Specified Gene.
	 *
	 * @param geneSymbol Gene Symbol.
	 * @return ArrayList of Interaction Objects.
	 */
	public static ArrayList<Interaction> getInteractions(String geneSymbol) {
		Session session = GlobalSession.getInstance().getSession();
		Query query = session.getNamedQuery("org.mskcc.netbox.getInteractionsByGeneSymbol");
		query.setString("geneSymbol", geneSymbol);
		Iterator<Interaction> iterator = query.iterate();
		ArrayList<Interaction> interactionList = new ArrayList<Interaction>();
		while (iterator.hasNext()) {
			Interaction interaction = iterator.next();
			interactionList.add(interaction);
		}
		return interactionList;
	}

	/**
	 * Init method to populate the global interaction and global graph.
	 *
	 * @throws GraphCreationException Graph Creation Error.
	 */
	private static void init() throws GraphCreationException {
		ProgressMonitor pMonitor = ProgressMonitor.getInstance();
		pMonitor.setCurrentMessage("\nLoading network from database.  "
				+ "This will take a few moments.");
		//globalInteractionList = new ArrayList<Interaction>();
		//Session session = GlobalSession.getInstance().getSession();
		//Query query = session.getNamedQuery("org.mskcc.netbox.getAllInteractions");
		//Iterator<Interaction> iterator = query.iterate();



		Iterator<Interaction> iterator = globalInteractionList.iterator();
		interactionMap = new HashMap<String, ArrayList<Interaction>>();
		while (iterator.hasNext()) {
			Interaction interaction = iterator.next();
			appendToInteractionList(interaction, interaction.getGeneA());
			appendToInteractionList(interaction, interaction.getGeneB());
		}
		globalGraph = InteractionToJung.createGraph(globalInteractionList);
	}

	/**
	 * Initializes the Graph Graph.
	 *
	 * @throws GraphCreationException Graph Creation Error.
	 */
	private static void initGraph() throws GraphCreationException {
		if (interactionMap == null) {
			init();
		}
		globalGraph = InteractionToJung.createGraph(globalInteractionList);
	}

	/**
	 * Appends to the Internal Interaction HashMap.
	 *
	 * @param interaction Current Interaction Object.
	 * @param gene        Current Target Gene.
	 */
	private static void appendToInteractionList(Interaction interaction, String gene) {
		ArrayList<Interaction> interactionList = interactionMap.get(gene);
		if (interactionList == null) {
			interactionList = new ArrayList<Interaction>();
			interactionList.add(interaction);
			interactionMap.put(gene, interactionList);
		} else {
			interactionList.add(interaction);
		}
	}
}
