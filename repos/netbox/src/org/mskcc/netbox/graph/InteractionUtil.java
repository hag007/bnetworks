package org.mskcc.netbox.graph;

import org.mskcc.netbox.model.Interaction;
import org.mskcc.netbox.query.InteractionQuery;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Interaction Utility Class.
 *
 * @author Ethan Cerami.
 */
public final class InteractionUtil {
    private static HashMap<String, ArrayList<Interaction>> interactionMap;

    /**
     * Private constructor to prevent instantiation.
     */
    private InteractionUtil() {
    }

    /**
     * Is this a direct interaction between input genes?
     *
     * @param interaction Interaction Object.
     * @param geneSymbol  Target Gene.
     * @param geneSet     Input Gene Set.
     * @return true or false.
     */
    public static boolean isDirectInteraction(Interaction interaction, String geneSymbol,
                                              ArrayList<String> geneSet) {
        String geneA = interaction.getGeneA();
        String geneB = interaction.getGeneB();
        String target = null;
        if (geneA.equalsIgnoreCase(geneSymbol)) {
            target = geneB;
        } else {
            target = geneA;
        }
        return geneSet.contains(target);
    }

    /**
     * Gets all Interactions Associated with Gene List.
     *
     * @param geneList Gene List.
     * @return ArrayList of Interactions Objects.
     * @throws GraphCreationException Graph Creation Error.
     */
    public static ArrayList<Interaction> getAllInteractions(ArrayList<String> geneList)
            throws GraphCreationException {
        if (interactionMap == null) {
            interactionMap = InteractionQuery.getAllInteractions();
        }
        ArrayList<Interaction> globalInteractionList = new ArrayList<Interaction>();
        for (String gene : geneList) {
            ArrayList<Interaction> interactionList = interactionMap.get(gene);
            if (interactionList != null) {
                for (Interaction interaction : interactionList) {
                    globalInteractionList.add(interaction);
                }
            }
        }
        return globalInteractionList;
    }

    /**
     * Connects the Genes into a Network, No Linkers.
     *
     * @param geneList Gene List.
     * @return ArrayList of Interactions Objects.
     * @throws GraphCreationException Graph Creation Error.
     */
    public static ArrayList<Interaction> connectGenesNoLinkers(ArrayList<String> geneList)
            throws GraphCreationException {
        if (interactionMap == null) {
            interactionMap = InteractionQuery.getAllInteractions();
        }
        ArrayList<Interaction> globalInteractionList = new ArrayList<Interaction>();
        for (String gene : geneList) {
            ArrayList<Interaction> interactionList = interactionMap.get(gene);
            if (interactionList != null) {
                for (Interaction interaction : interactionList) {
                    if (InteractionUtil.isDirectInteraction(interaction, gene, geneList)) {
                        globalInteractionList.add(interaction);
                    }
                }
            }
        }
        return globalInteractionList;
    }

    /**
     * Connects the Genes into a Network, With Linkers.
     *
     * @param geneList Gene List.
     * @return ArrayList of Interactions Objects.
     * @throws GraphCreationException Graph Creation Error.
     */
    public static ArrayList<Interaction> connectsGenesWithLinkers(ArrayList<String> geneList)
            throws GraphCreationException {
        if (interactionMap == null) {
            interactionMap = InteractionQuery.getAllInteractions();
        }

        HashSet<String> linkerSet = new HashSet<String>();

        //  First pass: iterate through all genes in list.
        //  For each gene X, get all interactions.
        ArrayList<Interaction> globalInteractionList = new ArrayList<Interaction>();
        for (String gene : geneList) {
            ArrayList<Interaction> interactionList = interactionMap.get(gene);
            if (interactionList != null) {
                for (Interaction interaction : interactionList) {
                    globalInteractionList.add(interaction);

                    //  Track the new linkers separately
                    String other = getOtherGene(gene, interaction);
                    if (!geneList.contains(other)) {
                        linkerSet.add(other);
                    }
                }
            }
        }

        //  Second pass: iterate through all the new linker.
        //  For each linker X, get all interactions, but only keep those interactions that link to
        //  genes in the original list or to genes in the linker set
        ArrayList<String> globalGeneList = new ArrayList<String>();
        globalGeneList.addAll(geneList);
        globalGeneList.addAll(linkerSet);

        for (String linker : linkerSet) {
            ArrayList<Interaction> interactionList = interactionMap.get(linker);
            if (interactionList != null) {
                for (Interaction interaction : interactionList) {
                    if (InteractionUtil.isDirectInteraction(interaction, linker, globalGeneList)) {
                        globalInteractionList.add(interaction);
                    }
                }
            }
        }

        return globalInteractionList;
    }

    /**
     * Gets the other gene in the interaction.
     *
     * @param gene        First gene.
     * @param interaction Interaction Object.
     * @return second gene.
     */
    public static String getOtherGene(String gene, Interaction interaction) {
        String geneA = interaction.getGeneA();
        String geneB = interaction.getGeneB();
        if (geneA.equalsIgnoreCase(gene)) {
            return geneB;
        } else {
            return geneA;
        }
    }
}
