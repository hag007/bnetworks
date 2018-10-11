package org.mskcc.netbox.util;

import edu.uci.ics.jung.algorithms.cluster.ClusterSet;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import org.mskcc.netbox.algorithm.LinkerGene;
import org.mskcc.netbox.algorithm.ModuleDetector;
import org.mskcc.netbox.genomic.GeneWithScore;
import org.mskcc.netbox.graph.GraphUtil;
import org.mskcc.netbox.graph.JungToSif;
import org.mskcc.netbox.graph.Module;
import org.mskcc.netbox.graph.NewmanGirvanModuleDetector;
import org.mskcc.netbox.model.Interaction;
import org.mskcc.netbox.netcarto.NetCartoAnnealing;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Utility Class for Outputting Networks to various simple text formats, including the
 * Cytoscape SIF format.
 *
 * @author Ethan Cerami.
 */
public final class NetworkWriter {

    /**
     * Private Constructor to Prevent Instantitation.
     */
    private NetworkWriter() {
    }

    /**
     * Outputs the Specified Graph to a File called "network.sif".
     *
     * @param g Graph Object.
     * @throws IOException IO Error.
     */
    public static void outputNetwork(Graph g) throws IOException {
        File out = new File("network.sif");
        ProgressMonitor.getInstance().setCurrentMessage("Network file suitable for loading "
                + "into Cytoscape is available at:  " + out.getAbsolutePath());

        FileWriter writer = new FileWriter(out);
        String sif = JungToSif.convertToSif(g);
        writer.write(sif);
        writer.close();
    }

    /**
     * Outputs the Specified Graph to a File called "network.sif".
     *
     * @param interactionList ArrayList of Interaction Objects.
     * @throws IOException IO Error.
     */
    public static void outputNetwork(ArrayList<Interaction> interactionList) throws IOException {
        File out = new File("network_full.sif");
        FileWriter writer = new FileWriter(out);
        for (Interaction interaction : interactionList) {
            writer.write(interaction.getGeneA() + "\t" + interaction.getInteractionType()
                    + "\t" + interaction.getGeneB() + "\n");
        }
        writer.close();
        out = new File("edge_data_source.txt");
        writer = new FileWriter(out);
        writer.write("Data_Source\n");
        for (Interaction interaction : interactionList) {
            writer.write(interaction.getGeneA() + " (" + interaction.getInteractionType()
                    + ") " + interaction.getGeneB() + " = " + interaction.getSource() + "\n");
        }
        writer.close();
    }

    /**
     * Outputs the Node Types to nodes.txt.
     *
     * @param geneList   ArrayList of Input Genes.
     * @param linkerList ArrayList of Linker Genes.
     * @throws IOException IO Error.
     */
    public static void outputNodeAttributes(ArrayList<String> geneList,
                                            ArrayList<LinkerGene> linkerList)
            throws IOException {
        File out = new File("node_type.txt");
        FileWriter writer = new FileWriter(out);

        ProgressMonitor.getInstance().setCurrentMessage("Attribute file containing "
                + "gene ALTERED / LINKER attributes, suitable for loading into Cytoscape is "
                + "available at:  " + out.getAbsolutePath());

        writer.write("NODE_TYPE\n");
        for (String gene : geneList) {
            writer.write(gene + " = ALTERED\n");
        }
        for (LinkerGene gene : linkerList) {
            writer.write(gene.getGene() + " = LINKER\n");
        }
        writer.close();
    }

    /**
     * Outputs the Modules to a set of Files, including modules.dat, modules.txt, and
     * module_x.sif.
     *
     * @param moduleDetector ModuleDetector Object.
     * @throws IOException IO Error.
     */
    public static void outputModules(NewmanGirvanModuleDetector moduleDetector) throws IOException {
        Graph g = moduleDetector.getOptimalGraph();
        StringLabeller labeller = StringLabeller.getLabeller(g);
        WeakComponentClusterer wcSearch = new WeakComponentClusterer();
        ClusterSet moduleSet = wcSearch.extract(g);

        File attrOut = new File("modules.txt");
        FileWriter attrWriter = new FileWriter(attrOut);
        attrWriter.write("MODULE\n");

        ProgressMonitor.getInstance().setCurrentMessage("Attribute file containing "
                + "gene to module assignments, suitable for loading into Cytoscape is "
                + "available at:  " + attrOut.getAbsolutePath());

        //  Iterate through each module
        for (int i = 0; i < moduleSet.size(); i++) {
            Graph gModule = moduleSet.getClusterAsNewSubGraph(i);
            Iterator iterator = gModule.getVertices().iterator();
            while (iterator.hasNext()) {
                Vertex vertex = (Vertex) iterator.next();
                String label = GraphUtil.getVertexLabel(labeller, vertex);
                if (gModule.numVertices() == 1) {
                    attrWriter.write(label + " = -1\n");
                } else {
                    attrWriter.write(label + " = " + i + "\n");
                }
            }
        }
        attrWriter.close();
    }

    /**
     * Outputs the Modules to a set of Files, including modules.dat, modules.txt, and
     * module_x.sif.
     *
     * @param detector NetCartoAnnealing Object.
     * @throws IOException IO Error.
     */
    public static void outputModules(NetCartoAnnealing detector) throws IOException {
        File modOut = new File("modules_sa.dat");
        FileWriter modWriter = new FileWriter(modOut);

        ProgressMonitor.getInstance().setCurrentMessage("Text file summarizing "
                + "all genes in each modules is available at: " + modOut.getAbsolutePath());

        File attrOut = new File("modules_sa.txt");
        FileWriter attrWriter = new FileWriter(attrOut);
        attrWriter.write("MODULE_NG\n");

        ProgressMonitor.getInstance().setCurrentMessage("Attribute file containing "
                + "gene to module assignments, suitable for loading into Cytoscape is "
                + "available at:  " + attrOut.getAbsolutePath());

        ArrayList<String> moduleList = detector.getFinalModuleList();
        HashMap<String, ArrayList<String>> moduleMap = detector.getGlobalModuleMap();

        //  Iterate through each module
        for (String moduleId : moduleList) {
            ArrayList<String> nodeList = moduleMap.get(moduleId);
            modWriter.write("Module_" + moduleId + "\t" + nodeList.size() + "\t");

            for (String node : nodeList) {
                modWriter.write(node + " ");
                if (nodeList.size() == 1) {
                    attrWriter.write(node + " = -1\n");
                } else {
                    attrWriter.write(node + " = " + moduleId + "\n");
                }
            }
            modWriter.write("\n");
        }
        modWriter.close();
        attrWriter.close();
    }

    /**
     * Outputs the Modules to a set of Files, including modules.dat, modules.txt, and
     * module_x.sif.
     *
     * @param detector ModuleDetector Object.
     * @throws IOException IO Error.
     */
    public static void outputModulesByPercentage(ModuleDetector detector)
            throws IOException {
        File modOut = new File("modules.dat");

        ProgressMonitor.getInstance().setCurrentMessage("Text file summarizing "
                + "all genes in each modules is available at: " + modOut.getAbsolutePath());

        FileWriter modWriter = new FileWriter(modOut);

        File attrOut = new File("modules.txt");

        ProgressMonitor.getInstance().setCurrentMessage("Attribute file containing "
                + "gene to module assignments, suitable for loading into Cytoscape is "
                + "available at:  " + attrOut.getAbsolutePath());

        FileWriter attrWriter = new FileWriter(attrOut);
        attrWriter.write("MODULE\n");

        //  Iterate through each module
        ArrayList<Module> moduleList = detector.getModules();
        for (Module module : moduleList) {
            ArrayList<GeneWithScore> geneList = module.getGeneList();
            modWriter.write("Module_" + module.getModuleId() + "\t" + geneList.size() + "\t");
            for (GeneWithScore gene : geneList) {
                modWriter.write(gene.getGene() + " ");
                if (moduleList.size() == 1) {
                    attrWriter.write(gene.getGene() + " = -1\n");
                } else {
                    attrWriter.write(gene.getGene() + " = " + module.getModuleId() + "\n");
                }
            }
            modWriter.write("\n");
        }
        modWriter.close();
        attrWriter.close();
    }

    /**
     * Outputs the Gene Frequencies as Cytoscape Attributes.
     *
     * @param geneList Gene List
     * @throws IOException IO Error.
     */
    public static void outputGeneAttributes(ArrayList<GeneWithScore> geneList)
            throws IOException {
        File modOut = new File("genes.txt");
        ProgressMonitor.getInstance().setCurrentMessage("Attribute file containing "
                + "gene alteration frequency values, suitable for loading into Cytoscape is "
                + "available at:  " + modOut.getAbsolutePath());
        FileWriter attrWriter = new FileWriter(modOut);

        attrWriter.write("ALTERATION_FREQUENCY\n");

        //  Iterate through each gene
        for (GeneWithScore gene : geneList) {
            if (gene.getScore() > 0.0001) {
                attrWriter.write(gene.getGene() + " = " + gene.getScore() + "\n");
            }
        }
        attrWriter.close();
    }

}
