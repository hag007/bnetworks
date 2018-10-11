##<< ##########################################################################################################################
##<< # ModuleDiscoverer - Tutorial                                                                                            #
##<< # Copyright (c) 2015 Leibniz-Institut f?r Naturstoff-Forschung und Infektionsbiologie e.V. - Hans-Knoell-Institut (HKI)  #
##<< # and Friedrich-Schiller-Universit?t Jena.
##<< # Contributors: Sebastian Vlaic <Sebastian.Vlaic@hki-jena.de>                                                            #
##<< # First version: November 2015                                                                                           #
##<< # Filename: Tutorial.r                                                                                                   #
##<< # Language: R                                                                                                            #
##<< #                                                                                                                        #
##<< # This program is free software; you can redistribute it and/or modify it under the terms of the                         #
##<< # GNU General Public License as published by the Free Software Foundation, Version 3.                                    #
##<< #                                                                                                                        #
##<< # This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied     #
##<< # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  #
##<< #                                                                                                                        #
##<< #                                                                                                                        #
##<< # This file contains the tutorial of ModuleDiscoverer. The web-version of this tutorial is available under:              #
##<< # http://www.hki-jena.de/index.php/0/2/490   Others/ModuleDiscoverer                                                     #
##<< ##########################################################################################################################
#
#
#ModuleDiscoverer [1] is an algorithm for the identification of regulatory modules based on large-scale, whole-genome
#protein-protein interaction networks (PPIN) in conjunction with gene expression data from high-throughput platforms
#such as microarrays or RNA-sequencing.
#
#Identification of regulatory modules using ModuleDiscoverer can be divided in three steps
#I to III. Given a PPIN and gene expression data from high-throughput experiments (Input) the algorithm first 
#approximates the PPIN's community structure by enumeration of protein cliques (I). The identified cliques are then 
#tested for their enrichment with proteins associated to differentially expressed genes (II). Finally, the 
#significantly enriched protein cliques are unified to form the regulatory module.
#
#The following tutorial explains the application of ModuleDiscoverer based on a minimal high-confidence PPI network
#of R. norvegicus obtained from the STRING database (version 10) [2]. In detail, tShe PPIN encloses all 277 first
#neighbors of (including) the peroxisome proliferator activated receptor alpha (PPAR-alpha) connected by 4902 edges
#with an score > 0.7. Additionally, we will use expression data of a rat model of diet induced non-alcoholic
#steatohepatitis (NASH) presented in [3]. The raw data can be accessed via Gene Omnibus Express (GEO) [4] using the
#ID GSE8253. Details on data pre-processing and identification of differentially expressed genes are provided in the
#publication. Since ModuleDiscoverer is implemented in the statistical language R [5], which can be obtained from the
#R-project homepage.

set.seed(42) 
library("igraph") 
library("parallel") 
library("doParallel") 
library("foreach")
library("Cairo") 
library("AnnotationDbi")
library('org.Rn.eg.db')
source('ModuleDiscoverer.R') 

#There are 7 objects stored in the RData file:
#
#pparaCenteredNetwork:	An igraph object storing the PPIN used in this tutorial. This network is a sub-network of the
#						PPIN of R. norvegicus as provided by the STRING database (v10). It is composed of all 277 first
#						neighbors of the transcription factor Ppara and all edges with an score > 0.7 (high-confidence
#						interctions) connecting them.
#A:						The adjacency matrix of the PPIN. (required by ModuleDiscoverer)
#proteins:				EnsemblProteinIds of proteins in the PPIN. The order corresponds to the order of rows and columns
#						in the adjacency matrix A.
#vlist:					Table of all 278 vertices in the PPIN. The order of rows corresponds to the order of rows and
#						columns in A. For each vertex in the graph we store the label (content), the number of proteins
#						represented by this node (weight) and the number of connecting edges (degree). (required by
#						ModuleDiscoverer)
#degs:					An array of EnsemblProteinIds corresponding to the 286 differentially expressed genes in the
#						dataset.
#background:			An array of EnsemblProteinIds corresponding to all 4590 genes for which expression was measured
#						in the dataset (the statistical background).
#degs.random.datasets:	A list of 10,000 EnsemblProteinId-sets corresponding to 286 randomly selected genes from the
#						statistical background.


#Step I: Approximation of the PPIN-underlying community structure
#
#In the first step we assess the community structure underlying the given PPIN by the enumeration of cliques in the
#network. ModuleDiscoverer enumerates cliques in the PPIN iteratively. Thus, calling the moduleDiscoverer.fragmentGraph
#function will return all cliques identified in a single iteration. Depending on the number of seed nodes used, the
#algorithm will return either one maximal clique (single-seed approach) or a table of cliques that do not have to be
#necessarily maximal (multi-seed approach). The difference between these two approaches is explained in detail in the
#supplement 1 of the publication. The following command runs a single iteration of ModuleDiscoverer using 1 seed node
#(the single-seed approach). Since the algorithm uses random numbers, the seed-parameter is used to specify the seed
#for the random number generator.

moduleDiscoverer.fragmentGraph(A=A, vlist=vlist, nbrOfSeeds=1, seed=1)
#		content													weight	degree
#[1,]	"74 151 100 68 214 117 248 82 263 40 278 253 246 157"	"14"	"1318"


#The function returns a table where each row corresponds to one identified clique. Since only one seed node was used,
#the identified clique is maximal. The table contains no rows if the selected seed protein is not part of a minimal
#clique of three proteins. In the above example, the algorithm identified a clique composed of 14 proteins (weight),
#which are connected to the remaining 264 proteins in the network via 1318 edges (degree). For each protein in the
#clique (content) the integer value corresponds to a position in the array proteins, which holds the EnsemblProteinID
#of the protein.

#To approximate the PPIN's underlying community structure we now have to run multiple iterations of ModuleDiscoverer
#and combine the enumerated clique in a single table. Since the iterations are independent of each other, this process
#can be easily parallelized. The code below runs 15000 iterations of ModuleDiscoverer in 15 work-packages at 1000
#iterations each. If you have more then one core, you can increase the integer value passed to the makeCluster function
#accordingly. 
cl = makeCluster(1) # creating a cluster with 1 workers (has to be adapted accordingly).
registerDoParallel(cl) # register doParallel package to be used with foreach %dopar%
packages = 15 # number of packages to run

db.results.singleSeed = foreach(j = 1:packages, .combine = 'append' ) %dopar% {
	cat(paste("processing package:",j,'\n')) # output to monitor process.
	set.seed(j) # setting a seed for working package prevents all workers to return identical results.
	times = 3 

	db.results = lapply(1:times, function(i){
		return(moduleDiscoverer.fragmentGraph(A=A, vlist=vlist, nbrOfSeeds=1))
	})
	
	return(db.results)
}
stopCluster(cl) # stop the cluster


#The above code produces a table that summarizes the results of each iteration of ModuleDiscoverer. This table is the
#input for the moduleDiscoverer.createDatabase function, which creates a database of cliques that is used by
#ModuleDiscoverer in all subsequent steps. Once created, this database can be also used for the analysis of other
#high-throughput expression data that is based on the same PPIN.
database.singleSeed = moduleDiscoverer.createDatabase(results=db.results.singleSeed, proteins=proteins)


#Step II: Identification of significantly enriched cliques
#
#In the second step of the algorithm all cliques in the created database are tested for their enrichment with proteins
#associated to differentially expressed genes obtained from experimental high-throughput gene expression data. Using
#the moduleDiscoverer.db.create_MD_object function, the algorithm first collects all required information and assembles
#them in a single input-object.
input.singleSeed = moduleDiscoverer.db.create_MD_object(database=database.singleSeed, foregrounds=list("NASH"=degs), cores=5, background=background, chunks=100, randomDataSets=list(degs.random.datasets))


#The moduleDiscoverer.db.create_MD_object has 9 parameters. While database, and foregrounds are mandatory, background
#and randomDataSets are recommended. The remaining parameters are obligatory. Their default values are based on our
#experience.
#
#database:					(default: NULL) The clique database as returned by the moduleDiscoverer.createDatabase function.
#foregrounds:				(default: NULL) A list of arrays containing all differentially expressed genes. The provided
#							IDs must match the IDs used in the PPIN, e.g., EnsemblProteinIDs for PPINS from the STRING-db.
#background:				(default: NULL) An array of genes present on the microarry platform. This defines the
#							statistical background.
#cores:						(default: 1) Number of cores used for parallelization. If cores > available cores the
#							algorithm will pick automatically available cores - 1.
#chunks:					(default: 5000) This defines the chunks that should be processed at once. The number of chunks
#							is also related to the resolution of the progress-bar that is shown in the terminal.
#randomDataSets:			(default: NULL) For each foreground you can provide a list of arrays with random sets of genes
#							from the statistical background. These lists are passed as a list.
#minBackgroundRequired:		(default: 3) Defines the number of proteins that need to be associated with a gene in the
#							statistical background.
#minBackgroundToCliqueRatio:(default: 0.5) Defines the ratio of the number of proteins in a clique that are associated
#							with a gene in the statistical background and the total number of proteins in that clique.
#minForegroundRequired:		(default: 1) Defines the number of proteins that need to be associated with a DEG.



#The input-object is then used in conjunction with the clique-database to calculate a p-value for each clique.
result.singleSeed = moduleDiscoverer.db.testForCliqueEnrichment(database=database.singleSeed, input=input.singleSeed)


#The moduleDiscoverer.db.testForCliqueEnrichment has 4 parameters. Database and input are mandatory. The remaining
#parameters are obligatory.
#
#database:	(default: NULL) The clique database as returned by the moduleDiscoverer.createDatabase function.
#input:		(default: NULL) The input object as returned by the moduleDiscoverer.db.create_MD_object function.
#cores:		(default: as stored in the input object) Number of cores used for parallelization. If cores > available
#			cores the algorithm will pick automatically available cores - 1.
#chunks:	(default: as stored in the input object) This defines the chunks that should be processed at once. The
#			number of chunks is also related to the resolution of the progress-bar that is shown in the terminal.


#For each clique, the P-values calculated based on random datasets can then be compared to p-values obtained using
#Fisher's-exact test. The method produces two plots, one that plots p-values for all cliques and one 'detailed' plot
#that is restricted to cliques with a user defined p-value (p.value). A filename can be specified to create a PDF
#including both plots.
moduleDiscoverer.db.plotComputedPValues(result=result.singleSeed, fileName=NULL, p.value=0.01)


#The moduleDiscoverer.db.plotComputedPValues has 3 parameters. The parameter result is mandatory. The remaining
#parameters are obligatory.
#
#result:	(default: NULL) The object returned by the moduleDiscoverer.db.testForCliqueEnrichment function.
#fileName:	(default: NULL) Name of the file that should be used to store the plots. If NULL, plots will be displayed
#			directly.
#p.value:	(default: 0.05) A p-value cutoff for the detailed plot.

#The resulting plot compares the permutation-based p-values with p-values of the one-sided Fisher's exact test. For
#this tutorial p-values calculated based on permutations are similar to the p-values calculated by Fisher's exact test
#with a over-optimistic tendency. Finally, the function moduleDiscoverer.db.extractEnrichedCliques is used to extract
#all significantly (defined by the parameter p.value) enriched cliques 

result.singleSeed.ec = moduleDiscoverer.db.extractEnrichedCliques(database=database.singleSeed, result=result.singleSeed, p.value=0.01)

#The moduleDiscoverer.db.extractEnrichedCliques has 5 parameters. The result result is mandatory. The remaining parameters
#are obligatory.
#
#result:					(default: NULL) The object returned by the moduleDiscoverer.db.testForCliqueEnrichment function.
#fileName:					(default: NULL) Name of the file that should be used to store the plots. If NULL, plots will
#							be displayed directly.
#p.value:					(default: 0.05) A p-value cutoff for the detailed plot.
#coef:						(default: 1) Extracts the significantly enriched cliques for the respective foreground.
#useFishersExactTestPvalue:	(default: FALSE) Fisher's exact test p-values are used instead of the permutation based
#							p-values.



#Step III: Assembly of the regulatory module
#
#In the last step of the algorithm the regulatory module is assembled by the unification of all significantly enriched
#cliques. In most cases these cliques will overlap in one or more proteins. Thus, the union of all significantly enriched
#cliques results in a large regulatory module. The function moduleDiscoverer.module.createModule assembles the regulatory
#module and returns an igraph object representing the identified regulatory module.

module.singleSeed = moduleDiscoverer.module.createModule(result=result.singleSeed.ec, module.name="diet induced NASH model")

#The moduleDiscoverer.module.createModule has 3 parameters. The parameter result is mandatory. The remaining parameters
#are obligatory.

#result:		(default: NULL) The object returned by the moduleDiscoverer.db.extractEnrichedCliques function.
#module.name:	(default: "ModuleDiscoverer - regulatory module") Title of the regulatory module.
#cores:			(default: as stored in the input object) Number of cores used for parallelization. If cores > available
#				cores the algorithm will pick automatically available cores - 1.

#The proteins in the regulatory module can then be annotated using additional annotation databases such the
#org.Rn.eg.db package.
module.singleSeed = moduleDiscoverer.module.annotateModule(module=module.singleSeed, annotateWith=c("SYMBOL","GENENAME"), annotation.db="org.Rn.eg.db", nodeIdentifier="ENSEMBLPROT")


#The moduleDiscoverer.module.annotateModule has 4 parameters. All parameters are mandatory.
#
#module:		(default: NULL) The object returned by the moduleDiscoverer.module.createModule function.
#annotation.db:	(default: NULL) The annotation database to be used.
#annotateWith:	(default: NULL) Type of annotation to use.
#nodeIdentifier:(default: NULL) Reference identifier for the proteins.


#Finally, the identified regulatory module can be stored in graphml format to be processed with tools such as
#Cytoscape [6]. To this end, we use the write.graph function provided by the igraph package for R.
write.graph(graph=module.singleSeed, file="ModuleSingleSeed.graphml", format="edgelist")


#Quality control of the regulatory module
#
#After the regulatory module has been successfully identified we can collect some information about its quality and
#biological relevance. First, we assess if the regulatory module is significantly enriched with proteins associated
#to DEGs by using the moduleDiscoverer.db.testNetworkEnrichment function. This test is especially important for the
#multi-seed approach (see additional file 1 of the publication). The function returns two p-values, one calculated
#based on the random data sets provided by the user and one that is based on the one-sided Fisher's exact test. For
#the network shown in figure 2, the calculated permutation-based p-value is < 1e-4 while the Fisher's exact test
#reported p-value is 3.037021e-11. Consequently, the identified regulatory module is significantly enriched with DEGs.

moduleDiscoverer.db.testNetworkEnrichment(database=database.singleSeed, result=result.singleSeed, p.value=0.01)

#The moduleDiscoverer.db.testNetworkEnrichment has 5 parameters. The parameters database and result are mandatory.
#The parameter p.value should be adjusted accordingly. The remaining parameters are obligatory.

#database:	(default: NULL) The clique-database.
#result:	(default: NULL) The object returned by the moduleDiscoverer.db.extractEnrichedCliques function.
#p.value:	(default: 0.05) The p-value cutoff used for identification of the regulatory module.
#cores:		(default: as stored in the input object) Number of cores used for parallelization. If cores > available
#			cores the algorithm will pick automatically available cores - 1.
#chunks:	(default: as stored in the input object) This defines the chunks that should be processed at once. The
#			number of chunks is also related to the resolution of the progress-bar that is shown in the terminal.

#Next we estimate the stability of the identified regulatory module, i.e., the extend to which the identified regulatory
#module would be recovered if we would perform another n iterations for clique discovery. This is important since our
#algorithm is a randomization-based heuristic. Since performing additional iterations is cost intense in terms of
#computational time we generate bootstrap samples (model-free re-sampling with replacement) of the clique-database
#simulating additional n iterations. Pairwise comparison of the bootstrap-sample-based regulatory modules is then used
#to estimate the stability of the identified regulatory module.

module.stability = moduleDiscoverer.db.computeNetworkStability(database=database.singleSeed, result=result.singleSeed, p.value=0.01)

#The moduleDiscoverer.db.computeNetworkStability has 8 parameters. The parameters database and result are mandatory. The parameter p.value should be adjusted accordingly. The remaining parameters are obligatory.

#database:			(default: NULL) The clique-database.
#result:			(default: NULL) The object returned by the moduleDiscoverer.db.extractEnrichedCliques function.
#p.value:			(default: 0.05) The p-value cutoff used for identification of the regulatory module.
#repeats:			(default: 100) The number of bootstrap samples.
#size:				(default: number of iterations used for the database) The number iterations to simulate for the
#					bootstrap samples.
#module.reference:	(default: NULL) A reference module to compare the identified regulatory module against.
#cores:				(default: as stored in the input object) Number of cores used for parallelization.
#					If cores > available cores the algorithm will pick automatically available cores - 1.
#chunks:			(default: as stored in the input object) This defines the chunks that should be processed at once.
#					The number of chunks is also related to the resolution of the progress-bar that is shown in the terminal.

module.stability$subSampling
#				median-nodeStability	CB-nodeStability.95%	median-edgeStability	CB-edgeStability.95%
#Foreground-1	1						0.9821429				0.9915254				0.9449153


#For the identified regulatory module shown in figure 2 the result of the stability analysis is shown in table 2. For nodes
#and edges the upper 95% bound is above 0.94 stating that 95% of all pairwise module comparisons resulted in a similarity
#of at least 94%.

#References:
#[1] Vlaic, S. et al. (----). ModuleDiscoverer: Identification of regulatory modules in protein-protein interaction networks. submitted, -(-), ---.
#[2] Szklarczyk, D. et al. (2015). String v10: protein-protein interaction networks, integrated over the tree of life. Nucleic Acids Res., 43(Database issue), D447-D452.
#[3] Baumgardner, J.N. et al. (2008). A new model for nonalcoholic steatohepatitis in the rat utilizing total enteral nutrition to overfeed a high-polyunsaturated fat diet. Am J Physiol Gastrointest Liver Physiol., 294(1), G27-G38.
#[4] Barrett, T. et al. (2013). Ncbi geo: archive for functional genomics data sets--update. Nucleic Acids Res., 41(Database issue), D991-D995.
#[5] R Core Team (2015). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.
#[6] Shannon, P. et al. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks Genome research Cold Spring Harbor Lab, 13, 2498-250

