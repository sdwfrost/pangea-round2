library("XML")
library("ape")

#Define main function
#Inputs:
# ~ seqs - sequences in matrix and character format
# ~ years - list of dates in same order as sequences
# ~ clusterList - list of clusters (of size >=3)
# ~ niters - number of MCMC iterations
# ~ samp - sampling frequency
# ~ outfile - name to give to output files (will be appended with .log and .trees)

getMultiLocusXML<-function(seqs,years,clusterList,niters,samp,outfile){
  doc<-newXMLDoc()
  beast<-newXMLNode("beast",doc=doc)
  
  nseqs<-dim(seqs)[1] #Get total number of sequences
  nam<-rownames(seqs) #Get names of sequences
  newXMLCommentNode(paste("The list of taxa to be analysed (can also include dates/ages).          -->
                          <!-- ntax=",nseqs,sep=""), parent = beast)
  
  taxa<-newXMLNode("taxa", attrs = c(id = "taxa"), parent = beast)
  
  for(i in 1:nseqs){
    taxon<-newXMLNode("taxon",attrs=c(id=nam[i]), parent = taxa)
    date<-newXMLNode("date",attrs=c(value=years[i],direction="forwards",units="years"),parent=taxon)
  }
  
  for(j in 1:length(clusterList)){
    newXMLCommentNode(paste("gene name = cluster",j," ntax= ",length(clusterList[[j]]),sep=""), parent = beast)
    alignTaxa<-newXMLNode("taxa",attrs=c(id=paste("cluster",j,".taxa",sep="")),parent=beast)
    
    for(k in 1:length(clusterList[[j]])){
      alignTaxon<-newXMLNode("taxon",attrs=c(idref=clusterList[[j]][k]),parent=alignTaxa)
    }
  }
  
  for(j in 1:length(clusterList)){
    thisSeqLength<-length(seqs[which(nam==clusterList[[j]][1]),])
    newXMLCommentNode("The sequence alignment (each sequence refers to a taxon above).", parent = beast)
    newXMLCommentNode(paste("ntax=",length(clusterList[[j]])," nchar=",thisSeqLength,sep=""))
    alignment<-newXMLNode("alignment",attrs=c(id=paste("alignment",j,sep=""),dataType="nucleotide"),parent=beast)
    
    for(k in 1:length(clusterList[[j]])){
      thisSeq<-paste(toupper(seqs[which(nam==clusterList[[j]][k]),]),collapse="")
      sequence<-newXMLNode("sequence",newXMLNode("taxon",attrs=c(idref=clusterList[[j]][k])),thisSeq,parent=alignment)
    }
  }
  
  for(j in 1:length(clusterList)){
    newXMLCommentNode("The unique patterns from 1 to end.", parent = beast)
    patterns<-newXMLNode("patterns",attrs=c(id=paste("cluster",j,".patterns",sep=""),from="1",strip="false"),parent=beast)
    newXMLNode("alignment",attrs=c(idref=paste("alignment",j,sep="")),parent=patterns)
  }
  
  newXMLCommentNode("A prior assumption that the population size has grown exponentially throughout the time spanned by the genealogy.")
  exponentialGrowth<-newXMLNode("exponentialGrowth",attrs=c(id="exponential", units="years"),parent=beast)
  populationSize<-newXMLNode("populationSize",newXMLNode("parameter",attrs=c(id="exponential.popSize", value="7000", lower="0.0")),parent=exponentialGrowth)
  growthRate<-newXMLNode("growthRate",newXMLNode("parameter",attrs=c(id="exponential.growthRate",value="0.02")),parent=exponentialGrowth)
  
  for(j in 1:length(clusterList)){
    newXMLCommentNode("Generate a random starting tree under the coalescent process",parent=beast)
    coalescentSimulator<-newXMLNode("coalescentSimulator",attrs=c(id=paste("cluster",j,".startingTree",sep="")),parent=beast)
    newXMLNode("taxa",attrs=c(idref=paste("cluster",j,".taxa",sep="")),parent=coalescentSimulator)
    newXMLNode("exponentialGrowth",attrs=c(idref="exponential"),parent=coalescentSimulator)
  }
  
  for(j in 1:length(clusterList)){
    newXMLCommentNode("Generate a tree model",parent=beast)
    treeModel<-newXMLNode("treeModel",attrs=c(id=paste("cluster",j,".treeModel",sep="")),parent=beast)
    newXMLNode("coalescentTree",attrs=c(idref=paste("cluster",j,".startingTree",sep="")),parent=treeModel)
    newXMLNode("rootHeight",newXMLNode("parameter",attrs=c(id=paste("cluster",j,".treeModel.rootHeight",sep=""))),parent=treeModel)
    newXMLNode("nodeHeights",attrs=c(internalNodes="true"),newXMLNode("parameter",attrs=c(id=paste("cluster",j,".treeModel.internalNodeHeights",sep=""))),parent=treeModel)
    newXMLNode("nodeHeights",attrs=c(internalNodes="true",rootNode="true"),newXMLNode("parameter",attrs=c(id=paste("cluster",j,".treeModel.allInternalNodeHeights",sep=""))),parent=treeModel)
  }
  
  for(j in 1:length(clusterList)){
    newXMLCommentNode("Taxon Sets",parent=beast)
    tmrcaStatistic<-newXMLNode("tmrcaStatistic",attrs=c(id=paste("tmrca(cluster",j,")",sep=""),includeStem="false"),parent=beast)
    newXMLNode("mrca",newXMLNode("taxa",attrs=c(idref=paste("cluster",j,".taxa",sep=""))),parent=tmrcaStatistic)
    newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep="")),parent=tmrcaStatistic)
  }
  
  for(j in 1:length(clusterList)){
    newXMLCommentNode("Generate a coalescent likelihood",parent=beast)
    coalescentLikelihood<-newXMLNode("coalescentLikelihood",attrs=c(id=paste("cluster",j,".coalescent",sep="")),parent=beast)
    newXMLNode("model",newXMLNode("exponentialGrowth",attrs=c(idref="exponential")),parent=coalescentLikelihood)
    newXMLNode("populationTree",newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep=""))),parent=coalescentLikelihood)
  }
  
  newXMLCommentNode("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut (2006) PLoS Biology 4, e88)",parent=beast)
  
  branchRates<-newXMLNode("discretizedBranchRates",attrs=c(id="cluster1.branchRates"),parent=beast)
  newXMLNode("treeModel",attrs=c(idref="cluster1.treeModel"),parent=branchRates)
  distribution<-newXMLNode("distribution",parent=branchRates)
  lnDistModel<-newXMLNode("logNormalDistributionModel",attrs=c(meanInRealSpace="true"),parent=distribution)
  newXMLNode("mean",newXMLNode("parameter",attrs=c(id="ucld.mean",value="0.003",lower="0.0")),parent=lnDistModel)
  newXMLNode("stdev",newXMLNode("parameter",attrs=c(id="ucld.stdev",value="0.4",lower="0.0")),parent=lnDistModel)
  newXMLNode("rateCategories",newXMLNode("parameter",attrs=c(id="cluster1.branchRates.categories")),parent=branchRates)
  
  meanRate<-newXMLNode("rateStatistic",attrs=c(id="cluster1.meanRate",name="cluster1.meanRate",mode="mean",internal="true",external="true"),parent=beast)
  newXMLNode("treeModel",attrs=c(idref="cluster1.treeModel"),parent=meanRate)
  newXMLNode("discretizedBranchRates",attrs=c(idref="cluster1.branchRates"),parent=meanRate)
  
  coefficientOfVariation<-newXMLNode("rateStatistic",attrs=c(id="cluster1.coefficientOfVariation",name=paste("cluster",j,".coefficientOfVariation",sep=""),mode="coefficientOfVariation",internal="true",external="true"),parent=beast)
  newXMLNode("treeModel",attrs=c(idref="cluster1.treeModel"),parent=coefficientOfVariation)
  newXMLNode("discretizedBranchRates",attrs=c(idref="cluster1.branchRates"),parent=coefficientOfVariation)
  
  covariance<-newXMLNode("rateCovarianceStatistic",attrs=c(id="cluster1.covariance",name=paste("cluster",j,".covariance",sep="")),parent=beast)
  newXMLNode("treeModel",attrs=c(idref="cluster1.treeModel"),parent=covariance)
  newXMLNode("discretizedBranchRates",attrs=c(idref="cluster1.branchRates"),parent=covariance)
  
  for(j in 2:length(clusterList)){
    branchRates<-newXMLNode("discretizedBranchRates",attrs=c(id=paste("cluster",j,".branchRates",sep="")),parent=beast)
    newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep="")),parent=branchRates)
    distribution<-newXMLNode("distribution",parent=branchRates)
    lnDistModel<-newXMLNode("logNormalDistributionModel",attrs=c(meanInRealSpace="true"),parent=distribution)
    newXMLNode("mean",newXMLNode("parameter",attrs=c(idref="ucld.mean")),parent=lnDistModel)
    newXMLNode("stdev",newXMLNode("parameter",attrs=c(idref="ucld.stdev")),parent=lnDistModel)
    newXMLNode("rateCategories",newXMLNode("parameter",attrs=c(id=paste("cluster",j,".branchRates.categories",sep=""))),parent=branchRates)
    
    meanRate<-newXMLNode("rateStatistic",attrs=c(id=paste("cluster",j,".meanRate",sep=""),name=paste("cluster",j,".meanRate",sep=""),mode="mean",internal="true",external="true"),parent=beast)
    newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep="")),parent=meanRate)
    newXMLNode("discretizedBranchRates",attrs=c(idref=paste("cluster",j,".branchRates",sep="")),parent=meanRate)
    
    coefficientOfVariation<-newXMLNode("rateStatistic",attrs=c(id=paste("cluster",j,".coefficientOfVariation",sep=""),name=paste("cluster",j,".coefficientOfVariation",sep=""),mode="coefficientOfVariation",internal="true",external="true"),parent=beast)
    newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep="")),parent=coefficientOfVariation)
    newXMLNode("discretizedBranchRates",attrs=c(idref=paste("cluster",j,".branchRates",sep="")),parent=coefficientOfVariation)
    
    covariance<-newXMLNode("rateCovarianceStatistic",attrs=c(id=paste("cluster",j,".covariance",sep=""),name=paste("cluster",j,".covariance",sep="")),parent=beast)
    newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep="")),parent=covariance)
    newXMLNode("discretizedBranchRates",attrs=c(idref=paste("cluster",j,".branchRates",sep="")),parent=covariance)
  }
  
  newXMLCommentNode("The general time reversible (GTR) substitution model",parent=beast)
  gtrModel<-newXMLNode("gtrModel",attrs=c(id="gtr"),parent=beast)
  frequencies<-newXMLNode("frequencies",parent=gtrModel)
  frequencyModel<-newXMLNode("frequencyModel",attrs=c(dataType="nucleotide"),parent=frequencies)
  newXMLNode("frequencies",newXMLNode("parameter",attrs=c(id="frequencies",value="0.378 0.178 0.228 0.216")),parent=frequencyModel)
  
  newXMLNode("rateAC",newXMLNode("parameter",attrs=c(id="ac",value="0.112",lower="0.0")),parent=gtrModel)
  newXMLNode("rateAG",newXMLNode("parameter",attrs=c(id="ag",value="0.927",lower="0.0")),parent=gtrModel)
  newXMLNode("rateAT",newXMLNode("parameter",attrs=c(id="at",value="0.114",lower="0.0")),parent=gtrModel)
  newXMLNode("rateCG",newXMLNode("parameter",attrs=c(id="cg",value="0.0968",lower="0.0")),parent=gtrModel)
  newXMLNode("rateGT",newXMLNode("parameter",attrs=c(id="gt",value="0.0948",lower="0.0")),parent=gtrModel)
  
  newXMLCommentNode("Site model",parent=beast)
  siteModel<-newXMLNode("siteModel",attrs=c(id="siteModel"),parent=beast)
  newXMLNode("substitutionModel",newXMLNode("gtrModel",attrs=c(idref="gtr")),parent=siteModel)
  newXMLNode("gammaShape",attrs=c(gammaCategories="4"),newXMLNode("parameter",attrs=c(id="alpha",value="0.5",lower="0.0")),parent=siteModel)
  newXMLNode("proportionInvariant",newXMLNode("parameter",attrs=c(id="pInv",value="0.3",lower="0.0",upper="1.0")),parent=siteModel)
  
  newXMLCommentNode("Likelihood for tree given sequence data",parent=beast)
  for(j in 1:length(clusterList)){
    treeLikelihood<-newXMLNode("treeLikelihood",attrs=c(id=paste("cluster",j,".treeLikelihood",sep=""),useAmbiguities="false"),parent=beast)
    newXMLNode("patterns",attrs=c(idref=paste("cluster",j,".patterns",sep="")),parent=treeLikelihood)
    newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep="")),parent=treeLikelihood)
    newXMLNode("siteModel",attrs=c(idref="siteModel"),parent=treeLikelihood)
    newXMLNode("discretizedBranchRates",attrs=c(idref=paste("cluster",j,".branchRates",sep="")),parent=treeLikelihood)
  }
  
  newXMLCommentNode("DefineOperators",parent=beast)
  operators<-newXMLNode("operators",attrs=c(id="operators",optimizationSchedule="default"),parent=beast)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="ac")),parent=operators)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="ag")),parent=operators)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="at")),parent=operators)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="cg")),parent=operators)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="gt")),parent=operators)
  
  newXMLNode("deltaExchange",attrs=c(delta="0.01",weight="0.1"),newXMLNode("parameter",attrs=c(idref="frequencies")),parent=operators)
  
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="alpha")),parent=operators)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="pInv")),parent=operators)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="ucld.mean")),parent=operators)
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="0.1"),newXMLNode("parameter",attrs=c(idref="ucld.stdev")),parent=operators)
  
  for(j in 1:length(clusterList)){
    newXMLNode("subtreeSlide",attrs=c(size="0.03",gaussian="true",weight="15"),newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep=""))),parent=operators)
    newXMLNode("narrowExchange",attrs=c(weight="15"),newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep=""))),parent=operators)
    newXMLNode("wideExchange",attrs=c(weight="3"),newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep=""))),parent=operators)
    newXMLNode("wilsonBalding",attrs=c(weight="3"),newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep=""))),parent=operators)
    newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="3"),newXMLNode("parameter",attrs=c(idref=paste("cluster",j,".treeModel.rootHeight",sep=""))),parent=operators)
    newXMLNode("uniformOperator",attrs=c(weight="30"),newXMLNode("parameter",attrs=c(idref=paste("cluster",j,".treeModel.internalNodeHeights",sep=""))),parent=operators)
  }
  
  newXMLNode("scaleOperator",attrs=c(scaleFactor="0.75",weight="3"),newXMLNode("parameter",attrs=c(idref="exponential.popSize")),parent=operators)
  newXMLNode("randomWalkOperator",attrs=c(windowSize="1.0",weight="3"),newXMLNode("parameter",attrs=c(idref="exponential.growthRate")),parent=operators)
  
  for(j in 1:length(clusterList)){
    newXMLNode("upDownOperator",attrs=c(scaleFactor="0.75",weight="3"),newXMLNode("up",newXMLNode("parameter",attrs=c(idref="ucld.mean"))),newXMLNode("down",newXMLNode("parameter",attrs=c(idref=paste("cluster",j,".treeModel.allInternalNodeHeights",sep="")))),parent=operators)
    newXMLNode("swapOperator",attrs=c(size="1",weight="10",autoOptimize="false"),newXMLNode("parameter",attrs=c(idref=paste("cluster",j,".branchRates.categories",sep=""))),parent=operators)
    newXMLNode("uniformIntegerOperator",attrs=c(weight="10"),newXMLNode("parameter",attrs=c(idref=paste("cluster",j,".branchRates.categories",sep=""))),parent=operators)
  }
  
  newXMLCommentNode("Define MCMC",parent=beast)
  mcmc<-newXMLNode("mcmc",attrs=c(id="mcmc",chainLength=niters, autoOptimize="true", operatorAnalysis=paste(outfile,".ops.txt",sep="")),parent=beast)
  posterior<-newXMLNode("posterior",attrs=c(id="posterior"),parent=mcmc)
  prior<-newXMLNode("prior",attrs=c(id="prior"),parent=posterior)
  
  newXMLNode("gammaPrior",attrs=c(shape="0.05",scale="10.0",offset="0.0"),newXMLNode("parameter",attrs=c(idref="ac")),parent=prior)
  newXMLNode("gammaPrior",attrs=c(shape="0.05",scale="10.0",offset="0.0"),newXMLNode("parameter",attrs=c(idref="ag")),parent=prior)
  newXMLNode("gammaPrior",attrs=c(shape="0.05",scale="10.0",offset="0.0"),newXMLNode("parameter",attrs=c(idref="at")),parent=prior)
  newXMLNode("gammaPrior",attrs=c(shape="0.05",scale="10.0",offset="0.0"),newXMLNode("parameter",attrs=c(idref="cg")),parent=prior)
  newXMLNode("gammaPrior",attrs=c(shape="0.05",scale="10.0",offset="0.0"),newXMLNode("parameter",attrs=c(idref="gt")),parent=prior)
  
  newXMLNode("uniformPrior",attrs=c(lower="0.0",upper="1.0"),newXMLNode("parameter",attrs=c(idref="frequencies")),parent=prior)
  newXMLNode("exponentialPrior",attrs=c(mean="0.5",offset="0.0"),newXMLNode("parameter",attrs=c(idref="alpha")),parent=prior)
  newXMLNode("uniformPrior",attrs=c(lower="0.0",upper="1.0"),newXMLNode("parameter",attrs=c(idref="pInv")),parent=prior)
  newXMLNode("exponentialPrior",attrs=c(mean="0.33333333333333",offset="0.0"),newXMLNode("parameter",attrs=c(idref="ucld.stdev")),parent=prior)
  newXMLNode("uniformPrior",attrs=c(lower="0.0001",upper="0.01"),newXMLNode("parameter",attrs=c(idref="ucld.mean")),parent=prior)
  newXMLNode("oneOnXPrior",newXMLNode("parameter",attrs=c(idref="exponential.popSize")),parent=prior)
  newXMLNode("laplacePrior",attrs=c(mean="0.0",scale="0.026315258205646237"),newXMLNode("parameter",attrs=c(idref="exponential.growthRate")),parent=prior)
  
  for(j in 1:length(clusterList)){
    newXMLNode("coalescentLikelihood",attrs=c(idref=paste("cluster",j,".coalescent",sep="")),parent=prior)
  }
  
  likelihood<-newXMLNode("likelihood",attrs=c(id="likelihood"),parent=posterior)
  for(j in 1:length(clusterList)){
    newXMLNode("treeLikelihood",attrs=c(idref=paste("cluster",j,".treeLikelihood",sep="")),parent=likelihood)
  }
  
  newXMLNode("operators",attrs=c(idref="operators"),parent=mcmc)
  
  newXMLCommentNode("Write log to screen",parent=mcmc)
  screenLog<-newXMLNode("log",attrs=c(id="screenLog", logEvery="1000"),parent=mcmc)
  newXMLNode("column",attrs=c(label="posterior",dp="4",width="12"),newXMLNode("posterior",attrs=c(idref="posterior")),parent=screenLog)
  newXMLNode("column",attrs=c(label="prior",dp="4",width="12"),newXMLNode("prior",attrs=c(idref="prior")),parent=screenLog)
  newXMLNode("column",attrs=c(label="likelihood",dp="4",width="12"),newXMLNode("likelihood",attrs=c(idref="likelihood")),parent=screenLog)
  
  for(j in 1:length(clusterList)){
        newXMLNode("column",attrs=c(label=paste("cluster",j,".treeModel.rootHeight",sep=""),sf="6",width="12"),newXMLNode("parameter",attrs=c(idref=paste("cluster",j,".treeModel.rootHeight",sep=""))),parent=screenLog)
  }
  newXMLNode("column",attrs=c(label="ucld.mean",sf="6",width="12"),newXMLNode("ucld.mean",attrs=c(idref="ucld.mean")),parent=screenLog)
  
  
  newXMLCommentNode("Write log to file",parent=mcmc)
  fileLog<-newXMLNode("log",attrs=c(id="fileLog",logEvery=samp,fileName=paste(outfile,".log",sep=""),overwrite="false"),parent=mcmc)
  newXMLNode("posterior",attrs=c(idref="posterior"),parent=fileLog)
  newXMLNode("prior",attrs=c(idref="prior"),parent=fileLog)
  newXMLNode("likelihood",attrs=c(idref="likelihood"),parent=fileLog)
  
  for(j in 1:length(clusterList)){
    newXMLNode("parameter",attrs=c(idref=paste("cluster",j,".treeModel.rootHeight",sep="")),parent=fileLog)
    newXMLNode("tmrcaStatistic",attrs=c(idref=paste("tmrca(cluster",j,")",sep="")),parent=fileLog)
  }
  
  newXMLNode("parameter",attrs=c(idref="exponential.popSize"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="exponential.growthRate"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="ac"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="ag"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="at"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="cg"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="gt"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="frequencies"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="alpha"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="pInv"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="ucld.mean"),parent=fileLog)
  newXMLNode("parameter",attrs=c(idref="ucld.stdev"),parent=fileLog)
    
  for(j in 1:length(clusterList)){
    newXMLNode("rateStatistic",attrs=c(idref=paste("cluster",j,".meanRate",sep="")),parent=fileLog)
    newXMLNode("rateStatistic",attrs=c(idref=paste("cluster",j,".coefficientOfVariation",sep="")),parent=fileLog)
    newXMLNode("rateCovarianceStatistic",attrs=c(idref=paste("cluster",j,".covariance",sep="")),parent=fileLog)
    
    newXMLNode("treeLikelihood",attrs=c(idref=paste("cluster",j,".treeLikelihood",sep="")),parent=fileLog)
    newXMLNode("coalescentLikelihood",attrs=c(idref=paste("cluster",j,".coalescent",sep="")),parent=fileLog)
  }

  newXMLCommentNode("Write tree log to file",parent=mcmc)
  for(j in 1:length(clusterList)){
    treeLog<-newXMLNode("logTree",attrs=c(id=paste("cluster",j,".treeFileLog",sep=""),logEvery=samp,fileName=paste(outfile,"_cluster",j,".trees",sep=""),nexusFormat="true",sortTranslationTable="true"),parent=mcmc)
    newXMLNode("treeModel",attrs=c(idref=paste("cluster",j,".treeModel",sep="")),parent=treeLog)
    newXMLNode("trait",attrs=c(name="rate",tag="rate"),newXMLNode("discretizedBranchRates",attrs=c(idref=paste("cluster",j,".branchRates",sep=""))),parent=treeLog)
    newXMLNode("posterior",attrs=c(idref="posterior"),parent=treeLog)
  }
  
  report<-newXMLNode("report",parent=beast)
  newXMLNode("property",attrs=c(name="timer"),newXMLNode("mcmc",attrs=c(idref="mcmc")),parent=report)
  
  
  
  saveXML(doc, file=paste(outfile,".xml",sep=""), prefix = "<?xml version=\"1.0\" standalone=\"yes\"?>\n
          <!-- Generated by BEAUTi v1.8.0                                              -->
          <!--       by Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard         -->
          <!--       Department of Computer Science, University of Auckland and        -->
          <!--       Institute of Evolutionary Biology, University of Edinburgh        -->
          <!--       David Geffen School of Medicine, University of California, Los Angeles-->
          <!--       http://beast.bio.ed.ac.uk/                                        -->
          ")
  
}

#Read in some test sequences (in matrix format)
#seqs<-read.dna("Align1.txt",format="fasta",as.character=T,as.matrix=T)

#Get some dates
#years<-1:dim(seqs)[1]

#Get list of clusters
#clusterList<-list(c("Seq1","Seq2","Seq3"),c("Seq4","Seq5","Seq6"))

#Example usage
#options(scipen=999)
#getMultiLocusXML(seqs,years,clusterList,10000000,1000,"test")
