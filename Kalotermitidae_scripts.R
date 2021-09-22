###CONTENTS
#package setup
#file import
#combine BEAST, IQTREE, RASP onto consensus tree
#stats for summary tree
#plot summary tree
#polytomise branches without consensus support
#Biogeographic realm map
#traverse Kalo-subtree to get tip-to-root node RASP probabilities
#fossil ages plot
#subsample Neotropical Kalos, (analyze with RASP), plot
#historical biogeography with Biogeobears

################
#package setup
###############
#install.packages("pacman")
pacman::p_load(sf,ggfortify,magrittr, pheatmap,tidyverse,
               ggplot2,data.table,gtable,grid,gridExtra,
               rlang,plotly,phangorn,ape,treeio,DECIPHER,
               dendextend,phylogram,dplyr,ggrepel,tidytree,
               ggimage,gtools,adephylo,ggtree,xlsx, 
               evobiR, #for beast tree subsampling
               phytools, #for getParent function
               gdata,
               combinat, #for combinations and permutations
               BioGeoBEARS
               ) 

###############
#file import
###############
#NOTE1: no spaces allowed in tree tip labels
#NOTE2: place "data" directory into the working directory

##directory with input files
inDir<-"data/"
##first sheet of excel table
KaloTable<-read.xlsx2(paste(inDir,"TableS1_v10.xlsx",sep=""), 1) #use read.xlsx2 so empty cells are imported as no values rather than NAs
#add "TipLabels" column by pasting: "Collection_code_updated"_"Genus"_"Species"
KaloTable<-add_column(KaloTable, TipLabels = paste(KaloTable$Collection_code_updated,KaloTable$Genus,KaloTable$Species,sep="_"), .after = 2) 
#add "TipLabels2" column by pasting: "Genus" "Species"
KaloTable<-add_column(KaloTable, TipLabels2 = paste(KaloTable$Genus,KaloTable$Species,sep=" "), .after = 2) 
##input geography for BioGeoBears
geogfn <- paste(inDir,"KaloOnly_geography.txt",sep="")
##list of Kalos to be removed for node age plots
RedundantKalos<-read_lines(paste(inDir,"Redundant-KaloTips_v2.txt",sep=""))
##trees from Treegraph imported as S3 phylo objects
#vector with tree files
PhyloTrees.files<-paste(inDir,
                        "TreeGraph/", 
                        c("Treegraph_BEAST6supports.tre", 
                        "Treegraph_BEAST6conflicts.tre",
                        "Treegraph_BEAST8supports.tre",
                        "Treegraph_IQTREE_with3rd_conflicts.tre",
                        "Treegraph_IQTREE_with3rd_supports.tre",
                        "Treegraph_IQTREE_without3rd_conflicts.tre",
                        "Treegraph_IQTREE_without3rd_supports.tre"),
                        sep="")
#vector with tree ids 
PhyloTrees.ids<-c("BEAST6supports", 
               "BEAST6conflicts",
               "BEAST8supports",
               "IQTREE_with3rd_conflicts",
               "IQTREE_with3rd_supports",
               "IQTREE_without3rd_conflicts",
               "IQTREE_without3rd_supports")
# newick trees as list, name list elements by tree IDs
PhyloTrees.lst<-lapply(PhyloTrees.files, function(x) treeio::read.newick(x))
names(PhyloTrees.lst)<-PhyloTrees.ids
##BEAST tree (with node age intervals) as backbone tree for all analyses and plotting
#tree0_treedata <- treeio::read.beast (file = paste(inDir,"BEAST6_1-5.tree",sep="")) #alternative bayesian backbone topology based on sequence data with 3rd codon positions, uncomment and comment above to plot alternative Figure1 with this topology
tree0_treedata <- treeio::read.beast (file = paste(inDir,"BEAST8_1-10_40kresampled.tree",sep=""))
##reformated RASP tables (manually fixed error in RASP table!)
RASP_indir<-paste(inDir,"reformated-RASP",sep="")
RASP_files<-list.files(RASP_indir,full.names=T)
RASP_files_refTopo<-RASP_files[ grepl("Beast8", RASP_files) ] #RASP files for only reference topology
#RASP_files_refTopo<-RASP_files[ grepl("Beast6", RASP_files) ] #RASP files for alternative reference topology
RASP_IDs<-paste(gsub(".txt","",list.files(RASP_indir)),"_prob",sep="") #get analysis IDs from file names by dropping ".txt" suffix and adding "_prob"
RASP_IDs_refTopo<-RASP_IDs[ grepl("Beast8", RASP_IDs) ] #RASP IDs for only reference topology
#RASP_IDs_refTopo<-RASP_IDs[ grepl("Beast6", RASP_IDs) ] #RASP IDs for alternative reference topology
RASP.lst<-lapply(RASP_files_refTopo, function(x) read_tsv(x))
  #=> list of 6 RASP analyses results
  #=> 173 internal node info
##list of tips with <2Myears distance to be dropped
tips2drop<-as.vector(read.csv(paste(inDir,"Below2MY-tips-to-drop.txt",sep=""),header = F)$V1) 
  #CHECK: how many tips?
  length(tips2drop)
    #=> 36
##fossil table (positions of fossils specified by node numbers => specific for each backbone tree topology)
fossils<-read_csv("Fossils_for_R_BEAST8.csv")
##kalotermitid fossil records from PaleoDB
PaleoDB<-read.xlsx(paste(inDir,"TableS7.xlsx",sep=""),sheetIndex=1, na = c("", "NA"),stringsAsFactors = FALSE) #import from excel, empty cells as NAs, strings not as factor to enable mutate later
  #CHECK
  nrow(PaleoDB)
  #=> 74 records
##reformatted RASP tables for sumbsampled Neotropics analysis
RASP_indir_subsample<-paste(inDir,"reformated-RASP-subsampled",sep="")
RASP_files_subsample<-list.files(RASP_indir_subsample,full.names=T)
RASP_IDs_subsample<-paste(gsub(".txt","",list.files(RASP_indir_subsample)),"_prob",sep="") #get analysis IDs from file names by dropping ".txt" suffix and adding "_prob"
RASP.lst_subsample<-lapply(RASP_files_subsample, function(x) read_tsv(x))
  
###############################################
#combine BEAST, IQTREE, RASP onto consensus tree
###############################################
tree0_phylo<-as.phylo(tree0_treedata) #convert tree0 to phylo object required for getMRCA(function)
tree_all<-as_tibble(tree0_treedata) #convert tree object to tibble for data manipulation (https://yulab-smu.top/treedata-book/chapter2.html)

#add supports and conflicts from Treegraph 
for (i in 1:length(PhyloTrees.ids)) {
  treeID<-PhyloTrees.ids[i]
  tree_all[treeID]<-as_tibble(PhyloTrees.lst[[i]])$label
  tree_all[!is.na(tree_all$label),treeID]<-NA #remove tip labels imported along with node labels (=supports,conflicts) from Treegraph trees
}

#add a YES/NO variable for nodes with above-threshold support from all trees, skip empty values
tree_all$support_by_all<-ifelse(as.numeric(tree_all$BEAST6supports)>=0.95 & tree_all$BEAST6supports!="" & 
                                as.numeric(tree_all$BEAST8supports)>=0.95 & tree_all$BEAST8supports!="" &
                                as.numeric(tree_all$IQTREE_with3rd_supports)>=95 & tree_all$IQTREE_with3rd_supports!="" &
                                as.numeric(tree_all$IQTREE_without3rd_supports)>=95 & tree_all$IQTREE_without3rd_supports!="", 
                                "YES", "NO")

#add a YES/NO variable for nodes with conflict
tree_all$conflict<-ifelse(tree_all$BEAST6conflicts!="" | tree_all$IQTREE_with3rd_conflicts!="" | tree_all$IQTREE_without3rd_conflicts!="", "YES", "NO")

### add fossil calibrations
#WARNING: revise node numbers for each backbone tree
#add fossil details to matching nodes of "tree_all" data
tree_all$Fossil<-fossils$Fossil[match(tree_all$node,fossils$Fossil_node)]
tree_all$Fossil_number<-fossils$Fossil_number[match(tree_all$node,fossils$Fossil_node)]

### add RASP probabilities
#1) in imported RASP tables rename probability columns to include: analysis ID (e.g. "Beast8_estimated_equal_prob") and realm ID (original colname)
for (i in 1:length(RASP.lst)) 
  {
  analysisID<-RASP_IDs_refTopo[i]
  newColNames<-paste(analysisID,colnames(RASP.lst[[i]][4:13]), sep = "_") 
  colnames(RASP.lst[[i]])[4:13]<-newColNames
}
#rename rasp node ("node") to "raspNode"
colnames(RASP.lst[[1]])[1]<-"raspNode" 

#2)add corresponding node numbers from tree0_phylo (~"BEAST8_1-10_40kresampled.tree") to RASP lists 
for (j in 1:length(RASP.lst)) 
  {
  for (i in 1:nrow(RASP.lst[[j]])) 
    {
    tip1<-as.character(RASP.lst[[j]][i,"tip1ID"]) #tip1 from RASP
    tip2<-as.character(RASP.lst[[j]][i,"tip2ID"]) #tip2 from RASP
    RASP.lst[[j]][i,"treeNode"]<-getMRCA(tree0_phylo,c(tip1,tip2)) #
    }
  }


#3) add rasp node numbers to tree_all
tree_all$raspNode<-RASP.lst[[1]]$raspNode[match(tree_all$node,RASP.lst[[1]]$treeNode)]

### add rasp probabilities to tree_all
for(j in 1:length(RASP.lst)) 
  {
  Realms<-colnames(RASP.lst[[j]][4:13])
  for(i in 1:length(Realms)) 
    {
    Realm<-Realms[i]
    tree_all[,Realm]<-RASP.lst[[j]][match(tree_all$node,RASP.lst[[j]]$treeNode),Realm]
    }
  }

#add realm information to Kalotermitidae tips from "KaloTable"
tree_all[,"Realm"]<-KaloTable[match(tree_all$label,KaloTable$TipLabels),"Biogeographic_realm"] #add column "Realm" with a value for each tip (i.e. row with tip label in "label" column)

#fill probability values for Kalotermitidae tips: 100 for the corresponding realm probability, 0 for all other realm probabilities
Realms<-as.character(subset(unique(tree_all$Realm), !is.na(unique(tree_all$Realm)) & unique(tree_all$Realm)!="EMPTY")) #get vector of realms as unique non-NA non-"EMPTY" values in tree_all$Realm column, enforce "character" type, otherwise automatically factor with levels including EMPTY => mess downstream
ProbColNames<-names(tree_all)[grepl("estimated|fixed",names(tree_all) )] #vector of all colnames with rasp probabilities (grepped as containing "fixed" or "estimated" string)
for(i in 1:length(Realms)) {
  Realm<-Realms[i]
  tip_and_realm.vec<-!is.na(tree_all$label) & tree_all$Realm==Realm 
  tree_all[tip_and_realm.vec & !is.na(tip_and_realm.vec),grepl(Realm, names(tree_all))]<-100 # assign probability "100" 
  tree_all[tip_and_realm.vec & !is.na(tip_and_realm.vec),!grepl(Realm, names(tree_all)) & names(tree_all) %in% ProbColNames]<-0 # assign probability "0" 
}

###add mean and SD RASP probabilities
for(i in 1:length(Realms)) {
  Realm<-Realms[i]  
  MeanColName<-paste("Means",Realm,sep="_") #name for the column of means for probs of a particular realm
  GreplPattern<-paste("prob_",Realm,sep="") #search pattern for column names with RASP probabilities
  Realm.cols<-grepl(GreplPattern, names(tree_all)) #TRUE/FALSE vector of columns with probs for the particular realm
  tree_all[,MeanColName]<-rowMeans(tree_all[,Realm.cols]) #calculate and add to tree_all rowmeans for each row of probs for a particular realm
  Realm.cols<-grepl(GreplPattern, names(tree_all)) #TRUE/FALSE vector of columns with probs for the particular realm, re-initialize because "tree_all" has an additional column with means
  SDColName<-paste("SDs",Realm,sep="_") #name for the column of SDs for probs of a particular realm
  tree_all[,SDColName]<-apply(tree_all[,Realm.cols],1,sd) #add SDs for the particular realm
}

### convert final annotated tree file (tibble) to treedata object class S4 (for ggtree plotting and tip dropping)
tree_all_treedata<-as.treedata(tree_all) 

### drop tips below 2MY distance and same realm/country
#calculatedistances between tips
tree_all_phylo<-as.phylo(tree_all_treedata) #convert treedata to phylo object
dists<-distTips(tree_all_phylo, tips = "all", method = "patristic", useC = TRUE) #calculate patristic distances
dists.mat<-as.matrix(dists) #convert to matrix
#TEST: check distance of two selected 
  dists.mat["MAL39_Cryptotermes_domesticus","CRYDOMESTI_Cryptotermes_domesticus"] 
#check tips which are less than 2M apart (sum of branch lengths - patristic distance)
below2My<-list()
for (i in 1:length(tree_all_phylo$tip.label)) {
  tipLabel<-tree_all_phylo$tip.label[[i]]
  column<-dists.mat[,i]
  aboveThresh<-column[column<2 & column!=0]
  if(length(aboveThresh)!=0){
    below2My[[tipLabel]]<-column[column<2 & column!=0] #named vector of tips within threshold distance
  }
}
#visually inspect, make list of 36 tips to be dropped (below 2MY distance and same realm/country)
  #=> Below2MY-tips-to-drop.txt
tree_all_treedata_drop<-drop.tip(tree_all_treedata, tips2drop)
  #=> dropped 36 tips (from 232 to 196 tips, total 391 nodes)
tree_all_tbl_drop<-as_tibble(tree_all_treedata_drop) #convert treedata to tbl for future data manipulations
tree_all_treedata_drop<-as.treedata(tree_all_tbl_drop)
tree_all_tbl<-as_tibble(tree_all_treedata) #convert treedata to tbl for future data manipulations
#generate alternative treedata with renamed tip labels excluding collection codes
tree_all_treedata_drop_simpleLab<-rename_taxa(tree_all_treedata_drop,data = KaloTable[,c("TipLabels","TipLabels2")],key = "TipLabels",value = "TipLabels2")

########################
#stats for summary tree
#########################
#drop non-Kalotermitids
NonKalos_lines<-KaloTable %>% filter(Family=="OUTGROUP") %>% select(TipLabels) %>% pull() #non-Kalo tip labels
tree_all_td_drop_KaloOnly <- drop.tip(tree_all_treedata_drop,NonKalos_lines)
tree_all_phylo_drop_KaloOnly <- as.phylo(tree_all_td_drop_KaloOnly) #convert to phylo
tree_all_tbl_drop_KaloOnly<-as_tibble(tree_all_td_drop_KaloOnly) #convert to tibble

#define input trees
inTree_tbl<-tree_all_tbl_drop
inTree_treedata<-tree_all_treedata_drop
inTree_tbl_kaloonly<-tree_all_tbl_drop_KaloOnly
#plot tree for visual inspection
ggtree(inTree_treedata,ladderize=TRUE, right=TRUE)+
  geom_tiplab() + #show all tip labels
  geom_text2(aes(label=node), hjust=-.3, size=7,color="red") + #show internal node number
  xlim_tree(250)  #add space for tip labels

#B)count highly-supported congruent, congruent and conflicting branches
nrow(subset(inTree_tbl, support_by_all=="YES"))
nrow(subset(inTree_tbl, subset=conflict=="NO" & support_by_all=="NO"))
nrow(subset(inTree_tbl, subset=conflict=="NO"))
  #=>168 nodes without conflict for tree with non-kalos
nrow(subset(inTree_tbl_kaloonly, subset=conflict=="NO"))
  #=>119 nodes without conflict for kalo-only tree
nrow(subset(inTree_tbl, subset=conflict=="YES" & support_by_all=="NO"))
  #=>27 nodes with conflict for tree with non-kalos
nrow(subset(inTree_tbl_kaloonly, subset=conflict=="YES" & support_by_all=="NO"))
  #=>18 nodes with conflict for tree with non-kalos

##BEAST6 supports vs. conflicts
nrow(subset(inTree_tbl,is.na(inTree_tbl$label) & inTree_tbl$BEAST6supports!="" & !is.na(inTree_tbl$raspNode)))
#=>124 nodes of Kalotermitidae subtree supported
nrow(subset(inTree_tbl,is.na(inTree_tbl$label) & inTree_tbl$BEAST6conflicts!="" & !is.na(inTree_tbl$raspNode)))
#=>13 nodes of Kalotermitidae subtree with conflict

###check mean+SD of RASP analysis for each tree
#define function
SDs_of_RASP_probs<-function(TreeID){ 
  #1)import reformatted RASP tables as named list
  RASP_analysis_types<-c("_estimated_equal","_estimated_gamma","_fixed_equal","_fixed_gamma")
  RASP_files<-paste(inDir,"reformated-RASP/",TreeID,RASP_analysis_types,".txt", sep="") 
  RASP_IDs<-paste(TreeID,RASP_analysis_types,sep="")
  RASP.lst<-lapply(RASP_files, function(x) read_tsv(x))
  #rename probability columns to include: analysis ID and realm ID
  for (i in 1:length(RASP.lst)) {
    analysisID<-RASP_IDs[i]
    newColNames<-paste(analysisID,colnames(RASP.lst[[i]][4:13]), sep = "_") 
    colnames(RASP.lst[[i]])[4:13]<-newColNames
  }
  #join lists into a dataframe
  RASP_joined<-RASP.lst %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="node"),.)
  #calculate mean+SD for each RASP realm probability
  for(i in 1:length(Realms)) {
    Realm<-Realms[i]  
    Realm.cols<-grepl(Realm, names(RASP_joined)) #TRUE/FALSE vector of columns with probs for the particular realm
    MeanColName<-paste("Means",Realm,sep="_") #name for the column of means for probs of a particular realm
    RASP_joined[,MeanColName]<-rowMeans(RASP_joined[,Realm.cols]) #calculate and add to tree_all rowmeans for each row of probs for a particular realm
    Realm.cols<-grepl(Realm, names(RASP_joined)) #TRUE/FALSE vector of columns with probs for the particular realm, re-initialize because "tree_all" has an additional column with means
    SDColName<-paste("SDs",Realm,sep="_") #name for the column of SDs for probs of a particular realm
    RASP_joined[,SDColName]<-apply(RASP_joined[,Realm.cols],1,sd) #add SDs for the particular realm
  }
  write.xlsx(RASP_joined, paste("RASP_joined_",TreeID,".xls",sep="")) #write merged table with probabilities and sumstats into excel table
  perRealm<-apply(RASP_joined[,grepl("SD",colnames(RASP_joined))],2, max, na.rm=TRUE) #max SD of probabilities per realm for each RASP analysis parameter
  amongRealms<-max(perRealm) #max SD of probabilities among all realms 
  output<-list(perRealm,amongRealms) #merge two dataframes into "output" list
  return(output)
}
#run function for each tree topology
SDs_of_RASP_probs("Beast8")
  #=>16%
SDs_of_RASP_probs("Beast6")
  #=>13%
SDs_of_RASP_probs("IQtree_no3rd")
  #=>27%
SDs_of_RASP_probs("IQtree3rd")
  #=>17%
  ##=> TableS2-5

##################
#plot summary tree
##################
summaryTreePlot<-function(RASPrunID, inData)
{
  ggtree(inData,ladderize=TRUE, right=TRUE)+ #ladderize with smallest clade on top
  geom_tiplab(offset = 2) + #show all tip labels, adjust position to left
  #TODO: shift tip labels right, keep left-alignment
  geom_range("height_0.95_HPD", color='blue', size=2, alpha=.6) + # add branch height (~age) credible intervals (95% HPD)
  geom_nodepoint(aes(shape=">=0.95 in all trees",subset=support_by_all=="YES"), size=12, color="black", stroke=1.5) + #draw node point when "tree_all$support_by_all" is "YES"  
  geom_nodepoint(aes(shape="<0.95 AND no conflict among trees",subset=conflict=="NO" & support_by_all=="NO"), size=12, color="black", stroke=1.5) +  #draw node point when "tree_all$conflict" is "NO"
  geom_nodepoint(aes(shape="conflict in at least one tree",subset=conflict=="YES" & support_by_all=="NO"), size=12, color="black", stroke=1.5) +   #draw node point when "tree_all$conflict" is "YES"
  geom_label2(aes(label=Fossil_number,fill=paste(Fossil_number,":",Fossil),subset=(Fossil_number!="")),vjust=-0.4, hjust=-0.4) + #add fossil labels (subset only nodes with fossil info, otherwise legend will include NA value), vjust+hjust for imperfect simulation of repel
  #geom_text2(aes(label=node, subset=!isTip), hjust=-.3) + #show internal node number
  xlim_tree(80) + #add space for tip labels
  theme_tree2(legend.position=c(0.2,0.3))+ #"theme_tree2" to add timescale, define legend parametres
  geom_vline(xintercept = seq(-140, -10, by = 10), linetype="dotted", alpha=.5)+ #add vertical lines
  scale_x_continuous(labels = abs, breaks=seq(-140, 0, by = 10))+ #define x-axis ticks every 10MYA
  scale_shape_manual(name="Branch support:",breaks=c(">=0.95 in all trees", "<0.95 AND no conflict among trees", "conflict in at least one tree"), values = c(">=0.95 in all trees" = 1, "<0.95 AND no conflict among trees" = 5, "conflict in at least one tree" = 4))+  #add custom legend name, sort maually values (via breaks) 
  scale_fill_discrete(name="Fossil calibrations:", breaks = mixedsort(paste(fossils$Fossil_number,":",fossils$Fossil))) + #manual legend name, sort legend items as mixed numerics and characters 
  ggtitle(paste ("RASP data mapped:",RASPrunID))
}

### plot tree without bootstraps, fossils
summaryTreePlotSimple<-function(RASPrunID, inData)
{
  ggtree(inData,ladderize=TRUE, right=TRUE)+ #ladderize with smallest clade on top
    geom_tiplab(offset = 2) + #show all tip labels, adjust position to left
    #TODO: shift tip labels right, keep left-alignment
    geom_range("height_0.95_HPD", color='blue', size=2, alpha=.6) + # add branch height (~age) credible intervals (95% HPD)
    xlim_tree(80) + #add space for tip labels
    theme_tree2(legend.position=c(0.2,0.3))+ #"theme_tree2" to add timescale, define legend parametres
    geom_vline(xintercept = seq(-140, -10, by = 10), linetype="dotted", alpha=.5)+ #add vertical lines
    scale_x_continuous(labels = abs, breaks=seq(-140, 0, by = 10))+ #define x-axis ticks every 10MYA
    ggtitle(paste ("RASP data mapped:",RASPrunID))
}

p1<-summaryTreePlot("Beast8_means_of_4xRASP",tree_all_treedata_drop_simpleLab)
  #=> Figure 1 (simple labels, dropped <2MY tips)

##uncomment below to plot instead Figure1 with collection codes in tip labels
#p1<-summaryTreePlot("Beast8_means_of_4xRASP",tree_all_treedata_drop)

##uncomment below to plot instead full tree with all samples including replicated species
#p1<-summaryTreePlot("Beast8_means_of_4xRASP",tree_all_treedata) 
#revts(p1) 
  #=>Figure S1

##uncomment below to plot instead a simple tree without bootstraps, fossils...
#p1<-summaryTreePlotSimple("Beast6_means_of_4xRASP",tree_all_treedata) 

#reverse order of time scale, make absolute numbers
p2 <- revts(p1)  

#add internal node probabilities as pie charts using ggtree::inset (see https://guangchuangyu.github.io/software/ggtree/vignettes/ggtree-inset.html#annotate-with-other-types-of-charts and http://www.randigriffin.com/2017/05/11/primate-phylogeny-ggtree.html)
tree_rasp_only<-tree_all_tbl_drop[!is.na(tree_all_tbl_drop$raspNode),] #rows with rasp data
#tree_rasp_only<-tree_all_tbl[!is.na(tree_all_tbl$raspNode),] #rows with rasp data
#CHECK: how many internal nodes with data
  nrow(tree_rasp_only)
  #=> 137 nodes
#extract colnames corresponding to mean RASP probabilities 
SelColNames<-colnames(tree_rasp_only[,grepl("Means_", colnames(tree_rasp_only))]) 
#manual colors for realms (corresponding to realmColors1 but without realm names)
realmColors<-c("khaki4", "lightsteelblue4", "maroon2" , "lightblue2", "firebrick1", "palevioletred2", "black", "lemonchiffon2", "blue1", "coral" ) 
#generate node pie charts for columns extracted above
pies <- nodepie(tree_rasp_only, cols=SelColNames,color=realmColors)
  #=> 137 pies
p3<-ggtree::inset(p2, pies, width=0.06,height=0.06,hjust=0.17,vjust=0.1)
#add tip "probabilities" as another set of pie charts 
tree_tips<-tree_all_tbl_drop[!is.na(tree_all_tbl_drop$label) & !is.na(tree_all_tbl_drop$Realm),] #tip rows with realm info
#tree_tips<-tree_all_tbl[!is.na(tree_all_tbl$label) & !is.na(tree_all_tbl$Realm),] #tip rows with realm info
pies2 <- nodepie(tree_tips, cols=SelColNames,color=realmColors) #generate node pie charts for columns extracted above
#plot everything
p4<-ggtree::inset(p3, pies2, width=0.16, height=0.16, hjust = -1.5)
#geom_text2(aes(label=node), hjust=-.3, size=7,color="red")  # uncomment and add to ggtree to show node nums
p4
###export pdf: 40x30inch
  #=> Figure 1

###############################################
#polytomise branches without consensus support
###############################################
#workaround: add node labels to phylo object, then polytomise nodes with conflicts
#nodes with conflict among analyses <-0, without conflict <- 100
tree_all_tbl_drop_relabel<-tree_all_tbl_drop %>% 
  mutate(support_by_all=replace(support_by_all,conflict=="YES",0)) %>%
  mutate(support_by_all=replace(support_by_all,conflict=="NO",100))  
#convert tbl to phylo
tree_all_phylo_drop_relabel <- as.phylo(tree_all_tbl_drop_relabel)
#node.label vector excluding NA values (~tips)
nodelabels_noNAs <- tree_all_tbl_drop_relabel %>% filter(!is.na(support_by_all)) %>% select(support_by_all) %>% pull()
#add non-tip labels as "node.label" to phylo tree 
tree_all_phylo_drop_relabel$node.label<-nodelabels_noNAs 
  #=> 195 node labels
write.tree(tree_all_phylo_drop_relabel,"tree_all_phylo_drop_relabel.newick")
#polytomise nodes with conflict
tree_all_phylo_drop_relabel_polytom<-as.polytomy(tree_all_phylo_drop_relabel, feature='node.label', fun=function(x) as.numeric(x)<100)
#CHECK: plot as cladogram
    ggtree(tree_all_phylo_drop_relabel_polytom,ladderize=TRUE, right=TRUE, branch.length = "none")+
      geom_tiplab()+
      xlim_tree(50)
#export newick
write.tree(tree_all_phylo_drop_relabel_polytom, "tree_all_phylo_drop_relabel_polytom.nwck")
#...colapse at genus level and relabel in figtree
  #=> Figure3

########################
#Biogeographic realm map
########################
#import (needs whole unzipped directory downloaded from https://macroecology.ku.dk/resources/wallace/cmec_regions___realms.zip)
cmec.regions <- sf::read_sf("HoltRealms/CMEC regions & realms/newRealms.shp")

#CHECK: default R color scheme (10 values)
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  n = 10
  cols = gg_color_hue(n)
  plot(1:n, pch = 16, cex = 2, col = cols)

#custom colors (default R colors and Panamian realm with the same color as Neotropic)
realmColors1 <- c("Afrotropical" = "khaki4", "Australian" = "lightsteelblue4", "Madagascan" = "maroon2" , "Nearctic" = "lightblue2", "Neotropical" = "firebrick1", "Oceanina" = "palevioletred2", "Oriental" = "black", "Palearctic" = "lemonchiffon2", "Panamanian" = "firebrick1", "Saharo-Arabian" = "blue1", "Sino-Japanese" = "coral" )  

# plot map
ggplot() +
  geom_sf(data = cmec.regions, aes(fill = Realm), size = .2, col = 0, alpha=1)+
  coord_sf(expand = FALSE)+
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0), panel.background = element_rect(fill = "white"))+
  scale_fill_manual(breaks=cmec.regions$Realm, values=realmColors1)  #custom colors (default R colors and Panamian realm with the same color as Neotropic)
  #=> Figure1

##################################################################
#traverse Kalo-subtree to get tip-to-root node RASP probabilities
##################################################################
#drop non-Kalo tips
NonKalos_lines<-KaloTable %>% filter(Family=="OUTGROUP") %>% select(TipLabels) %>% pull() #non-Kalot tip labels
tree_all_td_drop_KaloOnly <- drop.tip(tree_all_treedata_drop,NonKalos_lines)
tree_all_phylo_drop_KaloOnly <- as.phylo(tree_all_td_drop_KaloOnly) #convert to phylo
tree_all_tbl_drop_KaloOnly<-as_tibble(tree_all_td_drop_KaloOnly) #convert to tibble
#TEST
  nrow(tree_all_tbl_drop_KaloOnly)
  #=> 275 rows

##define intree variables for plotting below
#plotting of all >2MY Kalos
InTree_tbl<-tree_all_tbl_drop_KaloOnly
InTree_phy<-tree_all_phylo_drop_KaloOnly
InTree_td<-tree_all_td_drop_KaloOnly

#TEST: inspect Kalo subtree with node numbers
ggtree(InTree_td,ladderize=TRUE, right=TRUE)+ #ladderize with smallest clade on top
  geom_tiplab(offset = 2)+
  geom_text2(aes(label=node), hjust=-.3, size=7,color="red") + #show internal node number
  xlim_tree(250)  #add space for tip labels

#get numbers of tips in Kalo subtree
KaloTerminalNodeNums<-InTree_tbl %>% filter(!is.na(label)) %>% dplyr::select(node) %>% pull()
#CHECK 
  #length(KaloTerminalNodeNums)  
  #=> 138

#vector of realms in order matching the per-realm-RASP-mean columns by inverting order in Realms
Realms2<-rev(Realms) 

#function to generate vector of tip-to-root node numbers
GetParents<-function(NodeNum,inData)
  {
  HasParent<-1 #initialize object
  ParentsVector<-NodeNum
  while(HasParent!=0) 
    {
      #inData<-InTree_tbl #loop testing
      #NodeNum<-201 #loop testing
    ParentNodeInfo<-tidytree::parent(inData,NodeNum) 
    HasParent<-nrow(ParentNodeInfo)
    NodeNum<-ParentNodeInfo$node
    ParentsVector<-append(ParentsVector,NodeNum)
      #print(HasParent) #loop testing
    }
  ParentsVector
  } 

#TEST: GetParents with test terminal node "85"
  #GetParents(85,InTree_tbl)

#run GetParents with all terminal nodes
ParentsVectors<-lapply(KaloTerminalNodeNums,GetParents, inData=InTree_tbl)

#custom function to traverse a vector of nodes and generate node realm info for plotting
GetNodeInfo<-function(ParentsVector)
{
  ProbThresh<-90 #set RASP probability threshold for "AMBIGUOUS"
  TipsToRoot<-data.frame() #initialize empty output dataframe
  #loop1: through the vector of node numbers to generate TipsToRoot dataframe
  for (i in 1:length(ParentsVector)) 
  {
    Node<-ParentsVector[i]
    NodeInfo<-InTree_tbl%>%filter(node==Node)
    maxRASPnodeRealm<-Realms2[which.max(NodeInfo[,c(89,87,85,83,81,79,77,75,73,71)])] #realm name with max RASP probability averaged across 4 analyses
    Age<-InTree_tbl %>% filter(node==Node) %>% dplyr::select(height) %>% as.numeric() #is ~1e-14 for terminal nodes
    AgeHPD095<-InTree_tbl %>% filter(node==Node) %>% dplyr::select(height_0.95_HPD)
    AgeHPD095<-AgeHPD095$height_0.95_HPD[[1]] #convert from a multiclass object to numeric
    LowerAgeHPD095<-AgeHPD095[1]
    UpperAgeHPD095<-AgeHPD095[2]
    #build "TipsToRoot" dataframe
    TipsToRoot[i,"Node"]<-Node
    TipsToRoot[i,"maxRASPnodeValue"]<-max(NodeInfo[,c(89,87,85,83,81,79,77,75,73,71)]) #highest RASP probability among nodes
    TipsToRoot[i,"Oceanian"]<-max(NodeInfo[,89]) #RASP probability for realm
    TipsToRoot[i,"Madagascan"]<-max(NodeInfo[,87]) #RASP probability for realm
    TipsToRoot[i,"Palearctic"]<-max(NodeInfo[,85]) #RASP probability for realm
    TipsToRoot[i,"Saharo_Arabian"]<-max(NodeInfo[,83]) #RASP probability for realm
    TipsToRoot[i,"Sino_Japanese"]<-max(NodeInfo[,81]) #RASP probability for realm
    TipsToRoot[i,"Nearctic"]<-max(NodeInfo[,79]) #RASP probability for realm
    TipsToRoot[i,"Oriental"]<-max(NodeInfo[,77]) #RASP probability for realm
    TipsToRoot[i,"Neotropical"]<-max(NodeInfo[,75]) #RASP probability for realm
    TipsToRoot[i,"Afrotropical"]<-max(NodeInfo[,73]) #RASP probability for realm
    TipsToRoot[i,"Australian"]<-max(NodeInfo[,71]) #RASP probability for realm
    TipsToRoot[i,"maxRASPnodeRealm"]<-maxRASPnodeRealm
    if(i==1) #for first node (~tip), set by definition Age and intervals as 0
    {
      TipsToRoot[i,"Age"]<-0
      TipsToRoot[i,"LowerAgeHPD095"]<-0
      TipsToRoot[i,"UpperAgeHPD095"]<-0
    }
    else
    {
    TipsToRoot[i,"Age"]<-Age
    TipsToRoot[i,"LowerAgeHPD095"]<-LowerAgeHPD095
    TipsToRoot[i,"UpperAgeHPD095"]<-UpperAgeHPD095
    }
    #asign "maxRASPnodeRealm" if realm above probability threshold, otherwise "AMBIGUOUS"
    if(TipsToRoot[i,"maxRASPnodeValue"]>ProbThresh)
    {
      RealmAboveThreshold<-maxRASPnodeRealm
    }
    else
    {
      RealmAboveThreshold<-"AMBIGUOUS"
    }
    TipsToRoot[i,"RealmAboveThreshold"]<-RealmAboveThreshold
  }
  #loop2: through TipsToRoot for the age of last non-AMBIGUOUS node with realm different than tip realm (i.e. max estimate of dispersal to extant realm)
  #1) add tip info to "TipsToRoot"
  TipNodeInfo<-InTree_tbl%>%filter(node==ParentsVector[1])
  TipRealm<-Realms2[which.max(TipNodeInfo[,c(89,87,85,83,81,79,77,75,73,71)])]
  TipLabel<-TipNodeInfo %>% select(label) %>% as.character()
  TipsToRoot[,"TipLabel"]<-TipLabel
  TipsToRoot[,"TipRealm"]<-TipRealm
  #2) loop, starting with first internal node
  LastPreDispAge<-0 #set initial age of oldest predispersal node to 0 
  LowerLastPreDispAge<-0 #set initial age of oldest predispersal node to 0 
  UpperLastPreDispAge<-0 #set initial age of oldest predispersal node to 0 
  j<-2 #set loop counter to initial value
  for (i in 2:nrow(TipsToRoot))
  {
    if (TipsToRoot[i,TipRealm]<(100-ProbThresh)) # if current internal node realm is identical to the tip realm with below (100-ProbThresh) probability, save the internal node age as "LastDispAge" and exit loop
    {
      LastDispAge<-TipsToRoot[i,"Age"] 
      LowerLastDispAge<-TipsToRoot[i,"LowerAgeHPD095"]
      UpperLastDispAge<-TipsToRoot[i,"UpperAgeHPD095"]
      break()
    }
    else if (TipsToRoot[i,TipRealm]>ProbThresh & i==j) # ...else if current internal node realm is identical to the tip realm with above ProbThresh probability AND it is a part of an uninterrupted string of such nodes, save the internal node age as "LastPreDispAge"
      {
      LastPreDispAge<-TipsToRoot[i,"Age"] 
      LowerLastPreDispAge<-TipsToRoot[i,"LowerAgeHPD095"]
      UpperLastPreDispAge<-TipsToRoot[i,"UpperAgeHPD095"]
      j<-j+1
    }
    else if (identical(setdiff(TipsToRoot$RealmAboveThreshold,TipRealm),"AMBIGUOUS")) #...else if all node realms which are not identical to tip realm are  AMBIGUOUS, save the age of root (~81 MY) "LastDispAge" and exit loop (to be replaced with question marks in final figure)
    {
      LastDispAge<-max(InTree_tbl$height) 
      LowerLastDispAge<-max(InTree_tbl$height)
      UpperLastDispAge<-max(InTree_tbl$height)
      break()
    }
  }
  #assign tip info to each node descending from the tip
  TipsToRoot[,"LastDispAge"] <- LastDispAge 
  TipsToRoot[,"LastPreDispAge"] <- LastPreDispAge 
  TipsToRoot[,"LowerLastDispAge"] <- LowerLastDispAge 
  TipsToRoot[,"UpperLastDispAge"] <- UpperLastDispAge 
  TipsToRoot[,"LowerLastPreDispAge"] <- LowerLastPreDispAge 
  TipsToRoot[,"UpperLastPreDispAge"] <- UpperLastPreDispAge 
  TipsToRoot #output final dataframe
}   


#CHECK: with single vector of node numbers
  #ParentsVector<-ParentsVectors[[85]]
  #GetNodeInfo(ParentsVector)
  #=>OK
  #ParentsVector<-ParentsVectors[[1]]
  #GetNodeInfo(ParentsVector)
  #=>OK, LastDispAges ~ age of first AMBIGUOUS node
  
#make list of node infos for each tip %>% bind into one df  
TipsToRoots<-lapply(ParentsVectors,GetNodeInfo) %>% bind_rows()
#add columns from KaloTable: numerical groups for ordering y-axis in plot, "RevisedGenus", and "Genus"
TipsToRootsGroups<-left_join(x=TipsToRoots,y=KaloTable[,c("TipLabels","GroupingGenus","GroupingRealm","GenusRevised","Genus")], by=c("TipLabel"="TipLabels"))

#color palette as named vector
realmColors4 <- c("Neotropical" = "firebrick1", 
                  "Oriental" = "black", 
                  "Australian" = "lightsteelblue4",
                  "Afrotropical" = "khaki4",
                  "Madagascan" = "maroon2" , 
                  "Nearctic" = "lightblue2", 
                  "Oceanian" = "palevioletred2", 
                  "Saharo_Arabian" = "blue1", #with underscore
                  "Sino_Japanese" = "coral", #with underscore
                  "Palearctic" = "lemonchiffon2",
                  "AMBIGUOUS"="light grey" )  

##drop "biogeographically uninformative" tips: i.e. sister species belonging to the same genus, same revised genus, and same realm
#using a manually created list of "redundant" tips: i.e., when multiple tips belong to the same revised genus, same biogeographic realm, and their last common ancestor is inferred to share the realm with them, then select all but one representative tip
#CHECK: number of tips to drop
  #length(RedundantKalos)
#remove rows corresponding to "uninformative" tips
`%nin%` = Negate(`%in%`) #define not-in function
TipsToRootsGroups_dropped<-TipsToRootsGroups %>% 
  filter(TipLabel %nin% RedundantKalos) %>%
  drop.levels() #drop unused factor levels from before subsetting
  #CHECK: surviving tips
  #nrow(subset(TipsToRootsGroups,Age==0))
  #nrow(subset(TipsToRootsGroups_dropped,Age==0))

##plotting function for interval between dispersal node and predispersal non-foreign node
DisjunctPlot3<-function(inData) 
{ inData<-inData %>%
  filter(Age==LastDispAge | Age==LastPreDispAge) %>% #filter nodes with "LastDispAge" or "LasPreDispAge" age
  mutate(RealmAboveThreshold=factor(RealmAboveThreshold, levels = names(realmColors4))) %>% #set levels of "RealmAboveThreshold" according to "realmColors4" vector
  mutate(TipRealm=factor(TipRealm, levels = names(realmColors4))) #set levels of "TipRealm" according to "realmColors4" vector to order facets
#reorder factor levels of column "RealmAboveThreshold" (to have AMBIGUOUS as last in plot legend ~ same order of names in named vector "realmColors4")
inData<-inData[order(inData$GroupingRealm,inData$LastDispAge),] #order dataframe by first GroupingRealm, then LastDispAge
inData$TipLabel <- fct_inorder(inData$TipLabel) %>% fct_rev() #order factor levels by dataframe order
inData %>%
  ggplot(aes(x=Age,y=TipLabel)) + 
  geom_errorbar(aes(x=LastDispAge, y=TipLabel,xmin=LowerLastPreDispAge,xmax=UpperLastDispAge),size=2, width=0)+  #error bars as combined lower 95HPD range estimate for pre-dispersal node and higher 95HPD range estimate for dispersal node
  geom_point(data=subset(inData,RealmAboveThreshold!="AMBIGUOUS"),aes(color=RealmAboveThreshold),size=4)+ #draw age points for nodes with non-ambiguous realms
  geom_point(data=subset(inData,RealmAboveThreshold=="AMBIGUOUS"),fill="white",shape=21,size=4)+ #draw age points for nodes with non-ambiguous realms
  #geom_path(aes(x=LastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=1,alpha=0.9)+ #draw connecting line between "LastDispAge" nodes (use geom_path instead of geom_line, use group(1) to connect points when one variable (TipLabel) is factor)
  #geom_path(aes(x=LowerLastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=0.5,alpha=0.9)+ #line for lower 95%HPDage
  #geom_path(aes(x=UpperLastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=0.5,alpha=0.9)+ #line for upper 95%HPDage
  scale_x_reverse(breaks = seq(95,0,by=-5),limits=c(95,0))+ #reverse time scale and set manual limits
  scale_y_discrete(position = "right")+
  scale_colour_manual(values=realmColors4)+ #set manual color scale for realms to match other graphics
  ggthemes::theme_few()+ #frame around facets
  theme(plot.caption = element_text(hjust = 0), #set title alignment to rightmost, remove grey bckgroung and gridlines, remove facet background and add outline
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.05, color="grey" ) 
  )+
  ggtitle("Kalos\n>2MY distant\nSorted by realms of extant species\nWithin realms sorted by max age of last disjunction with >90% RASP probabl."  ) +
  facet_grid(rows=vars(TipRealm),scales = "free", space = "free") #each realm into separate row-facet
}

#plot, without uninformative groups
DisjunctPlot3(TipsToRootsGroups_dropped)
  #export 10x20inch => Figure2

##plotting function for Last dispersal age nodes only
DisjunctPlot<-function(inData) 
{ inData<-inData %>%
  filter(Age==LastDispAge) %>% #filter only nodes of last high-probability dispersals and all nodes for lineages with AMBIGUOUS last dispersals 
  mutate(RealmAboveThreshold=factor(RealmAboveThreshold, levels = names(realmColors4))) %>% #set levels of "RealmAboveThreshold" according to "realmColors4" vector
  mutate(TipRealm=factor(TipRealm, levels = names(realmColors4))) #set levels of "TipRealm" according to "realmColors4" vector to order facets
#reorder factor levels of column "RealmAboveThreshold" (to have AMBIGUOUS as last in plot legend ~ same order of names in named vector "realmColors4")
inData<-inData[order(inData$GroupingRealm,inData$LastDispAge),] #order dataframe by first GroupingRealm, then LastDispAge
inData$TipLabel <- fct_inorder(inData$TipLabel) %>% fct_rev() #order factor levels by dataframe order
inData %>%
  ggplot(aes(Age,TipLabel)) + 
  geom_point(aes(color=RealmAboveThreshold),size=4)+ #draw age points
  geom_errorbar(data=subset(inData,RealmAboveThreshold!="AMBIGUOUS"),aes(color=RealmAboveThreshold,x=LastDispAge, y=TipLabel,xmin=LowerLastDispAge,xmax=UpperLastDispAge),size=2,alpha=0.5, width=0)+  #error bars for non-ambiguous source realm
  geom_errorbar(data=subset(inData,RealmAboveThreshold=="AMBIGUOUS"),aes(color=RealmAboveThreshold,x=LastDispAge, y=TipLabel,xmin=LowerLastDispAge,xmax=UpperLastDispAge),size=2,alpha=0.5, width=0,linetype="dotted")+ #dotted error bars for ambiguous source realm
  #geom_path(aes(x=LastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=1,alpha=0.9)+ #draw connecting line between "LastDispAge" nodes (use geom_path instead of geom_line, use group(1) to connect points when one variable (TipLabel) is factor)
  #geom_path(aes(x=LowerLastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=0.5,alpha=0.9)+ #line for lower 95%HPDage
  #geom_path(aes(x=UpperLastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=0.5,alpha=0.9)+ #line for upper 95%HPDage
  scale_x_reverse(breaks = seq(95,0,by=-5),limits=c(95,0))+ #reverse time scale and set manual limits
  scale_colour_manual(values=realmColors4)+ #set manual color scale for realms to match other graphics
  theme(plot.caption = element_text(hjust = 0), #set title alignment to rightmost, remove grey bckgroung and gridlines, remove facet background and add outline
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_rect(fill="white"))+
  ggthemes::theme_few()+
  ggtitle("Kalos\n>2MY distant\nSorted by realms of extant species\nWithin realms sorted by max age of last disjunction with >90% RASP probabl."  ) +
  facet_grid(rows=vars(TipRealm),scales = "free", space = "free") #each realm into separate row-facet
}

#plot, without "uninformative" groups
DisjunctPlot(TipsToRootsGroups_dropped)
#export 10x20inch

##plotting function for all nodes
DisjunctPlot2<-function(inData) 
{
  #reorder factor levels of column "RealmAboveThreshold" (to have AMBIGUOUS as last in plot legend ~ same order of names in named vector "realmColors4")
  inData<- inData %>%
    mutate(RealmAboveThreshold=factor(RealmAboveThreshold, levels = names(realmColors4))) %>% #set levels of "RealmAboveThreshold" according to "realmColors4" vector
    mutate(TipRealm=factor(TipRealm, levels = names(realmColors4))) #set levels of "TipRealm" according to "realmColors4" vector to order facets
  ##ggplot, y-axis ordered by realms of extant species
  inData<-inData[order(inData$GroupingRealm,inData$LastDispAge),] #order dataframe by first GroupingRealm, then LastDispAge
  inData$TipLabel <- fct_inorder(inData$TipLabel) %>% fct_rev() #order factor levels by dataframe order
  inData %>%
    filter(Age!=0) %>% #filter 
    ggplot(aes(Age,TipLabel)) + 
    geom_point(aes(color=RealmAboveThreshold),size=4)+ #draw age points
    geom_errorbar(aes(color=RealmAboveThreshold,x=Age, y=TipLabel,xmin=LowerAgeHPD095,xmax=UpperAgeHPD095),size=2,alpha=0.5, width=0)+  #error bars for non-ambiguous source realm
    geom_path(aes(x=LastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=1,alpha=0.9)+ #draw connecting line between "LastDispAge" nodes (use geom_path instead of geom_line, use group(1) to connect points when one variable (TipLabel) is factor)
    #geom_path(aes(x=LowerLastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=0.5,alpha=0.9)+ #line for lower 95%HPDage
    #geom_path(aes(x=UpperLastDispAge, y=TipLabel,group = RealmAboveThreshold),color="black",subset(inData,Age==0),size=0.5,alpha=0.9)+ #line for upper 95%HPDage
    scale_x_reverse(breaks = seq(95,0,by=-5),limits=c(95,0))+ #reverse time scale and set manual limits
    scale_colour_manual(values=realmColors4)+ #set manual color scale for realms to match other graphics
    ggthemes::theme_few()+
    theme(plot.caption = element_text(hjust = 0), #set title alignment to rightmost, remove grey bckgroung and gridlines, remove facet background and add outline
          plot.title.position = "plot",
          plot.caption.position =  "plot",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid.major.x = element_blank() ,
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line( size=.1, color="black" ) 
          )+
    ggtitle("Kalos\n>2MY distant\nSorted by realms of extant species\nWithin realms sorted by max age of last disjunction with >90% RASP probabl."  ) +
    facet_grid(rows=vars(TipRealm),scales = "free", space = "free") #each realm into separate row-facet
}
#plot
DisjunctPlot2(TipsToRootsGroups_dropped)
#export 10x15inch

#CHECK: plot before dropping uninformative groups
  #DisjunctPlot(TipsToRootsGroups)
  #DisjunctPlot2(TipsToRootsGroups)

##genus labels
#order dataframe to match the y-axis label order
inData<-TipsToRootsGroups_dropped
inData<-inData[order(inData$GroupingRealm,inData$LastDispAge),] #order dataframe by first GroupingRealm, then LastDispAge
TipsToRootsGroups_dropped_reordered<-inData
#shorten Genus names by removing "termes", shorten Bifidi+Allo+ROisini+Epicalo to "B+A+R+E"
TipsToRootsGroups_dropped_reordered_selected<-TipsToRootsGroups_dropped_reordered %>% 
  filter(Age==0) %>% 
  select(TipLabel,Genus,GenusRevised,TipRealm) %>% 
  mutate(GenusRevised=str_replace(GenusRevised,"Bifiditermes\\+Allotermes\\+Roisinitermes\\+Epicalotermes", "B\\+A\\+R\\+E")) %>%
  mutate(Genus=str_replace(Genus,"termes", "")) %>% 
  mutate(GenusRevised=str_replace(GenusRevised,"termes", "")) 
write.csv(TipsToRootsGroups_dropped_reordered_selected,file = "TipsToRootsGroups_dropped_reordered_selected.csv")

###############
#fossil ages plot
###############
#retrieve Kalotermitidae records from PaleoDB
PaleoDB<- PaleoDB %>% filter(genus!="UNKNOWN") #remove fossils with unassigned genus, does not remove the oldest Kalotermitidae fossil in any realm
PaleoDB<- PaleoDB %>% filter(genus!="X_Cratokalotermes")  %>% filter(genus!="X_Kalotermes") #remove Cratokalo and Kalotermes piacentinii
PaleoDB<-PaleoDB %>% mutate(Realm=factor(Realm, levels = names(realmColors4))) #set factor levels according to "realmColors4"
#CHECK
  nrow(PaleoDB)
  #=> 57 records
#plot fossil ages from paleodatabase
PaleoDB$avgAge<-(PaleoDB$min_ma+PaleoDB$max_ma)/2 #calculate average age
pos <- position_jitter(seed = 1)
PaleoDB %>%
  ggplot(aes(x=avgAge,y=Realm))+
  geom_jitter(aes(color=Realm),size=4,alpha=1,position = pos)+ # "shape=GenusStatus" to differentiate affiliations to extinct vs. extant genera
  #geom_errorbar(aes(xmin = min_ma, xmax = max_ma,color=Realm),alpha=0.3,size = 3, width=0)+ #bars for age intervals, does not work properly with jitter
  geom_text_repel(aes(label = genus),position = pos,size=3)+
  #theme_classic()+
  facet_grid(rows=vars(Realm),scales = "free", space = "free")+ #each realm into separate row-facet
  xlab("Fossil age")+
  theme_linedraw()+
  scale_x_reverse(breaks = seq(125,0,by=-5),limits=c(125,0))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  scale_colour_manual(values=realmColors4) #set manual color scale for realms to match other graphics
#export:10x15inch
  #=> FigureS4

##########################################################
#subsample Neotropical Kalos, (analyze with RASP), plot
##########################################################
# 3 subsamplings, 2 subsampling depths (dropping 23 and 30 Neotropical kalotermitids), single model (JC=fixed equal)
#all neotropic tips (45)
Neotropicals<-subset(tree_all$label, tree_all$Realm == "Neotropical")
#drop non-Kalo tips
tree0_KaloOnly<-drop.tip(tree0_treedata, NonKalos_lines)
  #=>174 tips
#define seeds for 3 random samplings
seeds3<-c(1234,5678,9101)

#drop random samples of Neotropical tips from kalotermitid trees (define number of taxa to keep, number of sampling replicaets to run, basename of output files)
subSampling <-function(NumOfDropped,FileName)
{
  for (i in seq(1:3)) #run loop 3-times to generate 3 replicates, each time with different but pre-defined seed for random sampling
  {
    set.seed(seeds3[i])
    treeFileName<-paste(FileName,"_",NumOfDropped,"dropped-sample",i,".tre",sep="")
    tiplistFileName<-paste(FileName,"_",NumOfDropped,"dropped-sample",i,".txt",sep="")
    to_drop<-sample(Neotropicals, NumOfDropped)
    tree_dropped<-drop.tip(tree0_KaloOnly, to_drop) #drop "NumOfDropped" number of Neotropical tips
    tree_dropped_phylo<-as.phylo(tree_dropped)
    write.tree(tree_dropped_phylo, treeFileName) #write tree files
    writeLines(to_drop, tiplistFileName) #write list of dropped tips
  }
}

outDir<-paste(inDir,"subsampled_Neotropics_trees/",sep="")
dir.create(outDir)
outName<-paste(outDir,"BEAST8",sep="")
#run subsampling to drop 23, and 30 Neotropical tips, respectively
subSampling(23,outName)
  #=> 151 tips
subSampling(30,outName)
  #=>144 tips

##(OUTSIDE OF R)
#1) subsample distribution table of all specimens to match the replicated random sampling
#2) run RASP with subsampled input trees and the subsampled distribution tables
#3) reformat RASP tables prior importing into R

#list of input tree filenames
InTrees<-paste(outDir,c("BEAST8_23dropped-sample1.tre",
                     "BEAST8_23dropped-sample2.tre",
                     "BEAST8_23dropped-sample3.tre",
                     "BEAST8_30dropped-sample1.tre",
                     "BEAST8_30dropped-sample2.tre",
                     "BEAST8_30dropped-sample3.tre"
                   ),sep=""
               )

#B) rename probability columns to include analysis ID (e.g. "Beast8_estimated_equal_prob") and realm ID (original colname)
inRASP<-RASP.lst_subsample #use RASP data imported at the begining of script file
for (i in 1:length(inRASP)) {
  analysisID<-RASP_IDs_subsample[i]
  newColNames<-paste(analysisID,colnames(inRASP[[i]][4:13]), sep = "_") 
  colnames(inRASP[[i]])[4:13]<-newColNames
}
#C) initialize list holding all plots
p<-list()

##run tree plotting in loop, filling the plotlist "p"
for (i in 1:length(InTrees)) {
  #TEST: i<-5
  #1)import tree
  InTree <- treeio::read.newick (InTrees[i])
  #2)convert tree object to tibble for data manipulation
  InTree_tbl<-as_tibble(InTree) 
  #4)convert tree to phylo object required for getMRCA(function)
  InTree_phylo<-as.phylo(InTree)
  #5)rename rasp node ("node") to "raspNode"
  colnames(inRASP[[i]])[1]<-"raspNode" 
  #6)add node numbers from original phylogenetic trees to RASP lists 
  for (j in 1:nrow(inRASP[[i]])) 
    {
    #TEST: j<-33
    #TEST: inRASP[[i]][j,"tip1ID"]
    #TEST: inRASP[[i]][j,"tip2ID"]
    tip1<-as.character(inRASP[[i]][j,"tip1ID"])
    tip2<-as.character(inRASP[[i]][j,"tip2ID"])
    inRASP[[i]][j,"treeNode"]<-getMRCA(InTree_phylo,c(tip1,tip2))
    }
  #7)add rasp node number to tree
  InTree_tbl$raspNode<-inRASP[[i]]$raspNode[match(InTree_tbl$node,inRASP[[i]]$treeNode)]
  #8)add rasp probabilities to tree
  Realms<-colnames(inRASP[[i]][4:13])
  for(j in 1:length(Realms)) {
    Realm<-Realms[j]
    InTree_tbl[,Realm]<-inRASP[[i]][match(InTree_tbl$node,inRASP[[i]]$treeNode),Realm]
  }
  #9)add realm information to Kalotermitidae tips from "KaloTable"
  InTree_tbl[,"Realm"]<-KaloTable[match(InTree_tbl$label,KaloTable$TipLabels),"Biogeographic_realm"] #add column "Realm" with a value for each tip (i.e. row with tip label in "label" column)
  
  #10)fill-in probability values for tips: 100 for the corresponding realm probability, 0 for all other realm probabilities
  Realms<-subset(unique(InTree_tbl$Realm), !is.na(unique(InTree_tbl$Realm))) #get vector of realms as unique non-NA values in InTree_tbl$Realm column
  ProbColNames<-names(InTree_tbl)[grepl("_prob_",names(InTree_tbl) )] #vector of all colnames with rasp probabilities (grepped as containing "_prob_" string)
  for(j in 1:length(Realms)) {
    Realm<-Realms[j]
    tip_and_realm.vec<-!is.na(InTree_tbl$label) & InTree_tbl$Realm==Realm
    InTree_tbl[tip_and_realm.vec & !is.na(tip_and_realm.vec),grepl(Realm, names(InTree_tbl))]<-100 # assign probability "100" 
    InTree_tbl[tip_and_realm.vec & !is.na(tip_and_realm.vec),!grepl(Realm, names(InTree_tbl)) & names(InTree_tbl) %in% ProbColNames]<-0 # assign probability "0" 
  }
  #11) convert tibble to phylodata object class S4 for ggtree plotting
  InTree_treedata<-as.treedata(InTree_tbl) 
  #12) plot tree
  p1<-ggtree(InTree_treedata,ladderize=TRUE, right=TRUE)+ #ladderize with smallest clade on top
    geom_tiplab(offset = 3) + #show all tip labels, adjust position to left
    xlim_tree(130) + #add space for tip labels
    ggtitle(paste ("RASP data mapped:",RASP_IDs_subsample[i]))
  #13) add internal node probabilities as pie charts using ggtree::inset (see https://guangchuangyu.github.io/software/ggtree/vignettes/ggtree-inset.html#annotate-with-other-types-of-charts and http://www.randigriffin.com/2017/05/11/primate-phylogeny-ggtree.html)
  InTree_tbl_rasp_only<-InTree_tbl[!is.na(InTree_tbl$raspNode),] #subset rows with rasp data
  #export to check in excel
  #flatten the list to matrix, convert matrix to tibble
  #tree_rasp_only_flat <- as_tibble(apply(tree_rasp_only,2,as.character))
  #write as csv
  #write.xlsx(tree_rasp_only_flat, "tree_rasp_only_df.xls")
  #extract probability column names
  SelColNames<-colnames(InTree_tbl_rasp_only[,grepl("_prob_", colnames(InTree_tbl_rasp_only))]) 
  #manual colors for realms (corresponding to realmColors1 but without realm names)
  realmColors<-c("khaki4", "lightsteelblue4", "maroon2" , "lightblue2", "firebrick1", "palevioletred2", "black", "lemonchiffon2", "blue1", "coral" ) 
  pies<-nodepie(InTree_tbl_rasp_only, cols=SelColNames,color=realmColors) #generate node pie charts for columns extracted above
  p2<-ggtree::inset(p1, pies, width=0.1,height=0.1,hjust=0.17,vjust=0.1)
  #14)add tip "probabilities" as another set of pie charts 
  tree_tips<-InTree_tbl[!is.na(InTree_tbl$label) & !is.na(InTree_tbl$Realm),] #tip rows with realm info
  pies2 <- nodepie(tree_tips, cols=SelColNames,color=realmColors) #generate node pie charts for columns extracted above
  p3<-ggtree::inset(p2, pies2, width=0.16, height=0.16, hjust = -1.5)
  p[[i]]<-p3 #save final plot to plotlist "p"
}

dev.off() #reset the graphic device
#arrange into grid of three (larger grid typically crashes RStudio)
grid.arrange(p[[1]],p[[2]],p[[3]],ncol=3) # dropped 23 of 45 Neotropical tips
grid.arrange(p[[4]],p[[5]],p[[6]],ncol=3) # dropped 30 of 45 Neotropical tips
#export 30x40mm
  #=>FigureS6
  #=>FigureS7

#########################################
#Historical biogeography with Biogeobears
#########################################
#see http://phylo.wikidot.com/biogeobears#script for tutorial
tr<-drop.tip(tree0_treedata, NonKalos_lines)
tr<-as.phylo(tr)
tr<-ape::ladderize(tr,right="TRUE") #ladderize to have smallest clades on bottom
write.tree (tr, paste(inDir,"BEAST8_KaloOnly.tree",sep=""))
trfn<-paste(inDir,"BEAST8_KaloOnly.tree",sep="")
  #=> 174 tips

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
#INSPECT geography file
  #moref(geogfn)
#define maximum number of areas per species
max_range_size <- 3

##add fossils to constrain nodes (see http://phylo.wikidot.com/fossil-data-in-biogeographical-analysis-in-biogeobears#toc0)
#CHECK: node numbers
  tree<-read.tree(trfn)
#CHECK: plot in ggtree
  tree<-read.tree(trfn)
  ggtree(tree,ladderize=TRUE, right=TRUE)+
  geom_text2(aes(label=node), hjust=-.3)+  #show all node number
  geom_tiplab(offset = 2)+
  xlim_tree(150) 
#make vectors for constraining realms via "fixlikes" argument
#CHECK: order of ranges
  rangeListNums<-res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
#generate list of ranges in the same order as used by biogeobears
Areas<-c("Afr","Ori","Aus","Neot","Oce","Mad","Sah","Nea","Sin","Pal")
RangeGen<-function(NumAreas){
  Areas<-unname(as_tibble(combn(Areas, NumAreas)))
  apply(Areas,2,paste,collapse="")
}
Ranges2<-RangeGen(2)
Ranges3<-RangeGen(3)
Ranges_all<-c("null",Areas,Ranges2,Ranges3) #concatenate ranges with 2 areas and 3 areas, add "null" range

#CHECK lenght
length(Ranges_all)
  #=> 176 possible ranges for max 3

#generate "fixlikes" vectors
MakeFixlikeVecs<-function(MustInclude) 
{
  fixlikeVec<-vector()
  for (i in 1:length(Ranges_all)) 
  {
    Area<-Ranges_all[i]
    if (grepl(MustInclude,Area)) 
    {
      Area<-"1"
    } else {
      Area<-"0"
    }
    fixlikeVec<-append(fixlikeVec,Area)
  }
  as.numeric(fixlikeVec)
  #paste(fixlikeVec,collapse=" ")
}
#make vector of fixed probailities for "fixlikes" argument of Biogeobears run
FixOri<-MakeFixlikeVecs("Ori")
#FixPal<-MakeFixlikeVecs("Pal") #uncomment to include Palearct as another constraint
#FixLikes <- as.matrix(rbind(FixOri,FixPal)) #uncomment to bind multiple FixLikes into matrix

##set-up BioGeoBEARS_run_object
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn
# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
#load the dispersal multiplier matrix etc. from the text files into the model object
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
#add node fossil info
BioGeoBEARS_run_object$fixlikes<-FixOri
BioGeoBEARS_run_object$fixnode<-175 #node 175: all kalotermitids
  
# INSPECT the BioGeoBEARS_run_object
  #BioGeoBEARS_run_object
#define CPU number
BioGeoBEARS_run_object$num_cores_to_use = 8
## run Biogeobears
res = bears_optim_run(BioGeoBEARS_run_object)
#save results into file
#resFile<-"Biogeobears_BEAST8_3maxrange.Rdata"
resFile<-"Biogeobears_BEAST8_3maxrange_1fossil.Rdata"
save(res, file=resFile)
#LOAD the saved result
  #load(resFile)
  #=> "res" results

###Plot ancestral states
analysis_titletxt<-paste("analysis file:",resFile)
# Setup
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
# plot states
plot_BioGeoBEARS_results(results_object, 
                                analysis_titletxt, 
                                addl_params=list("j"), 
                                plotwhat="text", 
                                label.offset=0.1, 
                                tipcex=0.2, 
                                statecex=0.5, 
                                #splitcex=0.1, 
                                titlecex=0.8, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=TRUE, 
                                #tr=tr, 
                                #tipranges=tipranges,
                                plotsplits = FALSE, 
                                plotlegend = FALSE)
# plot pie chart
plot_BioGeoBEARS_results(results_object, 
                         analysis_titletxt, 
                         addl_params=list("j"), 
                         plotwhat="pie", 
                         label.offset=0.45, 
                         tipcex=0.7, 
                         statecex=0.7, 
                         splitcex=0.6, 
                         titlecex=0.8, 
                         cornercoords_loc=scriptdir, 
                         include_null_range=TRUE, 
                         #tr=tr, 
                         #tipranges=tipranges,
                         plotsplits = FALSE, 
                         plotlegend = FALSE)

#export:20x30inch
  #=> Figure S5