###CONTENTS
#package setup
#file import and pre-processing
#28sp time tree comparison
#Tanglegrams of UCE-based ML trees
#combine BEAST, IQTREE, RASP onto consensus tree
#stats for summary tree
#plot summary tree
#polytomise branches without consensus support
#Biogeographic realm map
#fossil ages plot
#subsample Neotropical Kalos, (analyze with RASP), plot
#historical biogeography with Biogeobears

################
#package setup
###############
#install.packages("pacman")
pacman::p_load(sf,
               tidyverse,
               ggplot2,
               gridExtra,
               plotly,
               phangorn,
               treeio,
               DECIPHER,
               dendextend,
               phylogram,
               dplyr,
               ggrepel,
               tidytree,
               ggimage,
               gtools,
               adephylo,
               ggtree,
               xlsx, 
               evobiR, #for beast tree subsampling
               phytools, #for getParent function
               BioGeoBEARS,
               phylobase, 
               ape) 

##############################
#file import and pre-processing
##############################
#NOTE1: no spaces allowed in tree tip labels
#NOTE2: place "data" directory into the working directory
inDir<-"data/"

#1) import sample data table (first sheet of excel table)
KaloTable<-read.xlsx2(paste(inDir,"TableS1_v16.xlsx",sep=""), 1) #use read.xlsx2 so empty cells are imported as no values rather than NAs
KaloTable[KaloTable==""]<-NA #replace empty cells with NAs
#add "TipLabels" column by pasting: "Collection_code_updated",_,"Genus",_,"Species"
KaloTable<-add_column(KaloTable, TipLabels = paste(KaloTable$Collection_code_updated,KaloTable$Genus,KaloTable$Species,sep="_"), .after = 2) 
#add "TipLabels2" column with species name by pasting: "Genus", ,"Species_revised"
KaloTable<-add_column(KaloTable, TipLabels2 = paste(KaloTable$Genus,KaloTable$Species_revised,sep=" "), .after = 2) 


#2)define geography file for BioGeoBears
geogfn <- paste(inDir,"KaloOnly_geography_forBiogeobears.txt",sep="")

#3) import trees with summarized support from Treegraph onto "BEAST_REV4_3a" topology
#vector with tree files in newick format
Treegraph.files<-paste(inDir,c("Treegraph_BEAST_REV9a_conflicts.tre", "Treegraph_BEAST_REV9a_supports.tre","Treegraph_IQTREE_REV2_conflicts.tre","Treegraph_IQTREE_REV2_supports.tre","Treegraph_IQTREE_REV3_conflicts.tre","Treegraph_IQTREE_REV3_supports.tre","Treegraph_BEAST_REV4_3a_supports.tre"),sep="")
#vector with tree ids (strip path and constant prefix and suffix)
Treegraph.ids<-sub(paste(inDir,"Treegraph_",sep=""),"",Treegraph.files) 
Treegraph.ids<-sub(".tre","",Treegraph.ids) 
# import trees as list, name list elements by tree IDs
Treegraph.lst<-lapply(Treegraph.files, function(x) treeio::read.newick(x))
names(Treegraph.lst)<-Treegraph.ids

#4) import backbone tree for all analyses and plotting (uncomment below to change backbone tree)
TreeTopo<-"BEAST_REV4_3a" #select backbone tree topology (NOTE: other than"BEAST_REV4_3a" topologies will not work for plots and analysis that use topology conflicts summarized in Treegraph)
#TreeTopo<-"BEAST_REV9a" #uncomment for RASP files for alternative topology #1
#TreeTopo<-"IQTREE_REV2" #uncomment for RASP files for alternative topology #2
#TreeTopo<-"IQTREE_REV3" #uncomment for RASP files for alternative topology #3
#define and run import function (trees with "BEAST" in name imported as beast trees, otherwise as newick)
ImportRefTree<-function(TreeTopo){
  if(grepl("BEAST",TreeTopo)) {
    treeio::read.beast (file = paste(inDir,TreeTopo,".tre",sep=""))
  }else {
    treeio::read.newick (file = paste(inDir,TreeTopo,"_Dendroscope_rooted.tree",sep=""))
  }
}
tree0_treedata<-ImportRefTree(TreeTopo)

#5) import trees for 28-species subtree comparison
#ToDrop28sp<-KaloTable %>% filter(UCE_subtree!=1 | Family=="OUTGROUP") %>% select(TipLabels) %>% pull() #non-UCE-subtree or non-Kalo tip labels
ToDrop28sp<-KaloTable %>% filter(UCE_subtree!=1 ) %>% select(TipLabels) %>% pull() #comment above and uncomment this to keep non-Kalos
#import function: import BEAST trees; drop tips; convert to tibble
SmallTreeImport<-function(treeFile) {
  treeio::read.beast(file = paste(inDir,treeFile,sep="")) %>%
    drop.tip(ToDrop28sp)
}
UCEtree28sp<-SmallTreeImport("BEAST_REV12.tre")
mtGtree230spProtModel<-SmallTreeImport("BEAST_REV4_3a.tre")
mtGtree230spNucModel<-SmallTreeImport("BEAST_REV9a.tre")
mtGtree28spProtModel<-SmallTreeImport("BEAST_REV10a.tre")

#6) define maximum-likelihood tree files for 82 vs. 691 loci comparison 
loci691Tree<-paste(inDir,"IQTREE_REV4.treefile_noALRT_FigTreed.tre",sep="")
loci82Tree<-paste(inDir,"IQTREE_REV5.treefile_noALRT_FigTreed.tre",sep="")
mtGtree230spProtModelFull<-paste(inDir,"IQTREE_REV2_Dendroscope_rooted.tree",sep="")
mtGtree230spNucModelFull<-paste(inDir,"IQTREE_REV3_Dendroscope_rooted.tree",sep="")

#7) import tables with tree node numbers for BEAST fossil calibrations (specific for each backbone tree topology)
fossils<-read_csv(paste(inDir,"Fossil_node_numbers_for_BEAST_REV4_3a.csv",sep=""))

#8) import reformated RASP tables ONLY for reference topology (NOTE: errors in RASP output - ocassionally missing zeroes - fixed manually)
RASP_indir<-paste(inDir,"RASP_REV_formated_manFixed/",sep="")
RASP_files<-list.files(RASP_indir,full.names=T)
RASP_files_refTopo<-RASP_files[ grepl(TreeTopo, RASP_files) ] #RASP files for only reference topology
RASP_IDs<-paste(gsub(".txt","",list.files(RASP_indir)),"_prob",sep="") #get analysis IDs from file names by dropping ".txt" suffix and adding "_prob"
RASP_IDs<-RASP_IDs[ grepl(TreeTopo, RASP_IDs) ] #RASP IDs for only reference topology
RASP.lst<-lapply(RASP_files_refTopo, function(x) read_tsv(x))
  #=> list of 4 RASP analyses results for reference tree topology
  #=> 171 internal node info

#9) import kalotermitid fossil records from PaleoDB
PaleoDB<-read.xlsx(paste(inDir,"TableS7.xlsx",sep=""),sheetIndex=1, na = c("", "NA"),stringsAsFactors = FALSE) #import from excel, empty cells as NAs, strings not as factor to enable mutate later
  #CHECK
  nrow(PaleoDB)
  #=> 74 records

#10) drop non-Kalotermitid tips from 280samples trees and export the trees in nwk format for analysis in RASP
NonKalos<-KaloTable %>% filter(Family=="OUTGROUP") %>% select(TipLabels) %>% pull() #non-Kalo tip labels
#import Treegraph tree files in newick format (use Treegraph trees - node values not important)
PhyloTrees.files<-paste(inDir,c("BEAST_REV9a.nwk","IQTREE_REV2.nwk","IQTREE_REV3.nwk","BEAST_REV4_3a.nwk"),sep="")  #vector with tree file paths
PhyloTrees.ids<-sub(inDir,"",PhyloTrees.files) # strip path prefix
PhyloTrees.ids<-sub(".nwk","",PhyloTrees.ids) #  strip suffix
# import trees as list, name list elements by tree IDs
PhyloTrees.lst<-lapply(PhyloTrees.files, function(x) treeio::read.newick(x))
names(PhyloTrees.lst)<-PhyloTrees.ids
for (i in 1:length(PhyloTrees.lst)) {
  write.tree(drop.tip(PhyloTrees.lst[[i]],NonKalos),paste(inDir,names(PhyloTrees.lst)[i],"_KaloOnly.tre",sep=""))
}
  
##########################
#28sp time tree comparison
##########################
#plotting function
SmallTreePlot<-function(TreeName, inData, BranchColor,BarTipColor,xAxisMax)
{
  ggtree(inData,ladderize=TRUE, right=TRUE,color=BranchColor)+ #ladderize with smallest clade on top
    geom_tiplab(offset = 2,color=BarTipColor,size=7) + #show all tip labels, adjust position to left
    geom_nodepoint(aes(shape="1",), size=3, color=BarTipColor) + 
    geom_range("height_0.95_HPD", color=BarTipColor, size=2, alpha=.6) + # add branch height (~age) credible intervals (95% HPD)
    xlim_tree(70) + #add space for tip labels
    theme_tree2(legend.position=c(0.2,0.3))+ #"theme_tree2" to add timescale, define legend parametres
    geom_vline(xintercept = seq(xAxisMax, -10, by = 10), linetype="dotted", alpha=.5)+ #add vertical lines
    scale_x_continuous(labels = abs, breaks=seq(xAxisMax, 0, by = 10))+ #define x-axis ticks every 10MYA
    ggtitle(paste (TreeName))+
    theme(legend.position = "none") #remove legend
}

#plot 4 tree with "mtGtree230spProtModel" as backbone tree with visible edges
p1<-SmallTreePlot("mtGtree230spProtModel",mtGtree230spProtModel,"black","black",-110) %>% revts() #backbone tree
p2<-SmallTreePlot("mtGtree28spProtModel",mtGtree28spProtModel,NA,"blue",-110) %>% revts()
p3<-SmallTreePlot("UCEtree28sp",UCEtree28sp,NA,"orange",-110) %>% revts()
p4<-SmallTreePlot("mtGtree230spNucModel",mtGtree230spNucModel,NA,"red",-110) %>% revts()
grid.arrange(p1,p2,p3,p4)
#export as pdf: 30x30inch

#plot all trees with visible edges, with extended x-axis (for plotting with non-Kalo outgroups after modifying the tree import)
p1<-SmallTreePlot("mtGtree230spProtModel",mtGtree230spProtModel,"black","black",-180) %>% revts() #backbone tree
p2<-SmallTreePlot("mtGtree28spProtModel",mtGtree28spProtModel,"blue","blue",-180) %>% revts()
p3<-SmallTreePlot("UCEtree28sp",UCEtree28sp,"orange","orange",-180) %>% revts()
p4<-SmallTreePlot("mtGtree230spNucModel",mtGtree230spNucModel,"red","red",-180) %>% revts()
grid.arrange(p1,p2,p3,p4,nrow=2)
#export as pdf: 30x40inch
  #=> Figure S3


############
#Tanglegrams
############
#define function for making tanglegrams
DoTanglegram <- function (leftTree, rightTree, leftname, rightname) {
  dendLeft <- read.dendrogram(leftTree) # import first newick as dendogram, ordering set in FigTree as increasing to have decreasing in the dendogram, re-rooted on Mastotermed
  dendRight <- read.dendrogram(rightTree)
  tanglegram(dendLeft, dendRight, 
             highlight_branches_lwd=FALSE, 
             common_subtrees_color_branches = FALSE, 
             common_subtrees_color_lines = FALSE, #horizontal lines' colors
             axes=FALSE, 
             highlight_distinct_edges=FALSE,
             main_left=leftname,
             main_right=rightname,
             #lab.cex = .5, 
             margin_inner = 15)
}

#plot tanglegram of IQtrees for 691 vs 82 loci
DoTanglegram(loci691Tree,loci82Tree, "691 loci","82 loci")

#plot tanglegram of IQtrees for nucleotide vs. nucleotide+protein model
DoTanglegram(mtGtree230spProtModelFull,mtGtree230spNucModelFull, "mtGtree230spProtModelFull","mtGtree230spNucModelFull")


###############################################
#combine BEAST, IQTREE, RASP onto consensus tree
###############################################
tree0_phylo<-as.phylo(tree0_treedata) #convert tree0 to phylo object required for getMRCA(function)
tree_all<-as_tibble(tree0_treedata) #convert tree object to tibble for data manipulation (https://yulab-smu.top/treedata-book/chapter2.html)

#add supports and conflicts from Treegraph 
for (i in 1:length(Treegraph.ids)) {
  treeID<-Treegraph.ids[i]
  tree_all[treeID]<-as_tibble(Treegraph.lst[[i]])$label
  tree_all[!is.na(tree_all$label),treeID]<-NA #remove tip labels imported along with node labels (=supports,conflicts) from Treegraph trees
}

#add a YES/NO variable for nodes with above-threshold support from all trees, skip empty values
tree_all$support_by_all<-ifelse(as.numeric(tree_all$BEAST_REV9a_supports)>=0.95 & tree_all$BEAST_REV9a_supports!="" & 
                                as.numeric(tree_all$BEAST_REV4_3a_supports)>=0.95 & tree_all$BEAST_REV4_3a_supports!="" &
                                as.numeric(tree_all$IQTREE_REV2_supports)>=95 & tree_all$IQTREE_REV2_supports!="" &
                                as.numeric(tree_all$IQTREE_REV3_supports)>=95 & tree_all$IQTREE_REV3_supports!="", 
                                "YES", "NO")

#add a YES/NO variable for nodes with conflict
tree_all$conflict<-ifelse(tree_all$BEAST_REV9a_conflicts!="" | tree_all$IQTREE_REV2_conflicts!="" | tree_all$IQTREE_REV3_conflicts!="", "YES", "NO")

### add fossil calibrations
#add fossil details to matching nodes of "tree_all" data
tree_all$Fossil<-fossils$Fossil[match(tree_all$node,fossils$Fossil_node)]
tree_all$Fossil_number<-fossils$Fossil_number[match(tree_all$node,fossils$Fossil_node)]

### add RASP probabilities calculated ONLY for reference topology
#1) pre-process RASP data
#1A) rename probability columns to include: analysis ID and realm ID (original colname)
for (i in 1:length(RASP.lst)) 
  {
  analysisID<-RASP_IDs[i]
  newColNames<-paste(analysisID,colnames(RASP.lst[[i]][4:13]), sep = "_") 
  colnames(RASP.lst[[i]])[4:13]<-newColNames
}
#1B) rename rasp node ("node") to "raspNode"
colnames(RASP.lst[[1]])[1]<-"raspNode" 
#1C) add corresponding node numbers from reference tree to RASP lists 
for (j in 1:length(RASP.lst)) 
  {
  for (i in 1:nrow(RASP.lst[[j]])) 
    {
    tip1<-as.character(RASP.lst[[j]][i,"tip1ID"]) #tip1 from RASP
    tip2<-as.character(RASP.lst[[j]][i,"tip2ID"]) #tip2 from RASP
    RASP.lst[[j]][i,"treeNode"]<-getMRCA(tree0_phylo,c(tip1,tip2)) #
    }
  }
#2) add rasp node numbers to tree_all
tree_all$raspNode<-RASP.lst[[1]]$raspNode[match(tree_all$node,RASP.lst[[1]]$treeNode)]
#3) add rasp probabilities to tree_all
for(j in 1:length(RASP.lst)) 
  {
  Realms<-colnames(RASP.lst[[j]][4:13])
  for(i in 1:length(Realms)) 
    {
    Realm<-Realms[i]
    tree_all[,Realm]<-RASP.lst[[j]][match(tree_all$node,RASP.lst[[j]]$treeNode),Realm]
    }
  }
#4) add realm information to Kalotermitidae tips from "KaloTable"
tree_all[,"Realm"]<-KaloTable[match(tree_all$label,KaloTable$TipLabels),"Biogeographic_realm"] #add column "Realm" with a value for each tip (i.e. row with tip label in "label" column)
#5) fill probability values for Kalotermitidae tips: 100 for the corresponding realm probability, 0 for all other realm probabilities
Realms<-KaloTable %>% filter(Family=="Kalotermitidae") %>% select(Biogeographic_realm) %>% unique() %>% pull() %>% as.character() #get vector of realms, convert to "character" type, otherwise automatically factor with levels causing mess downstream
ProbColNames<-names(tree_all)[grepl("_prob_",names(tree_all) )] #vector of all colnames with rasp probabilities (grepped as containing "fixed" or "estimated" string)
for(i in 1:length(Realms)) {
  Realm<-Realms[i]
  tip_and_realm.vec<-!is.na(tree_all$label) & tree_all$Realm==Realm 
  tree_all[tip_and_realm.vec & !is.na(tip_and_realm.vec),grepl(Realm, names(tree_all))]<-100 # assign probability "100" 
  tree_all[tip_and_realm.vec & !is.na(tip_and_realm.vec),!grepl(Realm, names(tree_all)) & names(tree_all) %in% ProbColNames]<-0 # assign probability "0" 
}
#6) add mean and SD RASP probabilities
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
#add column with max RASP %probability for each node
tree_all$maxRASP<-apply(tree_all[,grep("Means_",colnames(tree_all))], MARGIN=1, FUN=max)


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
#visually inspect and manually make a list of tips to be dropped (below 2MY distance AND same realm/country)
  #=> Below2MY-redundant-tips-to-drop.txt
tips2drop<-as.vector(read.csv(paste(inDir,"Below2MY-redundant-tips-to-drop.txt",sep=""),header = F)$V1) #list of tips with <2M years distance to be dropped
  #CHECK: how many tips?
  length(tips2drop)
  #=> 35
tree_all_treedata_drop<-drop.tip(tree_all_treedata, tips2drop)
  #=> from 230 to 195 tips
tree_all_tbl_drop<-as_tibble(tree_all_treedata_drop) #convert treedata to tbl for future data manipulations
tree_all_treedata_drop<-as.treedata(tree_all_tbl_drop)
tree_all_tbl<-as_tibble(tree_all_treedata) #convert treedata to tbl for future data manipulations
#generate alternative treedata with renamed tip labels excluding collection codes
tree_all_treedata_drop_simpleLab<-rename_taxa(tree_all_treedata_drop,data = KaloTable[,c("TipLabels","TipLabels2")],key = "TipLabels",value = "TipLabels2")


##############################
#stats for reduced summary tree
#############################
#tree only with Kalotermitids
#drop non-Kalotermitids
NonKalos_lines<-KaloTable %>% filter(Family=="OUTGROUP") %>% select(TipLabels) %>% pull() #non-Kalo tip labels
tree_all_td_drop_KaloOnly <- drop.tip(tree_all_treedata_drop,NonKalos_lines)
  #=> 137 tips
tree_all_phylo_drop_KaloOnly <- as.phylo(tree_all_td_drop_KaloOnly) #convert to phylo
tree_all_tbl_drop_KaloOnly<-as_tibble(tree_all_td_drop_KaloOnly) #convert to tibble

#define input trees
inTree_tbl<-tree_all_tbl_drop
inTree_treedata<-tree_all_treedata_drop
inTree_tbl_kaloonly<-tree_all_tbl_drop_KaloOnly
  #CHECK: uncomment to visually check tree
  #ggtree(inTree_treedata,ladderize=TRUE, right=TRUE)+
  #geom_tiplab() + #show all tip labels
  #geom_text2(aes(label=node), hjust=-.3, size=7,color="red") + #show internal node number
  #xlim_tree(250)  #add space for tip labels

#B)count highly-supported congruent, congruent and conflicting branches
nrow(subset(inTree_tbl_kaloonly, support_by_all=="YES"))
  #=>95 nodes with support from all trees
nrow(subset(inTree_tbl_kaloonly, subset=conflict=="NO" & support_by_all=="NO"))
  #=>21 nodes without conflict but also without support
  #=>95+21=116 nodes without conflict
nrow(subset(inTree_tbl_kaloonly, subset=conflict=="YES" & support_by_all=="NO"))
  #=>20 nodes with conflict

###import RASP data for all tree topologies; check mean+SD of RASP analysis for each tree
#define function
SDs_of_RASP_probs<-function(TreeID){ 
  #1)import reformatted RASP tables as named list
  RASP_analysis_types<-c("-JC","-JC+GAMMA","-F81","-F81+GAMMA")
  RASP_files<-paste(RASP_indir,TreeID,RASP_analysis_types,".txt", sep="") 
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
  NAcheck<-sum(is.na(RASP_joined)) #check that there are no missing values in the tables: 3rd item of output list should be zero
  output<-list(perRealm,amongRealms,NAcheck) #merge two dataframes into "output" list
  return(output)
}
#run function for each tree topology, output max SD per realm, check for "NAs" (~ errors in RASP output)
SDs_of_RASP_probs("BEAST_REV4_3a")
  #=>16%
SDs_of_RASP_probs("BEAST_REV9a")
  #=>25%
SDs_of_RASP_probs("IQTREE_REV2")
  #=>13%
SDs_of_RASP_probs("IQTREE_REV3")
  #=>12%
  ##=> TableS2-5

##################
#plot summary tree
##################
realmColors<-c("khaki4", "lightsteelblue4", "maroon2" , "lightblue2", "firebrick1", "palevioletred2", "black", "lemonchiffon2", "blue1", "coral" ) #manual colors for realms

#plotting function for summary tree
summaryTreePlot<-function(inTreeData,inTreeTbl,RASPanalysisSubset)
{ 
  tree_rasp_only<-inTreeTbl[!is.na(inTreeTbl$raspNode),] #rows with rasp data
  SelColNames<-colnames(tree_rasp_only[,grepl(RASPanalysisSubset, colnames(tree_rasp_only))]) #grepl for names of columns with desired subset of RASP probabilities 
  p1<-ggtree(inTreeData,ladderize=TRUE, right=TRUE)+ #ladderize with smallest clade on top
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
  ggtitle(paste("Backbone tree topology:",TreeTopo,"; RASP analysis subset:", RASPanalysisSubset))
  
  p2<-revts(p1) #reverse order of time scale, make absolute numbers
  pies<-nodepie(tree_rasp_only, cols=SelColNames,color=realmColors) #node pie charts for internal nodes
  p3<-ggtree::inset(p2, pies, width=0.06,height=0.06,hjust=0.17,vjust=0.1)
  tree_tips<-inTreeTbl[!is.na(inTreeTbl$label) & !is.na(inTreeTbl$Realm),] #non-outgroup tip rows 
  pies2 <- nodepie(tree_tips, cols=SelColNames,color=realmColors) #node pie charts for tips
  ggtree::inset(p3, pies2, width=0.16, height=0.16, hjust = -1.5)
}

#Figure 1 (labels without sample codes, dropped <2MY tips)
summaryTreePlot(tree_all_treedata_drop_simpleLab, tree_all_tbl_drop,"Means_")

#Figure S1 (labels with sample codes, all samples)
summaryTreePlot(tree_all_treedata,tree_all_tbl,"Means_")
  #export pdf: 40x30inch

#Plot separately probabilities for each RASP model
p1<-summaryTreePlot(tree_all_treedata_drop_simpleLab, tree_all_tbl_drop,"Means_")
p2<-summaryTreePlot(tree_all_treedata_drop_simpleLab, tree_all_tbl_drop,"F81_prob")
p3<-summaryTreePlot(tree_all_treedata_drop_simpleLab, tree_all_tbl_drop,"F81\\+GAMMA")
p4<-summaryTreePlot(tree_all_treedata_drop_simpleLab, tree_all_tbl_drop,"JC_prob")
p5<-summaryTreePlot(tree_all_treedata_drop_simpleLab, tree_all_tbl_drop,"JC\\+GAMMA")
grid.arrange(p1,p2,p3,p4,p5,ncol=3)

#simplified summary tree plotting (no supports, no fossils, no age intervals)
summaryTreePlotSimpler<-function(inTreeData,inTreeTbl,RASPanalysisSubset)
{ 
  tree_rasp_only<-inTreeTbl[!is.na(inTreeTbl$raspNode),] #rows with rasp data
  SelColNames<-colnames(tree_rasp_only[,grepl(RASPanalysisSubset, colnames(tree_rasp_only))]) #grepl for names of columns with desired subset of RASP probabilities 
  p1<-ggtree(inTreeData,ladderize=TRUE, right=TRUE)+ #ladderize with smallest clade on top
    geom_tiplab() + #show all tip labels
    ggtitle(paste("Backbone tree topology:",TreeTopo,"; RASP analysis subset:", RASPanalysisSubset))+
    geom_text2(aes(label=maxRASP,subset=!isTip),vjust=-0.5, hjust=-0.5)+ # add max RASP %prob as label to internal nodes
    coord_cartesian(clip = 'off')
  pies<-nodepie(tree_rasp_only, cols=SelColNames,color=realmColors) #node pie charts for internal nodes
  p3<-ggtree::inset(p1, pies, width=0.03,height=0.03)
  tree_tips<-inTreeTbl[!is.na(inTreeTbl$label) & !is.na(inTreeTbl$Realm),] #non-outgroup tip rows 
  pies2 <- nodepie(tree_tips, cols=SelColNames,color=realmColors) #node pie charts for tips
  ggtree::inset(p3, pies2, width=0.04, height=0.04)
}
#plot for previously selected reference topologies other than BEAST_REV4_3a
summaryTreePlotSimpler(tree_all_treedata_drop, tree_all_tbl_drop,"Means_")

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
#polytomise nodes with conflict
tree_all_phylo_drop_relabel_polytom<-as.polytomy(tree_all_phylo_drop_relabel, feature='node.label', fun=function(x) as.numeric(x)<100)
#CHECK: plot as cladogram
    ggtree(tree_all_phylo_drop_relabel_polytom,ladderize=TRUE, right=TRUE, branch.length = "none")+
      geom_tiplab()+
      xlim_tree(50)
#export newick
write.tree(tree_all_phylo_drop_relabel_polytom, "tree_all_phylo_drop_relabel_polytom.nwck")
#...colapse at genus level and relabel in figtree
  #=> tree for Figure3

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
  #=> map in Figure1

###############
#fossil ages plot
###############
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
  #=>172 tips
#define seeds for 3 random samplings
seeds3<-c(1234,5678,9101)

#drop random samples of Neotropical tips from kalotermitid trees (define number of taxa to keep, basename of output files)
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
outName<-paste(outDir,"BEAST_REV4_3a",sep="")
#run subsampling to drop 23, and 30 Neotropical tips, respectively
subSampling(23,outName)
  #=> 149 tips
subSampling(30,outName)
  #=>142 tips

##(OUTSIDE OF R)
#1) subsample distribution table of all specimens to match the replicated random sampling
#2) run RASP with JC model and the subsampled input trees and the subsampled distribution tables => 6 RASP analyses
#3) reformat RASP tables prior importing into R

#list of input tree filenames
InTrees<-paste(outName,c("_23dropped-sample1.tre",
                     "_23dropped-sample2.tre",
                     "_23dropped-sample3.tre",
                     "_30dropped-sample1.tre",
                     "_30dropped-sample2.tre",
                     "_30dropped-sample3.tre"
                   ),sep=""
               )

#A)import reformatted RASP tables for sumbsampled Neotropics analysis
RASP_indir_subsample<-paste(inDir,"RASP_REV_formated-subsampled",sep="")
RASP_files_subsample<-list.files(RASP_indir_subsample,full.names=T)
RASP_IDs_subsample<-paste(gsub(".txt","",list.files(RASP_indir_subsample)),"_prob",sep="") #get analysis IDs from file names by dropping ".txt" suffix and adding "_prob"
inRASP<-lapply(RASP_files_subsample, function(x) read_tsv(x))

#B) rename probability columns to include analysis ID and realm ID (original colname)
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
write.tree (tr, paste(inDir,"BEAST_REV4_3a_KaloOnly.tree",sep=""))
trfn<-paste(inDir,"BEAST_REV4_3a_KaloOnly.tree",sep="")
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
BioGeoBEARS_run_object$fixnode<-173 #node 173: root of all kalotermitids
  
# INSPECT the BioGeoBEARS_run_object
  #BioGeoBEARS_run_object
#define CPU number
BioGeoBEARS_run_object$num_cores_to_use = 8
## run Biogeobears
res = bears_optim_run(BioGeoBEARS_run_object)
#save results into file
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

#rerun and replot as above but without "#add node fossil info" 
#..and save
resFile<-"Biogeobears_BEAST8_3maxrange_NOfossil.Rdata"
save(res, file=resFile)

#export:20x30inch
  #=> Figure S5