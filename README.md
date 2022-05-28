# ###Under construction###
# Molecular phylogeny reveals the past transoceanic voyages of drywood termites (Isoptera, Kalotermitidae)
An online supplemental material containing input data and scripts used to generate most of the figures in our manuscript on phylogeny and biogeography of dry wood termites (Kalotermitidae).

If you use this resource, please cite it as follows: 
<b><br>Buček A., Wang M., Šobotník J., Sillam-Dussès D., Mizumoto N., Stiblík P., Clitheroe C., Lu T., González Plaza J. J., Mohagan A., Rafanomezantsoa J. J., Fisher B., Engel M. S., Roisin Y., Evans T. A., Scheffrahn R., Bourguignon T. 
<br><i>Transoceanic voyages of drywood termites (Isoptera: Kalotermitidae) inferred from extant and extinct species</i></b> 
<br>https://www.biorxiv.org/content/10.1101/2021.09.24.461667v2


## 1) R scripts
1) R scripts used to generate the article figures: <i>Kalotermitidae_scripts.R</i>
  
## 2) Bash helper scripts
All scripts are in Bash_scripts/ repository folder. IMPORTANT: See comments within scripts for file format and naming requirements. 
#TODO: add example data
1) for re-formating of RASP results prior to import into R: <i>reformat_RASP2.sh</i><br>
2) for re-formating of MitoZ feature tables for import into JalView as sequence annotation file: <i>MitoZ_table_to_JalView.sh</i>
* run with two arguments; use absolute paths to files and directories. Outputs 1: feature tables converted to JalView annotation tables, and 2: concatenated JalView annotation table withGFF header that can be imported as feature annotation file into JalView
```
MitoZ_table_to_JalView.sh <directory with feature tables names as "*_mitoscaf.fa.tbl"> <output name>
```
4) for semi-manual annotation of features: <i>Alignment_to_FeatureTable.sh</i> 
* the script asists with visual mitochondrial genome annotation (e.g. fixing of incorrect annotations produced by automatized annotation pipelines) 
* copy-pasted alignment block corresponding to a feature to be added to annotation is merged with feature tables generated by automatized annotation pipeline in MitoZ (or similar, but feature tables generated by other pipelines than MitoZ might need further formatting) 
* run script as follows with five arguments. Use absolute paths to files and directories.
```
Alignment_to_FeatureTable.sh <feature name> <file with feature alignment block> <file with feature template> <directory with feature tables to be updated> <directory with fasta files of annotated genomes>
```
  
## 3) Data
1) Assembled mitochondrial genomes and the feature annotations are available via GenBank under accession: OK506044 to OK506060 and OM991296 to OM991443
2) Evolutionary trees in Newick format: <i>Data S1/</i><br>
3) Input data for R scripts: <i>data/</i><br>
4) UCE loci: https://doi.org/10.5061/dryad.5mkkwh77v
