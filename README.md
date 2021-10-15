# Kalotermitidae_biogeography
An online supplemental material containing input data and scripts used to generate most of the figures in our manuscript on phylogeny and biogeography of dry wood termites (Kalotermitidae).

If you use this resource, please cite it as follows: 
<b><br>Buček A., Wang M., Šobotník J., Sillam-Dussès D., Mizumoto N., Stiblík P., Clitheroe C., Lu T., González Plaza J. J., Mohagan A., Rafanomezantsoa J. J., Fisher B., Engel M. S., Roisin Y., Evans T. A., Scheffrahn R., Bourguignon T. 
<br><i>Transoceanic voyages of drywood termites (Isoptera: Kalotermitidae) inferred from extant and extinct species</i></b> 
<br>https://www.biorxiv.org/content/10.1101/2021.09.24.461667v2

# Repository content  
## R scripts
1) R scripts used to generate the article figures: <i>Kalotermitidae_scripts.R</i>
  
## Bash helper scripts
See comments within scripts for file format and naming requirements. 
#TODO: add example data
1) for re-formating of RASP results prior to import into R: <i>Bash_scripts/reformat_RASP2.sh</i><br>
2) for re-formating of MitoZ feature tables for import into JalView as sequence annotation file: <i>Bash_scripts/MitoZ_table_to_JalView.sh</i>
3) for semi-manual annotation of features: <i>Alignment_to_FeatureTable.sh</i> 
* the script asists with visual mitochondrial genome annotation (e.g. fixing of incorrect annotations produced by automatized annotation pipelines) 
* copy-pasted alignment block corresponding to a feature to be added to annotation is merged with feature tables generated by automatized annotation pipeline in MitoZ (or similar, but feature tables generated by other pipelines than MitoZ might need further formatting) 
* run script as follows with five arguments. Use absolute paths to files and directories.
```
Alignment_to_FeatureTable.sh <feature name> <file with feature alignment block> <file with feature template> <directory with feature tables to be updated> <directory with fasta files of annotated genomes>
```
  
## Data
1) Assembled mitochondrial genomes and the feature annotations are available via GenBank under accession: \<GenBank accessions to be updated\><br>
2) Evolutionary trees in Newick format: <i>Data S1/</i><br>
3) Input data for R scripts: <i>data/</i><br>
