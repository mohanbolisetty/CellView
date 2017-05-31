# CellView

### ShinyApp to visualize and explore single cell datasets ###

> > Hosted at: https://mbolisetty.shinyapps.io/CellView/

  An .Rds file made up of three data frames is necessary for upload, with the following object names and column names.

  log2cpm – This is your genes vs cells expression matrix. The gene names need to be in ENSG (or ENSM) ids.

  tsne.data – This is your clustering result. This contains 4 columns, {V1,V2,V3} for the 3 dimensions and dbCluster containing numerical cluster assignments. The row names of this data frame correspond to the column names of your expression matrix

  featuredata- One of the attached files to convert ENSG ids to human useful gene names. Use human or mouse depending on the data. Or, make your own. The Chromosome.Name column should be integer and at the minimum should have this column along with Associated.Gene.Name column. Beware of duplicate gene names with unique ENSG ids.

#### Code to generate an .Rds file for upload to CellView: ####

```
log2cpm<-read.csv('Data/Expression.csv',row.names=1,stringsAsFactors = F, as.is=T, check.names=F)
featuredata<-read.csv('Databases/HG19_v74_FeatureData.csv',row.names=1,stringsAsFactors = F, as.is=T,sep=',’)
tsne.data<-read.csv('Data/TNSE_dbscan.csv',row.names=1,stringsAsFactors = F,as.is=T)

save(log2cpm,featuredata,tsne.data,file=‘Filename.Rds’)
```
