# 
# # These functions grab the genes and the gene summaries from genecards.org
# # and some calculate the mean and median fold change values across 
# # samples of treatment/control or diseased/healthy etc.
# 
# # To extract the genes:
# # - find25genes(protein) will grab the 25 genes associated with the protein from web
# # - getProteinGenes(protein) will print the genes associated with the protein
# # - getSummaries(gene, protein) will grab the gene protein summaries from web
# # -getGeneSummaries(protein) will print the gene summary of protein gene
# # 


# This script can be used by calling
# source('geneCards.R') in other script inside same folder. The table are specific
# to the functions to get the fold change values and in this folder.

library(rvest)
library(lubridate)
library(dplyr)

Gene_Path <- './gene scrapes'

if (dir.exists(Gene_Path)){
  unlink(Gene_Path, recursive=TRUE)
  dir.create(Gene_Path)
} else {
  dir.create(Gene_Path)
}

find25genes <- function(protein){
  
  url <- 'https://www.genecards.org/Search/Keyword?queryString=protein'
  
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','%20',protein)
  
  url <- as.character(url)
  url <- gsub('protein',protein, url)
  
  webpage <- read_html(url,encoding = "UTF-8")
  
  protein_html <- html_nodes(webpage,'.symbol-col a')
  protein1 <- html_text(protein_html)
  
  Protein <- as.data.frame(protein1)
  colnames(Protein) <- 'proteinType'
  Protein$proteinType <- as.character(paste(Protein$proteinType))
  Protein$proteinType <- gsub('\n','',Protein$proteinType)
  
  
  date <- as.data.frame(rep(date(),length(Protein$proteinType)))
  colnames(date) <- 'todaysDate'
  
  protein2 <- gsub('%20','-',protein)
  
  proteinName <- as.data.frame(rep(protein2,length(Protein$proteinType)))
  colnames(proteinName) <- 'proteinSearched'
  
  tableProtein <- cbind(Protein,proteinName,date)
  
  setwd(Gene_Path)
  
  
  write.table(tableProtein, 
              paste(protein2,".csv",sep=''), append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  names <- colnames(tableProtein)
  write.csv(names,paste('tableProteinHeader_',protein2,'.csv',sep=''),row.names=FALSE)
  
  setwd('../')
}


#find25genes('estrogen')


getProteinGenes <- function(protein){
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','-',protein)
  table <- read.csv(paste(Gene_Path,'/',protein,'.csv',sep=''),sep=',',
                    header=F,na.strings=c('',' ','NA'), stringsAsFactors = F)
  header <- read.csv(paste(Gene_Path,'/tableProteinHeader_',protein,'.csv',sep=''),
                     sep=',', header=T, na.strings=c('',' ','NA'), stringsAsFactors = F)
  names <- header$x
  colnames(table) <- names
  fileName <- paste('Top25',protein,'s.csv',sep='')
  write.csv(table, fileName, row.names=FALSE)
  return(list(table,fileName))
}


#getProteinGenes('estrogen')



getSummaries <- function(gene,protein){
  url <- 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE&keywords=protein'
  
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ',',',protein)
  gene <- as.character(gene)
  gene <- tolower(gene)
  
  url <- as.character(url)
  url <- gsub('GENE',gene,url)
  url <- gsub('protein',protein, url)
  
  webpage <- read_html(url,encoding = "UTF-8")
  
  Entrez_html <- html_nodes(webpage, '.gc-section-header+ .gc-subsection p')
  Entrez <- html_text(Entrez_html) 
  
  GeneCards_html <- html_nodes(webpage, '.gc-subsection-header+ p')
  GeneCards <- html_text(GeneCards_html) 
  
  UniProt_html <- html_nodes(webpage, '#summaries li:nth-child(1) div')
  UniProtKB <- html_text(UniProt_html) 
  
  Entrez0 <- ifelse(length(Entrez)==0, 'no summary',as.character(paste(Entrez)))
  Entrez1 <- as.data.frame(Entrez0)
  colnames(Entrez1) <- 'EntrezSummary'
  
  GeneCards0 <- ifelse(length(GeneCards)==0,'no summary',
                       as.character(paste(GeneCards)))
  GeneCards1 <- as.data.frame(GeneCards0)
  colnames(GeneCards1) <- 'GeneCardsSummary'
  
  UniProtKB0 <- ifelse(length(UniProtKB)==0,'no summary',
                       as.character(paste(UniProtKB)))
  UniProtKB1 <- as.data.frame(UniProtKB0)
  colnames(UniProtKB1) <- 'UniProtKB_Summary'
  
  Entrez1$EntrezSummary <- as.character(paste(Entrez1$EntrezSummary))
  Entrez1$EntrezSummary <- gsub('\n','',Entrez1$EntrezSummary)
  
  GeneCards1$GeneCardsSummary <- as.character(paste(GeneCards1$GeneCardsSummary))
  GeneCards1$GeneCardsSummary <- gsub('\n','',GeneCards1$GeneCardsSummary)
  
  UniProtKB1$UniProtKB_Summary <- as.character(paste(UniProtKB1$UniProtKB_Summary))
  UniProtKB1$UniProtKB_Summary <- gsub('\n','',UniProtKB1$UniProtKB_Summary)
  
  date <- as.data.frame(rep(date(),length(Entrez1$EntrezSummary)))
  colnames(date) <- 'todaysDate'
  
  protein2 <- gsub(',','-',protein)
  
  proteinName <- as.data.frame(rep(protein2,length(Entrez1$EntrezSummary)))
  colnames(proteinName) <- 'proteinSearched'
  
  gene <- as.data.frame(rep(toupper(gene),length(Entrez1$EntrezSummary)))
  colnames(gene) <- 'gene'
  
  tableProtein <- cbind(proteinName,gene,Entrez1,GeneCards1,UniProtKB1,date)
  
  setwd(Gene_Path)
  
  
  write.table(tableProtein, 
              paste(protein2,"summary.csv",sep=''), append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  names <- colnames(tableProtein)
  write.csv(names,paste('geneHeader_summary_',protein2,'.csv',sep=''),row.names=FALSE)
  
  setwd('../')
  return(gene)
}


#getSummaries('TP53','estrogen')


getGeneSummaries <- function(protein){
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','-',protein)
  
  table <- read.csv(paste(Gene_Path,'/',protein,'summary.csv',sep=''),
                    sep=',',header=F,na.strings=c('',' ','NA'), stringsAsFactors = F)
  
  header <- read.csv(paste(Gene_Path,'/geneHeader_summary_',protein,'.csv',sep=''),
                     sep=',', header=T, na.strings=c('',' ','NA'), stringsAsFactors = F)
  names <- header$x
  colnames(table) <- names
  
  fileName <- paste('proteinGeneSummaries_',protein,'.csv',sep='')
  write.csv(table, fileName, row.names=FALSE)
  return(list(table,fileName))
}


#getGeneSummaries('estrogen')

