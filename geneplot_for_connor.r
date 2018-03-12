#Install packages

packs<-c("genoPlotR","data.table","igraph","RColorBrewer", "Biostrings","DECIPHER","phangorn","ggtree","dplyr","tidyr")
#lapply(packs,install.packages)
lapply(packs, require, character.only=T)

#Set working directory

setwd("C:/Users/clewi/Google Drive/Serina")
#setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/beta-lactone-NP-discovery/r-scripts/BGC_clusts/")

get_dna_segs<-function(my_list){
  for(i in 1:length(my_list)) {
    
    tmp=my_list[[i]]
    name=unlist(lapply(1:length(tmp),function(x) { strsplit(names(tmp)," \\[")[[x]][1] }))
    start<-vector(length = length(tmp))
    end<-vector(length = length(tmp))
    for(j in 1:length(tmp)){
      if(j==1){
        start[j]<-10
        end[j]<-width(tmp)[j]
      }
      else{
        start[j]<-end[j-1] + 10
        end[j]<-(width(tmp)[j] + start[j] + 1)
      }
    }
    strand<-rep(1,times = length(tmp))
    df1<-data.frame(name=name,
                    start=start,
                    end=end,
                    strand=strand)
    dna_seg1<-dna_seg(df1)
    num<-length(dna_seg1$col)
    #dna_seg1$col<-colorRampPalette(brewer.pal(9,"Spectral"))(num)
    dna_seg1$col<-"gray"
    dna_seg1$gene_type<-"headless_arrows"
  }
  return(dna_seg1)
}


#Read in the lipstatin gene cluster
cys <- readAAStringSet("BGC0000382.fasta")
head(cys)
cys6 <- cys[1:6]

ll<-list(cys6)

dna<-get_dna_segs(ll)
dnal<-list(dna)
dnal[[1]]$col <- brewer.pal(length(dnal[[1]]$col), "Set2")
dnal[[1]]$fill  <- brewer.pal(length(dnal[[1]]$col), "Set2")


dnal[[1]]$col[5] <- "#E78AC3"
dnal[[1]]$col[4] <- "#A6D854"

dnal[[1]]$fill[5] <- "#E78AC3"
dnal[[1]]$fill[4] <- "#A6D854"

mid_pos <- middle(dnal[[1]])
mid_pos
lett <- paste0("Lst",LETTERS[1:6])

annot<-genoPlotR::annotation(x1=mid_pos,
                             text=lett,
                             rot=0, col="black")

pl <- plot_gene_map(dna_segs = dnal,gene_type="arrows",
                    annotation_height= 1,
                    annotations=annot, dna_seg_label_cex = 5)
str(pl)
