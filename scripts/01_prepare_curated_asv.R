## Stage 1: prepare the BLAST-curated ASV table and genus indicator results.
##
## This provenance script starts from DADA2/BLAST intermediate files that are
## not distributed in this GitHub repository. See data/README.md for the exact
## expected filenames and the public raw-read accession.
options(stringsAsFactors=FALSE); set.seed(20260609)
suppressPackageStartupMessages({library(tidyverse); library(vegan); library(cowplot)
  library(RColorBrewer); library(grid); library(labdsv); library(VennDiagram)})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
intermediate_dir <- file.path(root, "data", "intermediate")
processed_dir <- file.path(root, "data", "processed")
out <- file.path(root, "outputs", "reanalysis")
tb <- file.path(out, "tables")
fg <- file.path(out, "figures")
dir.create(tb, recursive = TRUE, showWarnings = FALSE)
dir.create(fg, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  fungal_set = file.path(intermediate_dir, "fungal_set_BLASTrebuilt.txt"),
  asv_counts = file.path(intermediate_dir, "asv_table_nonchim.tsv"),
  taxonomy = file.path(processed_dir, "asv_taxonomy_fungi.tsv"),
  metadata = file.path(processed_dir, "sample_metadata.tsv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop(
    "Stage 1 requires non-distributed DADA2/BLAST intermediates. Missing:\n- ",
    paste(missing_inputs, collapse = "\n- "),
    "\nStart from scripts/02_run_community_statistics.R to use the included processed data.",
    call. = FALSE
  )
}
go <- c("HD1","HD2","HD3","XJ1","XJ2","XJ3","DS1","DS2","DS3","DG1","DG2","DG3")

fung <- readLines(inputs[["fungal_set"]])
tabw <- read.delim(inputs[["asv_counts"]],check.names=FALSE)
M <- as.matrix(tabw[,-1]); rownames(M)<-tabw$sample; M <- M[,colnames(M)%in%fung,drop=FALSE]
tax <- read.delim(inputs[["taxonomy"]],check.names=FALSE); rownames(tax)<-tax$asv_id
meta <- read.delim(inputs[["metadata"]],check.names=FALSE)
meta <- meta[match(rownames(M),meta$sample),]; meta$site<-factor(meta$site,levels=c("HD","XJ","DS","DG"))
meta$depth_factor<-factor(meta$depth_order,levels=c(1,2,3))

D <- min(rowSums(M)); rare <- rrarefy(M, D); rare <- rare[,colSums(rare)>0,drop=FALSE]
cat("Rarefied to", D, "; ASVs after:", ncol(rare), "\n")
write.table(cbind(sample=rownames(rare),as.data.frame(rare)), file.path(tb,sprintf("asv_table_fungi_rarefied_%d.tsv",D)), sep="\t",quote=FALSE,row.names=FALSE)

## ---------- composition Figure 2 (honest, Unassigned on top) ----------
long <- as.data.frame(rare) %>% rownames_to_column("sample") %>%
  pivot_longer(-sample,names_to="asv_id",values_to="abun") %>% filter(abun>0) %>%
  left_join(tax[,c("asv_id","Phylum","Class","Genus")],by="asv_id") %>%
  left_join(meta[,c("sample","group1")],by="sample")
pcol <- function(n){b<-brewer.pal(12,"Paired"); if(n<=12)b[seq_len(n)] else colorRampPalette(b)(n)}
build <- function(rank,reqclass=FALSE,topn=NULL){
  d<-long; d$t<-d[[rank]]; d$t[is.na(d$t)|d$t==""]<-"Unassigned"
  if(reqclass) d$t[d$t!="Unassigned"&(is.na(d$Class)|d$Class=="")]<-"Other"
  d<-d%>%group_by(group1,t)%>%summarise(ab=sum(abun),.groups="drop")
  tot<-d%>%filter(!t%in%c("Unassigned","Other"))%>%group_by(t)%>%summarise(s=sum(ab),.groups="drop")%>%arrange(desc(s))
  nd<-tot$t
  if(!is.null(topn)&&length(nd)>topn){k<-nd[seq_len(topn)]; d<-d%>%mutate(t=ifelse(t%in%c(k,"Unassigned"),t,"Other"))%>%group_by(group1,t)%>%summarise(ab=sum(ab),.groups="drop"); nd<-k}
  ho<-"Other"%in%d$t; lv<-c("Unassigned",if(ho)"Other",rev(nd))
  d$t<-factor(d$t,levels=lv); d$group1<-factor(d$group1,levels=go); list(d=d,nd=nd,ho=ho)
}
panel<-function(b,ttl){cols<-setNames(pcol(length(b$nd)),b$nd); if(b$ho)cols["Other"]<-"grey55"; cols["Unassigned"]<-"grey80"
  ggplot(b$d,aes(group1,ab,fill=t))+geom_bar(stat="identity",position="fill",color="black",width=.7,size=.15)+
   scale_fill_manual(values=cols,breaks=c("Unassigned",if(b$ho)"Other",b$nd),name=ttl)+
   scale_y_continuous(label=function(y)paste0(100*y))+background_grid(major="y",minor="none")+
   labs(x="",y="Relative Abundance (%)")+theme(axis.text.x=element_text(size=8,angle=45,hjust=1),
   axis.text.y=element_text(size=8),legend.key.size=unit(.36,"cm"),legend.text=element_text(size=7),legend.title=element_text(size=8))}
bP<-build("Phylum"); bC<-build("Class",topn=10); bG<-build("Genus",reqclass=TRUE,topn=10)
cat("Phylum:",paste(bP$nd,collapse=", "),"\n"); cat("Genus:",paste(bG$nd,collapse=", "),"\n")

## Venn (4-set, by species presence)
pres <- rare>0
sets <- lapply(levels(meta$site),function(s) colnames(pres)[colSums(pres[meta$sample[meta$site==s],,drop=FALSE])>0])
names(sets)<-levels(meta$site)
futile.logger::flog.threshold(futile.logger::ERROR,name="VennDiagramLogger")
venn_grob <- venn.diagram(sets[c("HD","DS","DG","XJ")],filename=NULL,
  fill=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"),alpha=.55,cex=1,cat.cex=1.2,margin=.08)
pD <- ggdraw() + draw_grob(grid::grobTree(children = do.call(grid::gList, venn_grob)))
fig2<-plot_grid(panel(bP,"Phylum"),panel(bC,"Class"),panel(bG,"Genus"),pD,labels=c("A","B","C","D"),ncol=2,label_size=14)
ggsave(file.path(fg,"Figure2_final.pdf"),fig2,width=32,height=20,units="cm")

## ---------- genus biomarkers (KW+BH & IndVal) ----------
g<-tax[colnames(rare),"Genus"]; keep<-!is.na(g)&g!=""; mg<-rare[,keep]; g<-g[keep]
gc<-t(rowsum(t(mg),group=g)); rel<-gc/rowSums(rare); rel<-rel[,colMeans(rel)>=0.001,drop=FALSE]
site<-meta$site
kw<-lapply(colnames(rel),function(gn){v<-rel[,gn];kt<-kruskal.test(v~site);gm<-tapply(v,site,mean)
  data.frame(genus=gn,enriched=names(which.max(gm)),mean_RA=round(100*mean(v),3),max_RA=round(100*max(gm),3),KW_p=kt$p.value)})%>%bind_rows()
kw$KW_BH<-p.adjust(kw$KW_p,"BH")
iv<-indval(rel,site); kw$IndVal<-round(iv$indcls[match(kw$genus,colnames(rel))],3)
kw$IndVal_site<-c("HD","XJ","DS","DG")[iv$maxcls[match(kw$genus,colnames(rel))]]; kw$IndVal_p<-iv$pval[match(kw$genus,colnames(rel))]
kw<-kw%>%arrange(KW_BH)
write.table(kw,file.path(tb,"genus_biomarkers_final.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
sig<-kw%>%filter(KW_BH<0.05 | (IndVal_p<0.05 & IndVal>0.6))%>%arrange(enriched,KW_BH)
cat(sprintf("\nGenus biomarkers (n tested=%d): %d significant\n",nrow(kw),nrow(sig)))
print(sig[,c("genus","enriched","mean_RA","max_RA","KW_BH","IndVal_site","IndVal","IndVal_p")],row.names=FALSE)
cat("Per site:\n"); print(table(sig$enriched))
cat("\nDONE\n")
