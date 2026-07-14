## Stage 3: FungalTraits annotation, core genera, and PERMDISP.
##
## FungalTraits is a third-party resource and is not redistributed here. Place
## polme2020_genera.csv under data/external/ before running this optional stage.
options(stringsAsFactors=FALSE); set.seed(20260609)
suppressPackageStartupMessages({library(tidyverse); library(vegan); library(RColorBrewer); library(cowplot)})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
input_dir <- file.path(root, "data", "processed")
external_dir <- file.path(root, "data", "external")
out <- file.path(root, "outputs", "reanalysis")
tb <- file.path(out, "tables")
fg <- file.path(out, "figures")
dir.create(tb, recursive = TRUE, showWarnings = FALSE)
dir.create(fg, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  asv_counts = file.path(input_dir, "asv_table_fungi_rarefied_3819.tsv"),
  taxonomy = file.path(input_dir, "asv_taxonomy_fungi.tsv"),
  metadata = file.path(input_dir, "sample_metadata.tsv"),
  genus_biomarkers = file.path(input_dir, "genus_biomarkers_final.tsv"),
  fungaltraits = file.path(external_dir, "polme2020_genera.csv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop(
    "Missing input files:\n- ", paste(missing_inputs, collapse = "\n- "),
    "\nFungalTraits is not redistributed; see data/README.md.",
    call. = FALSE
  )
}

rr<-read.delim(inputs[["asv_counts"]],check.names=FALSE); rare<-as.matrix(rr[,-1]); rownames(rare)<-rr$sample
tax<-read.delim(inputs[["taxonomy"]],check.names=FALSE); rownames(tax)<-tax$asv_id
meta<-read.delim(inputs[["metadata"]],check.names=FALSE); meta<-meta[match(rownames(rare),meta$sample),]
meta$site<-factor(meta$site,levels=c("HD","XJ","DS","DG"))
FT<-read.csv(inputs[["fungaltraits"]],check.names=FALSE)
ft<-setNames(FT$primary_lifestyle, FT$GENUS); ft2<-setNames(FT$Secondary_lifestyle, FT$GENUS)
go<-c("HD1","HD2","HD3","XJ1","XJ2","XJ3","DS1","DS2","DS3","DG1","DG2","DG3")

g <- tax[colnames(rare),"Genus"]; life <- unname(ft[g]); life[is.na(g)|g==""]<-NA
asvlife <- ifelse(is.na(life)|life=="","Unassigned/no-trait",life)

## (1) community guild composition (read-weighted) per group1 + per site AMF/ECM
df <- as.data.frame(rare) %>% rownames_to_column("sample") %>% pivot_longer(-sample,names_to="asv",values_to="ab") %>% filter(ab>0)
df$life <- asvlife[match(df$asv, colnames(rare))]
df <- df %>% left_join(meta[,c("sample","group1","site")],by="sample")
gl <- df %>% group_by(group1,life) %>% summarise(ab=sum(ab),.groups="drop") %>% group_by(group1) %>% mutate(rel=ab/sum(ab)) %>% ungroup()
write.table(gl, file.path(tb,"fungaltraits_guild_by_group.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
top<-gl%>%group_by(life)%>%summarise(s=sum(ab))%>%arrange(desc(s))%>%pull(life); top<-setdiff(top,"Unassigned/no-trait")
keep<-head(top,9); gl$lab<-ifelse(gl$life%in%keep,gl$life,ifelse(gl$life=="Unassigned/no-trait","Unassigned/no-trait","Other guild"))
gl$lab<-factor(gl$lab,levels=c(keep,"Other guild","Unassigned/no-trait")); gl$group1<-factor(gl$group1,levels=go)
pal<-c(setNames(brewer.pal(max(3,length(keep)),"Set3")[seq_along(keep)],keep),"Other guild"="grey60","Unassigned/no-trait"="grey85")
pg<-ggplot(gl,aes(group1,rel,fill=lab))+geom_bar(stat="identity",position="fill",color="black",size=.12,width=.75)+
  scale_fill_manual(values=pal,name="Primary lifestyle\n(FungalTraits)")+scale_y_continuous(label=function(y)100*y)+
  labs(x="",y="Relative abundance (%)")+theme_bw(base_size=8)+theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(fg,"FigS_fungaltraits_guilds.pdf"),pg,width=20,height=10,units="cm")
cat("% fungal reads with a FungalTraits guild:",round(100*sum(df$ab[df$life!="Unassigned/no-trait"])/sum(df$ab),1),"\n")
# AMF & ECM read fraction per site
amf<-df%>%group_by(site)%>%summarise(AMF=round(100*sum(ab[life=="arbuscular_mycorrhizal"])/sum(ab),2),
  ECM=round(100*sum(ab[life=="ectomycorrhizal"])/sum(ab),2),
  saprotroph=round(100*sum(ab[grepl("saprotroph",life)])/sum(ab),1),
  plant_pathogen=round(100*sum(ab[life=="plant_pathogen"])/sum(ab),1))
cat("AMF / ECM / saprotroph / plant_pathogen read% by site:\n"); print(as.data.frame(amf))

## (2) biomarker guilds
bm<-read.delim(inputs[["genus_biomarkers"]],check.names=FALSE)
sig<-bm[bm$KW_BH<0.05 | (bm$IndVal_p<0.05 & bm$IndVal>0.6),]
sig$primary_lifestyle<-ft[sig$genus]; sig$secondary_lifestyle<-ft2[sig$genus]
write.table(sig[,c("genus","enriched","mean_RA","KW_BH","IndVal","primary_lifestyle","secondary_lifestyle")],
  file.path(tb,"biomarkers_with_fungaltraits.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("\nBiomarker guilds:\n"); print(sig[,c("genus","enriched","primary_lifestyle")],row.names=FALSE)

## (3) core genera — STRICT prevalence criterion: present in >=50% of samples at EVERY site.
## (Distinct from the Venn overlap, which counts ASVs merely detected >=1x in all four sites.)
pres<-rare>0
prev <- sapply(levels(meta$site), function(s){ idx<-meta$site==s; colSums(pres[idx,,drop=FALSE])/sum(idx) })
core <- rownames(prev)[apply(prev>=0.5, 1, all)]
venn_detected <- length(Reduce(intersect, lapply(levels(meta$site),
  function(s) colnames(pres)[colSums(pres[meta$site==s,,drop=FALSE])>0])))
cg<-tax[core,"Genus"]; cg<-cg[!is.na(cg)&cg!=""]
n_unassigned <- length(core) - length(cg)
coretab<-as.data.frame(sort(table(cg),decreasing=TRUE)); names(coretab)<-c("genus","n_core_ASVs"); coretab$primary_lifestyle<-ft[as.character(coretab$genus)]
write.table(coretab,file.path(tb,"core_genera.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat(sprintf("\nStrict core (>=50%% prevalence at every site): %d ASVs; %d named genera (+%d unassigned).\n  (For comparison, %d ASVs are merely detected in all four sites = Venn overlap.)\nCore genera:\n",
  length(core),nrow(coretab),n_unassigned,venn_detected))
print(coretab,row.names=FALSE)

## (4) PERMDISP
b<-vegdist(rare,"bray"); pd<-lapply(c("site","depth_order"),function(gv){bd<-betadisper(b,factor(meta[[gv]]));a<-anova(bd);data.frame(grouping=gv,F=round(a$`F value`[1],3),p=round(a$`Pr(>F)`[1],3))})
pd<-do.call(rbind,pd); write.table(pd,file.path(tb,"permdisp_final.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("\nPERMDISP:\n"); print(pd,row.names=FALSE); cat("\nDONE\n")
