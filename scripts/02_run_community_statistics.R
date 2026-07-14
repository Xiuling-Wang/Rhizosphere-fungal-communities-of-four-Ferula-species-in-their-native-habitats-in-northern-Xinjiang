## Stage 2: reproduce the community statistics from the frozen processed data.
## Outputs are written under outputs/reanalysis; published tables in
## data/processed are never overwritten.
options(stringsAsFactors=FALSE); set.seed(20260609)
suppressPackageStartupMessages({library(vegan); library(ape); library(ggplot2); library(ggrepel)
  library(pheatmap); library(rdacca.hp); library(RColorBrewer); library(reshape2); library(cowplot)})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
input_dir <- file.path(root, "data", "processed")
out <- file.path(root, "outputs", "reanalysis")
tb <- file.path(out, "tables")
fg <- file.path(out, "figures")
dir.create(tb, recursive = TRUE, showWarnings = FALSE)
dir.create(fg, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  asv_counts = file.path(input_dir, "asv_table_fungi_rarefied_3819.tsv"),
  taxonomy = file.path(input_dir, "asv_taxonomy_fungi.tsv"),
  metadata = file.path(input_dir, "sample_metadata.tsv"),
  genus_biomarkers = file.path(input_dir, "genus_biomarkers_final.tsv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing processed input files:\n- ", paste(missing_inputs, collapse = "\n- "), call. = FALSE)
}

rr <- read.delim(inputs[["asv_counts"]], check.names=FALSE)
rare <- as.matrix(rr[,-1]); rownames(rare) <- rr$sample
tax <- read.delim(inputs[["taxonomy"]], check.names=FALSE); rownames(tax)<-tax$asv_id
meta <- read.delim(inputs[["metadata"]], check.names=FALSE); meta <- meta[match(rownames(rare),meta$sample),]
meta$species_code <- factor(meta$site, levels=c("HD","XJ","DS","DG"))
meta$site <- meta$species_code; meta$depth_order <- factor(meta$depth_order, levels=c(1,2,3))
env_vars <- c("SM","OM","TN","TP","TK","Nitrate_N","Ammonium_N","Olsen_P","AK","pH","EC","TDS")
envz <- as.data.frame(scale(meta[,env_vars]))

############ FIG 3 — alpha diversity + soil correlation ############
{
  est <- t(estimateR(rare))
  a <- data.frame(sample=rownames(rare), Observed=specnumber(rare), Chao1=est[,"S.chao1"], ACE=est[,"S.ACE"],
                  Shannon=diversity(rare,"shannon"), Simpson=diversity(rare,"simpson"))
  a$Pielou <- a$Shannon/log(a$Observed)
  A <- merge(meta[,c("sample","species_code","depth_order")], a, by="sample")
  write.table(A, file.path(tb,"alpha_diversity_final.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  mets <- c("Observed","Chao1","ACE","Shannon","Simpson","Pielou")
  # KW tests overall (species, depth) + within-species depth, BH
  rows <- list()
  for(m in mets){
    for(gv in c("species_code","depth_order")){ kt<-kruskal.test(A[[m]]~A[[gv]]); rows[[length(rows)+1]]<-data.frame(metric=m,grouping=gv,p=kt$p.value)}
    for(sp in levels(A$species_code)){ d<-A[A$species_code==sp,]; if(length(unique(d$depth_order))>1){kt<-kruskal.test(d[[m]]~d$depth_order); rows[[length(rows)+1]]<-data.frame(metric=m,grouping=paste0("depth_within_",sp),p=kt$p.value)}}
  }
  kw<-do.call(rbind,rows); kw$p_BH<-p.adjust(kw$p,"BH")
  write.table(kw, file.path(tb,"alpha_kruskal_final.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  cat("Fig3 alpha: any BH<0.05?", any(kw$p_BH<0.05), "\n")
  ## Friedman repeated-measures test for depth (plant as block; respects nested design)
  meta$plant <- paste0(as.character(meta$site), meta$replicate)
  Af <- merge(A, meta[,c("sample","plant")], by="sample")
  ctab <- table(Af$plant, Af$depth_order); compl <- rownames(ctab)[rowSums(ctab>0)==3]
  fr <- sapply(mets, function(mm){
    M <- tapply(Af[[mm]], list(Af$plant, Af$depth_order), function(z) z[1])[compl, , drop=FALSE]
    friedman.test(M)$p.value })
  frd <- data.frame(metric=names(fr), Friedman_depth_p=round(unname(fr),4))
  write.table(frd, file.path(tb,"alpha_friedman_depth_final.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  cat("Friedman (depth, plant as block) — any p<0.05?", any(fr<0.05),
      "| range", sprintf("%.3f-%.3f", min(fr), max(fr)), "\n")
  # plot Shannon + Chao1 by depth x species
  long <- melt(A[,c("sample","species_code","depth_order","Shannon","Chao1")], id=c("sample","species_code","depth_order"))
  p3 <- ggplot(long, aes(depth_order,value,color=species_code))+geom_boxplot(outlier.shape=NA,width=.68)+
    geom_point(position=position_jitterdodge(.12,dodge.width=.68),size=1.3,alpha=.8)+
    facet_wrap(~variable,scales="free_y")+scale_color_brewer(palette="Dark2")+
    labs(x="Soil depth class",y=NULL,color="Species/site")+theme_bw(base_size=9)
  ggsave(file.path(fg,"Fig3AB_alpha.pdf"),p3,width=18,height=8,units="cm")
  # Fig3C correlation matrix: alpha indices + soil + depth
  cm_in <- cbind(A[,mets], meta[,env_vars], depth=as.numeric(as.character(meta$depth_order)))
  cc <- cor(cm_in, method="spearman", use="pairwise")
  pv <- matrix(NA,ncol(cm_in),ncol(cm_in),dimnames=dimnames(cc))
  for(i in 1:ncol(cm_in)) for(j in 1:ncol(cm_in)) pv[i,j]<-suppressWarnings(cor.test(cm_in[,i],cm_in[,j],method="spearman",exact=FALSE)$p.value)
  star <- ifelse(pv<0.001,"***",ifelse(pv<0.01,"**",ifelse(pv<0.05,"*","")))
  pdf(file.path(fg,"Fig3C_alpha_soil_corr.pdf"),width=9,height=8)
  pheatmap(cc, clustering_method="ward.D2", color=colorRampPalette(c("#2166ac","white","#b2182b"))(60),
           breaks=seq(-1,1,length.out=61), display_numbers=star, fontsize=8, main="Alpha indices, soil & depth (Spearman) - final")
  dev.off(); cat("Fig3 done\n")
}

############ helper: aggregate rank ############
agg <- function(rank){ v<-tax[colnames(rare),rank]; k<-!is.na(v)&v!=""; cnt<-t(rowsum(t(rare[,k]),group=v[k])); cnt/rowSums(rare) }

############ FIG 4A — class z-score heatmap ############
{
  cl <- agg("Class"); cls_site <- apply(cl,2,function(x)tapply(x,meta$site,mean))
  topc <- names(sort(colMeans(cl),decreasing=TRUE))[seq_len(min(14,ncol(cl)))]
  z <- scale(cls_site[,topc])
  pdf(file.path(fg,"Fig4A_class_heatmap.pdf"),width=7,height=6)
  pheatmap(t(z),clustering_method="ward.D2",color=colorRampPalette(c("#2166ac","white","#b2182b"))(60),
           main="Class z-score by Ferula site (final)",fontsize=8); dev.off(); cat("Fig4A done\n")
}

############ FIG 4B — genus biomarker dotplot (LEfSe substitute) ############
{
  bm <- read.delim(inputs[["genus_biomarkers"]], check.names=FALSE)
  sig <- bm[bm$KW_BH<0.05 | (bm$IndVal_p<0.05 & bm$IndVal>0.6),]
  sig$enriched <- factor(sig$enriched, levels=c("HD","XJ","DS","DG"))
  sig <- sig[order(sig$enriched, sig$IndVal),]; sig$genus <- factor(sig$genus, levels=sig$genus)
  p4b <- ggplot(sig, aes(IndVal, genus, color=enriched, size=mean_RA))+geom_point()+
    scale_color_brewer(palette="Dark2")+labs(x="IndVal (fidelity)",y=NULL,color="Enriched site",size="Mean RA (%)",
      title="Genus biomarkers (KW+BH & IndVal)")+theme_bw(base_size=8)
  ggsave(file.path(fg,"Fig4B_biomarkers.pdf"),p4b,width=14,height=12,units="cm"); cat("Fig4B done\n")
}

############ FIG 5A — dbRDA biplot + PERMANOVA ; 5B — rdacca.hp ############
bray <- vegdist(rare,"bray")
{
  ## Nested design: the three depths are sampled WITHIN each plant.
  ##  - Species/site -> tested at the PLANT level (depths aggregated; plant = unit of replication, n=12).
  ##  - Depth        -> tested as a WITHIN-PLANT factor (permutations restricted within plant).
  library(permute)
  meta$plant <- factor(paste0(as.character(meta$site), meta$replicate))
  plant_tab  <- rowsum(as.matrix(rare), group = as.character(meta$plant))
  meta_plant <- unique(meta[, c("plant","site")])
  meta_plant <- meta_plant[match(rownames(plant_tab), as.character(meta_plant$plant)), ]
  ## Reset seeds immediately before each permutation test so plot rendering or
  ## package internals cannot change the sampled permutation sequence.
  set.seed(20260609)
  ad_site  <- adonis2(vegdist(plant_tab, "bray") ~ site, data = meta_plant, permutations = 9999)
  set.seed(20260610)
  ad_depth <- adonis2(bray ~ depth_order, data = meta,
                      permutations = how(blocks = meta$plant, nperm = 9999))
  permtab <- data.frame(
    Df       = c(ad_site$Df[1],       ad_depth$Df[1]),
    SumOfSqs = c(ad_site$SumOfSqs[1], ad_depth$SumOfSqs[1]),
    R2       = c(ad_site$R2[1],       ad_depth$R2[1]),
    F        = c(ad_site$F[1],        ad_depth$F[1]),
    `Pr(>F)` = c(ad_site$`Pr(>F)`[1], ad_depth$`Pr(>F)`[1]),
    row.names = c("site","depth_order"), check.names = FALSE)
  write.table(permtab, file.path(tb,"permanova_final.tsv"), sep="\t", quote=FALSE, col.names=NA)
  # VIF selection
  sel<-env_vars; repeat{ db<-capscale(as.formula(paste("bray~",paste(sel,collapse="+"))),data=envz); v<-vif.cca(db); v<-v[is.finite(v)]
    if(length(v)==0||max(v)<=10||length(sel)<=2)break; sel<-setdiff(sel,names(which.max(v))) }
  db <- capscale(as.formula(paste("bray~",paste(sel,collapse="+"))),data=envz)
  r2 <- RsquareAdj(db)
  set.seed(20260611)
  trm <- anova.cca(db,by="term",permutations=999)
  write.table(as.data.frame(trm), file.path(tb,"dbrda_terms_final.tsv"), sep="\t", quote=FALSE, col.names=NA)
  format_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  cat(sprintf("Fig5A dbRDA: vars=%s ; R2=%.3f adjR2=%.3f ; PERMANOVA site R2=%.3f p=%s depth R2=%.3f p=%s\n",
      paste(sel,collapse=","), r2$r.squared, r2$adj.r.squared,
      permtab["site","R2"], format_p(permtab["site","Pr(>F)"]),
      permtab["depth_order","R2"], format_p(permtab["depth_order","Pr(>F)"])))
  st<-as.data.frame(scores(db,display="sites",choices=1:2)); st$sample<-rownames(st); st<-merge(meta[,c("sample","species_code","depth_order")],st,by="sample")
  bp<-as.data.frame(scores(db,display="bp",choices=1:2)); bp$var<-rownames(bp); sc<-.8*max(abs(st[,c("CAP1","CAP2")]))/max(abs(bp[,1:2]))
  p5a<-ggplot(st,aes(CAP1,CAP2,color=species_code,shape=depth_order))+geom_point(size=2.3)+
    geom_segment(data=bp,aes(0,0,xend=CAP1*sc,yend=CAP2*sc),inherit.aes=FALSE,arrow=arrow(length=unit(.14,"cm")),color="grey25")+
    geom_text_repel(data=bp,aes(CAP1*sc,CAP2*sc,label=var),inherit.aes=FALSE,size=2.5,color="grey15")+
    scale_color_brewer(palette="Dark2")+labs(color="Species/site",shape="Depth")+theme_bw(base_size=9)
  ggsave(file.path(fg,"Fig5A_dbRDA.pdf"),p5a,width=15,height=12,units="cm")
  # Fig5B rdacca.hp
  hp<-rdacca.hp(bray, envz[,sel], method="dbRDA", type="adjR2", scale=FALSE)
  hpd<-data.frame(var=rownames(hp$Hier.part), ind=hp$Hier.part[,"Individual"])
  write.table(hpd, file.path(tb,"rdaccahp_final.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  p5b<-ggplot(hpd,aes(reorder(var,ind),ind*100))+geom_col(fill="#4575b4",width=.7)+coord_flip()+
    labs(x=NULL,y="Individual effect (% adj R2)")+theme_bw(base_size=9)
  ggsave(file.path(fg,"Fig5B_varpart.pdf"),p5b,width=11,height=9,units="cm"); cat("Fig5A/B done\n")
}

############ FIG 5C — genus vs soil Spearman heatmap ############
{
  g<-agg("Genus"); topg<-names(sort(colMeans(g),decreasing=TRUE))[seq_len(min(30,ncol(g)))]
  G<-g[,topg]; rho<-matrix(NA,length(topg),length(env_vars),dimnames=list(topg,env_vars)); pv<-rho
  for(i in topg)for(j in env_vars){ct<-suppressWarnings(cor.test(G[,i],meta[[j]],method="spearman",exact=FALSE));rho[i,j]<-ct$estimate;pv[i,j]<-ct$p.value}
  star<-ifelse(pv<0.001,"***",ifelse(pv<0.01,"**",ifelse(pv<0.05,"*","")))
  write.table(data.frame(genus=rownames(rho),rho,check.names=FALSE), file.path(tb,"Fig5C_genus_soil_final.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  pdf(file.path(fg,"Fig5C_genus_soil.pdf"),width=9,height=10)
  pheatmap(rho,clustering_method="ward.D2",color=colorRampPalette(c("#2166ac","white","#b2182b"))(60),breaks=seq(-1,1,length.out=61),
           display_numbers=star,fontsize=8,fontsize_number=7,main="Top genera vs soil (Spearman, Ward.D2) - final"); dev.off()
  cat("Fig5C done\n")
}
cat("\nALL DONE\n")
