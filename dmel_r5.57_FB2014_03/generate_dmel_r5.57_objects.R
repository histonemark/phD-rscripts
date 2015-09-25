require(vtrackR)
require(GenomicRanges)

GFF_names = c("seqnames", "source", "feature",
   "start", "end", "score", "strand", "frame", "attr")

# Read exons and genes from GFF files.
exons = read.table("dmel-all-exons-r5.57.gff", sep="\t", quote="")
genes = read.table("dmel-all-genes-r5.57.gff", sep="\t", quote="")

colnames(exons) = GFF_names
colnames(genes) = GFF_names

# Get expression from RNA seq obtained in the lab.
Kc_expression = read.table("Kc_exp_color_MC.txt")
active_FBgn = with(Kc_expression, V7[V5 > median(V5)])

# Use it to separate active versus inactive genes.
active = sub('ID=(FBgn[0-9]{7});.*', '\\1', genes$attr) %in% active_FBgn
act_genes = subset(genes, active)
inact_genes = subset(genes, !active)

active_names = sub('.*?Name=([^;]*);.*', '\\1', act_genes$attr)
active = sub('Name=([^:]*):.*', '\\1', exons$attr) %in% active_names
inactive = !sub('Name=([^:]*):.*', '\\1', exons$attr) %in% active_names

act_exons = subset(exons, active)
inact_exons = subset(exons, inactive)

exons_r5.57 = vtag(makeGRangesFromDataFrame(exons, keep=TRUE))
act_exons_r5.57 = vtag(makeGRangesFromDataFrame(act_exons, keep=TRUE))
inact_exons_r5.57 = vtag(makeGRangesFromDataFrame(inact_exons, keep=TRUE))

genes_r5.57 = vtag(makeGRangesFromDataFrame(genes, keep=TRUE))
act_genes_r5.57 = vtag(makeGRangesFromDataFrame(act_genes, keep=TRUE))
inact_genes_r5.57 = vtag(makeGRangesFromDataFrame(inact_genes, keep=TRUE))

introns_r5.57 = vtag(setdiff(genes_r5.57, exons_r5.57))
act_introns_r5.57 = vtag(setdiff(act_genes_r5.57, exons_r5.57))
inact_introns_r5.57 = vtag(setdiff(inact_genes_r5.57, exons_r5.57))

promoters_r5.57 = vtag(promoters(genes_r5.57), up=1000, down=0)
act_promoters_r5.57 = vtag(promoters(act_genes_r5.57), up=1000, down=0)
inact_promoters_r5.57 = vtag(promoters(inact_genes_r5.57), up=1000, down=0)

# Save it all as R objects.
save(exons_r5.57, file="exons_r5.57.rda")
save(act_exons_r5.57, file="act_exons_r5.57.rda")
save(inact_exons_r5.57, file="inact_exons_r5.57.rda")

save(genes_r5.57, file="genes_r5.57.rda")
save(act_genes_r5.57, file="act_genes_r5.57.rda")
save(inact_genes_r5.57, file="inact_genes_r5.57.rda")

save(introns_r5.57, file="introns_r5.57.rda")
save(act_introns_r5.57, file="act_introns_r5.57.rda")
save(inact_introns_r5.57, file="inact_introns_r5.57.rda")

save(promoters_r5.57, file="promoters_r5.57.rda")
save(act_promoters_r5.57, file="act_promoters_r5.57.rda")
save(inact_promoters_r5.57, file="inact_promoters_r5.57.rda")
