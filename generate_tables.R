require(vtrackR)
require(GenomicRanges)

base = c("barcode","chr","strand","pos","ipcr")
rep22 = c(base, "gdna1", "gdna2", "rna1", "rna2")
rep12 = c(base, "gdna1", "rna1", "rna2")
rep21 = c(base, "gdna1", "gdna2", "rna1")
rep11 = c(base, "gdna1", "rna1")

fnames = list(
p5 = "p5/p5_RB_7081_CGATGT_iPCR_insertions.txt",
p8 = c("p8/p8_RA_7077_ACTTGA_iPCR_insertions.txt",
       "p8/p8_RB_7082_TTAGGC_iPCR_insertions.txt"),
p38 = c("p38/p38_10K_7233_TTAGGC_iPCR_insertions.txt",
        "p38/p38_20K_7239_GATCAG_iPCR_insertions.txt"),
p39 = c("p39/p40_10K_7235_ACAGTG_iPCR_insertions.txt",
        "p39/p39_20K_7240_TAGCTT_iPCR_insertions.txt"),
p40 = c("p40/p39_10K_7234_TGACCA_iPCR_insertions.txt",
        "p40/p40_20K_7241_GGCTAC_iPCR_insertions.txt"),
p41 = c("p41/p41_10K_7236_GCCAAT_iPCR_insertions.txt",
        "p41/p41_20K_7242_CTTGTA_iPCR_insertions.txt"),
p14 = c("p14/p14_RA_7080_GGCTAC_iPCR_insertions.txt",
        "p14/p14_RB_7085_GCCAAT_iPCR_insertions.txt")
)

avail = list(
   p5  = list(rep22),
   p8  = list(rep22, rep22),
   p38 = list(rep11, rep11),
   p39 = list(rep11, rep11),
   p40 = list(rep11, rep11),
   p41 = list(rep11, rep11),
   p14 = list(rep21, rep21)
)

dflist = list()
total_insertions = 0
for (i in 1:length(fnames)) {

   prom = names(fnames)[i]
   reps = avail[[i]]

   for (j in 1:length(fnames[[i]])) {

      fname = fnames[[i]][j]
      reps = avail[[i]][[j]]

      print (fname)

      dat = read.table(fname, comment.char="#")
      colnames(dat) = reps

      total_insertions = total_insertions +
         sum(!dat$chr %in% c("spike", "pT2"))

      # Get the pike unit, add 1/2 unit to RNA.
      rna1Spk = with(dat, mean(rna1[chr == "spike" & rna1 > 0]))
      gdna1Spk = with(dat, mean(gdna1[chr == "spike" & gdna1 > 0]))
      dat$rna1 = dat$rna1 + rna1Spk/2

      # Keep only integrations in Drosophila genome,
      # and with positive gDNA counts.
      if (identical(reps, rep22)) {
         # If there is a second gDNA replicate, get the spike unit.
         gdna2Spk = with(dat, mean(gdna2[chr == "spike" & gdna2 > 0]))
         # If there is an RNA replicate, add 1/2 spike unit.
         rna2Spk = with(dat, mean(rna2[chr == "spike" & rna2 > 0]))
         dat$rna2 = dat$rna2 + rna2Spk/2
      }
      if (identical(reps, rep12)) {
         # If there is an RNA replicate, add 1/2 spike unit.
         rna2Spk = with(dat, mean(rna2[chr == "spike" & rna2 > 0]))
         dat$rna2 = dat$rna2 + rna2Spk/2
      }
      if (identical(reps, rep21)) {
         # If there is a second gDNA replicate, get the spike unit.
         gdna2Spk = with(dat, mean(gdna2[chr == "spike" & gdna2 > 0]))
      }


      if (identical(reps, rep12) || identical(reps, rep11)) {
         # If there is only one gDNA replicate, keep only the
         # rows with positive count.
         dat = subset(dat,
                !(chr %in% c("spike", "pT2")) & gdna1 > 0)
         dat$gdna1 = dat$gdna1 + gdna1Spk/2
      }
      else {
         # If there are two gDNA replicates, keep only the rows where
         # both are positive.
         dat = subset(dat,
                !(chr %in% c("spike", "pT2")) & gdna1+gdna2 > 0)
         dat$gdna1 = dat$gdna1 + gdna1Spk/2
         dat$gdna2 = dat$gdna2 + gdna2Spk/2
      }


      # Normalize expression.
      if (identical(reps, rep22)) {
         nexp = with(dat,
            log2(rna1/sum(rna1) + rna2/sum(rna2)) -
            log2(gdna1/sum(gdna1) + gdna2/sum(gdna2))
         )
      }
      else if (identical(reps, rep12)) {
         # The bottom term should be adjusted, but this
         # is not necessary due to mean-standardization.
         nexp = with(dat,
            log2(rna1/sum(rna1) + rna2/sum(rna2)) -
            log2(gdna1/sum(gdna1))
         )
      }
      else if (identical(reps, rep21)) {
         # The top term should be adjusted, but this
         # is not necessary due to mean-standardization.
         nexp = with(dat,
            log2(rna1/sum(rna1)) -
            log2(gdna1/sum(gdna1) + gdna2/sum(gdna2))
         )
      }
      else if (identical(reps, rep11)) {
         nexp = with(dat,
            log2(rna1/sum(rna1)) -
            log2(gdna1/sum(gdna1))
         )
      }

      dat$nexp = round(scale(nexp, center=TRUE, scale=FALSE), 4)
      dat$prom = prom
      dat$rep = j

      dflist[[fname]] = dat[,c("barcode", "chr", "strand",
         "pos", "nexp", "prom", "rep")]

   }

}

print(total_insertions)

ins = Reduce(rbind, dflist)
ins = ins[order(ins$chr, ins$pos),]

write.table(ins, file="allprom_nochromP.txt",
   sep = "\t", quote = FALSE, row.names=FALSE)

chromP = read.delim("GSE36175_norm_aggregated_tiling_arrays.txt.gz")
chromP$chromosome = sub("^chr", "", chromP$chromosome)
colnames(chromP)[2] = "chr"

GRins = GRanges(Rle(ins$chr), IRanges(start=ins$pos, width=1))
GRchromP = makeGRangesFromDataFrame(chromP[,2:4])

# GATC fragments used for chromatin proteins overlap by 4
# nucleotides (always GATC), but are otherwise disjoint.
ov = as.matrix(findOverlaps(GRins, GRchromP))
ins = data.frame(ins[ov[,1],], chromP[ov[,2],5:ncol(chromP)])

write.table(ins, file="allprom.txt",
   sep = "\t", quote = FALSE, row.names=FALSE)

