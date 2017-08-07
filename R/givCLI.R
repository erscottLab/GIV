### --qualityscore 35 (DEFAULT 30)

suppressMessages(library(devtools))
suppressMessages(library(Rsamtools))
suppressMessages(library(seqinr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(argparse))
suppressMessages(library(car))

GIV_Calculator = function()

{

 ## Function 1:
importseq <- function(bam_path, fasta_path) {

  BAM <- BamFile(bam_path)
  p_param <- PileupParam( min_mapq= args$mmap, min_base_quality= args$mbq, distinguish_strands=TRUE, distinguish_nucleotides=TRUE, include_deletions=FALSE, include_insertions=FALSE, max_depth = args$md )
  phred2ASCIIOffset(phred = integer(args$mbq), scheme= c("Illumina 1.8+"))
  sbp_first <- ScanBamParam(flag=scanBamFlag(isFirstMate=TRUE),what=c("rname", "strand", "pos", "seq", "qname"))
  sbp_second <- ScanBamParam(flag=scanBamFlag(isSecondMate=TRUE),what=c("rname", "strand", "pos", "seq", "qname"))

  refseq = read.fasta(file=fasta_path, seqtype="DNA", as.string=FALSE, set.attributes= TRUE)
  refseq_name = names(refseq)
  ref_length = length(refseq[[refseq_name]])
  ref = refseq[[refseq_name]][1:ref_length]
  refseq = data_frame(ref)
  refseq$ref = toupper(refseq$ref)
  refseq$pos = seq(1,ref_length)

  #res = pileup(BAM, pileupParam=p_param, scanBamParam = param)
  #pileup_df = merge(res, refseq, by="pos") #### ASK- why did we get ride of, now coverage_df will not run

  pileup_fwd = pileup(BAM, pileupParam=p_param, scanBamParam=sbp_first)
  pileup_fwd = merge(pileup_fwd, refseq, by="pos")
  pileup_fwd$nucleotide[pileup_fwd$strand=="-"] = recode(pileup_fwd$nucleotide[pileup_fwd$strand=="-"],"'A'='T';'T'='A';'G'='C';'C'='G'")
  pileup_fwd$ref[pileup_fwd$strand=="-"] = recode(pileup_fwd$ref[pileup_fwd$strand=="-"],"'A'='T';'T'='A';'G'='C';'C'='G'")
  #print(sum(pileup_fwd$count))
  #print(head(pileup_fwd))
  #print(dim(pileup_fwd))
  #print(' ')

  pileup_rev = pileup(BAM, pileupParam=p_param, scanBamParam=sbp_second)
  pileup_rev = merge(pileup_rev, refseq, by="pos")
  pileup_rev$nucleotide[pileup_rev$strand=="-"] = recode(pileup_rev$nucleotide[pileup_rev$strand=="-"],"'A'='T';'T'='A';'G'='C';'C'='G'")
  pileup_rev$ref = recode(pileup_rev$ref,"'A'='T';'T'='A';'G'='C';'C'='G'")
  #print(sum(pileup_rev$count))
  #print(head(pileup_rev))
  #print(dim(pileup_rev))


  #coverage_df = dcast(pileup_df, seqnames+pos ~ nucleotide, value.var="count", fun.aggregate=sum)
  #coverage_df = coverage_df[rowSums(coverage_df[,c("A", "T", "G", "C")]) >= args$minc , ]
  #coverage_df = coverage_df[rowSums(coverage_df[,c("A", "T", "G", "C")]) <= args$maxc , ]
  #pileup_df = pileup_df[rownames(coverage_df), ]

  #pileup_fwd = subset(pileup_df, res$strand=="+")
  #pileup_rev = subset(pileup_df, res$strand=="-")

  pileupList = list("pileup_fwd" = pileup_fwd, "pileup_rev" = pileup_rev)
  return(pileupList)
  }


#Function 2
nuc_colsums = function(pileup_df) {
  nuc_df = dcast(pileup_df, seqnames+pos ~ nucleotide, value.var="count", fun.aggregate=sum)
  nuc_colsums = colSums(nuc_df[,c("A", "T", "G", "C")])
  return(nuc_colsums)
  }

bc_colsums = function(pileup_df) {
  pileup_df$bc = paste(pileup_df$ref, pileup_df$nucleotide, sep="_")
  bc_df = dcast(pileup_df, seqnames+pos ~ bc, value.var="count", fun.aggregate=sum)

  basechanges =     c("G_T",
                     "C_A",
                     "A_T",
                     "A_G",
                     "A_C",
                     "G_A",
                     "G_C",
                     "C_G",
                     "C_T",
                     "T_A",
                     "T_G",
                     "T_C")

  bc_colsums = colSums(bc_df[, basechanges])
  return(bc_colsums)
  }


###Function 3:
  GIVscore = function(fwd_nuc_colsums, fwd_bc_colsums,
           rev_nuc_colsums, rev_bc_colsums,
           vt) {

  vt_rc = vartypes_rc[[vt]]
  vt1 = unlist(strsplit(vt, '_'))[1]
  vt1_rc = unlist(strsplit(vt_rc, '_'))[1]

  # Number of G to T variants in R1
  C1v = fwd_bc_colsums[vt]  #vt
  # Number of C to A variants in R2
  C2v = rev_bc_colsums[vt_rc]  #vt_rc

  #Total number of G in R1
  C1 = fwd_nuc_colsums[vt1]  #vt[1]
  #Total number of C in R2
  C2 = rev_nuc_colsums[vt1_rc]  #vt_rc[1]

  #Total number of C in R1
  C1_RC = fwd_nuc_colsums[vt1_rc]  #vt_rc[1]
  #Total number of G in R2
  C2_RC = rev_nuc_colsums[vt1]  #vt[1]

  # Number of C to A variants in R1
  C1v_RC = fwd_bc_colsums[vt_rc]  #vt_rc
  #Number of G to T variants in R2
  C2v_RC = rev_bc_colsums[vt]  #vt

  GIV = ((C1v + C2v) / (C1 + C2)) / ((C1v_RC + C2v_RC)/(C1_RC + C2_RC))

  giv_results = paste( c(vt, GIV, bam_path, paste(fasta_path, '\n', sep='')), sep=',')
  cat(giv_results, sep=',')
  #return(GIV)
  return()
  }


### Function 4: command line options

parser <- ArgumentParser()
parser$add_argument("--b")
parser$add_argument("--f")
parser$add_argument("--md", type="integer", default=8000)
parser$add_argument("--mbq", type="integer", default =0)
parser$add_argument("--mmap", type="integer", default=0)
parser$add_argument("--mcd", type="integer")
parser$add_argument("--mad", type="integer")
parser$add_argument("--iq")
parser$add_argument("--lb")
parser$add_argument("--qb")
parser$add_argument("--cb")
parser$add_argument("--minc", type="integer", default = 1)
parser$add_argument("--maxc", type="integer", default = 1000)
parser$add_argument("--o")


'usage: giv_data.R  [-h] [-md] [-mbq] [-mmap] [-mcd] [-mad] [-iq] [-i] [-lb] [-qb] [-cb] [-minc] [-max] [-o]

optional arguments:
    - h       help
    - md     max_depth
    - mbq    min_base_quality  [default 30]
    - mmap   min_mapq
    - mcd    min_nucleotide_depth
    - mad    min_minor_allele_depth
    - iq     ignore_query_Ns
    - lb     left_bins
    - qb     query_bins
    - cb     cycle_bins
    - minc   min_depth [DEFAULT 1]
    - maxc   max_depth [DEFAULT 100]
    - o      output file [default "giv_data.csv"]
'


  #args = commandArgs(trailingOnly = TRUE)
  args = parser$parse_args()
  bam_path = args$b
  fasta_path = args$f

  #mbq = args$mbq

  pileups = importseq(bam_path, fasta_path)

  vartypes_rc = list("G_T"="C_A",
                     "C_A"="G_T",
                     "A_T"="T_A",
                     "A_G"="T_C",
                     "A_C"="T_G",
                     "G_A"="C_T",
                     "G_C"="C_G",
                     "C_G"="G_C",
                     "C_T"="G_A",
                     "T_A"="A_T",
                     "T_G"="A_C",
                     "T_C"="A_G")

  rc =  list("A"="T",
             "G"="C",
             "T"="A",
             "C"="G")


 fwd_nuc_colsums = nuc_colsums(pileups$pileup_fwd)
 #print(fwd_nuc_colsums)
 fwd_bc_colsums = bc_colsums(pileups$pileup_fwd)
 #print(fwd_bc_colsums)

 rev_nuc_colsums = nuc_colsums(pileups$pileup_rev)
 #print(rev_nuc_colsums)
 rev_bc_colsums = bc_colsums(pileups$pileup_rev)
 #print(rev_bc_colsums)
  for (vt in names(vartypes_rc)) {

       GIVscore(fwd_nuc_colsums, fwd_bc_colsums,
                rev_nuc_colsums, rev_bc_colsums,
                vt)
       }
### add command line out (save output -o)


}

giv = GIV_Calculator()
