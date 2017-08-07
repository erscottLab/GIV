### SHINY APP 


getwd()
setwd("/Users/rebeccamolinsky/Google Drive/viral_dna_damage/code/ERS_LAB/runApp/")
setwd("/Users/rebeccamolinsky/Google\ Drive/viral_dna_damage/ebola_data/samples/SRR1553419")

library(shiny)






ui <- fluidPage(
  fileInput(inputId = "bam",  label="Bam File"), ## restrict the bam file to the desired region
  fileInput(inputId = "fasta",  label="Fasta File"),
  fileInput(inputId = "Meta_Data",label= "Additional Data"),
  numericInput(inputId = "depth", label="Max Depth", value = 8000, min = 0, max = 8000),  
  sliderInput(inputId = "num", label="Min Base Quality", min = 0 , max = 50 , value = 0, ticks = 5 ),
  sliderInput(inputId = "map", label="Min Map Quality", min = 0 , max = 100 , value = 10, ticks = 5),
  
  plotOutput(outputId = "mainplot")
  #tableOutput(outputId = "g_df")
  #plotOutput(outputId = "logGIV")
  #plotOutput(outputId = "GIV")
  
)


server <- function(input,output) {
  options(shiny.maxRequestSize=600*1024^2) 
  
  GIV_Calculator = function()
    
  {
    
    bam_path = input$bam$datapath #"SRR1553419_fake_a.bam"
    fasta_path = input$fasta$datapath #"EBOV_2014_EM096.fa"
    
    #bam_path = "Aligned.sortedByCoord.out.bam"
    #fasta_path = "KM034562v1.fa"
    
    ## Function 1:
    importseq <- function(bam_path, fasta_path) {
      
      BAM <- BamFile(bam_path)
      p_param <- PileupParam( min_mapq= input$map, min_base_quality= input$num, distinguish_strands=TRUE, distinguish_nucleotides=TRUE, include_deletions=FALSE, include_insertions=FALSE, max_depth = input$depth )
      phred2ASCIIOffset(phred = integer(20), scheme= c("Illumina 1.8+"))
      sbp_first <- ScanBamParam(flag=scanBamFlag(isFirstMate=TRUE),what=c("rname", "strand", "pos", "seq", "qname"))
      sbp_second <- ScanBamParam(flag=scanBamFlag(isSecondMate=TRUE),what=c("rname", "strand", "pos", "seq", "qname"))
      
      refseq = read.fasta(file=fasta_path, seqtype="DNA", as.string=FALSE, set.attributes= TRUE)
      refseq_name = names(refseq)
      ref_length = length(refseq[[refseq_name]])
      ref = refseq[[refseq_name]][1:ref_length]
      refseq = data_frame(ref)
      refseq$ref = toupper(refseq$ref)
      refseq$pos = seq(1,ref_length)
      
      
      
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
      
      giv_results = cbind(vt, GIV, bam_path, fasta_path)
      
      return(giv_results)
    }
    
    
    
    
    
    
    ### Function 4: command line options
    
    #parser <- ArgumentParser()
    #parser$add_argument("--b")
    #parser$add_argument("--f")
    #parser$add_argument("--md", type="integer", default= 8000)
    #parser$add_argument("--mbq", type="integer", default = 0)
    #parser$add_argument("--mmap", type="integer", default= 10)
    #parser$add_argument("--mcd", type="integer")
    #parser$add_argument("--mad", type="integer")
    #parser$add_argument("--iq")
    #parser$add_argument("--lb")
    #parser$add_argument("--qb")
    #parser$add_argument("--cb")
    #parser$add_argument("--minc", type="integer", default = 1)
    #parser$add_argument("--maxc", type="integer", default = 100)
    #parser$add_argument("--o")
    
    
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
    #args = parser$parse_args()
    #bam_path = args$b
    #fasta_path = args$f
    
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
    
    giv_l = list()
    for (vt in names(vartypes_rc)) {
      
      temp_giv <- GIVscore(fwd_nuc_colsums, fwd_bc_colsums,
                           rev_nuc_colsums, rev_bc_colsums,
                           vt)
      giv_l[[vt]] <- temp_giv
    }
    giv_df = do.call("rbind", giv_l)
    giv_df = as.data.frame(giv_df)
    giv_df$GIV = as.numeric(as.character(giv_df$GIV))
    return(giv_df)
  }
  
  
  
  
  
  output$mainplot <- renderPlot( {
  output$g_df <- GIV_Calculator()
  plot(giv$vt, giv$GIV)
    })  #renderplot
  
  
  #output$nuc_colsums <- renderDataTable(table(nuc_colsums))
  
  
  #output$bc_colsums <- renderDataTable(table(bc_colsums))
  
  
  #output$logGIV <- renderPlot(plot(log(GIV_Calculator)))       
  
  #output$GIV <- renderPlot(plot(GIV_Calculator))
  
  #giv$logGIV <- log(giv$GIV)  
  
  #output$logGIV <- renderPlot ({  
  #  plot(giv$vt, giv$logGIV)
    
  #})
    
    
  }  #server

  
   
  
  


shinyApp(ui = ui, server = server)




