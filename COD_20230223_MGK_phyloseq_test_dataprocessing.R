library(dada2); packageVersion("dada2") # 1.16 or later
library(DECIPHER); packageVersion("DECIPHER")
library(ShortRead)
library(Biostrings)
library(tidyverse)
library(ggplot2)
library(ggtext)
path <- "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20170815_SICAS2_filter_dust_phyloseq_test_of_skills/data/fastq_raw"
list.files(path)


# It seems it is generaged by cellranger-arc mkfastq software. In that case,outputs are
#I1: Dual index i7 read 
#R1: Read 1
#R2: Dual index i5 read
#R3: Read 2 
#refer
#https://www.biostars.org/p/9540021/#9540044

#or 
#https://resources.qiagenbioinformatics.com/manuals/biomedicalgenomicsanalysis/2110/index.php?manual=Illumina_Custom_Reads.html
#For example, to import paired reads from R1, R2, R3 fastq files where R1 has forward reads, R3 has reverse reads, and R2 has molecular indices, set the custom reads options to 'R2 R1, R3'. The imported paired sequence list can now be used by tools that support symbols at the begining of read 1 of paired reads.


#Demultiplexing conducted for the raw files. note demultiplexing at 
#/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20170815_SICAS2_filter_dust_phyloseq_test_of_skills/data/demultiplexed/readme.rtf



# Files original and decoded 

barcodeFile= "/Users/minsikkim/idemp/MGK/barcode.txt"#, sep = "\t", check.names = F)%>% select(c("BarcodeSequence", "#SampleID"))
I1File="/Users/minsikkim/idemp/MGK/Raw_Read2_Barcodes.fq"
R1File="/Users/minsikkim/idemp/MGK/Raw_Read1.fq" 
R2File="/Users/minsikkim/idemp/MGK/Raw_Read3.fq"
decodedFile=list.files("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20170815_SICAS2_filter_dust_phyloseq_test_of_skills/data/demultiplexed", pattern="*.fastq.gz$")

# Check decoded files -----------------------------------------------------
#### Read in barcode table, sequence names and barcodes ####
read.csv(barcodeFile, sep = "\t", check.names = F)
firstLine=readLines(barcodeFile,1)
barcodeTable=read.table(barcodeFile, header=grepl("arcode", firstLine))
### Read in seqquence names and barcodes ###
barcodeReads=readFastq(I1File)
I1Names=unname(sapply(as.character(id(barcodeReads)), function(x) (strsplit(x," ")[[1]][1] ) ) )
I1Barcodes=as.character(sread(barcodeReads))

#### Loop through decoded files
setwd("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20170815_SICAS2_filter_dust_phyloseq_test_of_skills/data/demultiplexed")

for( i in 1:length(decodedFile) ) {
        cat(i, " of ", length(decodedFile), "\n")
        cat(decodedFile[i],"\n")
        decodedReads=readFastq(decodedFile[i])
        seqNames=unname(sapply(as.character(id(decodedReads)), function(x) (strsplit(x," ")[[1]][1] ) ) )
        seqBarcodes=as.character(sread(decodedReads))
        
        idx=match(seqNames, I1Names)
        #table(I1Barcodes[idx])
        #sort(table(I1Barcodes[idx]), decreasing=T)
        bctable=sort(table(I1Barcodes[idx]), decreasing=T)
        if ( length(bctable)> 1000 ) {
                message(length(bctable), " barcodes found out of ", length(seqNames), " reads.") 
        }
        if ( length( grep(names(bctable)[1], barcodeTable[,1]) )==0 ) {
                message("most frequent barcode ",names(bctable)[1], " not matched.") 
                next
        }
        if ( length(bctable)> 1000 ) next
        
        cat(bctable,"\n")
        codeMax=names(bctable)[1]
        cat(codeMax,"\n")
        if( sort(table(I1Barcodes[idx]), decreasing=T)[1] < length(idx)*0.9 ) {
                message("suspicious decoding\n")
                message(sort(table(I1Barcodes[idx]), decreasing=T),"\n","Total reads:",length(idx),"\n")
        }
        sampleid = as.character(barcodeTable[barcodeTable[,1]==codeMax,2][1])
        print( sampleid )
        print( grepl(sampleid, decodedFile[i]) )  
        if ( ! grepl(sampleid, decodedFile[i]) ) 
                message("SampleID not in file name\n", sampleid, "\n", decodedFile[i], "\n")
        ##if ( i>10 ) break; 
}

q("no")


# DADA2 pipeline ----------------------------------------------------------

######################################################################
#####  GOAL: Choose and validate appropriate dada2 parameters    #####
#####         for processing  paired-end Iluumina 16S data       #####
######################################################################
#####   DATA is taken from the following paper                   #####
#
# Pyrethroid exposure alters internal and cuticle surface bacterial communities in Anopheles albimanus, ISMEJ, 2019.
# https://doi.org/10.1038/s41396-019-0445-5
# Sequencing data: First 18 samples from https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA512122
#


# Primer check ------------------------------------------------------------
#Get file path
path <- getwd()
abolute_path <- paste(getwd(), "data_dada2", sep = "/")

# Read in the forward and reverse fastq filenames
fnFs <- list.files(pattern="Read1.fq", full.names=TRUE)
fnRs <- list.files(pattern="Read3.fq", full.names=TRUE)
head(fnFs)
        # Check - Do those filenames look like what we expect?

plotQualityProfile(fnFs[1:2]) #seqeucnign set of forward is similar
plotQualityProfile(fnRs[1:2]) #seqeucnign set of forward is similar
#Read length: 250
#Read quality: good for forward. Bad for reverse

#What is amplicon length?
#16S V4

#Read3: GTGCCAGCMGCCGCGGTAA = 515 F
#Read1: GGACTACHVGGGTWTCTAAT - 806 R
# Expected mplicon length: 252
# --> We can throw out many overlapping reads.

# Define the paths to filtered files we are going to create
filtFs <- file.path(path, "filtered", basename(fnFs)) 
filtRs <- file.path(path, "filtered", basename(fnRs))
        # The filtered files will be in the `filtered/` subdirectory within `path`

allOrients <- function(primer) {
# Create all orientations of the input sequence
        require(Biostrings)
        dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                     RevComp = reverseComplement(dna))
        return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients("GGACTACHVGGGTWTCTAAT")
REV.orients <- allOrients("GTGCCAGCMGCCGCGGTAA")
FWD.orients
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
# Counts number of reads in which the primer is found
        nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
        return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#The linker primer is sequenced
                                        #However, the  16S pirmer target sequence is the part of the actual seqeuence.
                                        #Mihai Pop - we don't need to remove the region.
# Here linker primer sequence is not the part of biological sequence. need to be removed


#refer https://github.com/benjjneb/dada2/issues/1422


# DADA2 -------------------------------------------------------------------



# filter length -----------------------------------------------------------


#Sequence length
#16s v3-v4 regions were amplified
# Perform filtering and trimming
out <- filterAndTrim(fnFs.filtN, filtFs, fnRs.filtN, filtRs, maxEE=2, 
                     trimLeft=c(0, 0), # REPLACE XXX/YYY with proper parameter choices
                     truncLen=c(200, 150)) # REPLACE XXX/YYY with proper parameter choices
out
#out <- list.files(path = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20170815_SICAS2_filter_dust_phyloseq_test_of_skills/data/demultiplexed/filtered", pattern="fastq.gz", full.names=TRUE)



# Error calculation -------------------------------------------------------


# Were most reads retained during filtering? 
# How might the depth of sequencing in this data affect the questions that can be addressed?

# Learn the error model from the filtered data.
errF <- learnErrors(filtFs, multi=TRUE) # `multi=TRUE` activated multithreading to reduce computation time
errR <- learnErrors(filtRs, multi=TRUE)

# Visualize the error model. Points are observations, black line is fitted error model.`
plotErrors(errF, nominalQ = T)
plotErrors(errR, nominalQ = T)
                # Do the fitted error models look reasonable?
#Black line (estimated error rates) are alining with the actual error (black dots).
#Plus error rate gets tends to get smaller when QC score gets higher.
#Fitting looks reasonable


# Run the DADA2 method using the fitted error model.
# dada2 will result in denoised reads based on the predicted error

#ddF <- dada(filtFs, errF, pool=F, multi=TRUE) # for pseudo-pooling
ddF <- dada(filtFs, errF, pool=FALSE, multi=TRUE)
#ddR <- dada(filtRs, errF, pool=F, multi=TRUE) # for pseudo-pooling
ddR <- dada(filtRs, errR, pool=FALSE, multi=TRUE)
                # What pooling option makes sense? FALSE (default), "pseudo", or TRUE? See ?dada for more info
#pseudo-pooling can fill rare taxa undetected. Here, We do not proceed with pseudo-pooling
                #Pooling results in slightly more unique sequences relative to the no pooling.
                # For more pseudo-pooling detail: https://benjjneb.github.io/dada2/pseudo.html#pseudo-pooling

# Merge the denoised forward and reverse reads together.
mm <- mergePairs(ddF, filtFs, ddR, filtRs, verbose=TRUE)
head(mm[[1]])

        # Were most reads retained during merging. If not, why not?
#As some of the reads are overlapped; forward and reverse reads should be complementary 

# Construct a sequence table: rows are samples, columns are ASVs, values are abundances.
sta <- makeSequenceTable(mm)
dim(sta)
# How many samples and ASVs are in this table?
# Inspect distribution of sequence lengths
table(nchar(getSequences(sta)))


# Remove chimeric ASVs and construct a new chimera-free sequence table.
st <- removeBimeraDenovo(sta, multi=TRUE, verbose=TRUE)
sum(st)/sum(sta)

        # Were most reads retained during chimera removal? How about ASVs? 
# 3.3 % of chimera were removed
        # How does the fraction of chimeric ASVs compare to the fraction of chimeric reads?
        # Why is that?

####################################################################
######  Inspect the number of reads passing through each step ######
######  THIS IS THE SINGLE MOST IMPORTANT SANITY CHECK!!      ######
####################################################################

# Code derived from the dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(ddF, getN), sapply(ddR, getN), sapply(mm, getN), rowSums(st))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- basename(fnFs)

write_csv(track %>% as.data.frame(), "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/1_Daily_Log/MGK/DAT_20230224_MGK_phyloseq_test_dataprocessing_reads_QC.csv")

        # In this case, most reads should make it through the entire pipeline!
        # Most importantly, a large majority (>80% of reads) should merge successfully,
        #  and almost all (>95%)  reads should pass chimera filtering.
        # IF THAT ISN'T THE CASE, you have a problem, and need to revisit your truncation lengths
        #   (merging problem) or primer removal (trimLeft, chimera problem).
#This read is having low quality (about 50-70% reads survived). However, there was high overlap (amplicon: 291. F+R = 500), lower merged reads seems reasonable.



# The dada2 taxonomic reference page https://benjjneb.github.io/dada2/training.html has links to a 
#   number of reference databases formatted to work with `assignTaxonomy`.Ben Callahan recommends Silva for 16S.


#silva database: https://benjjneb.github.io/dada2/training.html
tax <- assignTaxonomy(st, "/Users/minsikkim/Downloads/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
unname(head(tax))
str(tax)
tax

#Alternatives: The recently developed IdTaxa taxonomic classification method is also available via the DECIPHER Bioconductor package. The paper introducing the IDTAXA algorithm reports classification performance that is better than the long-time standard set by the naive Bayesian classifier. Here we include a code block that allows you to use IdTaxa as a drop-in replacement for assignTaxonomy (and it’s faster as well!). Trained classifiers are available from http://DECIPHER.codes/Downloads.html. Download the SILVA SSU r132 (modified) file to follow along.


#dna <- DNAStringSet(getSequences(st)) # Create a DNAStringSet from the ASVs
#load("/Users/minsikkim/Downloads/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
#ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
#ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
#taxid <- t(sapply(ids, function(x) {
#        m <- match(ranks, x$rank)
#        taxa <- x$taxon[m]
#        taxa[startsWith(taxa, "unclassified_")] <- NA
#        taxa
#}))
#colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)







        # Challenge 2: Use `addSpecies` or `assignSpecies` to do species-level assignemtn, where appropriate.
        #tax <- addSpecies(tax2, "/Users/MK1059/Dropbox (Personal)/Macbook Pro 16 2021/Downloads/silva_species_assignment_v138.1.fa")
        #species matchign - exact match
        ##sometimes 99% is used in publications
#No match for exact sequences...no species assigned
# should change the variables






# Generating phyloseq object ----------------------------------------------

# Phyloseq is a package for the manipulation and analysis of microbiome data.
# Here we use it briefly to produce an ordination of our sequenced communities.
library(phyloseq); library(ggplot2)

# We define a very simple data.frame that records the 3 experimental groups these samples 
#   came from (see paper for more info)
samdf <- 
        data.frame(row.names = rownames(st), 
                    sample_id = rownames(st) %>% gsub("Raw_Read1.fq_Phipatanakul.|.fastq.gz", "", .))

#changing sample names
ps <- phyloseq(sample_data(samdf), otu_table(st, taxa_are_rows=FALSE), tax_table(tax))
sample_names(ps) <- ps %>% sample_data() %>% .$sample_id

#changing taxa name
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/1_Daily_Log/MGK/DAT_20230224_MGK_phyloseq_test_dataprocessing.R")
ps <- readRDS("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/1_Daily_Log/MGK/DAT_20230224_MGK_phyloseq_test_dataprocessing.R")

#mock community
ps_mock <- subset_samples(ps, sample_id == "Mock") %>% prune_taxa(taxa_sums(.) != 0, .)
otu_table_mock <- ps_mock %>% otu_table() %>% data.frame()



#plot mock community
plot_bar(transform_sample_counts(ps_mock,function(x){x/sum(x)}), fill="Genus") +
        ylab("Relative abundnace")


#Comparing with old phyloseq 

ps_old <- import_biom("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Analysis/phyloseq_test_of_skills/data_test/OTU_Table_16S.biom")

ps_old_mock <- merge_phyloseq(otu_table(ps_old) %>% data.frame() %>% select("Phipatanakul.Mock") %>% otu_table(taxa_are_rows = T),
                              gsub(tax_table(ps_old)[, colnames(tax_table(ps_old))], pattern = "__", replacement = "") %>% 
                                      data.frame() %>%
                                      rename(Kingdom = Rank1, Phylum = Rank2, Class = Rank3, Order = Rank4, Family = Rank5, Genus = Rank6) %>% as.matrix() %>% tax_table) %>% prune_taxa(taxa_sums(.) != 0, .)
sample_names(ps_old_mock) <- ps_old_mock %>% sample_names() %>% gsub("Phipatanakul.", "", .) 

#separate plots
plot_bar(transform_sample_counts(ps_old_mock,function(x){x/sum(x)}), fill="Genus") +
        ggtitle("Old phyloseq object") +
        ylab("Relative abundnace") +
        theme(plot.title = element_text(size = 30), axis.text.x = element_blank())

plot_bar(transform_sample_counts(ps_mock,function(x){x/sum(x)}), fill="Genus") +
        ylab("Relative abundnace") +
        ggtitle("New phyloseq object") +
        theme(plot.title = element_text(size = 30), axis.text.x = element_blank())

#Merged plots
transform_sample_counts(ps_mock,function(x){x/sum(x)})



compare_bar <- cbind("phyloseq - OTU clustered",
                     transform_sample_counts(ps_old_mock, function(x){x/sum(x)}) %>% tax_table() %>% data.frame() %>% .$Genus,
                     transform_sample_counts(ps_old_mock, function(x){x/sum(x)}) %>% otu_table() %>% data.frame() %>% .$Mock)
compare_bar <- rbind(compare_bar, cbind("phyloseq - ASV analyzed",
                           transform_sample_counts(ps_mock, function(x){x/sum(x)}) %>% tax_table() %>% data.frame() %>% .$Genus,
                           transform_sample_counts(ps_mock, function(x){x/sum(x)}) %>% otu_table() %>% t() %>% data.frame() %>% .$Mock))

compare_bar <- rbind(compare_bar, cbind("Theoretical",
                                        c("Pseudomonas", "Escherichia", "Salmonella", "Lactobacillus", "Enterococcus", "Staphylococcus", "Listeria", "Bacillus"),
                                        c(0.042, 0.101, 0.104, 0.184, 0.099, 0.155, 0.141, 0.174))) %>% data.frame() %>% rename(data = X1, Genus = X2, abundance = X3)
compare_bar$abundance <- compare_bar$abundance %>% as.numeric()
compare_bar %>% subset(compare_bar$data == "Old phyloseq") %>% select(abundance) %>% colSums()
compare_bar$Genus <- paste("*", compare_bar$Genus, "*", sep = "")
compare_bar$Genus
# Stacked
compare_bar$data=factor(compare_bar$data, levels = c("Theoretical", "phyloseq - OTU clustered", "phyloseq - ASV analyzed")) 
        
ggplot(compare_bar, aes(fill = Genus, y = abundance, x = data)) + 
        theme_classic(base_size = 15, base_family = "serif") +
        geom_bar(position="stack", stat="identity") +
        facet_wrap(~data, scales = "free", ) +
        theme(legend.text = element_markdown()) +
#        scale_x_discrete(limits = c("Theoretical", "Old phyloseq", "New phyloseq")) +
        ylab("Relative abundance") +
        xlab("Type of data")


compare_bar_unified_taxa <- compare_bar
compare_bar_unified_taxa$Genus <- ifelse(compare_bar_unified_taxa$Genus == "*Escherichia-Shigella*", "*Escherichia*",
                                         ifelse(compare_bar_unified_taxa$Genus == "*Limosilactobacillus*", "*Lactobacillus*",
                                                compare_bar_unified_taxa$Genus)
                                         )


ggplot(compare_bar_unified_taxa, aes(fill = Genus, y = abundance, x = data)) + 
        theme_classic(base_size = 12, base_family = "serif") +
        geom_bar(position="stack", stat="identity") +
        facet_wrap(~data, scales = "free", ) +
        theme(legend.text = element_markdown()) +
        #        scale_x_discrete(limits = c("Theoretical", "Old phyloseq", "New phyloseq")) +
        ylab("Relative abundance") +
        scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
        xlab("Type of data")

ggsave("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Analysis/phyloseq_test_of_skills/data_test/REP_20230224_MGK_phyloseq_testskill_dada2_compare.pdf",
        plot = last_plot(),
        scale = 1,
        width = 8,
        height = 6,
        units = c("in"),
        dpi = 300)
#Cross check samples

cbind({ps %>% sample_names()}[order({ps %>% sample_names()}, decreasing = T)],
      {ps_old %>% sample_names()}[order({ps_old %>% sample_names()}, decreasing = T)])
ps_newone <- subset_samples(ps, sampleID == "SICAS2.021T") %>% prune_taxa(taxa_sums(.) != 0, .)


#This sample was missing in the previous set.... Double-check the total read counts

sample_data(ps)$sequencing_depth <- otu_table %>% rowSums()
sample_data(ps)

cbind(order({ps %>% sample_names()}, decreasing = T), order({ps_old %>% sample_names()}, decreasing = T))
# Use phyloseq to plot a bray-curtis NMDS odrination of our samples, colored by treatment.
plot_ordination(ps, ordinate(ps, method="CAP", distance="bray"), color="treatment") + 
        aes(size=4) + theme_bw()




