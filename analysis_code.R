## set up packages

library(readr)
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(lubridate)
library(tidyr)
library(stringr)

write.df <- function(dataframe, filename , sep="\t", quote = F, row.names = F, ...){
  write.table(dataframe, filename, sep = sep, quote = quote, row.names =F)
}

######################
## Set up data 
## The code assumes you have a data frame composed of:
# * MRN :unique identifier for the patient 
# * Normal : identifier for the sample sequenced
# * Major and minor contamination : these values are derived from the MSK-IMPACT assay and inidicate the level of contamination a sample might have. If your data doesn not have this, then ignore the part of the code that filters away mutations
# * VF_N : Variant allele fraction of the mutation in the normal
# * DP_N : Total number of reads at the mutation position in the normal
# * AD_N : Total number of reads supporting the mutation in the normal
# * VF_T : Variant allele fraction of the mutation in the tumor
# * DP_T : Total number of reads at the mutation position in the tumor
# * AD_T : Total number of reads supporting the mutation in the tumor
# * CosmicCount : Number of occurences of a mutation in the Cosmic database
# * COSMIC : description of the tumor types where Cosmic database has an entry for the mutation
# * PrcntOccuranceInNormals : For each mutation, what percent of a set of standard blood normals harbor the mutation. This is used for eliminating sequencing artifacts
# * OccuranceCountInCohort : Count of how many times a mutation is observed within the cohort


all_mutations <- all_mutations %>% filter(MinorCont <= 0.01) %>% filter(MajorCont < 0.5) 

all_mutations$Key <- with(all_mutations, paste(Chr, Pos, Ref, Alt, sep=":"))

genes.leukemia = c("ASXL1", "CBL", "DNMT3A", "GNAS", "JAK2", "NRAS", "SF3B1", "TP53", "U2AF1", "SETBP1", "IDH2", "BCOR", "PPM1D", "TET2", "IDH1", "IDH2", "SRSF2", "RUNX1", "SH2B3", "ZRSR2", "PHF6", "STAT3", "KRAS", "MYD88", "ATM", "CALR", "CEBPA", "ETV6", "EZH2", "FLT3", "KIT", "MPL", "NPM1", "SRSF2", "STAG2", "WT1", "SETD2", "CREBBP")

#these blacklisted mutations were identified as false positives based on IGV review of the mapped sequencing data
blacklist_key <- c("20:31023821:G:T", "7:151945054:T:TC", "20:30946591:CAGA:C", "7:151875096:C:CT", "22:29091840:T:C", "1:120612002:CGG:C", "1:27105930:TG:T", "1:120611960:C:T", "X:66765158:TGCAGCAGCAGCA:T", "X:66765155:T:TGCA", "X:66765158:T:TGCAGCA", "6:29911227:G:C", "16:11349151:GCGTGCGAACGGAATGTGCGGAAGTGCGTGTCGCCGGGGGCCGGGGCCGGGACCGCGGGGCACGGCCGC:G", "16:11349205:GCGGGGCACGGCCGCGGGCGCGCGGGGGC:G")
all_mutations <- all_mutations %>% filter(!Key %in% blacklist_key)

########
## Determine CH-PD mutations

# threshold for presence in cosmic
cosmic <- 5

#initial filering of variants based on variant genotypes
mutations <- all_mutations %>% 
  filter(VF_N >= 0.1) %>%
  filter(VF_N < 0.35) %>%
  filter(VF_T < 0.35) %>%
  filter(DP_N > 20) %>%
  filter(DP_T > 20) %>%
  filter(AD_N >=8) %>%
  filter(VF_N/VF_T > 2)

# high VAF variants that are hotspots
high_freq_mutations <- all_mutations %>% 
  filter(Gene %in% genes.leukemia) %>%
  filter(VF_N >= 0.35) %>%
  filter(DP_N > 20) %>%
  filter(DP_T > 20) %>% 
  filter(AD_N >= 8) %>%
  filter(VF_N/VF_T > 2) %>%
  filter(CosmicCount > 10) 

# select mutations that appear in COSMIC 
Cosmic.hotspot <- mutations[grepl("haematopoietic", mutations$COSMIC), ]
Cosmic.hotspot <- Cosmic.hotspot %>% 
  filter(CosmicCount >= cosmic)

# select mutations that have roles in hematopoietic disease development based on literature
exons <- paste("exon", seq(7, 23, 1), sep="")
DNMT3A.disruptive <- mutations %>% filter(Gene == "DNMT3A") %>% filter(Exon %in% exons) %>% filter(VariantClass != "nonsynonymous_SNV")  %>% ungroup()
ASXL1.disruptive <- mutations %>% filter(Gene == "ASXL1") %>% filter(VariantClass != "nonsynonymous_SNV") %>% ungroup()
TET2.disruptive <- mutations %>% filter(Gene == "TET2") %>% filter(VariantClass != "nonsynonymous_SNV") %>% ungroup()
PPM1D.disruptive <- mutations %>% filter(Gene == "PPM1D") %>% filter(VariantClass != "nonsynonymous_SNV") %>% ungroup()
JAK2.disruptive <- mutations %>% filter(Gene == "JAK2") %>% filter(AA == "p.V617F")%>% ungroup()
TP53.disruptive <- mutations %>% filter(Gene == "TP53") %>% filter(VariantClass != "nonsynonymous_SNV")%>% ungroup()
CBL.mutations <-  mutations %>% filter(Gene == "CBL") %>% filter(AA == "p.C404S" | AA == "p.C416S" | AA == "p.E366K" | AA == "p.C384Y")%>% ungroup()
RAD21.disruptive <- mutations %>% filter(Gene == "RAD21") %>% filter(VariantClass != "nonsynonymous_SNV") %>% ungroup()
STAG2.disruptive <- mutations %>% filter(Gene == "STAG2") %>% filter(VariantClass != "nonsynonymous_SNV") %>% ungroup()
CALR.mutations <- mutations %>% filter(Gene == "CALR") %>% filter(Exon == "exon9")
SETD2.distruptive <- mutations %>% filter(Gene == "SETD2")  %>% filter(AA == "p.R1625C")
ATM.disruptive <- mutations %>% filter(Gene == "ATM") %>% filter(VariantClass != "nonsynonymous_SNV") %>% ungroup()
NF1.disruptive <- mutations %>% filter(Gene == "NF1") %>% filter(VariantClass != "nonsynonymous_SNV") %>% ungroup()
MPL.distruptive <- mutations %>% filter(Gene == "MPL")  %>% filter(AA == "p.W515S")

CH-PD <- rbind(Cosmic.hotspot, 
            ASXL1.disruptive, 
            TET2.disruptive, 
            PPM1D.disruptive, 
            JAK2.disruptive, 
            TP53.disruptive, 
            high_freq_mutations, 
            DNMT3A.disruptive, 
            CBL.mutations, 
            RAD21.disruptive, 
            STAG2.disruptive, 
            CALR.mutations, 
            SETD2.distruptive, 
            ATM.disruptive, 
            NF1.disruptive, 
            MPL.distruptive)

CH-PD <- CH-PD %>% filter(PrcntOccuranceInNormals < 0.2) 

## Compare to mutations identified in Papaemmanuil ey al. 2016  NEJM study
elli <- read.delim("Elli_mutations.txt")
elli <- elli %>% filter(Result == "ONCOGENIC")

elli$Key <- with(elli, paste(CHR, POSITION, WT, MT, sep=":"))

additional.mutations <- mutations %>% 
  filter(Key %in% elli$Key) %>% 
  filter(!Key %in% CH-PD$Key) %>% 
  mutate(Source = "Elli-only")

CH-PD <- rbind(CH-PD, additional.mutations)
write.df(PD, "MSK-IMPACT_CH_Mutation_List_Potential_Drivers.txt")

##########################
### CH ALL GENES
####

first_pass_mutations <- rbind(all_mutations %>%
                                filter(Gene %in% genes.leukemia) %>%
                                filter(VF_N >= 0.02) %>%
                                filter(VF_N < 0.35) %>%
                                filter(VF_T < 0.35) %>%
                                filter(DP_N > 20) %>%
                                filter(DP_T > 20) %>%
                                filter(AD_N >=8) %>%
                                filter(VF_N/VF_T > 2),
                              all_mutations %>%
                                filter(!Gene %in% genes.leukemia) %>%
                                filter(VF_N >= 0.05) %>%
                                filter(VF_N < 0.35) %>%
                                filter(VF_T < 0.35)%>%
                                filter(DP_N > 20) %>%
                                filter(DP_T > 20) %>%
                                filter(AD_N >=8) %>%
                                filter(VF_N/VF_T > 2))


second_pass_mutations <- rbind(all_mutations  %>%
                                 filter(Gene %in% genes.leukemia) %>%
                                 filter(VF_N >= 0.01) %>%
                                 filter(VF_N < 0.02) %>%
                                 filter(VF_T < 0.35) %>%
                                 filter(DP_N > 20) %>%
                                 filter(DP_T > 20)%>% 
                                 filter(AD_N >=8) %>%
                                 filter(VF_N/VF_T > 2) %>%
                                 filter(Normal %in% first_pass_mutations$Normal),
                               all_mutations  %>%
                                 filter(!Gene %in% genes.leukemia) %>%
                                 filter(VF_N >= 0.03) %>%
                                 filter(VF_N < 0.05) %>%
                                 filter(VF_T < 0.35) %>%
                                 filter(DP_N > 20) %>%
                                 filter(DP_T > 20)%>%
                                 filter(AD_N >=8) %>%
                                 filter(VF_N/VF_T > 2) %>%
                                 filter(Normal %in% first_pass_mutations$Normal))


CH <- rbind(first_pass_mutations, second_pass_mutations, high_freq_mutations)


# remove variants based on the occurence counts in the cohort and in COSMIC databse
n = 5
CH.count5 <- rbind(CH %>% filter(OccuranceCountInCohort < n), 
                   CH %>% filter(OccuranceCountInCohort >= n & OccuranceCountInCohort < 100 ) %>% filter(CosmicCount > 1 ),
                   CH %>% filter(OccuranceCountInCohort >= 100 ) %>% filter(CosmicCount > 15 )) 
CH.count5 <- CH.count5 %>% 
  filter(PrcntOccuranceInNormals < 0.2)

write.df(CH.count5, "MSK-IMPACT_CH_Mutation_List_All_Genes.txt")
