library(vegan)
library(GUniFrac)
library(plyr)

#####16S loading and calculations#####
taxa.16s.all <- t(read.table("total_4m_16s_tax.txt", sep = "\t", fill = TRUE, header = TRUE, row.names = 1))
m.all <- read.table("microbiome_data_Probiotics.txt", sep = "\t", header = TRUE, fill = TRUE, row.names = 1)
otu.16s.all <- t(read.table("total_4m_16s_ref.txt", sep = "\t", fill = TRUE, header = TRUE, row.names = 1))

# cleans up the map
m.all[m.all == "UNKNOWN" ] <- NA
m.all$deliverymode <- factor(m.all$deliverymode, levels = c("CS", "VD"))
m.all$Sample_collection_time_M <- factor(m.all$Sample_collection_time_M, levels = c("1","6"))
m.all$EBF_6mo <- factor(m.all$EBF_6mo, levels = c("0","1"))

m.all$pre_probiotics_levels <- factor(m.all$pre_probiotics_levels, levels = c("None", "Low", "High"))
m.all$pre_probiotics_no_high <- factor(m.all$pre_probiotics_no_high, levels = c("None", "High"))
m.all$pre_probiotics_low_high <- factor(m.all$pre_probiotics_low_high, levels = c("Low", "High"))
m.all$pre_probiotics_no_low <- factor(m.all$pre_probiotics_no_low, levels = c("None", "Low"))

m.all$m6_mother_probiotic <- factor(m.all$m6_mother_probiotic, levels = c("No", "Yes"))
m.all$m1_mother_probiotic <- factor(m.all$m1_mother_probiotic, levels = c("No", "Yes"))

m.all$pre_probiotics <- factor(m.all$pre_probiotics, levels = c("No", "Yes"))

m.all$m1_mother_probiotic_levels  <- factor(m.all$m1_mother_probiotic_levels, levels = c("None", "Low", "High"))
m.all$m1_mother_probiotic_no_high <- factor(m.all$m1_mother_probiotic_levels, levels = c("None", "High"))
m.all$m1_mother_probiotic_low_high <- factor(m.all$m1_mother_probiotic_levels, levels = c("Low", "High"))
m.all$m1_mother_probiotic_no_low <- factor(m.all$m1_mother_probiotic_levels, levels = c("None", "Low"))

m.all$New_1mo_Probiotic_Group  <- factor(m.all$New_1mo_Probiotic_Group, levels = c("None", "Low", "High"))
m.all$New_1mo_Probiotic_None_Low <- factor(m.all$New_1mo_Probiotic_Group, levels = c("None", "Low"))
m.all$New_1mo_Probiotic_High_Low <- factor(m.all$New_1mo_Probiotic_Group, levels = c("High", "Low"))
m.all$Combined_1mo_Probiotic_NoLow_High <- ifelse((m.all$New_1mo_Probiotic_Group==c("Low")|(m.all$New_1mo_Probiotic_Group==c("None"))), yes="< 1 month", ifelse(m.all$New_1mo_Probiotic_Group==c("High"), yes="> 1 month", no=NA))
m.all$Combined_1mo_Probiotic_NoLow_High <- factor(m.all$Combined_1mo_Probiotic_NoLow_High, levels = c("< 1 month", "> 1 month"))

m.all$Combined_1mo_Probiotic_High <- factor(m.all$Combined_1mo_Probiotic_NoLow_High, levels = c("> 1 month"))


m.all$New_6mo_Probiotic_Group  <- factor(m.all$New_6mo_Probiotic_Group, levels = c("None", "Low", "High"))
m.all$New_6mo_Probiotic_None_Low <- factor(m.all$New_6mo_Probiotic_Group, levels = c("None", "Low"))
m.all$New_6mo_Probiotic_High_Low <- factor(m.all$New_6mo_Probiotic_Group, levels = c("High", "Low"))
m.all$Combined_6mo_Probiotic_NoLow_High <- ifelse((m.all$New_6mo_Probiotic_Group==c("Low")|(m.all$New_6mo_Probiotic_Group==c("None"))), yes="< 1 month", ifelse(m.all$New_6mo_Probiotic_Group==c("High"), yes="> 1 month", no=NA))
m.all$Combined_6mo_Probiotic_NoLow_High <- factor(m.all$Combined_6mo_Probiotic_NoLow_High, levels = c("< 1 month", "> 1 month"))

m.all$Combined_6mo_Probiotic_High <- factor(m.all$Combined_6mo_Probiotic_NoLow_High, levels = c("> 1 month"))

m.all$ANY_PROBIOTICS_1to6MO <- ifelse((m.all$Combined_1mo_Probiotic_NoLow_High==c("> 1 month")|(m.all$Combined_6mo_Probiotic_NoLow_High==c("> 1 month"))), yes = "> 1 month", no = "< 1 month")

m.all$PROBIOTICS_1and6MO <- ifelse((m.all$Combined_1mo_Probiotic_NoLow_High==c("> 1 month")&(m.all$Combined_6mo_Probiotic_NoLow_High==c("> 1 month"))), yes = "> 1 month", no = "< 1 month")

  

m.all$New_1mo_Mat_Antibiotic_Group <- factor(m.all$New_1mo_Mat_Antibiotic_Group, levels = c("No", "Yes"))
m.all$New_6mo_Mat_Antibiotic_Group <- factor(m.all$New_6mo_Mat_Antibiotic_Group, levels = c("No", "Yes"))
m.all$New_1mo_Infant_Antibiotic_Group <- factor(m.all$New_1mo_Infant_Antibiotic_Group, levels = c("No", "Yes"))
m.all$New_6mo_Infant_Antibiotic_Group <- factor(m.all$New_6mo_Infant_Antibiotic_Group, levels = c("No", "Yes"))
m.all$Any_Probiotics_6mo <- factor(m.all$Any_Probiotics_6mo, levels = c("None", "Yes"))
m.all$New_1mo_Probiotic_Group_High_None <- factor(m.all$New_1mo_Probiotic_Group, levels = c("None", "High"))
m.all$New_6mo_Probiotic_Group_High_None <- factor(m.all$New_6mo_Probiotic_Group, levels = c("None", "High"))
m.all$New_1mo_Infant_Probiotic_Group <- factor(m.all$m1_infant_supp___3, levels = c("Unchecked", "Checked"))
m.all$New_6mo_Infant_Probiotic_Group <- factor(m.all$m6_infant_supp___3, levels = c("Unchecked", "Checked"))

m.all$BMI_High_Normal <- ifelse(m.all$mat_bmi_cat==c(">= 30")|m.all$mat_bmi_cat==c("25.00 - 29.99"), yes="High", no="Healthy")

m.all$total_m1_pro_duration <- factor(m.all$total_m1_pro_duration, levels = c("None", "less than 1 week", "1-3 weeks","1-2 months","2-4 months", "4-6 months", "more than 6 months"))
m.all$m6_mother_probio_duration <- factor(m.all$m6_mother_probio_duration, levels = c("None", "less than 1 week", "1-3 weeks","1-2 months","2-4 months", "4-6 months", "more than 6 months"))

### New probiotics categories###

m.all$pre_probiotics_Combined <- ifelse((m.all$pre_probiotics_levels==c("Low")|(m.all$pre_probiotics_levels==c("None"))), yes="< 1 month", ifelse(m.all$pre_probiotics_levels==c("High"), yes="> 1 month", no=NA))
m.all$Probiotics_Birth_to_1month_Combined <- ifelse((m.all$m1_mother_probiotic_levels==c("Low")|(m.all$m1_mother_probiotic_levels==c("None"))), yes="< 1 month", ifelse(m.all$m1_mother_probiotic_levels==c("High"), yes="> 1 month", no=NA))

m.all$Probiotics_Birth_to_6month_Combined <- ifelse((m.all$m6_mother_probio_duration==c("1-2 months")|
                                                       (m.all$m6_mother_probio_duration==c("2-4 months"))|
                                                       (m.all$m6_mother_probio_duration==c("4-6 months"))|
                                                       (m.all$m6_mother_probio_duration==c("more than 6 months"))|
                                                       (m.all$Probiotics_Birth_to_1month_Combined==c("> 1 month"))), yes="> 1 month", no="< 1 month")
m.all$Any_Probiotics_prenatal_to6month_combined <- ifelse((m.all$Probiotics_Birth_to_6month_Combined==c("> 1 month")|
                                                             (m.all$pre_probiotics_Combined==c("> 1 month"))|
                                                             (m.all$Probiotics_Birth_to_1month_Combined==c("> 1 month"))), yes="> 1 month", no="< 1 month")

m.all$Any_mat_abx_birthto6 <- ifelse((m.all$m1_mother_rx==c("Yes")|
                                                             (m.all$m6_mother_rx==c("Yes"))), yes="Yes", no="No")

m.all$Any_infant_abx_birthto6 <- ifelse((m.all$m1_infant_infect_rx==c("Yes")|
                                        (m.all$m6_infant_infect_rx==c("Yes"))), yes="Yes", no="No")


##getting overall n for probiotics yes/no###
m.all$Yesmatpro <- ifelse((m.all$pre_probiotics_Combined==c("> 1 month")|
                                                       (m.all$Probiotics_Birth_to_1month_Combined==c("> 1 month"))|
                                                       (m.all$Probiotics_Birth_to_6month_Combined==c("> 1 month"))|
                                                      (m.all$Combined_1mo_Probiotic_NoLow_High==c("> 1 month"))|
                                                       (m.all$Combined_6mo_Probiotic_NoLow_High==c("> 1 month"))|
                                                       (m.all$ANY_PROBIOTICS_1to6MO==c("> 1 month"))), yes="Yes", no="No")


##adding metadata from the additional metadata file##
m.add <- read.csv("microbiome_data_ADD.txt", sep = "\t", header = TRUE, fill = TRUE, row.names = 1)
m.add <- m.add[rownames(m.all),]
m.all$bm_il_6_1 <- m.add$bm_il_6_1
m.all$log_bm_il_6_1 <- log(m.add$bm_il_6_1)

##adding more metadata from the additional metadata file##
m.add.more <- read.csv("4M_microbiome_data_MASTER_SAMPLE_ID.txt", sep = "\t", header = TRUE, fill = TRUE, row.names = 1)
m.add.more <- m.add.more[rownames(m.all),]
m.all$GBS_Positive <- m.add.more$GBS_Positive
m.all$GBS_Positive <- factor(m.all$GBS_Positive, levels = c("Yes", "No"))
m.all$GBS_CS <- ifelse((m.all$GBS_Positive==c("Yes")|(m.all$deliverymode==c("CS"))), yes = "Yes", no = "No")

m.all$m1_mother_rx <- factor(m.all$m1_mother_rx, levels = c("Yes", "No"))
m.all$m1_mother_rx_GBS <- ifelse((m.all$GBS_Positive==c("Yes")|(m.all$m1_mother_rx==c("Yes"))), yes = "Yes", no = "No")
m.all$m1_mother_rx_GBS_CS <- ifelse((m.all$m1_mother_rx_GBS==c("Yes")|(m.all$deliverymode==c("CS"))), yes = "Yes", no = "No")


##adding ERP###
setwd("~/Desktop/R Projects/ERP") #Set to directory where results are desired

m.erp.add <- read.csv("ERP_NSW.txt", sep = "\t", header = TRUE, fill = TRUE, row.names = 1)
m.erp.add <- m.erp.add[rownames(m.all),]
m.all$LC_SW_Familiar_6mo <- m.erp.add$LC_SW_Familiar_6mo
m.all$LF_SW_Familiar_6mo <- m.erp.add$LF_SW_Familiar_6mo
m.all$RC_SW_Familiar_6mo <- m.erp.add$RC_SW_Familiar_6mo
m.all$RF_SW_Familiar_6mo <- m.erp.add$RF_SW_Familiar_6mo

m.all$LC_SW_Novel_6mo <- m.erp.add$LC_SW_Novel_6mo
m.all$LF_SW_Novel_6mo <- m.erp.add$LF_SW_Novel_6mo
m.all$RC_SW_Novel_6mo <- m.erp.add$RC_SW_Novel_6mo
m.all$RF_SW_Novel_6mo <- m.erp.add$RF_SW_Novel_6mo

m.all$LC_SW_Fam_Novel_6mo <- m.erp.add$LC_SW_Fam_Novel_6mo
m.all$LF_SW_Fam_Novel_6mo <- m.erp.add$LF_SW_Fam_Novel_6mo
m.all$RC_SW_Fam_Novel_6mo <- m.erp.add$RC_SW_Fam_Novel_6mo
m.all$RF_SW_Fam_Novel_6mo <- m.erp.add$RF_SW_Fam_Novel_6mo


##do this for infants too##
setwd("~/Desktop/R Projects/4M_Probiotics") #Set to directory where results are desired

m.all$m1_infant_infect_rx <- factor(m.all$m1_infant_infect_rx, levels = c("Yes", "No"))
m.all$m1_infant_infect_rx_GBS <- ifelse((m.all$GBS_Positive==c("Yes")|(m.all$m1_infant_infect_rx==c("Yes"))), yes = "Yes", no = "No")
m.all$m1_infant_infect_rx_GBS_CS <- ifelse((m.all$m1_infant_infect_rx_GBS==c("Yes")|(m.all$deliverymode==c("CS"))), yes = "Yes", no = "No")


taxa.16s.all <- as.data.frame(taxa.16s.all)
taxa.16s.table <- taxa.16s.all
taxa.16s.table[is.na(taxa.16s.table)] <- 0
m.taxa.16s <- m.all[(rownames(taxa.16s.table)),] # map for the taxa table


####Taxa table QC####
taxa.16s.norm <- sweep(taxa.16s.table, 1, rowSums(taxa.16s.table), '/') * 100
taxa.16s.filt <- taxa.16s.table[rowSums(taxa.16s.table) > 1000,]
taxa.16s.filt <- taxa.16s.filt[,colSums(taxa.16s.filt) > 0]
m.filt <- m.taxa.16s[rownames(taxa.16s.filt),]

m.filt.infant <- m.filt[rownames(m.filt[m.filt$Sample_type == c("Infant_fecal"),]),]
m.filt.milk <- m.filt[rownames(m.filt[m.filt$Sample_type == c("breastmilk"),]),]

####Compress taxa table to genus/lowest taxonomy level####
taxa.16s.filt.t <- t(taxa.16s.filt)
splittaxa.16s <- strsplit(rownames(taxa.16s.filt.t), ";")
for(i in 1:length(splittaxa.16s)){
  if(length(splittaxa.16s[[i]]) > 5){
    rownames(taxa.16s.filt.t)[i] <- splittaxa.16s[[i]][6] #Genus-level taxa
  } else {
    rownames(taxa.16s.filt.t)[i] <- paste(splittaxa.16s[[i]][(length(splittaxa.16s[[i]]))],collapse = "_") #Deepest taxa if not genus
  }
}
taxa.16s.filt.t <- aggregate(taxa.16s.filt.t, by=list(rownames(taxa.16s.filt.t)),sum)
rownames(taxa.16s.filt.t) <- taxa.16s.filt.t$Group.1
taxa.16s.filt.t <- taxa.16s.filt.t[,2:length(colnames(taxa.16s.filt.t))]
taxa.16s.filt.genus <- as.data.frame(t(taxa.16s.filt.t))

####Taxonomy table subsets####
taxa.16s.filt.norm.genus <- taxa.16s.filt.genus
taxa.16s.infants <- taxa.16s.filt.norm.genus[rownames(m.filt.infant),]
taxa.16s.1month <- taxa.16s.filt.norm.genus[rownames(m.filt.infant[m.filt.infant$Sample_collection_time_M == 1,]),]
m.16s.1month <- m.filt.infant[rownames(taxa.16s.1month),]
taxa.16s.6month <- taxa.16s.filt.norm.genus[rownames(m.filt.infant[m.filt.infant$Sample_collection_time_M == 6,]),]
m.16s.6month <- m.filt.infant[rownames(taxa.16s.6month),]
taxa.16s.milk <- taxa.16s.filt.norm.genus[rownames(m.filt.milk),]
m.16s.1monthandmilk <- m.filt[m.filt$Sample_collection_time_M == 1,]
taxa.16s.1monthandmilk <- taxa.16s.filt.norm.genus[rownames(m.16s.1monthandmilk),]
taxa.16s.6month.bf <- taxa.16s.6month[rownames(m.16s.6month[m.16s.6month$EBF_6mo == 1,]),]
m.16s.6month.bf <- m.16s.6month[rownames(taxa.16s.6month.bf),]
taxa.16s.1month.6EBF <- taxa.16s.1month[rownames(m.16s.1month[m.16s.6month$EBF_6mo == 1,]),]
m.16s.1month.6EBF <- m.16s.1month[rownames(taxa.16s.1month.6EBF),]

m.16s.1month.CS <- m.16s.1month[rownames(m.16s.1month[m.16s.1month$deliverymode == c("CS"),]),]
taxa.16s.1month.CS <- taxa.16s.1month[rownames(m.16s.1month.CS),]

taxa.16s.1month.VD <- taxa.16s.1month[rownames(m.16s.1month[m.16s.1month$deliverymode == c("VD"),]),]
m.16s.1month.VD <- m.16s.1month[rownames(taxa.16s.1month.VD),]


####CLR transformations####
taxa.16s.genus.clr <- t(taxa.16s.filt.genus)
eps <- 0.5
taxa.16s.genus.clr <- taxa.16s.genus.clr*(1-rowSums(taxa.16s.genus.clr==0)*eps/rowSums(taxa.16s.genus.clr))
taxa.16s.genus.clr[taxa.16s.genus.clr==0] <- eps
taxa.16s.genus.clr <- sweep(taxa.16s.genus.clr,1,rowSums(taxa.16s.genus.clr),'/')
ltaxa.16s <- log(taxa.16s.genus.clr)
taxa.16s.genus.clr <- t(ltaxa.16s - rowMeans(ltaxa.16s))
taxa.16s.genus.clr <- taxa.16s.genus.clr[,!is.nan(colSums(taxa.16s.genus.clr))]
m.clr <- m.all[rownames(taxa.16s.genus.clr),]
taxa.16s.1month.clr <- taxa.16s.genus.clr[rownames(taxa.16s.1month),]
taxa.16s.6month.clr <- taxa.16s.genus.clr[rownames(taxa.16s.6month),]
taxa.16s.milk.clr <- taxa.16s.genus.clr[rownames(taxa.16s.milk),]
taxa.16s.infants.clr <- taxa.16s.genus.clr[rownames(taxa.16s.infants),]
taxa.16s.1monthandmilk.clr <- taxa.16s.genus.clr[rownames(taxa.16s.1monthandmilk),]
taxa.16s.6month.bf.clr <- taxa.16s.genus.clr[rownames(taxa.16s.6month.bf),]
taxa.16s.1month.6EBF.clr <- taxa.16s.genus.clr[rownames(taxa.16s.1month.6EBF),]
