##' A script for comparing differential expression between HER2+ and HER2- breast cancers and comparing the results to those observed when comparing HER2+ and HER2- samples in GDSC.

#' Load libraries
library("GenomicRanges")
library("ggplot2")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("org.Hs.eg.db")
library("GenomicFeatures")

#' Set the root directory.
theRootDir <- "/media/user/My Passport/UofC_backup/bionimbus/data_scratch/finalData/"

#' Load the TCGA breast cancer expression data.
brcaDataLoc <- paste(theRootDir, "dataIn/rnaSeq/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2015082100.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", sep="")
tpmDatMat_bc <- read.delim(brcaDataLoc, as.is=T)
tpmDatMat_bc_tpm <- apply(tpmDatMat_bc[-1,which(tpmDatMat_bc[1,] == "scaled_estimate")], 2, as.numeric)
tpmDatMat_bc_tpm <- tpmDatMat_bc[-1,which(tpmDatMat_bc[1,] == "scaled_estimate")]
tpmDatMat_bc_tpm <- apply(tpmDatMat_bc_tpm, 2, as.numeric)
geneNames <- do.call(cbind, strsplit(tpmDatMat_bc[, "Hybridization.REF"], "|", fixed=TRUE))[1,][-1]
rownames(tpmDatMat_bc_tpm) <- geneNames
colnames(tpmDatMat_bc_tpm) <- substr(colnames(tpmDatMat_bc_tpm), 1, 28)
tpmDatMat_bc_tpm_logged <- log((tpmDatMat_bc_tpm*1000000)+1)
tpmDatMat_bc_tpm_logged_tumor <- tpmDatMat_bc_tpm_logged[, do.call(rbind, strsplit(colnames(tpmDatMat_bc_tpm_logged), ".", fixed=T))[,4] == "01A"]
colnames(tpmDatMat_bc_tpm_logged_tumor) <- substring(colnames(tpmDatMat_bc_tpm_logged_tumor), 1, 12)
colnames(tpmDatMat_bc_tpm_logged_tumor) <- gsub(".", "-", colnames(tpmDatMat_bc_tpm_logged_tumor),  fixed=T)

#' Load the matched clinical data for BRCA. This file contains the HER2 status, as measured by immunohistochemistry.
clinicalDataLocation <- paste(theRootDir, "dataIn/clinical/nationwidechildrens.org_clinical_patient_brca.txt", sep="")
clinDataBrca <- read.delim(clinicalDataLocation, as.is=T)
her2status <- clinDataBrca[, "her2_status_by_ihc"]
names(her2status) <- clinDataBrca[, "bcr_patient_barcode"]

#' Load the tumor purity estimates.
nComsProps <- read.csv("/media/user/My Passport/UofC_backup/bionimbus/data_scratch/prediXcanProj/ncomms9971-s2.csv", as.is=T)
rownames(nComsProps) <- nComsProps[,1]
nComsProps_brca <- nComsProps[nComsProps[,2] == "BRCA" & substring(nComsProps[,1], 14, 16) == "01A", "CPE"] # exctract 01A (primary site) breast cancer samples.
names(nComsProps_brca) <- nComsProps[nComsProps[,2] == "BRCA" & substring(nComsProps[,1], 14, 16) == "01A", 1]
nComsProps_brca_noNa <- nComsProps_brca[!is.na(nComsProps_brca)]
names(nComsProps_brca_noNa) <- substring(names(nComsProps_brca_noNa), 1, 12)
samplesInAllThree <- intersect(intersect(colnames(tpmDatMat_bc_tpm_logged_tumor), names(her2status)), names(nComsProps_brca_noNa))

#' Test association of gene expression and purity as a sanity check (we expect lots of associations - this looks as expected).
options(warn=-1)
pTumPurExpr <- numeric()
for(i in 1:nrow(tpmDatMat_bc_tpm_logged_tumor))
{
    pTumPurExpr[i] <- cor.test(tpmDatMat_bc_tpm_logged_tumor[i,samplesInAllThree], nComsProps_brca_noNa[samplesInAllThree])$p.value
}
hist(pTumPurExpr)

#' Define her2 + to her2- TCGA patient samples.
her2Neg <- intersect(names(her2status[her2status == "Negative"]), samplesInAllThree)
her2Pos <- intersect(names(her2status[her2status == "Positive"]), samplesInAllThree)

#' Perform a simple differential expression analysis 
exprData <- tpmDatMat_bc_tpm_logged_tumor[, c(her2Neg, her2Pos)]
caseControlStatus <- c(rep(0, length(her2Neg)), rep(1, length(her2Pos)))
propsFin <- 1- nComsProps_brca_noNa[c(her2Neg, her2Pos)]
pValsNaiveTcga <- numeric()
pValsIntTcga <- numeric()
betaInt <- numeric()
betaNaive <- numeric()
for(i in 1:nrow(exprData))
{
    pValsNaiveTcga[i] <- coef(summary(lm(exprData[i,]~caseControlStatus)))[2,4] # Bulk tumor analysis.
    pValsIntTcga[i] <- coef(summary(lm(exprData[i,]~caseControlStatus*propsFin)))[2,4] # Cancer specific analysis.
    betaInt[i]  <- coef(summary(lm(exprData[i,]~caseControlStatus*propsFin)))[2,1]
    betaNaive[i] <- coef(summary(lm(exprData[i,]~caseControlStatus)))[2,1]
}
names(pValsIntTcga) <- rownames(exprData)
names(pValsNaiveTcga) <- rownames(exprData)
names(betaInt) <- rownames(exprData)
names(betaNaive) <- rownames(exprData)
fdrNaiveTcga <- p.adjust(pValsNaiveTcga, method="BH")
fdrIntTcga <- p.adjust(pValsIntTcga, method="BH")
sigInt <- names(sort(fdrIntTcga[which(fdrIntTcga < 0.05)]))
sigNaive <- names(sort(fdrNaiveTcga[which(fdrNaiveTcga < 0.05)]))

#' Now investigate the association between ERBB2 amplification and lapatinib sensitivty in the GDSC cell lines.
#' Load the CNV data from GDSC This was acquired from cancerrxgene.org. Generated using Affymetrix SNP 6.0 data. This data was not mapped to genes by CGP, so we do that below. 
allCnvs <- read.csv(paste(theRootDir, "dataIn/cell_lines_copy_number.csv", sep=""), as.is=T, header=T)
cellLines_cnv_list <- split(allCnvs, allCnvs[,1])

#' For each cell line in CGP, create a GenomicRanges object, which contains the locations and magnitutes of all Copy Number measurements. This data is GRCh37/hg19.
grCnvsList <- list()
for(i in 1:length(cellLines_cnv_list))
{
  grCnvsList[[i]] <- GRanges(seqnames=Rle(paste("chr", cellLines_cnv_list[[i]][, "chr_37"], sep="")), ranges=IRanges(cellLines_cnv_list[[i]][, "startpos_37"], cellLines_cnv_list[[i]][, "endpos_37"]), segMeans=cellLines_cnv_list[[i]][, "totalCN"])
}
names(grCnvsList) <- names(cellLines_cnv_list)

#' Load the gene ranges for HG19 using.
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
geneRanges <- genes(txdb)
e2s = toTable(org.Hs.egSYMBOL)
syms <- e2s[, "symbol"]
names(syms) <- e2s[, "gene_id"]
theGeneSymsOrd <- syms[as.character(geneRanges$gene_id)]

#' We will now intersect the gene ranges with the CNV data in order to establish the copy number for each gene.
numGenesQuantifid <- numeric()
theCnvQuantVecList <- list()
for(i in 1:length(grCnvsList))
{
    grCnvs <- grCnvsList[[i]]

    # Use count overlaps to find genes that unambiguously overlap a single peak. Give it an NA it it doesn't overlap a single peak. Assign it the value of the peak if it unambiguously overlaps a peak. PC.
    numOverlaps <- countOverlaps(geneRanges, grCnvs)
    numGenesQuantifid[i] <- sum(numOverlaps == 1)
    inCnv <- which(numOverlaps == 1) # take only gene unambiguously overlaping a peak, this is usually most genes.
    
    theCnvQuantVec <- rep(NA, length(geneRanges))
    olaps <- findOverlaps(geneRanges, grCnvs, type="within", ignore.strand=TRUE)
    theCnvQuantVec[queryHits(olaps)] <- grCnvs$segMeans[subjectHits(olaps)]
    theCnvQuantVecList[[i]] <- theCnvQuantVec
    names(theCnvQuantVecList[[i]]) <- theGeneSymsOrd
}
names(theCnvQuantVecList) <- names(grCnvsList)
theCnvQuantVecList_mat <- do.call(cbind, theCnvQuantVecList)
erbb2CnvVec <- theCnvQuantVecList_mat["ERBB2",]

#' We will now do some analysis testing the gene CVNs against lapatinib sensitivty here in GDSC
#' Load the CGP expression and phenotype data.
load("/media/user/My Passport/UofC_backup/bionimbus/data_scratch/finalData/dataIn/drugAndPhenoCgp.RData") # "drugSensitivityDataCgp" "drugToCellLineDataCgp" "gdsc_brainarray_syms
load("/media/user/My Passport/UofC_backup/bionimbus/data_scratch/finalData/dataIn/cgp2016ExprRma.RData") # "cgp2016ExprRma" 
breastCancerCellLines <- drugSensitivityDataCgp[, "Cell.Line"][drugSensitivityDataCgp[, "Tissue"] == "breast"]

#' Intersect the expression and CNV data.
isBrca_hasExpr_hasCnv <- intersect(intersect(breastCancerCellLines, colnames(cgp2016ExprRma)), names(erbb2CnvVec))
brcaErbb2CnvOrd <- erbb2CnvVec[isBrca_hasExpr_hasCnv]
brcaExprOrd <- cgp2016ExprRma[,isBrca_hasExpr_hasCnv]

#' Identify DE genes associated with ERBB2 amplification in GDSC. We will treat nominally significant genes with concordant directionality as a gold standard.
pVals <- numeric()
betaCellLines <- numeric()
for(i in 1:nrow(brcaExprOrd))
{
    pVals[i] <- coef(summary(lm(brcaExprOrd[i,]~brcaErbb2CnvOrd)))[2,4]
    betaCellLines[i] <- coef(summary(lm(brcaExprOrd[i,]~brcaErbb2CnvOrd)))[2,1]
}
names(pVals) <- rownames(brcaExprOrd)
names(betaCellLines) <- rownames(brcaExprOrd)
nomSigCellLines <- names(pVals)[pVals < 0.05]
pAdj <- p.adjust(pVals, method="BH")
sum(pAdj < 0.05)
pAdj[pAdj < 0.05]
theDeGenes <- names(sort(pAdj[pAdj < 0.05]))
genesBothDatasets <- intersect(rownames(exprData), rownames(cgp2016ExprRma))

#' Subset to genes in both TCGA and GDSC datasets and compare the overlaps.
theDeGenes_both <- theDeGenes[theDeGenes %in% genesBothDatasets]
sigInt_both <- sigInt[sigInt %in% genesBothDatasets]
sigNaive_both <- sigNaive[sigNaive %in% genesBothDatasets]
nomSigCellLines_both <- nomSigCellLines[nomSigCellLines %in% genesBothDatasets]

#' How many differentially expressed genes overlap between the two datasets (compare FDR < 0.05 from both datasets)
sum(theDeGenes_both %in% sigInt_both) # cancer specific model
sum(theDeGenes_both %in% sigNaive_both) # naive bulk model

#' What proportion overlap (overlap is 5 times higher for cancer-specific model)
sum(theDeGenes_both %in% sigInt_both) / length(sigInt_both)
sum(theDeGenes_both %in% sigNaive_both) / length(sigNaive_both)

#' Do the above analysis, treat nominally significant genes as the gold standard.
sum(nomSigCellLines_both %in% sigInt_both) # cancer specific model
sum(nomSigCellLines_both %in% sigNaive_both) # naive bulk model

#' What proportions overlap? Has improved from 18% to 26%
sum(nomSigCellLines_both %in% sigInt_both) / length(sigInt_both)
sum(nomSigCellLines_both %in% sigNaive_both) / length(sigNaive_both)


#' calculate fishers's exact test
numOlapInt <- sum(sigInt_both %in% nomSigCellLines_both)
numNoOlapInt <- sum(!sigInt_both %in% nomSigCellLines_both)
numOlapNaive <- sum(sigNaive_both %in% nomSigCellLines_both)
numNoOlapNaive <- sum(!sigNaive_both %in% nomSigCellLines_both)
print(fisher.test(t(matrix(c(numOlapInt, numNoOlapInt, numOlapNaive, numNoOlapNaive), nrow=2))))

#' Calculate the proportion of times directionality is concordant, is that improved? (ans: yes)
concDirInt <- betaInt[sigInt_both] * betaCellLines[sigInt_both] # how often is directionality concordant 
concDirNaive <- betaNaive[sigNaive_both] * betaCellLines[sigNaive_both]
sum(concDirInt > 0) / length(concDirInt)
sum(concDirNaive > 0) / length(concDirNaive)













