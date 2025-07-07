#!/usr/bin/env python
# coding: utf-8


import sys
import warnings
warnings.filterwarnings("ignore")
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
#sc.settings.logfile = sys.stdout
## fix palantir breaking down some plots
import seaborn 
seaborn.reset_orig()

sc.set_figure_params(figsize=(4, 4))


import random 
import scanpy as sc
import os
import numpy as np
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
adata_ref = sc.read(os.path.join(directory,'ref100k.h5ad'), compression='gzip')
adata_ref=adata_ref[np.random.choice(adata_ref.obs.index,5000,replace=False),:]
adata_ref.X = adata_ref.layers['counts']

#load query data and normalize
adata_query = sc.read(os.path.join(directory,'query_HM2Dall_guides.h5ad'), compression='gzip')
adata_query = adata_query[adata_query.obs['gene_NKS'] == 'non-targeting']
adata_query = adata_query[adata_query.obs['batch'] == 'HD7']
#adata_query=adata_query[np.random.choice(adata_query.obs.index,5000,replace=False),:]


import scanpy as sc
adata3 = sc.read('/wynton/group/pollen/jding/brainchromatin/iPSC/He2024/hnoca_cleanedmeta.50k.h5ad')
adata3=adata3[np.random.choice(adata3.obs.index,5000,replace=False),:]


adata4 = sc.read('/wynton/group/pollen/jding/brainchromatin/iPSC/Lin2021/matrices_timecourse/adata.h5ad')
adata4=adata4[np.random.choice(adata4.obs.index,5000,replace=False),:]

adata5 = sc.read('/wynton/group/pollen/jding/brainchromatin/perturb/Cellbender/Slice_NKS.h5ad')
adata5 = adata5.raw.to_adata()
adata5 = adata5[adata5.obs['gene_NKS'] == 'non-targeting']
#adata5=adata5[np.random.choice(adata5.obs.index,5000,replace=False),:]
print(adata5.X)

adata6 = sc.read('/wynton/group/pollen/jding/brainchromatin/iPSC/Hendriks2024/query.h5ad')
adata6 = adata6.raw.to_adata()
adata6=adata6[np.random.choice(adata6.obs.index,5000,replace=False),:]


import anndata
adata = anndata.AnnData.concatenate(adata_ref,adata5,adata_query, adata3,adata4,adata6, join='outer',batch_categories=['In Vivo Atlas - Wang 2025','Primary Slice - This Study','Primary 2D - This Study','Organoids - He 2024','iNeurons - Lin 2021', 'FEBO - Hendriks 2024']) 

sc.pp.normalize_total(adata, target_sum=1e4,exclude_highly_expressed=True)
sc.pp.log1p(adata)
HALLMARK_GLYCOLYSIS = ["ABCB6","ADORA2B","AGL","AGRN","AK3","AK4","AKR1A1","ALDH7A1","ALDH9A1","ALDOA","ALDOB","ALG1","ANG","ANGPTL4","ANKZF1","ARPP19","ARTN","AURKA","B3GALT6","B3GAT1","B3GAT3","B3GNT3","B4GALT1","B4GALT2","B4GALT4","B4GALT7","BIK","BPNT1","CACNA1H","CAPN5","CASP6","CD44","CDK1","CENPA","CHPF","CHPF2","CHST1","CHST12","CHST2","CHST4","CHST6","CITED2","CLDN3","CLDN9","CLN6","COG2","COL5A1","COPB2","CTH","CXCR4","CYB5A","DCN","DDIT4","DEPDC1","DLD","DPYSL4","DSC2","ECD","EFNA3","EGFR","EGLN3","ELF3","ENO1","ENO2","ERO1A","EXT1","EXT2","FAM162A","FBP2","FKBP4","FUT8","G6PD","GAL3ST1","GALE","GALK1","GALK2","GAPDHS","GCLC","GFPT1","GFUS","GLCE","GLRX","GMPPA","GMPPB","GNE","GNPDA1","GOT1","GOT2","GPC1","GPC3","GPC4","GPR87","GUSB","GYS1","GYS2","HAX1","HDLBP","HK2","HMMR","HOMER1","HS2ST1","HS6ST2","HSPA5","IDH1","IDUA","IER3","IGFBP3","IL13RA1","IRS2","ISG20","KDELR3","KIF20A","KIF2A","LCT","LDHA","LDHC","LHPP","LHX9","MDH1","MDH2","ME1","ME2","MED24","MERTK","MET","MIF","MIOX","MPI","MXI1","NANP","NASP","NDST3","NDUFV3","NOL3","NSDHL","NT5E","P4HA1","P4HA2","PAM","PAXIP1","PC","PDK3","PFKFB1","PFKP","PGAM1","PGAM2","PGK1","PGLS","PGM2","PHKA2","PKM","PKP2","PLOD1","PLOD2","PMM2","POLR3K","PPFIA4","PPIA","PPP2CB","PRPS1","PSMC4","PYGB","PYGL","QSOX1","RARS1","RBCK1","RPE","RRAGD","SAP30","SDC1","SDC2","SDC3","SDHC","SLC16A3","SLC25A10","SLC25A13","SLC35A3","SLC37A4","SOD1","SOX9","SPAG4","SRD5A3","STC1","STC2","STMN1","TALDO1","TFF3","TGFA","TGFBI","TKTL1","TPBG","TPI1","TPST1","TXN","UGP2","VCAN","VEGFA","VLDLR","XYLT2","ZNF292"]
sc.tl.score_genes(adata, HALLMARK_GLYCOLYSIS,score_name='Glycolysis')

RESPONSE_TO_ER_STRESS = ["ABCA7","ABL1","AFF4","AGR2","AIFM1","ALOX15","ALOX5","AMFR","ANKS4B","ANKZF1","APAF1","AQP11","ATF3","ATF4","ATF6","ATF6B","ATG10","ATP2A1","ATP2A2","ATP2A3","ATXN3","AUP1","BAG6","BAK1","BAX","BBC3","BCAP31","BCL2","BCL2L1","BCL2L11","BFAR","BHLHA15","BOK","BRSK2","CALR","CALR3","CANX","CASP4","CAV1","CCDC47","CCND1","CDK5RAP3","CEBPB","CERT1","CFTR","CHAC1","CLGN","CLU","COPS5","CREB3","CREB3L1","CREB3L2","CREB3L3","CREBZF","CTH","CXCL8","DAB2IP","DDIT3","DDRGK1","DDX3X","DERL1","DERL2","DERL3","DNAJB12","DNAJB2","DNAJB9","DNAJC10","DNAJC3","ECPAS","EDEM1","EDEM2","EDEM3","EIF2AK2","EIF2AK3","EIF2AK4","EIF2B5","EIF2S1","EIF4G1","ELAVL4","ERLEC1","ERLIN1","ERLIN2","ERMP1","ERN1","ERN2","ERO1A","ERP27","ERP29","ERP44","FAF1","FAF2","FAM8A1","FBXO17","FBXO2","FBXO27","FBXO44","FBXO6","FCGR2B","FGF21","FICD","FLOT1","FOXRED2","GET4","GORASP2","GRINA","GSK3B","HERPUD1","HERPUD2","HM13","HSP90B1","HSPA1A","HSPA5","HYOU1","ITPR1","JKAMP","JUN","KCNJ8","LPCAT3","LRRK2","MAGEA3","MAN1A1","MAN1A2","MAN1B1","MAN1C1","MANF","MAP3K5","MARCHF6","MARCKS","MBTPS1","MBTPS2","MIR199A1","MIR200C","NCCRP1","NCK1","NCK2","NFE2L2","NHLRC1","NIBAN1","NOD1","NOD2","NPLOC4","NR1H2","NR1H3","NRBF2","NUPR1","OPA1","OS9","P4HB","PARK7","PARP16","PARP6","PARP8","PDIA2","PDIA3","PDIA4","PDIA6","PDX1","PIGBOS1","PIK3R1","PIK3R2","PMAIP1","PML","PPP1R15A","PPP1R15B","PPP2CB","PRKN","PSMC6","PTPN1","PTPN2","QRICH1","RACK1","RASGRF1","RASGRF2","RCN3","RHBDD1","RHBDD2","RNF103","RNF121","RNF139","RNF175","RNF183","RNF185","RNF186","RNF5","RNFT1","RNFT2","RPAP2","SCAMP5","SDF2L1","SEC16A","SEC61B","SEL1L","SEL1L2","SELENOK","SELENOS","SERINC3","SERP1","SERP2","SESN2","SGF29","SGTA","SIRT1","SRPX","STC2","STT3B","STUB1","SVIP","SYVN1","TARDBP","TBL2","THBS1","THBS4","TMBIM6","TMCO1","TMED2","TMEM117","TMEM129","TMEM238L","TMEM258","TMEM259","TMEM33","TMEM67","TMTC3","TMTC4","TMUB1","TMUB2","TMX1","TNFRSF10B","TOR1A","TP53","TRAF2","TRIB3","TRIM13","TRIM25","TTC23L","TXNDC12","UBA5","UBAC2","UBE2G2","UBE2J1","UBE2J2","UBE4A","UBE4B","UBQLN1","UBQLN2","UBXN1","UBXN10","UBXN2A","UBXN4","UBXN6","UBXN8","UFC1","UFD1","UFL1","UFM1","UGGT1","UGGT2","UMOD","USP13","USP14","USP19","USP25","VAPB","VCP","WFS1","XBP1","YOD1"]
sc.tl.score_genes(adata, RESPONSE_TO_ER_STRESS,score_name='ER_Stress')

HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = ["ABCC1","ATOX1","CAT","CDKN2D","EGLN2","ERCC2","FES","FTL","G6PD","GCLC","GCLM","GLRX","GLRX2","GPX3","GPX4","GSR","HHEX","HMOX2","IPCEF1","JUNB","LAMTOR5","LSP1","MBP","MGST1","MPO","MSRA","NDUFA6","NDUFB4","NDUFS2","NQO1","OXSR1","PDLIM1","PFKP","PRDX1","PRDX2","PRDX4","PRDX6","PRNP","PTPA","SBNO2","SCAF4","SELENOS","SOD1","SOD2","SRXN1","STK25","TXN","TXNRD1","TXNRD2"]
sc.tl.score_genes(adata, HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, score_name='ROS')

HALLMARK_APOPTOSIS = ["ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1","BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1","BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1","CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP","EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D","H1-0","HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6","IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","LUM","MADD","MCL1","MGMT","MMP2","NEDD9","NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PLPPR4","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA","RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1","TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP"]
sc.tl.score_genes(adata, HALLMARK_APOPTOSIS,score_name='Apoptosis')


sc.pl.violin(adata, keys=['Glycolysis','ER_Stress','ROS','Apoptosis'], groupby='batch', rotation=90, order = ['In Vivo Atlas - Wang 2025','Primary 2D - This Study','Primary Slice - This Study','iNeurons - Lin 2021','FEBO - Hendriks 2024', 'Organoids - He 2024' ],save='stress')
