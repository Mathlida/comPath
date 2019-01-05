SPIAcom <- function(DEexpression){

    #overlap of all proteins and all phosphoproteins
    all_com=unique(c(DEexpression$allpro,DEexpression$allphos))

    #overlap of DEpro and DEphos and the sum of their fold change
    DE=DEexpression$DEpro
    DE_phos=DEexpression$DEphos
    fc_pro=abs(DEexpression$fcpro)
    fc_phos=abs(DEexpression$fcphos)
    DE_com=c(DE,DE_phos)
    fc_DE_com=c(fc_pro[DE],fc_phos[DE_phos])
    fc_DE_com=tapply(fc_DE_com,names(fc_DE_com),sum)[intersect(DE,DE_phos)]

    #do SPIA pathway analysis
    result <- spia(de= fc_DE_com,all=all_com,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher")
    return(result)
}
