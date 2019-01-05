GeneSetcom <- function(DEexpression){
	
	
	#get the pathway imformation
	kg.hsa=kegg.gsets("hsa")
    PW_KEGG=kg.hsa$kg.sets
	ListGSC <- list(PW_KEGG=PW_KEGG)
	
	DE=DEexpression$DEpro
    DE_phos=DEexpression$DEphos
    fc_pro=abs(DEexpression$fcpro)
    fc_phos=abs(DEexpression$fcphos)
	fc_all=c(fc_pro, fc_phos)
	fc_com_name = intersect(names(fc_pro),names(fc_phos))
	fc_com = tapply(fc_all,names(fc_all),sum)
	
	data4enrich_com = sort(fc_com[fc_com_name])
	data4enrich_com_name = names(data4enrich_com)
	data4enrich_com = as.numeric(data4enrich_com)
	names(data4enrich_com) = data4enrich_com_name

	#hits_com <- unique(c(DE,DE_phos))
	hits_com <- intersect(DE,DE_phos)
	hits_com = as.character(unlist(hits_com))

	gsca_com <- new("GSCA", listOfGeneSetCollections=ListGSC,geneList=data4enrich_com, hits=hits_com)
	gsca_com <- preprocess(gsca_com, species="Hs", initialIDs="Entrez.gene",
                        keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                        orderAbsValue=FALSE)

	results<-analyze(gsca_com, para=list(pValueCutoff=0.05, nPermutations=2000, exponent=1,pAdjustMethod="BH"))
	return(results)
}