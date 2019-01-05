
Dexpress <- function(exprs, phosexprs,case,ctrl){
    
	#The distribution of data
	exprsdis <- readline("Have you normalized exprs dataset? Yes(y) or NO(n)?")
	phosexprsdis <- readline("Have you normalized phosexprs dataset? Yes(y) or NO(n)?")
	
    #all proteins & all phosphoproteins
    allpro=rownames(exprs)
    allphos=rownames(phosexprs)
    
    #Differentially expressed proteins
    fcall=vector()
    pall=vector()
	if(exprsdis == 'y'){
		for(i in 1:nrow(exprs)){
			x=as.numeric(exprs[i,case])
			y=as.numeric(exprs[i,ctrl])
			fc=mean(x)- mean(y)
			fcall=c(fcall,fc)
			p=as.numeric(t.test(x,y)[3])
			pall=c(pall,p)
		}
	}else if(exprsdis == 'n'){
		for(i in 1:nrow(exprs)){
			x=as.numeric(exprs[i,case])
			y=as.numeric(exprs[i,ctrl])
			fc=mean(x)- mean(y)
			fcall=c(fcall,fc)
			p=as.numeric(wilcox.test(x,y)[3])
			pall=c(pall,p)
		}
	}else{
		stop("Please print 'y' or 'n' to the question about exprs")
		
	}
	names(fcall)=allpro
    pall_ad=p.adjust(pall)
    names(pall_ad)=allpro
    DE=names(pall_ad[pall_ad<0.05])
    
    #Differentially expressed phosphoproteins
    fcall_phos=vector()
    pall_phos=vector()
    if(phosexprsdis == 'y'){
		for(i in 1:nrow(phosexprs)){
			x=as.numeric(phosexprs[i,case])
			y=as.numeric(phosexprs[i,ctrl])
			fc=mean(x)- mean(y)
			fcall_phos=c(fcall_phos,fc)
			p=as.numeric(t.test(x,y)[3])
			pall_phos=c(pall_phos,p)
		}
	}else if(phosexprsdis == 'n'){
		for(i in 1:nrow(phosexprs)){
			x=as.numeric(phosexprs[i,case])
			y=as.numeric(phosexprs[i,ctrl])
			fc=mean(x)- mean(y)
			fcall_phos=c(fcall_phos,fc)
			p=as.numeric(wilcox.test(x,y)[3])
			pall_phos=c(pall_phos,p)
		}
	}else{
		stop("Please print 'y' or 'n' to the question about Phosexprs")
	}
    names(fcall_phos)=allphos
    pall_ad_phos=p.adjust(pall_phos)
    names(pall_ad_phos)=allphos
    DE_phos=names(pall_ad_phos[pall_ad_phos<0.05])
 
    #results
    DEexpression <- list(DEpro=DE,DEphos=DE_phos,fcpro=fcall,fcphos=fcall_phos,allpro=allpro,allphos=allphos)
    return(DEexpression)
}