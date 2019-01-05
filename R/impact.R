impact <- function(DEexpression){
    
    #create plot file
	mainDir=getwd()
	subDir="plots"
	dir.create(file.path(mainDir,subDir),showWarnings = FALSE)
	setwd(file.path(mainDir,subDir))
	
    #get the pathway imformation
    kg.hsa=kegg.gsets("hsa")
    PW_KEGG=kg.hsa$kg.sets
    
    #if there exists difference between
    allpro=DEexpression$allpro
    allphos=DEexpression$allphos
    fc_pro=zscore(DEexpression$fcpro)
    fc_phos=zscore(DEexpression$fcphos)
	minfc=min(c(fc_pro,fc_phos))
	maxfc=max(c(fc_pro,fc_phos))
	limit=max(abs(minfc),abs(maxfc))
    n=0
    m=0
    Dlist=list()
    Slist=list()
    for(i in 1:length(PW_KEGG)){
        pathmen=as.vector(unlist(PW_KEGG[i]))
        path_fc_pro=vector()
        path_fc_phos=vector()
        for(j in 1:length(pathmen)){
            if(pathmen[j] %in% allpro){
                path_fc_pro=c(path_fc_pro,fc_pro[pathmen[j]])
            }
            if(pathmen[j] %in% allphos){
                path_fc_phos=c(path_fc_phos,fc_phos[pathmen[j]])
            }
        }
        if(length(path_fc_pro)>1 && length(path_fc_phos)>1){
            p=wilcox.test(path_fc_pro,path_fc_phos)[3]
			temp=c(path_fc_pro,path_fc_phos)
            type=c(rep('pro',length(path_fc_pro)),rep('phos',length(path_fc_phos)))
            pathdata=cbind.data.frame(genes=names(temp),fc=unname(temp),type=type)
			pathname=names(PW_KEGG[i])
            if(p<0.05){
                n=n+1
                file=paste(substr(pathname,1,8),"_boxplot",".png",sep="")
                #Boxplot
				ggboxplot(pathdata,pathname)
                ggsave(file,width = 4, height = 4)
				#KEGG pahtway map
				#download xml & png file
				filedir=system.file("extdata",package="KEGGprofile")
				download.kegg(pathway.id = substr(pathname,4,8), species = "hsa", kegg.dir = filedir,file.type=c("xml", "png"))
				mapfile=paste(filedir,"/","map",substr(pathname,4,8),".png",sep="")
				hsafile=paste(filedir,"/",substr(pathname,1,8),".png",sep="")
				file.copy(from=mapfile,to=hsafile,overwrite=TRUE)
				temp <- rbind.fill.matrix(t(path_fc_pro), t(path_fc_phos))
				temp[is.na(temp)] <- 0
				temp <- t(temp)
				colnames(temp)=c("Proteome","Phosphoproteome")	
				col<-col_by_value(temp,col=colorRampPalette(c('yellow','white','red'))(1024),range=c(-limit,limit))
				invisible(capture.output(keggpath <- plot_pathway(temp,type="bg",bg_col=col,text_col="black",magnify=1.2,species='hsa',database_dir=filedir,pathway_id=substr(pathname,4,8),groups=c("Proteome ","Phosphoproteome"),pathway_min=3)))
				img=readPNG(paste(substr(pathname,1,8),"_profile_bg",".png",sep=""))
				width=ncol(img)
				height=nrow(img)
				png(file=paste(substr(pathname,1,8),"_profile_bg",".png",sep=""),width=width,height=height)
				par(yaxs="i")
				par(xaxs="i")
				par(mar=c(0,0,0,0))
				plot(c(0,width),c(0,height),  type='n',xlab="",ylab="")
				rasterImage(img,  0,  0,  width,  height,interpolate=F)
				image.plot(legend.only=T,zlim=c(-limit,limit),col=colorRampPalette(c('yellow','white','red'))(1024), horizontal =TRUE,legend.shrink=0.1,smallplot = c(.9,.99, .92,.95),axis.args = list(cex.axis = 2))
				dev.off()
                #file.remove("hsa",paste(substr(pathname,4,8),".xml",sep=""))
				Dlist[[pathname]]=pathdata[order(pathdata$fc,decreasing = TRUE),]
            }else{
                m=m+1
                Slist[[pathname]]=pathdata[order(pathdata$fc,decreasing = TRUE),]
            }
        }
    }
	print('Proteins and Phosphoproteins have different impacts on the pathways:')
	print(names(Dlist))
    results <- list(DiffPath=n,NODiffPath=m,Dlist=Dlist,Slist=Slist)
    setwd('..')
	return(results)
}