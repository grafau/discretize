#!/usr/bin/Rscript

##### COMMAND LINE OPTIONS

args = commandArgs(trailingOnly=TRUE)


if (length(args)<3) 
{
  stop("At least three arguments have to be supplied:\n
\tARG1\tObligatory\t[INTEGER] upper bound number of k subpopulations
\tARG2\tObligatory\t[FILE] meta information with ID, latitude, longitude, country
\tARG3\tObligatory\t[FILE] pairwise distance matrix, without headers or row names
\tARG4\tOptional\t[FILE] color scheme for graphics, columns correspond to each k, rows represent population"
, call.=FALSE)
}

max_k = as.numeric(args[1])


##### LOAD LIBRARIES
library(shape)
library(cluster)
library(geosphere)
library(RColorBrewer)
library(maps)


##### LOAD DATA

meta <- read.table(args[2])
nams <- c("Mauritania", "Senegal", "Guinea-Bissau", "Gambia", "Guinea", "Sierra Leone", "Liberia", "Ivory Coast", "Mali", "Ghana", "Togo", "Benin", "Nigeria", "Burkina Faso", "Niger", "Cameroon")
#nams <- map("world", namesonly=TRUE, plot=FALSE)
 

gen_dist <- read.table(args[3])
geo_dist <- as.matrix(distm(meta[,c(3,2)], fun=distHaversine))

gen_null <- sum(unlist(gen_dist))
geo_null <- sum(unlist(geo_dist))


if (length(args)==3) 
{
  col_pal <- matrix(rep(brewer.pal(max_k, "Paired"), max_k),nrow=max_k)
} else {
col_pal <- read.table(args[4])
}


###### INIT VARIABLES

# get global MDS
glob <- cmdscale(gen_dist, eig=TRUE,k=3)
xl  <- round(glob$eig[1]/sum(glob$eig[which(glob$eig > 0)]),3)*100
yl <- round(glob$eig[2]/sum(glob$eig[which(glob$eig > 0)]),3)*100

## Prep evaluation tools 

dgen_Sil <- c(NA)	# total silhouette in genetic space for total clustering
dgeo_Sil <- c(NA)	# total silhouette in geographic space for total clustering

sgen_Sil <- c(NA)	# total silhouette in genetic space for filtered clustering
sgeo_Sil <- c(NA)	# total silhouette in geographic space for filtered clustering
 
dgen_TSS <- c(NA)		# proportion of genetic sum of squares explained by total clustering
dgeo_TSS <- c(NA)		# proportion of geographic sum of squares explained by total clustering

sgen_TSS <- c(NA)		# proportion of genetic sum of squares explained by filtered clustering
sgeo_TSS <- c(NA)		# proportion of geographic sum of squares explained by filtered clustering

sind_Prop <- c(NA)	# proportion of individuals that are kept after silhoutte filtering


##### INTERATE THROUGH DIFFERENT K's 
for (k in 2:max_k)
{

### 1 ### CALCULATE EMPIRICAL SILHOUTTE SCORES FROM SIMULATIONS OF EVOLVING 50:50 POPULATIONS
# simulate admixture
	clu <- pam(gen_dist, k, diss =T, keep.diss = F)

# pre-filter data for bad fit indviduals

	list <- as.numeric(substring(row.names(clu$silinfo$widths[which(clu$silinfo$widths[,3] > 0),]),2))

## Iteratively estimate silhouette thresholds and filter out individuals
	old <- 0

	while ( length(list) != old )
	{	
		sil <- c()
		old <- length(list)
		flt_dist <- gen_dist[list,list]
		flu <- pam(flt_dist,  k, diss =T, keep.diss = F)
	
	# choose which two populations (medoids) are closest 
		meds <- as.matrix(flt_dist[flu$id.med,flu$id.med])
		diag(meds) <- 1
	
		sist <- apply(meds,1,which.min)	
		
		sil <- c()
	
		for (I in 1:k)
		{
			pop1 <- I
			pop2 <- sist[I] 
	
	
			ind1 <- which(flu$clustering == pop1)
			ind2 <- which(flu$clustering == pop2)
	## calculate mean of empirical silhouette scores for simulated admixtures
		
		
			sim <- apply(rbind(apply(flt_dist[which(flu$clustering == pop1),], 2, mean)*1, apply(flt_dist[which(flu$clustering == pop2),], 2, mean)*1), 2, mean)
			sim_dist <- cbind(rbind(flt_dist, sim), c(sim,0))
			rownames(sim_dist)[nrow(sim_dist)] <- 5050
	
			sim <- apply(rbind(apply(flt_dist[which(flu$clustering == pop1),], 2, mean)*1.2, apply(flt_dist[which(flu$clustering == pop2),], 2, mean)*0.8), 2, mean)
	                sim_dist <- cbind(rbind(sim_dist, c(sim, 0 )), c(sim,0,0))
	                rownames(sim_dist)[nrow(sim_dist)] <- 6040
	
			sim <- apply(rbind(apply(flt_dist[which(flu$clustering == pop1),], 2, mean)*0.8, apply(flt_dist[which(flu$clustering == pop2),], 2, mean)*1.2), 2, mean)
	                sim_dist <- cbind(rbind(sim_dist, c(sim, 0, 0)), c(sim, 0, 0, 0))
	                rownames(sim_dist)[nrow(sim_dist)] <- 4060
		
	# get empirical Silhouette scores
		
				
			sdis1 <- apply(sim_dist[(nrow(sim_dist)-2):nrow(sim_dist), c(ind1,ind1)],1,mean)
			sdis2 <- apply(sim_dist[(nrow(sim_dist)-2):nrow(sim_dist), c(ind2,ind2)],1,mean)		
		
			silmax <- apply(rbind(sdis1,sdis2),2,max)
		
			sils <- abs(sdis1-sdis2)/silmax
	
				
			sil <- c(sil,max(sils))
			
		}
		print(paste("k: ",k,", Kept indviduals: ",old, sep =""))
	### 2 ### RUN KMEDOID CLUSTERING AND FILTER IT FOR SILHOUETTE  SCORES
	
	# run clustering algorithm
		tlu <- pam(gen_dist[list,list], k, diss =T, keep.diss = F)
		
		list <- c()
		for (o in 1:k)
		{
			list <- c(list, as.numeric(row.names(tlu$silinfo$widths[which(tlu$silinfo$widths[,1] == o & tlu$silinfo$widths[,3] >= sil[o]),])) )
		}
		if (length(list) <= k)
		{
			break
		}
	}

### END of filtering iterations

### GET COLOUR PALETTES FOR THIS K


	pal <- as.character(col_pal[,k-1])
	dim <- paste(pal, 80, sep ="")
	
# write tables for PLINK files filtering 
	
	pops <- c()
	for (K in 1:k)
	{
		poop <- names(which.max(table(meta[list[which(clu$clustering[list] == K)],4])))
		if (poop %in% pops)
		{
			poop <- paste(poop,"A",sep="")
		}
		pops <- c(pops, poop )
	}

	write.table(cbind(sort(as.character(meta[list,1])), sort(as.character(meta[list,1])) ),paste("list2keepK",k,".txt",sep=""),quote = F, sep = "\t",row.names = F, col.names = F)
	opps <- rep(NA, nrow(meta))
	opps[list] <- pops[clu$clustering[list]]
	write.table(cbind(opps,as.character(meta[,1]))[order(meta[,1]),], paste("list2annoK",k,".txt",sep=""),quote = F, sep = " ",row.names = F, col.names = F)


	pdf(paste("maps_",k,".pdf",sep=""))
# plot MAP
	par(mfrow = c(1,1))
	map("world", regions = nams)
	points(meta[,3],meta[,2],col = rgb(0.8,0.8,0.8,0.5),pch =19,cex=0.6)
	points(meta[list,3],meta[list,2],col = dim[clu$clustering[list]],pch =19,cex=0.6)
	dev.off()

	pdf(paste("clust_",k,".pdf",sep=""))
# plot MDS

        plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)
	text(-glob$points[,2],-glob$points[,1], label = meta[,4], col = rgb(0.8,0.8,0.8,0.5), cex = 0.5)
	if (length(list) > 0 )
	{
        	text(-glob$points[list,2],-glob$points[list,1], label = meta[list,4],col = pal[clu$clustering[list]], cex = 0.5)
	}
#	points(glob$points[clu$id.med,1],glob$points[clu$id.med,2], col = pal, cex = 1, pch =8)

# plot PIES
	if (length(list) > 0 )
	{
		ftab <- table(as.data.frame(cbind(as.character(meta[,4]),clu$clustering)))
		stab <- table(as.data.frame(cbind(as.character(meta[list,4]),clu$clustering[list])))
		
		fsum <- apply(ftab,2,sum)
		ssum <- apply(stab,2,sum)
	
		fscal <- fsum/max(fsum)
		sscal <- ssum/max(fsum)

		allc <- c(brewer.pal(11,"RdYlBu"), brewer.pal(11,"PRGn"))
		pal2 <- allc[1:length(unique(meta[list,4]))]
	
		par(mfrow = c(3,2))
		for (i in 1:k)
		{
			to <- grepl(paste("^",i,"$",sep=""), colnames(stab))
			if (any(to))
			{	
	        		pie(stab[,which(to)], radius = sqrt(sscal[which(to)]))
	        		text(x=1,y=1,labels = paste(ssum[which(to)],"/",fsum[i],sep=""))
			} else {
				pie(c(0,1), col = NA, radius = 0.001)
	        		text(x=1,y=1,labels = paste("0/",fsum[i],sep=""))
			}
			plotcircle(sqrt(fscal[which(to)]))
		}
	}

# sort and order individuals
	cent <- matrix(NA, k,2)
	for (o in 1:k)
	{
		cent[o,1] <- mean(meta[which(clu$clustering ==o),2])
		cent[o,2] <- mean(meta[which(clu$clustering ==o),3])
	}
	ind <- hclust(dist(cent))
	ord0 <- order(ind$order)[clu$silinfo$widths[,1]]
	ord1 <- as.numeric(substring(row.names(clu$silinfo$widths[order(ord0),]),2))

# plot HEAT
	par(mfrow = c(1,1))
	image(as.matrix(gen_dist[ord1,ord1]),ylab = NA ,xlab = NA, xaxt = 'n', yaxt = 'n')
	par(mfrow = c(4,1))
	barplot(height = rep(10,dim(gen_dist)[1]), width = rep(1,dim(gen_dist)[1]), rep(1,dim(gen_dist)[1]), col = pal[clu$clustering[ord1]], border = NA, space = 0, yaxt = 'n', xaxt = 'n' )
	barplot(height = rep(10,dim(gen_dist)[1]), width = rep(1,dim(gen_dist)[1]), rep(1,dim(gen_dist)[1]), col = pal2[meta[ord1,4]], border = NA, space = 0, yaxt = 'n', xaxt = 'n')

	dev.off()

### Calculate evaluation statistics
	dgeo_Pil <- c()
	sgeo_Pil <- c()

	dgen_WSS <- c()
	dgeo_WSS <- c()
	sgen_WSS <- c()
	sgeo_WSS <- c()

	for (J in 1:k)
        {
                pop1 <- J
                pop2 <- sist[J]

                ind1 <- which(clu$clustering == pop1)
                ind2 <- which(clu$clustering == pop2)	

		for (S in ind1)
		{
			d1 <- mean( unlist(geo_dist[S,ind1])) 
			d2 <- mean( unlist(geo_dist[S,ind2]))
			dgeo_Pil <- c(dgeo_Pil, abs(d1-d2)/max(c(d1,d2)))
		}

		dgen_WSS <- c(dgen_WSS, sum(unlist(gen_dist[ind1,ind1])) )
	        dgeo_WSS <- c(dgeo_WSS, sum(unlist(geo_dist[ind1,ind1])) )

		ind1 <- intersect(ind1,list)						# only filtered inds
		ind2 <- intersect(ind2,list)						# only filtered inds

                for (S in ind1)
                {
                        d1 <- mean( unlist(geo_dist[S,ind1]))   
                        d2 <- mean( unlist(geo_dist[S,ind2]))
                        sgeo_Pil <- c(sgeo_Pil, abs(d1-d2)/max(c(d1,d2)))
                }

	        sgen_WSS <- c(sgen_WSS, sum(unlist(gen_dist[ind1,ind1])) )
       		sgeo_WSS <- c(sgeo_WSS, sum(unlist(geo_dist[ind1,ind1])) )		
	}

	dgen_Sil <- c(dgen_Sil, mean(clu$silinfo$widths[,3]))
	dgeo_Sil <- c(dgeo_Sil, mean(dgeo_Pil))

	dgen_TSS <- c(dgen_TSS, sum(dgen_WSS)/gen_null)
	dgeo_TSS <- c(dgeo_TSS, sum(dgeo_WSS)/geo_null)

	sgen_Sil <- c(sgen_Sil, mean(c(tlu$silinfo$widths[,3],rep(0,length(clu$clustering)-length(list))) ))
	sgeo_Sil <- c(sgeo_Sil, mean(c(sgeo_Pil,rep(0,length(clu$clustering)-length(list))), na.rm = T ))

	sgen_null <- sum(unlist(gen_dist[list,list]))
	sgeo_null <- sum(unlist(geo_dist[list,list]))

	sgen_TSS <- c(sgen_TSS, sum(sgen_WSS)/sgen_null)
	sgeo_TSS <- c(sgeo_TSS, sum(sgeo_WSS)/sgeo_null)

		
	sind_Prop <- c(sind_Prop,length(list)/nrow(meta))

}

pdf("WholeEval.pdf", useDingbats = F)
plot(dgen_Sil, type = 'l', col = rgb(0.1,0.8,0.3), ylim = c(0,1), xlab  = "K", ylab = "Fraction")
lines(dgeo_Sil, col = rgb(0.8,0.1,0.3))
lines(1-dgen_TSS, col = rgb(0.4,0.8,0.5), lty = 2)
lines(1-dgeo_TSS, col = rgb(0.8,0.5,0.4), lty = 2)
dev.off()

pdf("FilterEval.pdf", useDingbats = F)
plot(sgen_Sil, type = 'l', col = rgb(0.1,0.8,0.3), ylim = c(0,1), xlab  = "K", ylab = "Fraction")
lines(sgeo_Sil, col = rgb(0.8,0.1,0.3))
lines(1-sgen_TSS, col = rgb(0.4,0.8,0.5), lty = 2)
lines(1-sgeo_TSS, col = rgb(0.8,0.5,0.4), lty = 2)
lines(sind_Prop, col = rgb(0.2,0.2,0.8), lty = 3)
dev.off()


