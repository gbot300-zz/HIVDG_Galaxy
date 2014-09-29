###LOAD PACKAGES
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(seqinr))

#INPUT and specify outputfiles
nucfile <- commandArgs()[6]			#e.g., 'euclidean_nuc.fas'
namesfile <- commandArgs()[7]			#e.g., 'control_names.txt'
pvalf <- commandArgs()[8]
plotf <- commandArgs()[9]

#fxn to convert nucs to aa 

amino2idx <- function(amino)
{	
	switch(EXPR=amino, F=1, L=2, I=3, M=4, V=5, S=6, P=7, T=8, A=9, Y=10, 
				     H=11, Q=12, N=13, K=14, D=15, E=16, C=17, W=18, R=19, G=20, "*"=21, X=22, "-"=23)
}

tabulate <- function(x)
{
	table(factor(x, levels = c("F","L","I","M","V","S","P","T","A","Y","H","Q","N","K","D","E","C","W","R","G","*","X","-")))
}


#read the input data

align <- read.dna(nucfile, format="fasta", as.character=T) 
align <- t(apply(align, 1, translate))  #to aa

nseqs <- nrow(align)
rnames <- rownames(align)
vacnames <- as.character(read.csv(namesfile, header=F)[[1]])
vacseqs <- c()

for(i in 1:length(rnames))
	{ if(rnames[i] %in% vacnames) 
	vacseqs <- c(vacseqs, TRUE) else vacseqs <- c(vacseqs, FALSE) } 

control <- apply(align[!vacseqs,], 2, tabulate)				#TABULATE BY SITE not seq, i.e. row = codon site,rows == aa, col=count for each residue at that sitenobu
vaccinated <- apply(align[vacseqs,], 2, tabulate)

permdata <- matrix(mapply(amino2idx, align), nrow=nseqs, byrow=F)		#REPLACES EACH AA WITH ITS NUM CODE,nobu

n_con=sum(control[,1])		#Number of control sequences
n_vac=sum(vaccinated[,1])
n_tot=n_con+n_vac

####### check that all columns sum to group size
for(i in 1:ncol(control)){
	if(sum(control[,i])!=n_con){
		print(paste("Error: Site ", i, " of control group does not sum to ", n_con, "!", sep=""))
	}
}

for(i in 1:ncol(vaccinated)){
	if(sum(vaccinated[,i])!=n_vac){
		print(paste("Error: Site ", i, " of vaccinated group does not sum to ", n_vac, "!", sep=""))
	}
}
#######

# proportions i.e for each residue at each site
p_con=control/n_con
p_vac=vaccinated/n_vac

# variances
v_con=p_con*(1-p_con)/n_con
v_vac=p_vac*(1-p_vac)/n_vac

#differences
p_dif=p_con-p_vac
v_dif=(n_con-1)*v_con/(n_tot-2)+(n_vac-1)*v_vac/(n_tot-2)


#### check for sites where the proportion of an AA = 1 in control and different AA = 1 vaccinated
for(site in 1:ncol(p_con)){
	if(!all(p_vac[,site]==p_con[,site])){   #look at sites with diff only
		aa_con=-999
		aa_vac=-999
		for(aa in 1:nrow(p_con)){
			if(p_con[,site][aa]==1){	#if aa has prop =1 then conserved site
				aa_con=aa
			}			
			if(p_vac[,site][aa]==1){
				aa_vac=aa
			}
		}
		if(aa_con!=aa_vac && aa_con!=-999 && aa_vac!=-999){
			print(paste("Completely different site: ", site, sep=""))
		}
	}
}

#######################################################
########## STANDARDIZED EUCLIDEAN DISTANCE ############
#######################################################

Z=rep(NA, ncol(p_con))
for(site in 1:ncol(p_con)){
	z=p_dif[,site]^2
	for(aa in 1:nrow(p_con)){
		if(v_dif[,site][aa]>0 || p_dif[,site][aa]!=0){
			z[aa]=z[aa]/v_dif[,site][aa]
		}else{
			z[aa]=0
		}
	}
	Z[site]=sum(z)
}

### check for infinitely large test stats associated with "completely different" sites

if(!all((Z==Inf)==FALSE)){
	print("At least one test statistic is infinite")
}

########### PERMUTATION TEST #########
permnum=100
Zperm=matrix(NA, ncol=ncol(p_con), nrow=permnum)
for(rep in 1:permnum){
	random=sample.int(n_tot, n_con)
	permidx_con=permdata[random,]
	permidx_vac=permdata[-random,]
	permfreq_con=permfreq_vac=c()
	for(site in 1:ncol(permidx_con))
	{
	  permfreq_con=cbind(permfreq_con, table(factor(permidx_con[,site], levels=1:23)))
	  permfreq_vac=cbind(permfreq_vac, table(factor(permidx_vac[,site], levels=1:23)))
	}

	# proportions
	p_permcon=permfreq_con/n_con
	p_permvac=permfreq_vac/n_vac

	# check that p-hats sum to one
	for(site in 1:ncol(p_permcon)){
		if(round(sum(p_permcon[,site]),15)!=1 || round(sum(p_permvac[,site]),15)!=1){
			print(paste("Error: p-hats at site ", site, " do not sum to 1!", sep=""))
			print(paste(sum(p_permcon[,site]), sum(p_permcon[,site]), sep="  "))
		}
	}

	# variances
	v_permcon=p_permcon*(1-p_permcon)/n_con
	v_permvac=p_permvac*(1-p_permvac)/n_vac

	#differences
	p_permdif=p_permcon-p_permvac
	v_permdif=(n_con-1)*v_permcon/(n_tot-2)+(n_vac-1)*v_permvac/(n_tot-2)

	#compute test stat
	for(site in 1:ncol(p_permcon)){
		zperm=p_permdif[,site]^2
		for(aa in 1:nrow(p_permcon)){
			if(v_permdif[,site][aa]>0 || p_permdif[,site][aa]!=0){
				zperm[aa]=zperm[aa]/v_permdif[,site][aa]
			}else{
				zperm[aa]=0
			}
		}
		Zperm[,site][rep]=sum(zperm)
		rm(zperm)
	}
	rm(permfreq_con, permfreq_vac, p_permcon, p_permvac, v_permcon, v_permvac, p_permdif, v_permdif)
} # end permutations

# Calculate p-values based on permutation distribution
Zpvalue=rep(NA, length(Z))
for(site in 1:length(Z)){
	Zpvalue[site]=length(which(Zperm[,site]>=Z[site]))/permnum
}

#save.image(paste("Euclidean.RData", sep=""))


#### GRAPH #### 

numsites=length(Zpvalue)
sig_sites=rep("", numsites)
for(i in 1:numsites){
   if(Zpvalue[i]<0.05)
   {
  	sig_sites[i]=paste(i)
   }
 }
sig90_sites=rep("", numsites)  #sites with pvalues btwn 0.05 and 0.1
for(i in 1:numsites){
   if(Zpvalue[i]<0.1)
   {
	if(Zpvalue[i]>0.049){
 	  sig90_sites[i]=paste(i)}
    }
} 

#write permutation pvalue for each site

#write.table(Zpvalue,"EuclideanPvalues.txt")
write.table(as.matrix(Zpvalue),file=pvalf,sep="\t")

#plot and save 1-pvalues for each site
# pdf("EuclideanDistancesPLOT.pdf")
pdf(plotf)
 plot(1-Zpvalue, type="h", ylim=c(0,1.05), xlim=c(0,numsites+10),xaxs="i",yaxs="i", cex.axis=0.9,cex.lab=0.9, xlab="", main=paste("Comparison of Standardized Euclidean Distance"), cex.main=0.9, ylab="1-pvalue")
 abline(h=0.95, col="red3")
 abline(h=0.90, col="blue3")
 text(1-Zpvalue+0.02, labels=sig_sites,font=2, cex=0.35)
 text(1-Zpvalue+0.01, labels=sig90_sites, cex=0.35)
 text(numsites-5, 0.915, label="0.1 sig", cex=0.55)
 text(numsites-5, 0.965, label="0.05 sig", cex=0.6, font=2)
 dev.off()

