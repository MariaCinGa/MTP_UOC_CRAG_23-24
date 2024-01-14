args = commandArgs(trailingOnly=TRUE)

if(length(args)<3) {
  args <- array(dim=c(3))
  args[1] <- "/Users/sramos/Desktop/TFM-UOC-Proposals-2023/NEW_simulations_EHHS-GWAS_2023/Step1-SLiM_plus_Stats/"
  args[2] <- "20"
  args[3] <- "500" 
}
setwd(args[1])
scenario <- c(1:as.numeric(args[2]))
nsam <- as.numeric(args[3])

sum_funct_dom_ph <- function(genotype,effsize) {
  pos <- genotype[,1]
  nsnps <- length(pos)
  count.row.1 <- which(genotype[,2]==1)
  count.pos.1 <- genotype[count.row.1,1]
  count.row.2 <- which(genotype[,2]==2)
  count.pos.2 <- genotype[count.row.2,1]
  
  eff.size.1 <- NULL
  for(i in c(1:length(count.pos.1))) {
    eff.size.1 <- c(eff.size.1,which(effsize[,1]==count.pos.1[i])) #heterozygotes
  }
  eff.size.2 <- NULL
  for(i in c(1:length(count.pos.2))) {
    eff.size.2 <- c(eff.size.2,which(effsize[,1]==count.pos.2[i])) #homozygotes derived
  }
  
  effect.total <- 0
  effect.total <- sum(effsize[eff.size.2,3]) + sum(effsize[eff.size.1,2]*effsize[eff.size.1,3])
  return(effect.total)
}
h <- 0.7 #relationship genotype vs phenotype 
s <- 6 
for(s in scenario) {
  pdf(file=sprintf("Plots_Phenotype_fitness_stats_s%.0f.pdf",s))
  # par(mfrow=c(3,2))
  par(mfrow=c(2,2))
  d <-0 #iteration
  ph <- read.table(file=sprintf("./slim.output_s%.0f_it%.0f.p1p2_table.txt",s,d),skip=27,header=T)
  #this geno includes the genotype for 1000 p1 ind and 1000 p2 ind for each position
  geno <- read.table(file=sprintf("./slim.output_s%.0f_it%.0f.p1p2.txt.geno.txt",s,d),header=T)
  
  #remove all positions that have no effect
  pos.ph <- which(ph$EFFSIZE!=0)
  ph <- ph[pos.ph,]
  #keep only geno rows where positions are in ph 
  pos.geno <- NULL
  j <- 1;
  for(i in ph[,1]) {
    while(j<length(geno[,1]) && i!=geno[j,2]) {
      j <- j + 1
    }
    pos.geno <- c(pos.geno,j);
    j <- j + 1
  }
  geno<-geno[pos.geno,]
  
  geno1 <- geno[,c(1,2,3:(3+nsam-1))]
  geno2 <- geno[,c(1,2,(3+nsam):(3+2*nsam-1))]

  #plot variants effects and selective coefficient
  plot(ph[,c(2,5)],main=sprintf("S%.0f P1: SelCoef vs EffSize size",s))
  plot(ph[,c(3,5)],main=sprintf("S%.0f P2: SelCoef vs EffSize size",s),col="blue",type="p")

  #calculate the phenotype of each individual, 
  # considering the number of mutations and the dominance that each individual contains
  # The dominance gives a proportion on the phenotype.
  N1 <- length(geno1[1,])-2
  M1 <- length(geno1[,1])-1
  pheno1 <- array(dim=c(N1))
  for(ind in c(1:N1)) {
    pheno1[ind] <- sum_funct_dom_ph(genotype=geno1[,c(2,ind+2)],effsize=ph[,c(1,4,5)])
  }
  N2 <- length(geno2[1,])-2
  M2 <- length(geno2[,1])-1
  pheno2 <- array(dim=c(N2))
  for(ind in c(1:N2)) {
    pheno2[ind] <- sum_funct_dom_ph(genotype=geno2[,c(2,ind+2)],effsize=ph[,c(1,4,5)])
  }

  #To do fair comparison in relation to the given models, all scenarios should have equal heritability (say 80%)
  #include environmental variance, put the same quantity than genotypic variance to achieve h2 ~80%
  # variance.snps <- function(pf,dom,a) { #estimate genetic variance d=dominance, a=effsize, p=frequency
  #   alph <- a + dom*(1-2*pf)
  #   v <- 2*pf*(1-pf)*alph*alph
  #   return(v)
  # }
  # varF.snp1 <- variance.snps(pf=ph[,7],dom=ph[,4],a=ph[,5])
  # sum.varF.snp1 <- sum(varF.snp1)
  # sdF.pheno1 <- sqrt(sum.varF.snp1) 
  # 
  # varF.snp2 <- variance.snps(pf=ph[,8],dom=ph[,4],a=ph[,5])
  # sum.varF.snp2 <- sum(varF.snp2)
  # sdF.pheno2 <- sqrt(sum.varF.snp2) 
  # 
  # envF.eff1<- rnorm(n=N1,mean=0,sd=(sdF.pheno1/h-sdF.pheno1)) #for 50%, gen.eff=env.eff
  # envF.eff2<- rnorm(n=N1,mean=0,sd=(sdF.pheno2/h-sdF.pheno2))
  # pheno1Ff <- pheno1 + envF.eff1
  # pheno2Ff <- pheno2 + envF.eff2
  
  #OR include environmental variance using a normal distribution (h=80%)
  sd.pheno1 <- sd(pheno1)
  env.eff1 <- rnorm(n=N1,mean=0,sd=(sd.pheno1/((exp(h)-1)/(exp(1)-1))-sd.pheno1))
  sd.pheno2 <- sd(pheno2)
  env.eff2 <- rnorm(n=N2,mean=0,sd=(sd.pheno2/((exp(h)-1)/(exp(1)-1))-sd.pheno2))

  pheno1f <- pheno1 + env.eff1
  pheno2f <- pheno2 + env.eff2
  
  #write files for phenotypes + env.eff (Normal) per individual
  write.table(x=pheno1f,file=sprintf("./scenario_%s_phenotype_individual_Normal_pop1.txt",s),quote=F,col.names=F,row.names=F)
  write.table(x=pheno2f,file=sprintf("./scenario_%s_phenotype_individual_Normal_pop2.txt",s),quote=F,col.names=F,row.names=F)
  
  #write files for phenotypes + env.eff (Ve) per individual
  # write.table(x=pheno1Ff,file=sprintf("./scenario_%s_phenotype_individual_FV_pop1.txt",s),quote=F,col.names=F,row.names=F)
  # write.table(x=pheno2Ff,file=sprintf("./scenario_%s_phenotype_individual_FV_pop2.txt",s),quote=F,col.names=F,row.names=F)

  #write files for phenotpes (no env) per individual
  write.table(x=pheno1,file=sprintf("./scenario_%s_phenotype_noenv_individual_pop1.txt",s),quote=F,col.names=F,row.names=F)
  write.table(x=pheno2,file=sprintf("./scenario_%s_phenotype_noenv_individual_pop2.txt",s),quote=F,col.names=F,row.names=F)

  #final phenotype
  plot(x=pheno1,y=pheno1 + env.eff1, main=sprintf("S%.0f P1 \n corr=%.3f",s,cor(pheno1,pheno1 + env.eff1)),xlab="Genotype",ylab="Phenotype")
  plot(x=pheno2,y=pheno2 + env.eff2, main=sprintf("S%.0f P2 \n corr=%.3f",s,cor(pheno2,pheno2 + env.eff2)),col="blue",xlab="Genotype",ylab="Phenotype")

  # plot(x=pheno1,y=pheno1 + envF.eff1, main=sprintf("S%.0f P1 (VF)\n corr=%.3f",s,cor(pheno1,pheno1 + env.eff1)),xlab="Genotype",ylab="Phenotype")
  # plot(x=pheno2,y=pheno2 + envF.eff2, main=sprintf("S%.0f P2 (VF)\n corr=%.3f",s,cor(pheno2,pheno2 + env.eff2)),col="blue",xlab="Genotype",ylab="Phenotype")
  dev.off()
}

