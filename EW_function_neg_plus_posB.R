############################################################
#Eyre-Walker function PNAS 2010 (function relationship s,a)
# Relation fitness versus a QTL.
# include positive s values correlated with a QTL.
# include skew par to mimmic directional selection
# Transport selective values to a second scenario for a fixed z (phenotype)
# Calculate new mutations for a second scenario (SEPARATELY)
# Calculate h for each s
############################################################

z_function <- function(delta,Ne,s,eps,tau,tau2,max_s) {
  z <- delta * (4*Ne*s)^(tau) * (1 + eps) / (4*Ne*max_s)^(tau2) 
  return(z)
}
s.function <- function(z,tau,tau2,eps,Ne,disp,max_s) { #do for positive and negative separately and weight for skew
  s <- (4*Ne*max_s)^(tau2/tau) / (4*Ne) * (abs(z-disp)/(1+eps))^(1/tau)
  return(s)
}
#normalize.Zdist <- function(z,disp) {
#  mx <- max(z,na.rm=T); 
#  mn <- min(z,na.rm=T)
#  #a2 <- (max(a,na.rm=T)-min(a,na.rm=T))/(max(ai2,na.rm=T)-min(ai2,na.rm=T)) * (ai2-min(ai2,na.rm=T)) + min(a,na.rm=T)
#  z.new <- z
#  z.new[which(z>=disp)] <- (mx-0)/(mx-disp) *  (z[which(z>=disp)]-disp) + 0;
#  z.new[which(z< disp)] <- (0-mn)/(disp-mn) *  (z[which(z< disp)]-mn  ) + mn;
#  z.new
#}

#R --vanilla --args 250000 25000 0.01 0.01 1000 0.2 -120 1 50 50 0.36 0.00 0.00 0.60 muts_EW_tau00s1_mutations
#R --vanilla --args 250000 25000 0.01 0.10 1000 0.2 -120 1 50 5  0.36 0.00 0.00 0.60 muts_EW_tau00s2_mutations
#R --vanilla --args 250000 25000 0.01 0.01 1000 0.2 -120 1 50 50 0.36 0.25 0.10 0.50 muts_EW_tau00s3_mutations
#R --vanilla --args 250000 25000 0.01 0.10 1000 0.2 -120 1 50 5  0.36 0.25 0.10 0.50 muts_EW_tau00s4_mutations
#R --vanilla --args 250000 25000 0.01 0.01 1000 0.2 -120 1 50 50 0.36 0.10 0.25 0.40 muts_EW_tau00s5_mutations
#R --vanilla --args 250000 25000 0.01 0.10 1000 0.2 -120 1 50 5  0.36 0.10 0.25 0.40 muts_EW_tau00s6_mutations
#R --vanilla --args 250000 25000 0.01 0.01 1000 0.2 -120 1 50 50 0.36 0.05 0.50 0.30 muts_EW_tau00s7_mutations
#R --vanilla --args 250000 25000 0.01 0.10 1000 0.2 -120 1 50 5  0.36 0.05 0.50 0.30 muts_EW_tau00s8_mutations
#R --vanilla --args 250000 25000 0.01 0.01 1000 0.2 -120 1 50 50 0.36 0.02 0.80 0.10 muts_EW_tau00s9_mutations
#R --vanilla --args 250000 25000 0.01 0.10 1000 0.2 -120 1 50 5  0.36 0.02 0.80 0.10 muts_EW_tau00s10_mutations

#R --vanilla --args " + 
# totalm2_muts_p1 + " " + 
# totalm2_muts_p2 + " " + 
# prop_beneficial_p1 + " " + 
# prop_beneficial_p2 + " " + 
# Ne_pop + " " + 
# shape_deleterious + " " +  
# s_mean_deleterious + " " +  
# shape_beneficial + " " + 
# s_mean_beneficial_p1 + " " + 
# s_mean_beneficial_p2 + " " + 
# h_mean + " " + 
# disp + " " + 
# tau + " " + 
# sigma + " " + 
# filePathData + file_output1 + "_mutations" + " " + 
# refn + " " + "< " + 
# file_input1
args = commandArgs(trailingOnly=TRUE)

if(length(args)<16) {
  args <- array(dim=c(14))
  args[1] <- 2.5e5
  args[2] <- 2.5e4
  args[3] <- 0.01
  args[4] <- 0.1
  args[5] <- 1000
  args[6] <- 0.2
  args[7] <- -0.03
  args[8] <- 1
  args[9] <- 0.0125
  args[10] <- 0.00125
  args[11] <- 0.36
  args[12] <- 0.20
  args[13] <- 0.50
  args[14] <- 0.30
  args[15] <- "EW_function_neg_plus_pos"
  args[16] <- 1
}

#mutations
n.mutations <- as.numeric(args[1])# + as.numeric(args[2]) #1e4
n.mutations.p2 <- as.numeric(args[2]) #1e4
prop.positive <- as.numeric(args[3]) #0.01
prop.positive.p2 <- as.numeric(args[4]) #0.1
#number and proportion of negative and positive mutations and Ne
nit_neg <- round(n.mutations * (1 - prop.positive),0)
nit_pos <- round(n.mutations * (    prop.positive),0)
Ne <- as.numeric(args[5]) #1000
#the fitness has two distributions:  positive (exponential) and  negative (gamma).
shape_s_neg <- as.numeric(args[6]) #0.2
mean_s_neg  <- -as.numeric(args[7]) #-4000
shape_s_pos <- as.numeric(args[8]) #1
mean_s_pos1  <- as.numeric(args[9]) #100 #exponential
mean_s_pos2  <- as.numeric(args[10]) #100 #exponential
#mean dominance
h_mean <- as.numeric(args[11]) #0.36
#environmental change displacement
disp_values <- as.numeric(args[12]) #0.1
#parameters for the Eyre-Walker function for the phenotype:
tau_values   <- as.numeric(args[13]) #c(0.000,0.25,0.50,0.80)
sigma_values <- as.numeric(args[14]) #c(0.600,0.40,0.30,0.10)
output.file <- args[15]
counter.sc <- as.numeric(args[16])
show(sprintf("R --vanilla --args %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s < EW_function_neg_plus_posB.R",args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],args[10],args[11],args[12],args[13],args[14],args[15],args[16] ));

#test parameters for the Eyre-Waalker function:
#tau_values   <- c(0.00,0.10,0.25,0.50,0.80)
#sigma_values <- c(0.60,0.50,0.40,0.30,0.10)
#disp_values  <- c(0.00,0.10,0.05,0.02,0.01) #corr 0 means no displacement!

pdf(sprintf("%s.pdf",output.file),width=15,height=5)
#par(mfrow=c(2,2))
par(mfrow=c(1,3))
ts <- 1
#for(ts in c(1:length(tau_values))) {
  ####  MUTATIONS IN FIRST SCENARIO AND STANDING AND NEW MUTATIONS IN SECOND: #####
  ####  MUTATIONS IN FIRST SCENARIO AND STANDING AND NEW MUTATIONS IN SECOND: #####
  
  diff.pos.2vs1 <- mean_s_pos2/mean_s_pos1
  #obtain s values from distributions:
  s_neg <- -rgamma(nit_neg,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg) 
  #s_neg <- -rgamma(nit_neg,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg)*1e-9
  s_pos <- +rgamma(nit_pos,shape=shape_s_pos,rate=shape_s_pos/mean_s_pos1) 
  #distribution of positive and negative values in phenotype
  delta <- sample(x=c(-1,1),size=nit_neg+nit_pos,replace=T,prob=c(1,1))
  sigma <- sigma_values[ts]
  tau <- tau_values[ts]
  eps <- rnorm(nit_neg+nit_pos,mean=0,sd=sigma)
  #for negative values
  zn <- c(z_function(delta=delta[1:nit_neg],Ne=Ne,s=-s_neg[1:nit_neg],eps=eps[1:nit_neg],tau=tau,tau2=tau,max_s=max(abs(s_neg))))
  zn.no.eps <- c(z_function(delta=delta[1:nit_neg],Ne=Ne,s=-s_neg[1:nit_neg],eps=0,tau=tau,tau2=tau,max_s=max(abs(s_neg))))
  #for positive values
  zp <- c(z_function(delta=delta[(nit_neg+1):(nit_neg+nit_pos)],Ne=Ne,s=s_pos[1:nit_pos],eps=eps[(nit_neg+1):(nit_neg+nit_pos)],tau=-tau,tau2=tau,max_s=max(1/(s_pos))))
  zp.no.eps <- c(z_function(delta=delta[(nit_neg+1):(nit_neg+nit_pos)],Ne=Ne,s=s_pos[1:nit_pos]*diff.pos.2vs1,eps=0,tau=-tau,tau2=tau,max_s=max(1/(s_pos))))
  #z.no.eps <- c(zn.no.eps,zp.no.eps)
  s <- c(s_neg,s_pos)
  z <- c(zn,zp)
  #disp estimated from proportion of distribution
  disp.norm <- (max(z)-0) * disp_values[ts]
  #disp.norm <- disp_values[ts]
  #normalization for next environment
  zn.no.eps.disp <- zn.no.eps - disp.norm #c(normalize.Zdist(z=zn.no.eps,disp=disp.norm)) #to have same limits from initial
  zp.no.eps.disp <- zp.no.eps - disp.norm #c(normalize.Zdist(z=zp.no.eps,disp=disp.norm)) #to have same limits from initial
  
  #STANDING VARIANTS GOING TO SECOND SCENARIO
  #calculate the negative s values for a new scenario (environmental change) Assume first that all (+/-) are negative
  if(tau==0) {
    s2 <- s
  } else {
    s2 <- c(-s.function(z=zn.no.eps.disp[1:nit_neg],tau=tau,tau2=tau,eps=0,Ne=Ne,disp=0,max_s=max(abs(s_neg))))
    s2 <- c(s2,-s.function(z=zp.no.eps.disp[1:nit_pos],tau=tau,tau2=tau,eps=0,Ne=Ne,disp=0,max_s=max(abs(s_neg))))
  }
  #A proportion of s2 values are being positive (n.mutations*prop.positive.p2). 
  #Subtract from negative values. 
  #Find the z values that are closer to zp.no.eps+disp and change the s2 value to positive s
  #obtain the density values in several arbitrary sections for z
  ndiv <- 5
  qq2.p <- density(x=zp+disp.norm,n=ndiv)
  #mean.zp.no.eps <- disp.norm
  #sd.zp.no.eps <- sd(zp.no.eps)
  #qq2.p <- density(x=rnorm(n=1000,mean=mean.zp.no.eps,sd=sd.zp.no.eps),n=ndiv)# use a normal distrib. instead
  sum.dens.y <- sum(qq2.p$y)
  #locate the positions in where z is inside the possible positive values
  # Divide in segments and change it for the z value + disp
  n.pos.thr <- which(z >= min(zp+disp.norm) & z <= max(zp+disp.norm))
  count.pos <- array(0,dim=c(ndiv-1))
  for(i in n.pos.thr) {
    thr <- min(which(qq2.p$x > z[i])) - 1
    if(s[i]>0) {
      s2[i] <- s[i]*diff.pos.2vs1
      count.pos[thr] <- count.pos[thr] + 1
    } else {
      if(count.pos[thr] < (prop.positive.p2 * n.mutations * qq2.p$y[thr+1]/sum.dens.y)) {
          s2[i] <- s_pos[sample(c(1:(prop.positive.p2 * length(s_pos))),size=1)]*diff.pos.2vs1
          count.pos[thr] <- count.pos[thr] + 1
      }
    }
  }
  
  #DOMINANCE: from Wang, Caballero, Keightley and Hill
  alpha_n <- shape_s_neg / (mean_s_neg)
  alpha_p <- shape_s_pos / (mean_s_pos1)
  Kn = alpha_n * ((2.0*h_mean)^((-1.0)/shape_s_neg)-1.0);
  Kp = alpha_p * ((2.0*h_mean)^((-1.0)/shape_s_pos)-1.0);
  hn <- runif(n=length(s_neg),0,exp(+s_neg*Kn))
  hp <- runif(n=length(s_pos),0,exp(-s_pos*Kp))
  h <- c(hn,hp)
  
  #threshold to s=-1    
  s[s < -1] <- -1
  s2[s2 < -1] <- -1
  
  #comment <- function() {
    #### NEW MUTATIONS IN SECOND SCENARIO: #####
    #### NEW MUTATIONS IN SECOND SCENARIO: #####
    #mutations for the second scenario: same than first, inlcluding new prop.positive.new.env and displ
    #number and proportion of negative and positive mutations for scenario 2
    nit_neg.p2 <- round(n.mutations.p2 * (1 - prop.positive.p2),0)
    nit_pos.p2 <- round(n.mutations.p2 * (    prop.positive.p2),0)
    #obtain s values from distributions:
    s_neg.p2 <- -rgamma(nit_neg.p2,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg)
    #s_neg.p2 <- -rgamma(nit_neg.p2,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg)*1e-9
    s_pos.p2 <- +rgamma(nit_pos.p2,shape=shape_s_pos,rate=shape_s_pos/mean_s_pos2)
    #distribution of positive and negative values in phenotype
    #normalize by same extremes and calculate delta according displacement, and the apply disp 
    delta.pos.bias <- 1 - disp.norm/(max(z)-0);
    #delta.pos.bias <- 1
    delta.p2 <- sample(x=c(-1-(1-delta.pos.bias),delta.pos.bias),size=nit_neg.p2+nit_pos.p2,replace=T,prob=c(1+(1-delta.pos.bias),delta.pos.bias))
    sigma.p2 <- sigma_values[ts]
    tau.p2 <- tau_values[ts]
    eps.p2 <- rnorm(nit_neg.p2+nit_pos.p2,mean=0,sd=sigma.p2)
    #for negative values
    zn.p2 <- c(z_function(delta=delta.p2[1:nit_neg.p2],Ne=Ne,s=-s_neg.p2[1:nit_neg.p2],eps=eps.p2[1:nit_neg.p2],tau=tau,tau2=tau,max_s=max(abs(s_neg.p2))) + disp.norm)
    #for positive values
    zp.p2 <- c(z_function(delta=delta.p2[(nit_neg.p2+1):(nit_neg.p2+nit_pos.p2)],Ne=Ne,s=s_pos.p2[1:nit_pos.p2],eps=eps.p2[(nit_neg.p2+1):(nit_neg.p2+nit_pos.p2)],tau=-tau,tau2=tau,max_s=max(1/(s_pos.p2))) +disp.norm)
    #s_neg.p2[s_neg.p2 < -1] <- -1
    s.p2 <- c(s_neg.p2,s_pos.p2)
    z.p2 <- c(zn.p2,zp.p2)
    
    #DOMINANCE: from Wang, Caballero, Keightley and Hill
    alpha_n <- shape_s_neg / (mean_s_neg)
    alpha_p <- shape_s_pos / (mean_s_pos2)
    Kn.p2 = alpha_n * ((2.0*h_mean)^((-1.0)/shape_s_neg)-1.0);
    Kp.p2 = alpha_p * ((2.0*h_mean)^((-1.0)/shape_s_pos)-1.0);
    hn.p2 <- runif(n=length(s_neg.p2),0,exp(+s_neg.p2*Kn.p2))
    hp.p2 <- runif(n=length(s_pos.p2),0,exp(-s_pos.p2*Kp.p2))
    h.p2 <- c(hn.p2,hp.p2)
    
    #threshold to s=-1    
    s.p2[s.p2 < -1] <- -1
  #}
  #In case of NAs, eliminate all
  wna <- which(!is.na(s) & !is.na(s2) & !is.na(z) & !is.na(h))
  wna2 <- which(!is.na(s.p2) & !is.na(z.p2) & !is.na(h.p2))
  
  #### PLOTS #####
  r1 <- sample(x=wna,size=1e4,replace=T)
  r2 <- sample(x=wna2,size=1e4,replace=T)

  #plots for scenario 1 
  plot(s[r1],z[r1],main=sprintf("SC: %.0f tau=%.3f sd=%.3f disp=%.2f\ns1=%.4f prop.s1=%.2f \n ENV1 corr(|s|,|z|)=%.3f",counter.sc,tau,sigma,disp_values[ts],mean_s_pos1,prop.positive,cor(abs(s[r1]),abs(z[r1]),method="spearman")),xlim=c(-1,0.1),xlab="s",ylab="Phenotype",ylim=c(-3,3))#,ylim=c(min(z[r1],z.p2[r2]),max(z[r1],z.p2[r2])))
  abline(v = 0, h = 0,col="blacK")
  ##lines(s2,z.no.eps,col="red",pch=20,type="p")
  ##plot(s[r],h[r],pch=20,main=sprintf("Relation s versus dominance.\n kn: %.3e kp: %.3e",Kn,Kp),xlab="Sel.Coeff",ylab="dominance",xlim=c(-1,0.1),ylim=c(0,1)) 
  
  #plots for scenario 1 and standing mutations
#  plot(s[r1],z[r1],main=sprintf("SC: %.0f tau=%.3f sd=%.3f disp=%.2f\ns1=%.4f prop.s1=%.2f s2=%.4f prop.s2=%.3f\n ENV1 (black) Standing on ENV2 (purple) corr(|s|,|z|)=%.3f",counter.sc,tau,sigma,disp_values[ts],mean_s_pos1,prop.positive,mean_s_pos2,prop.positive.p2,cor(abs(s[r1]),abs(z[r1]),method="spearman")),xlim=c(-1,0.1),xlab="s",ylab="Phenotype",ylim=c(-3,3))#,ylim=c(min(z[r1],z.p2[r2]),max(z[r1],z.p2[r2])))
#  abline(v = 0, h = 0,col="blacK")
#  lines(s2[r1],z[r1],type="p",col="purple",pch=20); abline(h = disp.norm,col="blue")
  ##lines(s2,z.no.eps,col="red",pch=20,type="p")
  ##plot(s[r],h[r],pch=20,main=sprintf("Relation s versus dominance.\n kn: %.3e kp: %.3e",Kn,Kp),xlab="Sel.Coeff",ylab="dominance",xlim=c(-1,0.1),ylim=c(0,1)) 
  
  #plots for scenario 2 of new and standing mutations
  #plot(s2[r1],z[r1],col="purple",main=sprintf("SC: %.0f tau=%.3f sd=%.3f disp=%.04f \n s2=%.4f prop.s2=%.3f \n Standing from ENV1 on ENV2 (purple). New on ENV2 (blue) corr(|s.p2|,|z.p2|)=%.3f",counter.sc,tau,sigma,disp_values[ts],mean_s_pos2,prop.positive.p2,cor(abs(s.p2[r2]),abs(z.p2[r2]-disp.norm),method="spearman")),xlim=c(-1,0.1),ylim=c(-3,3),xlab="s.p2",ylab="Phenotype")
  #lines(s.p2[r2],z.p2[r2],type="p",col="blue",pch=20); abline(h = disp.norm,col="blue")
#  plot(s.p2[r2],z.p2[r2],col="blue",main=sprintf("SC: %.0f tau=%.3f sd=%.3f disp=%.04f \n s2=%.4f prop.s2=%.3f \n Standing from on ENV2 (purple). New on ENV2 (blue) \n corr(|s.p2|,|z.p2|)=%.3f",counter.sc,tau,sigma,disp_values[ts],mean_s_pos2,prop.positive.p2,cor(abs(s.p2[r2]),abs(z.p2[r2]-disp.norm),method="spearman")),xlim=c(-1,0.1),ylim=c(-3,3),xlab="s.p2",ylab="Phenotype")
#  lines(s2[r1],z[r1],type="p",col="purple",pch=20); abline(h = disp.norm,col="blue")
#  abline(v = 0, h = 0,col="blacK"); abline(h = disp.norm,col="blue")
  ## plot(s.p2[r],h.p2[r],pch=20,main=sprintf("Relation s versus dominance.\n kn.p2: %.3e kp.p2: %.3e",Kn.p2,Kp.p2),xlab="Sel.Coeff",ylab="dominance",xlim=c(-1,0.1),ylim=c(0,1))
  
  #plots for standing mutations
  plot(s2[r1],z[r1],type="p",col="purple",pch=20,main=sprintf("SC: %.0f tau=%.3f sd=%.3f disp=%.04f \n s2=%.4f prop.s2=%.3f \n Standing ENV1 on ENV2. corr(|s2|,|z|)=%.3f",counter.sc,tau,sigma,disp_values[ts],mean_s_pos2,prop.positive.p2,cor(abs(s2[r1]),abs(z[r1]-disp.norm),method="spearman")),xlim=c(-1,0.1),ylim=c(-3,3),xlab="s2",ylab="Phenotype")
  abline(v = 0, h = 0,col="blacK"); abline(h = disp.norm,col="blue")
  
  #plots for scenario 2 of new mutations
  plot(s.p2[r2],z.p2[r2],col="blue",main=sprintf("SC: %.0f tau=%.3f sd=%.3f disp=%.04f \n s2=%.4f prop.s2=%.3f \n ENV2 corr(|s.p2|,|z.p2|)=%.3f",counter.sc,tau,sigma,disp_values[ts],mean_s_pos2,prop.positive.p2,cor(abs(s.p2[r2]),abs(z.p2[r2]-disp.norm),method="spearman")),xlim=c(-1,0.1),ylim=c(-3,3),xlab="s.p2",ylab="Phenotype")
  abline(v = 0, h = 0,col="blacK"); abline(h = disp.norm,col="blue")
  ## plot(s.p2[r],h.p2[r],pch=20,main=sprintf("Relation s versus dominance.\n kn.p2: %.3e kp.p2: %.3e",Kn.p2,Kp.p2),xlab="Sel.Coeff",ylab="dominance",xlim=c(-1,0.1),ylim=c(0,1))
  
  #WRITE TABLES
  write.table(x=s[wna],file=sprintf("%s_SelCoeffs.txt",output.file),quote=F,sep="\t",row.names=F,col.names=F,eol="\n")
  write.table(x=s2[wna],file=sprintf("%s_SelCoeffs2.txt",output.file),quote=F,sep="\t",row.names=F,col.names=F,eol="\n")
  write.table(x=z[wna],file=sprintf("%s_PhenEffs.txt",output.file),quote=F,sep="\t",row.names=F,col.names=F,eol="\n")
  write.table(x=h[wna],file=sprintf("%s_dominance.txt",output.file),quote=F,sep="\t",row.names=F,col.names=F,eol="\n")
  
  write.table(x=s.p2[wna2],file=sprintf("%s_SelCoeffs.p2.txt",output.file),quote=F,sep="\t",row.names=F,col.names=F,eol="\n")
  write.table(x=z.p2[wna2],file=sprintf("%s_PhenEffs.p2.txt",output.file),quote=F,sep="\t",row.names=F,col.names=F,eol="\n")
  write.table(x=h.p2[wna2],file=sprintf("%s_dominance.p2.txt",output.file),quote=F,sep="\t",row.names=F,col.names=F,eol="\n")
#}
dev.off()
