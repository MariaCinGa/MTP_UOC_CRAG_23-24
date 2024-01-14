setwd("/Users/sramos/Desktop/TFM-UOC-Proposals-2023/NEW_simulations_EHHS-GWAS_2023/Step1-SLiM_plus_Stats/")

scenarios <- c(1:20)
scenarios_sweep <- c(6:10,16:20)
files <- sprintf("slim.output_s%.0f_it0.p1p2.txt.geno.txt_Results_Tang.txt",scenarios)
file_sweep <- sprintf("slim.output_s%.0f_it0.p2_sweeps.txt",scenarios);
file_sel.effects <- sprintf("slim.output_s%.0f_it0.p1p2_table.txt",scenarios)  

pdf("Rsb_Results.pdf")
par(mfrow=c(1,1))
scen <- 5
for(scen in scenarios) {
  ff <- read.table(file=files[scen],header=T)
  sel.pos <- read.table(file=file_sel.effects[scen],skip=27,header=T)
  #eliminate those columns that are not coincident
  snps <- intersect(ff[,1],sel.pos[,1])
  j <- 1;k<-1
  pos.no.snps <- NULL
  while(j<=length(sel.pos[,1]) && k<=length(snps)) {
    if(sel.pos[j,1]==snps[k]) {j<-j+1;k<-k+1;}
    else {pos.no.snps  <- c(pos.no.snps,j); j<-j+1;}
  }
  sel.pos<-sel.pos[-c(pos.no.snps),]
  
  #Plot Rsb 
  #plot(x=ff$Position,y=ff$log_Rsb.POP1.POP2.,type="p",xlab="Position",ylab="Rsb",main=sprintf("Rsb for S%.0f",scen),col="blue")
  #Plot RsbN
  plot(x=ff$Position,y=ff$log_RsbN.POP1.POP2.,type="p",pch=20,xlab="Position",ylab="RsbN",main=sprintf("RsbN for S%.0f",scen))
  abline(h=c(-4,0,4));text(x=1e6,y=4.25,"SIG.POP1");text(x=1e6,y=-4.25,"SIG.POP2")
  #plot located selective sweeps, (if they are)
  if(sum(scen==scenarios_sweep)>0) {
    fs <- read.csv(file=file_sweep[scen],header=F,sep=" ")
    abline(v=fs[,4],col="red")
  }
  #plot possible select sweeps from pop diffs
  # diff.sel <- sel.pos$SCOEF2-sel.pos$SCOEF1
  # ss <- sel.pos[diff.sel<(-0.1) & (sel.pos$FREQP2-sel.pos$FREQP1)>0.5,]
  # abline(v=ss$POS,col="magenta")
}
dev.off()
