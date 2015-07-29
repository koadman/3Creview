#!/usr/bin/env Rscript
# script to plot rates of proximity ligation in a few datasets
# not written to support application to arbitrary datasets
#
cutsite_densities <- c(
"gi|110645304|ref|NC_002516.2|"=0.0113629325312,
"gi|685631213|gb|CP009361.1|"=0.000480054007875,
"gi|685633859|gb|CP009362.1|"=0.000145502164345,
"gi|556503834|ref|NC_000913.3|"=0.0052375748979,
"gi|255767013|ref|NC_000964.3|"=0.00344647958087)

names <- c("gi|255767013|ref|NC_000964.3|"="B. subtilis 168", "gi|556503834|ref|NC_000913.3|"="E. coli K12 MG1655",
"gi|110645304|ref|NC_002516.2|"="P. aeruginosa PAO1", "gi|685631213|gb|CP009361.1|"="S. aureus ATCC 25923", "gi|685633859|gb|CP009362.1|"="S. aureus ATCC 25923 plasmid")


plotseries <- function(dat,ylabel){
	plot(0,0,xlim = c(1.5,5.5),ylim = c(0,max(dat)*1.3),type = "n",xlab="Formalin concentration (%)", ylab=ylabel)

	chromos <- unique(as.character(plig$V2))
	for (i in 1:length(chromos)){
	    lines(plig$V1[plig$V2==chromos[i]],dat[plig$V2==chromos[i]],col = i,type = 'b',lwd=2)
	}

	legend("topright",legend=names[chromos],col=seq(1,length(chromos)),lwd=3)

}

plig <- read.table(commandArgs(TRUE)[1])
chromos <- unique(as.character(plig$V2))

# normalize by cut site density
#for (i in 1:length(chromos)){
#	plig$V3[plig$V2 == chromos[i]] <- plig$V3[plig$V2 == chromos[i]] / (cutsite_densities[chromos[i]]^2)
#}

pdf(commandArgs(TRUE)[2])
plotseries(plig$V3,"Raw proximity ligation read rate")
dev.off()

pdf("abundances.pdf")
plotseries(plig$V4,"Fraction of reads mapping to organism")
dev.off()



