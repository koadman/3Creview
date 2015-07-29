#!/usr/bin/env Rscript
# script to plot cumulative insertion distances in a few datasets
# not intended for general purpose use
#

bin_size <- 50000
xmax <- 2000000
names <- c("gi|255767013|ref|NC_000964.3|"="B. subtilis 168", "gi|556503834|ref|NC_000913.3|"="E. coli K12 MG1655",
"gi|110645304|ref|NC_002516.2|"="P. aeruginosa PAO1", "gi|685631213|gb|CP009361.1|"="S. aureus ATCC 25923", "gi|685633859|gb|CP009362.1|"="S. aureus ATCC 25923 plasmid")
lengths <- c("gi|110645304|ref|NC_002516.2|"=6264404,"gi|685631213|gb|CP009361.1|"=2778854,"gi|685633859|gb|CP009362.1|"=27491,"gi|556503834|ref|NC_000913.3|"=4641652,"gi|255767013|ref|NC_000964.3|"=4215606)

burton_names <- c("gi|255767013|ref|NC_000964.3|"="Bacillus subtilis subsp. subtilis str. 168",
"gi|70728250|ref|NC_004129.6|"="Pseudomonas protegens Pf-5",
"gi|39933080|ref|NC_005296.1|"="Rhodopseudomonas palustris CGA009",
"gi|39840937|ref|NC_005297.1|"="Rhodopseudomonas palustris CGA009 plasmid",
"gi|45357563|ref|NC_005791.1|"="Methanococcus maripaludis S2",
"gi|172087630|ref|NC_006840.2|"="Vibrio fischeri ES114 chromosome I",
"gi|172087787|ref|NC_006841.2|"="Vibrio fischeri ES114 chromosome II",
"gi|59714356|ref|NC_006842.1|"="Vibrio fischeri ES114 plasmid pES100",
"gi|83716035|ref|NC_007650.1|"="Burkholderia thailandensis E264 chromosome II",
"gi|83718394|ref|NC_007651.1|"="Burkholderia thailandensis E264 chromosome I",
"gi|146297766|ref|NC_009441.1|"="Flavobacterium johnsoniae UW101")

burton_lengths <- c("gi|255767013|ref|NC_000964.3|"=4215606,
"gi|70728250|ref|NC_004129.6|"=7074893,
"gi|39933080|ref|NC_005296.1|"=5459213,
"gi|39840937|ref|NC_005297.1|"=8427,
"gi|45357563|ref|NC_005791.1|"=1661137,
"gi|172087630|ref|NC_006840.2|"=2897536,
"gi|172087787|ref|NC_006841.2|"=1330333,
"gi|59714356|ref|NC_006842.1|"=45849,
"gi|83716035|ref|NC_007650.1|"=2914771,
"gi|83718394|ref|NC_007651.1|"=3809201,
"gi|146297766|ref|NC_009441.1|"=6096872)

#names <- burton_names
#lengths <- burton_lengths


insd <- read.table(commandArgs(TRUE)[1])

pdf(commandArgs(TRUE)[2])
#plot(0,0,xlim = c(0,xmax),ylim = c(5,15),type = "n",xlab="Insert size (nt)", ylab="log number of read pairs (100kbp bins)")
plot(0,0,xlim = c(0,xmax),ylim = c(0,1),type = "n",xlab="Insert size (nt)", ylab="quantile", main="Cumulative distribution of read pairs > 1000nt")
exps <- unique(as.character(insd$V1))
ltype <- vector()
#exps <- c("3")
legnames <- vector()
for(j in 1:length(exps)){
	chromos <- unique(as.character(insd$V2[insd$V1 == exps[j]]))
	bb = seq(from=0,to=10000000,by=bin_size)
	ltype <- c(ltype,rep(j, length(names)))
	for (i in 1:length(names)){
		legnames <- c(legnames, paste(exps[j], names[i]))
		if(names(names)[i] == "gi|685633859|gb|CP009362.1|")
			next
		ins <- insd$V3[insd$V1 == exps[j] & insd$V2 == names(names)[i]]
		# fold the length distribution back over the circular chromosome
		ins[ins > lengths[names(names)[i]]/2] <- lengths[names(names)[i]] - ins[ins > lengths[names(names)[i]]/2]
		ins <- ins[ins > 1000]
		hh <- hist(ins,breaks=bb,plot=F)
#		hh$counts <- log(hh$counts)
#		hh$counts[is.infinite(hh$counts)] <- 0
		hh$counts <- hh$counts / sum(hh$counts) # normalise
		cs <- cumsum(hh$counts)
		last <- floor(lengths[names(names)[i]] / (2*bin_size))
		lines(bb[1:last] + bin_size/2,cs[1:last],col = i,lwd=2 + floor(i/8), lty=j)
	}
}

legend("bottomright",legend=legnames, col=rep(seq(1,length(names)),length(exps)),lwd=2,lty=ltype)

dev.off()
