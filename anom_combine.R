##########
# combine BAF or LOH anomaly files
# Usage: R --args config.file type end by < anom_combine.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]

# read end and by scan numbers
if (length(args) < 4) stop("missing type, end, by")
type <- args[2] # "BAF" or "LOH"
END <- as.integer(args[3])         #total number of samples
by <- as.integer(args[4])         # number of samples in a block of samples run on cluster


stt<-floor(END/by)*by
st<-seq(1,stt,by)
end<-st+by-1
if (end[length(end)] < END) {
  st<-c(st,stt+1)
  end<-c(end,END)
}

tmp <- list()
for(i in 1:length(st)){
  file <- file.path(config["out_anom_dir"], paste(config["project"], type, paste("filtered",st[i],end[i],sep="_"), "RData",sep="."))
  tmp[[i]]<-getobj(file)
}
filtered <- do.call("rbind", tmp)

saveas(filtered, paste(config["project"], type, "filtered.all", sep="."), config["out_anom_dir"])
