library(rjson)
library(data.table)

make_df <- function(trace,pgas_samples,p_file,spikefile_init=NULL,spikefile=NULL,jsonfile=NULL,keep=20){

    MCMC.samples = fread(pgas_samples)
    nsamples = as.numeric(MCMC.samples[,.(max(index))])+1

    par = as.matrix(read.csv(p_file))
    
    MCMC.samples[,time:=trace[,1],by=index]
    MCMC.samples[,f:=B+C]
    MCMC.summary = MCMC.samples[index>nsamples-keep+1, .(prob=mean(S), 
                                                         B=mean(B),
                                                         burst=mean(burst),
                                                         filter_mean=mean(f),
                                                         filter_sd=sd(f)
                                                         ), by=time]

    time = trace[,1]
    data=trace[,2]

    gtbase=c()

    if(!is.null(spikefile)){
        gtdata = read.table(spikefile)
        spiketimes=gtdata$V1/1000
    } else {
        print("no spikefile")
        spiketimes=c()
    }

    if(!is.null(spikefile_init)){
        spiketimes_init=read.table(spikefile_init)$V1
    } else {
        spiketimes_init=c()
    }


    L=list(time=MCMC.summary[,time],
           prob=MCMC.summary[,prob],
           B=MCMC.summary[,B],
           burst=MCMC.summary[,burst],
           gtbase=gtbase,
           Y=data,
           filter_mean=MCMC.summary[,filter_mean],
           filter_sd=MCMC.summary[,filter_sd],
           spiketimes=spiketimes,
           spiketimes_init=spiketimes_init
    )

    jsondata=toJSON(L)
    if(is.null(jsonfile)) jsonfile="summary.json"
    
    write(jsondata,jsonfile)
}

args = commandArgs(trailingOnly=TRUE)
#  wd args[1]
#  spikefile
#  label
#  keep

wd=args[1]
datafilename=args[2]
slice = as.numeric(args[3])
tag = args[4]
stimfile = args[5]

data          = read.table(datafilename)
pgas_samples  = paste0(wd,"/traj_samples_",tag,".dat")
paramfile     = paste0(wd,"/param_samples_",tag,".dat")

# set stimulations as spiketimes
if (file.exists(stimfile)){
    spikefile=stimfile
} else {
    spikefile=NULL
}

make_df(data[,c(1,as.numeric(slice)+1)],
        pgas_samples,
        paramfile,
        spikefile=spikefile,
        jsonfile=paste0(wd,"/summary_",tag,".json"),
        keep=as.numeric(args[6])) 


