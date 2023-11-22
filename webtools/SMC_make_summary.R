library(MASS)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(plotly)
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

validation_param <- function(paramfile,figurefile){
    par=read.csv(paramfile)
    h_Amax = ggplot(par)+ geom_histogram(aes(x=Amax)) + ylab("")
    h_c0 = ggplot(par)+ geom_histogram(aes(x=c0)) + ylab("")
    h_r1 = ggplot(par)+ geom_histogram(aes(x=r1)) + ylab("")
    h_sigma2 = ggplot(par)+ geom_histogram(aes(x=sigma2)) + ylab("")
    h_decay = ggplot(par)+geom_histogram(aes(x=decay_time)) + ylab("")
    h_rise = ggplot(par)+geom_histogram(aes(x=rise_time)) + ylab("")

    G=plot_grid(h_Amax,h_c0,h_sigma2,h_r1,h_decay,h_rise)
    ggsave(figurefile,G)
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

setwd(wd)
ssplit        = strsplit(wd,"/")
prefix        = ssplit[[1]][length(ssplit[[1]])]
data          = read.table(datafilename)
pgas_samples  = paste0("traj_samples_",tag,".dat")
paramfile     = paste0("param_samples_",tag,".dat")

# set stimulations as spiketimes
if (file.exists(stimfile)){
    spikefile=stimfile
} else {
    spikefile=NULL
}

make_df(data[,c(1,as.numeric(slice)+2)],
        pgas_samples,
        paramfile,
        spikefile=spikefile,
        jsonfile=paste0("summary_",tag,".json"),
        keep=as.numeric(args[6])) 


