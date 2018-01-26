###FUNCTIONS AND PACKAGES####
#install.packages(c("reshape","plyr","ggplot2","Hmisc","car","xlsx"))
library(reshape); library(plyr); library(ggplot2) ;library(Hmisc);  library(car);  library(xlsx)
#install.packages("QuantPsy") ;library(QuantPsyc); -----> says QuantPsy not available for version 3.2.2. May be good to get source code for packages so can run offline.
#install.packages("nlme")
library(nlme)

gainscorefxn<-function(data){
  #calculate gainscores from data (to be used in ddply with different tasks)
  #not used for change detection
  #sorts data by session
  data <- data[order(data$session),] 
  #calculates avg game score for initial sessions (first and second sessions) and final game score (last and second to last sessions completed)
  initialsess<-c(1,2)
  finalsess<-c(nrow(data)-1,nrow(data)) #usually 19 and 20 but a few subjects did not completed all 20 sessions
  initial<-data[initialsess,"diff"]
  final<-data[finalsess,"diff"]
  numsess<-length(data$diff)
  initMean<-mean(initial,na.rm=TRUE)
  finalMean<-mean(final,na.rm=TRUE)
  gainscore<-finalMean-initMean
  datasumm<-data.frame(numsess=numsess,initial=initMean,final=finalMean,gainscore=gainscore)
  return(datasumm)
}
compscoresfxn<-function(tasklist,fdata){
#calculates composite scores based on data frame with a list of variables from that data frame to be used in calculation
col<-0 
tempdata<-fdata[,tasklist] 
for(task in tasklist){col<-col+1
attach(tempdata)
tempdata[,col]<-scale(fdata[,task])
detach(tempdata)}
comp<-rowMeans(tempdata,na.rm=TRUE)
return(comp)}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

corstarsmelt<-function(var1,var2){
  #creates correlation coefficient with signficance star (to be used in ddply for melted variables)
  require(Hmisc) 
  x<-cbind(var1,var2)
  R <- round(rcorr(x)$r,2)[2]
  p <- rcorr(x)$P  [2]
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  corr_withstars<-paste(R,mystars,sep="")
  return(corr_withstars)
}
bootcorr<-function(idgainscores,data,i){
  #bootstrap function that uses same samples for multiple correlation coefficients
  d<-data[i,]
  data_melt<-melt(d,id=idgainscores) 
  data_melt<-rename(data_melt,c(variable="transfermeas",value="transferscore"))
  data_melt2<-melt(data_melt,id=c("transfermeas","transferscore"))
  data_melt2<-rename(data_melt2,c(variable="trainingmeas",value="trainingscore"))
  data_mat<-ddply(data_melt2,.(trainingmeas,transfermeas),summarise,corr=cor(trainingscore,transferscore,use="pairwise.complete.obs"))
  data_corr<-data_mat$corr
  return(data_corr)}
bootcifxn<-function(data,boottest){
  
  for(var in 1:nrow(data)){
    bootci<-boot.ci(boottest,type="bca",index=var)
    cichar<-paste("[",round(bootci$bca[4],2),",",round(bootci$bca[5],2),"]",sep="")
    data$BCa[var]<-cichar
  }
  return(data)}
#regfxn gives dataframe in vertical format (rather than horizontal which texreg and stargazer does)
#-->below are commented out variations of stats reported
regfxn<-function(lmobject){boot<-Boot(lmobject,R=1000)
summboot<-summary(boot) #summary of boot object
confboot<-confint(lmobject) #conf interval from boot object
summ<-summary(lmobject) #summary obj (summary on lm object)
rsq<-round(summ$r.squared,2) #num (R squared value)                 
adjrsq<-round(summ$adj.r.squared,2) #num (Adjusted R^2 from summary of lm object)
fstat<-round(summ$fstatistic,2) #num (F stat from summary of lm object)
se_norm<-round(summ$coefficients[,2],2) #standard error based on normal distribution
tval_norm<-round(summ$coefficients[,3],2) #tvalue based on normality assump
pval_norm<-round(summ$coefficients[,4],2) #pval based on norm assumption
AIC<-round(AIC(lmobject),2) #num (AIC stat)
BIC<-round(BIC(lmobject),2) #num (BIC stat)
coeff<-round(summboot$original,2) #num (coefficient values from Boot and lmobject are the same-- use summ$coefficients[,1] for same thing)
#bootstrap stats
bootlower<-round(confboot[,1],2) #lower confidence interval (2.5% is default)
bootupper<-round(confboot[,2],2) #upper confidence interval (97.5% is default)
bootbias<-summboot[,3] #bias of bootstrap (original estimate- the average of bootstrap sample estimates)
bootSE<-summboot[,4] #standard deviation of the bootstrap replicates

coeffnames<-row.names(summ$coefficients)  #fixed effects/predictor/coefficient names
model<-as.character(lmobject$call)[2] #character (type in lm code)
#standardized beta coefficients
stdcoeff<-round(lm.beta(lmobject),2) #standardized beta coefficient from lm.beta package (need NA filler as it doesn't have intercept parameter estimate)
filler<-rep(NA,length(coeff)-1) #filler for variables with only one value              
#beta coefficient with confidence intervals (bootstrapped)
bootconf<-paste("[",bootlower,",",bootupper,"]",sep="")
betawbootconf<-paste(coeff,bootconf,sep="\n")
#model out: all stats in table format for one model 
#                  modelout<-data.frame(model=c(model,filler),coefficients=coeffnames,
#                             betaestimates=coeff, normSE=se_norm,tvalue=tval_norm,pvalue=pval_norm,
#                             stdbetaestimates=c(NA,stdcoeff),
#                             bootstrap_lowerconf=bootlower,bootstrap_upperconf=bootupper,
#                             bootBias=bootbias, bootSE=bootSE,vif=c(NA,filler),
#                             rsq=c(rsq,filler),adjrsq=c(adjrsq,filler),                         
#                            AIC=c(AIC,filler),BIC=c(BIC,filler),Fchange=c(NA,filler))

# #reduced stats
#                   modelout<-data.frame(model=c(model,filler),coefficients=coeffnames,
#                                           betaestimates=betawbootconf, stdbetaestimates=c(NA,stdcoeff),tvalue=tval_norm,pvalue=pval_norm,
#                                           rsq=c(rsq,filler),adjrsq=c(adjrsq,filler),                         
#                                        AIC=c(AIC,filler),BIC=c(BIC,filler),Fchange=c(NA,filler)) 


#reduced stats for summary in paper
redoxmodelout<-data.frame(model=model,coefficients=paste("Step ",length(coeffnames)-1,": ",
                                                         coeffnames[length(coeffnames)],sep=""),
                          betaestimates=coeff[length(coeff)],bootstrapci=bootconf[length(bootconf)], 
                          normSE=se_norm[length(se_norm)],tvalue=tval_norm[length(tval_norm)],
                          pvalue=pval_norm[length(pval_norm)],
                          stdbetaestimates=stdcoeff[length(stdcoeff)],
                          rsq=rsq,adjrsq=adjrsq,                         
                          AIC=AIC,BIC=BIC,Fchange=NA) 

modelout<-redoxmodelout
return(modelout)}

#hregfxn combines hierarchical models and also includes model comparison statistics (e.g., F change) 
#and stats used for specifically multiple predictors(e.g.,variance inflation factor)
#the dvarind loop places a row before first step in each model and writes title of dependent var predicting
hregfxn<-function(lmlist,dvarind){
  i=0
  for(lmobj in lmlist){
    i=i+1
    if(i==1){
      lmtable<-regfxn(lmobj)
      if(dvarind==1){
        firstrow<-regfxn(lmobj)[1,]
        firstrow[1,]<-as.character("NA")
        predictvar<-paste("Predicting",names(lmobj$model[1]),sep=" ")
        firstrow$coefficients<-as.character(firstrow$coefficients)
        firstrow$coefficients[1]<-predictvar
        firstrow$coefficients<-as.factor(firstrow$coefficients)
        lmtable<-rbind(firstrow,lmtable)
      }
    }else{
      lmtable2<-regfxn(lmobj)
      comp1<-anova(lmlist[[i-1]],lmlist[[i]])  #lmlist[[i]] same as lmobj (compares F stat of this model to previous nested model)
      vif1<-vif(lmobj) #variance inflation factor for each predictor--> p. 293 of Andy Fields book: average above 10=cause for concern, average above 1= substantially > than 1--> some bias
      #   lmtable2$vif<-c(NA,vif1) 
      Fchange<-comp1[2,"F"]
      Fsigchange<-comp1[2,"Pr(>F)"]
      lmtable2$Fchange[1]<-round(Fsigchange,3)
      lmtable<-rbind(lmtable,lmtable2) 
    }  
  }
  return(lmtable)}
####IMPORT COMBINED DATA (created from the PACTcontroldata_dataimport code)####
datawd<-"/Volumes/CognitiveTraining/PACT/GameData/ActiveControl/DataCombined"
setwd(datawd)
mCDdata<-read.csv(file="CDdatacombined_nodrops_manualadds.csv",header=TRUE)
mVSdata<-read.csv(file="VSdatacombined_nodrops.csv",header=TRUE)

mf_taskwd<-"/Volumes/CognitiveTraining/PACT/GameData/MindFrontiers"
setwd(mf_taskwd)
mfdata<-read.csv(file="MFgamedata_casted.csv",header=TRUE) 
mfdata<-plyr::rename(mfdata,c("User"="subs","Session"="session"))

#####CLEAN DATA######

#change detection
mCDdata_redox_noerrorsubs<-subset(mCDdata,subset=!(subs %in% c(NA,1102,1104,1109,1122,1138,1144))) #exclude drops and subs who had old maxduration version
mCDdata_redox<-mCDdata #only exclude drops
#mCDdata_redox<-ddply(mCDdata_redox,.(task,setsize),outlierfxn,variable="avgduration") #exclude cases 3.29 stds from mean

#visual search
mVSdata_redox<-mVSdata #exclude drops only 
#mVSdata_redox<-ddply(mVSdata_redox,.(task,session),outlierfxn,variable="avglvl")

#mind frontiers

#mfdata_avg_melt_redox<-ddply(mfdata_avg_melt,.(Game,session),outlierfxn,variable="Average_Difficulty")
mfdata_avg<-subset(mfdata,select=grep("AvgDiff|subs|session",names(mfdata)))
names(mfdata_avg)<-gsub("AvgDiff","",names(mfdata_avg))
mfdata_avg_melt<-melt(mfdata_avg,id=c("subs","session"))
mfdata_avg_melt<-plyr::rename(mfdata_avg_melt,c(variable="Game",value="Average_Difficulty"))
mfdata_avg_melt_redox<-mfdata_avg_melt

####PRACTICE EFFECTS: LINEAR MIXED MODELS#####
cdlmefxn<-function(data){
  #anova with lme function for change detection task to assess training improvement across sessions
  model = lme(avgduration~session,random=~session|subs,data = data,na.action=na.omit,control=list(maxIter=2000,opt="optim"))
  anovamodel<-anova(model)
  return(anovamodel)}
vslmefxn<-function(data){
  #anova with lme function for visual search task to assess training improvement across sessions
  model = lme(avglvl~session,random=~session|subs,data = data,na.action=na.omit,control=list(maxIter=2000,opt="optim"))
  anovamodel<-anova(model)
  return(anovamodel)}
mflmefxn<-function(data){
  #anova with lme function for each mind frontiers task to assess training improvement across sessions
  model = lme(Average_Difficulty~session,random=~session|subs,data = data,na.action=na.omit,control=list(maxIter=2000,opt="optim"))
  anovamodel<-anova(model)
  return(anovamodel)}
cdanalysis_noerrorsubs<-dlply(mCDdata_redox_noerrorsubs,.(task,setsize),cdlmefxn)
cdanalysis_WITHerrorsubs<-dlply(mCDdata_redox,.(task,setsize),cdlmefxn)
vsanalysis<-dlply(mVSdata_redox,.(task),vslmefxn)
mfanalysis<-dlply(mfdata_avg_melt_redox,.(Game),mflmefxn)

###MAKE GRAPHS OF PRACTICE EFFECTS#####
####change detection plot in paper (without error subs with old max durations)####
mCDdata_graph<-mCDdata_redox_noerrorsubs
mCDdata_graph$task<-revalue(mCDdata_graph$task,c("cars"="Cars","letters"="Letters","shapes"="Shapes"))
mCDdata_graph$setsize<-revalue(mCDdata_graph$setsize,c("setsize3"="Set Size 3","setsize5"="Set Size 5"))
mCDdata_graph<-rename(mCDdata_graph,c(task="Task"))
#CD: stat summary separately and then ggplot
mCDdata_graph.summary<-ddply(subset(mCDdata_graph), .(session,Task,setsize), summarise,
                              avgduration.mean=mean(avgduration, na.rm=TRUE),
                              avgduration.se=sd(avgduration,na.rm=TRUE)/sqrt(length(avgduration[!is.na(avgduration)])))

CDggplot<-ggplot(mCDdata_graph.summary,aes(session,avgduration.mean,shape=Task))+
  geom_point(stat="identity",aes(group=Task),size=3.5)+
  geom_line(stat="identity",aes(group=Task))+
  geom_errorbar(aes(ymax = avgduration.mean + avgduration.se, ymin = avgduration.mean - avgduration.se,group=Task),width=.15,size=.2)+
  facet_grid(setsize~.,scales="free_y")+labs(x="Session",y="Average Duration")+
  scale_shape_manual(values=c(16,17,15),name="Change Detection")+
  theme_bw()+
  theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=12),
        strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),
        legend.text=element_text(size=10,face="plain"),legend.position=c(.6,.91),axis.title=element_text(size=12),
        legend.title=element_text(size=10,face="plain"))


mCDdata_graph<-mCDdata_redox_noerrorsubs
mCDdata_graph$task<-revalue(mCDdata_graph$task,c("cars"="Cars","letters"="Letters","shapes"="Shapes"))
mCDdata_graph$setsize<-revalue(mCDdata_graph$setsize,c("setsize3"="Set Size 3","setsize5"="Set Size 5"))
mCDdata_graph<-rename(mCDdata_graph,c(task="Task"))
#CD: stat summary separately and then ggplot
mCDdata_graph.summary<-ddply(subset(mCDdata_graph), .(session,Task,setsize), summarise,
                             avgduration.mean=mean(avgduration, na.rm=TRUE),
                             avgduration.se=sd(avgduration,na.rm=TRUE)/sqrt(length(avgduration[!is.na(avgduration)])))


####change detection plot (WITH ERRORS and in supplemental)####
mCDdata_graph.witherrors<-mCDdata_redox
mCDdata_graph.witherrors$task<-revalue(mCDdata_graph.witherrors$task,c("cars"="Cars","letters"="Letters","shapes"="Shapes"))
mCDdata_graph.witherrors$setsize<-revalue(mCDdata_graph.witherrors$setsize,c("setsize3"="Set Size 3","setsize5"="Set Size 5"))
mCDdata_graph.witherrors<-rename(mCDdata_graph.witherrors,c(task="Task"))
mCDdata_graph.summary.witherrors<-ddply(subset(mCDdata_graph.witherrors), .(session,Task,setsize), summarise,
                             avgduration.mean=mean(avgduration, na.rm=TRUE),
                             avgduration.se=sd(avgduration,na.rm=TRUE)/sqrt(length(avgduration[!is.na(avgduration)])))
CDggplot.witherrors<-ggplot(mCDdata_graph.summary.witherrors,aes(session,avgduration.mean,shape=Task))+
  geom_point(stat="identity",aes(group=Task),size=2.3)+
  geom_line(stat="identity",aes(group=Task))+
  geom_errorbar(aes(ymax = avgduration.mean + avgduration.se, ymin = avgduration.mean - avgduration.se,group=Task),width=.15,size=.2)+
  facet_grid(setsize~.,scales="free_y")+labs(x="Session",y="Average Duration")+
  scale_shape_manual(values=c(16,17,15),name="Change Detection")+
  theme_bw()+
  theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=12),
        strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),
        legend.text=element_text(size=10,face="plain"),legend.position=c(.6,.91),axis.title=element_text(size=12),
        legend.title=element_text(size=10,face="plain"))

##CD: ggplot calculating stat summaries (same as above but directly from raw data)

# CDggplot<-ggplot(mCDdata_graph,aes(session,avgduration,shape=Task))+
#   stat_summary(fun.y=mean,geom="point",aes(group=Task),size=3.5)+
#   stat_summary(fun.y=mean,geom="line",aes(group=Task))+
#   stat_summary(fun.data=mean_cl_normal,geom="errorbar",aes(group=Task),width=.15,mult=1,size=.2)+#mult=1 reps 1 std error from mean (defaults to 95% confidence interval)
#   facet_grid(setsize~.,scales="free_y")+labs(x="Session",y="Average Duration")+
#   theme_bw()+
#   theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=14),
#         strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),
#         legend.text=element_text(size=10,face="bold"),legend.position=c(.6,.91),axis.title=element_text(size=12))

####visual search plot in paper####
VSdata_graph<-mVSdata_redox
#VS:rename for graphs
VSdata_graph$task<-revalue(VSdata_graph$task,c("coloredps"="Colored Ps", "ls"="Ls","original"="Original"))
VSdata_graph<-rename(VSdata_graph,c(task="Task"))
#VS: stat summary separately and then ggplot
VSdata_graph.summary<-ddply(subset(VSdata_graph), .(session,Task), summarise,
                             avglvl.mean=mean(avglvl, na.rm=TRUE), 
                            avglvl.se=sd(avglvl,na.rm=TRUE)/sqrt(length(avglvl[!is.na(avglvl)])))

VSggplot<-ggplot(VSdata_graph.summary,aes(session,avglvl.mean,shape=Task))+
  geom_point(stat="identity",aes(group=Task),size=3.5)+
  geom_line(stat="identity",aes(group=Task))+
  geom_errorbar(aes(ymax = avglvl.mean + avglvl.se, ymin = avglvl.mean - avglvl.se,group=Task),width=.15,size=.2)+
  scale_shape_manual(values=c(0,1,2),name="Visual Search")+
  labs(x="Session",y="Average Level")+
  theme_bw()+
  theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=14),
        strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),
        legend.text=element_text(size=10,face="plain"),legend.position=c(.25,.91),axis.title=element_text(size=12),
        legend.title=element_text(size=10,face="plain"))


#VS: ggplot calculating stat summaries (same as above but directly from raw data)
# VSggplot<-ggplot(VSdata_graph,aes(session,avglvl,shape=Task))+
#   stat_summary(fun.y=mean,geom="point",aes(group=Task),size=3.5)+
#   scale_shape_manual(values=c(0,1,2))+
#   stat_summary(fun.y=mean,geom="line",aes(group=Task))+
#   stat_summary(fun.data=mean_cl_normal,geom="errorbar",aes(group=Task),width=.15,mult=1,size=.2)+#mult=1 reps 1 std error from mean (defaults to 95% confidence interval)
#   labs(x="Session",y="Average Difficulty Level")+
#   theme_bw()+
#   theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=14),
#         strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),
#         legend.text=element_text(size=10,face="bold"),legend.position=c(.35,.91),axis.title=element_text(size=12))


#####mind frontiers plot in paper####
#MF: rename for graphs
mfdata_graph<-mfdata_avg_melt_redox
mfdata_graph$Game<-revalue(mfdata_graph$Game,
                                    c("SentryDuty"="Sentry Duty","SafeCracker"="Safe Cracker",
                                      "RidingShotgun"="Riding Shotgun","PenEmUp"="Pen 'Em Up",
                                      "SupplyRun"="Supply Run"))
mfdata_graph$Task<-mfdata_graph$Game
#MF: stat summary separately and then ggplot
mfdata_graph.summary<-ddply(subset(mfdata_graph), .(session,Task), summarise,
                            Average_Difficulty.mean=mean(Average_Difficulty, na.rm=TRUE),
                            Average_Difficulty.se=sd(Average_Difficulty,na.rm=TRUE)/sqrt(length(Average_Difficulty[!is.na(Average_Difficulty)])))

MFggplot<-ggplot(mfdata_graph.summary,aes(session,Average_Difficulty.mean))+
  geom_point(stat="identity",aes(shape=Task),size=3.5)+
  geom_line(stat="identity",aes(group=Task))+
  geom_errorbar(aes(ymax = Average_Difficulty.mean + Average_Difficulty.se, ymin = Average_Difficulty.mean - Average_Difficulty.se,group=Task),width=.15,size=.2)+
  scale_shape_manual(values=c(16,17,15,0,1,2),name="Mind Frontiers")+
  labs(x="Session",y="Average Level")+
  theme_bw()+
  theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=14),
        strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),
        legend.text=element_text(size=10,face="plain"),legend.position=c(.25,.845),axis.title=element_text(size=12),
        legend.title=element_text(size=10,face="plain"))



#MF: ggplot calculating stat summaries
# MFggplot<-ggplot(mfdata_graph,aes(session,Average_Difficulty))+
#   stat_summary(fun.y=mean,geom="point",aes(shape=Task),size=2.3)+
#   scale_shape_manual(values=c(16,17,15,0,1,2))+
#   stat_summary(fun.y=mean,geom="line",aes(group=Task))+
#   stat_summary(fun.data=mean_cl_normal,geom="errorbar",aes(group=Task),width=.15,mult=1,size=.2)+#mult=1 reps 1 std error from mean (defaults to 95% confidence interval)
#   labs(x="Session",y="Average Difficulty Level")+
#   theme_bw()+
#   theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=14),
#         strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),
#         legend.text=element_text(size=10,face="bold"),legend.position=c(.28,.85),axis.title=element_text(size=12))

#####figure 3 with all plots####
multiplot(MFggplot,CDggplot,VSggplot,cols=3)

####TRAINING GAIN CALC####
cdSetsize<-mCDdata_redox
cdSetsize<-rename(cdSetsize,c("avgduration"="diff"))
#ddply function below is same as gainscorefxn except initialMean-finalMean as smaller durations= better performance
#change detection
cdfxn<-function(data){
  #sorts data by session
  data <- data[order(data$session),] 
  #calculate gainscores from data 
  #calculates avg game score for initial sessions (first and second sessions) and final game score (last and second to last sessions completed)
  initialsess<-c(1,2)
  finalsess<-c(nrow(data)-1,nrow(data)) #usually 19 and 20 but a few subjects did not completed all 20 sessions
  initial<-data[initialsess,"diff"]
  final<-data[finalsess,"diff"]
  numsess<-length(data$diff)
  initMean<-mean(initial,na.rm=TRUE)
  finalMean<-mean(final,na.rm=TRUE)
  gainscore<-initMean-finalMean
  datasumm<-data.frame(numsess=numsess,initial=initMean,final=finalMean,gainscore=gainscore)
  return(datasumm)
}
gainscores_cd<-ddply(cdSetsize,.(subs,task,setsize),cdfxn)

#visual search
vsgain<-subset(mVSdata_redox,,c(subs,session,avglvl,task))
vsgain<-rename(vsgain,c("avglvl"="diff"))
gainscores_vs<-ddply(vsgain,.(subs,task),gainscorefxn)
#mind frontiers
mfgain<-mfdata_avg_melt
mfgain<-plyr::rename(mfgain,c("Average_Difficulty"="diff","Game"="task"))
gainscores_mf<-ddply(mfgain,.(subs,task),gainscorefxn)

#MIND FRONTIERS
mfgain_melt<-melt(subset(gainscores_mf,,c(subs,task,initial:gainscore)),id=c("subs","task"))
gaineachgame<-cast(mfgain_melt,subs~task+variable)
#average gain scores across all games for each sub
avg_gain_overall<-compscoresfxn(grep("_gainscore",names(gaineachgame),value=TRUE),gaineachgame)
#average initial score all games
avg_initial_overall<-compscoresfxn(grep("_initial",names(gaineachgame),value=TRUE),gaineachgame)
#average final score all games
avg_final_overall<-compscoresfxn(grep("_final",names(gaineachgame),value=TRUE),gaineachgame)
#combined average composite scores
mf_compgainscores<-data.frame(subs=gaineachgame$subs,overall_gainscore=avg_gain_overall,overall_initial=avg_initial_overall)
mf_gs_all<-join(mf_compgainscores,gaineachgame,by="subs",type="full",match="all")
#ACTIVE CONTROL
#vs
vsgain_melt<-melt(subset(gainscores_vs,,c(subs,task,initial:gainscore)),id=c("subs","task"))
vsgaineachgame<-cast(vsgain_melt,subs~task+variable)

avg_gain_vs<-data.frame(subs=vsgaineachgame$subs,avg_gain_vs=compscoresfxn(grep("_gainscore",names(vsgaineachgame),value=TRUE),vsgaineachgame))
avg_initial_vs<-data.frame(subs=vsgaineachgame$subs,avg_initial_vs=compscoresfxn(grep("_initial",names(vsgaineachgame),value=TRUE),vsgaineachgame))
#cd
cdgain_melt<-melt(subset(gainscores_cd,,c(subs,task,setsize,initial:gainscore)),id=c("subs","task","setsize"))
cdgaineachgame<-cast(cdgain_melt,subs~task+variable+setsize)

cdtaskavgscores<-data.frame(subs=cdgaineachgame$subs,
                            letters_initial=compscoresfxn(grep("letters_initial",names(cdgaineachgame),value=TRUE),cdgaineachgame),
                            shapes_initial=compscoresfxn(grep("shapes_initial",names(cdgaineachgame),value=TRUE),cdgaineachgame),
                            cars_initial=compscoresfxn(grep("cars_initial",names(cdgaineachgame),value=TRUE),cdgaineachgame),
                            letters_gainscore=compscoresfxn(grep("letters_gainscore",names(cdgaineachgame),value=TRUE),cdgaineachgame),
                            shapes_gainscore=compscoresfxn(grep("shapes_gainscore",names(cdgaineachgame),value=TRUE),cdgaineachgame),
                            cars_gainscore=compscoresfxn(grep("cars_gainscore",names(cdgaineachgame),value=TRUE),cdgaineachgame))

#avg_gain_cd<-compscoresfxn(grep("_gainscore",names(cdgaineachgame),value=TRUE),cdgaineachgame)
#avg_initial_cd<-compscoresfxn(grep("_initial",names(cdgaineachgame),value=TRUE),cdgaineachgame)

avg_gain_cd<-data.frame(subs=cdtaskavgscores$subs,avg_gain_cd=rowMeans(cdtaskavgscores[,grep("_gainscore",names(cdtaskavgscores),value=TRUE)]))
avg_initial_cd<-data.frame(subs=cdtaskavgscores$subs,avg_initial_cd=rowMeans(cdtaskavgscores[,grep("_initial",names(cdtaskavgscores),value=TRUE)]))

overall_gaindataframe<-join(avg_gain_vs,avg_gain_cd,by="subs",type="full")
overall_initialdataframe<-join(avg_initial_vs,avg_initial_cd,by="subs",type="full")

overall_gain<-rowMeans(cbind(overall_gaindataframe$avg_gain_vs,overall_gaindataframe$avg_gain_cd),na.rm=TRUE)
overall_initial<-rowMeans(cbind(overall_initialdataframe$avg_initial_vs,overall_initialdataframe$avg_initial_cd),na.rm=TRUE)

ac_gs_all<-data.frame(subs=overall_gaindataframe$subs,overall_gainscore=overall_gain,overall_initial=overall_initial,cd_gainscore=overall_gaindataframe$avg_gain_cd,
                      cd_initial=overall_initialdataframe$avg_initial_cd,vs_gainscore=overall_gaindataframe$avg_gain_vs,vs_initial=overall_initialdataframe$avg_initial_vs)


#data frame to be used for correlation analyses (1 for active control and 1 for mind frontiers)-->created above
#mf_gs_all
#ac_gs_all
#import NP composite gain scores
npcorrected<-xlsx::read.xlsx(file="/Users/michaelkranz/Documents/Analyses/PACT_behavpaper_Pauline/USED_PACTBehTaskGainScores_081315.xlsx",sheetName="GainScores_forRscript")
npdata<-npcorrected[,c("Subject","GainWMback","GainWMspan",
                       "GainGFacc","GainGFrt","GainSelectiveAtt",
                       "GainDivAtt","GainPerceptualSpeed","GainEpisMem")]
npdata<-rename(npdata,c(Subject="subs"))

ac_gs_np<-join(ac_gs_all,npdata,type="left",match="all")
ac_gs_np$group<-as.factor("ActiveControl")
mf_gs_np<-join(mf_gs_all,npdata,type="left",match="all")
mf_gs_np$group<-as.factor("MindFrontiers")

#training gain~transfer gain analyses

#correlations and table for paper
#mf_corr_redox_tooutput<-subset(mf_gs_np,select=c(subs,overall_gainscore,updating_gainscore,reason_gainscore,Wmnback:EpisodicMemory))
#ac_corr_redox_tooutput<-subset(ac_gs_np,select=c(subs,overall_gainscore,cd_gainscore,vs_gainscore,Wmnback:EpisodicMemory))

#new spreadsheet has different names
mf_corr_redox<-mf_gs_np[,grep("gainscore|Gain",value=TRUE,names(mf_gs_np))]
ac_corr_redox<-ac_gs_np[,grep("gainscore|Gain",value=TRUE,names(ac_gs_np))]

mf_corr_redox_tooutput<-mf_gs_np[,grep("subs|gainscore|Gain",value=TRUE,names(mf_gs_np))]
ac_corr_redox_tooutput<-ac_gs_np[,grep("subs|gainscore|Gain",value=TRUE,names(ac_gs_np))]
#######mind frontiers correlation table######
#mf_idgainscores=c("overall_gainscore","updating_gainscore","reason_gainscore")
mf_idgainscores<-grep("_",names(mf_corr_redox),value=TRUE) #take out reasoning and updating composite scores (3-6)
mf_measuregainscores<-grep("Gain",names(mf_corr_redox),value=TRUE) #new names
#melt for ddply (so each id can be summarized as one correlation coefficient)
mf_corr_melt<-melt(mf_corr_redox,id=mf_idgainscores,measure=mf_measuregainscores)
mf_corr_melt<-rename(mf_corr_melt,c(variable="transfermeas",value="transferscore"))
mf_corr_melt2<-melt(mf_corr_melt,id=c("transfermeas","transferscore"))
mf_corr_melt2<-rename(mf_corr_melt2,c(variable="trainingmeas",value="trainingscore"))
#get correlation coefficient with probablity value of significance stars for data frame
mf_corr_mat<-ddply(mf_corr_melt2,.(trainingmeas,transfermeas),summarise,corr=corstarsmelt(trainingscore,transferscore))
mf_corr_cast<-cast(mf_corr_mat,trainingmeas~transfermeas,value="corr")
#######active control correlation table######
ac_idgainscores<-grep("_",names(ac_corr_redox),value=TRUE)
#melt for ddply (so each id can be summarized as one correlation coefficient)
ac_corr_melt<-melt(ac_corr_redox,id=ac_idgainscores)
ac_corr_melt<-rename(ac_corr_melt,c(variable="transfermeas",value="transferscore"))
ac_corr_melt2<-melt(ac_corr_melt,id=c("transfermeas","transferscore"))
ac_corr_melt2<-rename(ac_corr_melt2,c(variable="trainingmeas",value="trainingscore"))
#get correlation coefficient with probablity value of significance stars for data frame
ac_corr_mat<-ddply(ac_corr_melt2,.(trainingmeas,transfermeas),summarise,corr=corstarsmelt(trainingscore,transferscore))
ac_corr_cast<-cast(ac_corr_mat,trainingmeas~transfermeas,value="corr")

#correlations with bootstrapped confidence intervals
#mf
mf_boottest<-boot(data=mf_corr_redox,statistic=bootcorr,R=2000,idgainscores=mf_idgainscores)
mf_corr_mat.bootstrap<-bootcifxn(mf_corr_mat,mf_boottest)
#ac
ac_boottest<-boot(data=ac_corr_redox,statistic=bootcorr,R=2000,idgainscores=ac_idgainscores)
ac_corr_mat.bootstrap<-bootcifxn(ac_corr_mat,ac_boottest)

#write.xlsx(ac_corr_mat.bootstrap,file="trainingcorrs.xlsx",sheetName="activecontrol_withbootstrap",row.names = FALSE)
#write.xlsx(mf_corr_mat.bootstrap,file="trainingcorrs.xlsx",sheetName="MF_withbootstrap",row.names = FALSE,append = TRUE)
#write.xlsx(ac_corr_cast,file="trainingcorrs.xlsx",sheetName="activecontrol_corrmatrix",row.names = FALSE,append = TRUE)
#write.xlsx(mf_corr_cast,file="trainingcorrs.xlsx",sheetName="MF_corrmatrix",row.names = FALSE,append = TRUE)



#####regressions/scatterplots######

mf_scatter<-mf_gs_np[,grep("subs|overall_gainscore|Gain|overall_initial",value=TRUE,names(mf_gs_np))]
mf_scatter$group<-as.factor("MIND FRONTIERS")
ac_scatter<-ac_gs_np[,grep("subs|overall_gainscore|Gain|overall_initial",value=TRUE,names(ac_gs_np))]
ac_scatter$group<-as.factor("ACTIVE CONTROL")

scatterall<-rbind(ac_scatter,mf_scatter)
scatmelt<-melt(scatterall,id=c("subs","overall_gainscore","overall_initial","group"))
scatmelt$variable<-revalue(scatmelt$variable,
                           c("GainWMback"="Working\nMemory\n(nback)","GainWMspan"="Working\nMemory\n(span)",
                             "GainGFacc"="Reasoning\n(accuracy)","GainGFrt"="Reasoning\n(RT)",
                             "GainSelectiveAtt"="Selective\nAttention",
                             "GainDivAtt"="Divided\nAttention",
                             "GainPerceptualSpeed"="Perceptual\nSpeed",
                             "GainEpisMem"="Episodic\nMemory"))

scatmelt$group= factor(scatmelt$group,levels(scatmelt$group)[c(2,1)])
ggplot(scatmelt,aes(overall_gainscore,value))+theme_bw()+geom_point(aes(color=group))+geom_smooth(method="lm",aes(fill=group,color=group),alpha=.1)+
  labs(x="Standardized Training Gain Score",y="Standardized Improvement Score (Transfer) ")+
  facet_grid(.~variable)+
  scale_colour_manual(values = c("MIND FRONTIERS" = "red","ACTIVE CONTROL" = "black"))+
  scale_fill_manual(values = c("MIND FRONTIERS" = "red","ACTIVE CONTROL" = "black"))+
  theme_bw()+scale_x_continuous(limits=c(-2, 2))+
  theme(text=element_text(family="Helvetica"),panel.grid=element_blank(),strip.text= element_text(size=14,face="bold"),
        strip.background = element_rect(colour="black",fill="white"),axis.text=element_text(size=12),legend.title=element_blank(),
        legend.text=element_text(size=12,face="bold"),legend.direction="horizontal",legend.position=c(.6,.9),axis.title=element_text(size=12))

ggplot(scatmelt,aes(overall_gainscore))+geom_histogram()+facet_grid(.~variable)
overallregmelt<-melt(scatterall,id=c("subs","overall_gainscore","overall_initial","group"))
gaingroupcast<-cast(data=overallregmelt,subs+overall_gainscore+overall_initial+group~variable)

og1_WMn=lm(GainWMback~group,data=gaingroupcast)
og2_WMn=update(og1_WMn,.~.+overall_gainscore)
og3_WMn=update(og2_WMn,.~overall_gainscore*group)

og1_WMspan=lm(GainWMspan~group,data=gaingroupcast)
og2_WMspan=update(og1_WMspan,.~.+overall_gainscore)
og3_WMspan=update(og2_WMspan,.~overall_gainscore*group)

og1_ReasACC=lm(GainGFacc~group,data=gaingroupcast)
og2_ReasACC=update(og1_ReasACC,.~.+overall_gainscore)
og3_ReasACC=update(og2_ReasACC,.~overall_gainscore*group)

og1_ReasRT=lm(GainGFrt~group,data=gaingroupcast)
og2_ReasRT=update(og1_ReasRT,.~.+overall_gainscore)
og3_ReasRT=update(og2_ReasRT,.~overall_gainscore*group)

og1_SelectiveAttention=lm(GainSelectiveAtt~group,data=gaingroupcast)
og2_SelectiveAttention=update(og1_SelectiveAttention,.~.+overall_gainscore)
og3_SelectiveAttention=update(og2_SelectiveAttention,.~overall_gainscore*group)

og1_DividedAttention=lm(GainDivAtt~group,data=gaingroupcast)
og2_DividedAttention=update(og1_DividedAttention,.~.+overall_gainscore)
og3_DividedAttention=update(og2_DividedAttention,.~overall_gainscore*group)

og1_PerceptualSpeed=lm(GainPerceptualSpeed~group,data=gaingroupcast)
og2_PerceptualSpeed=update(og1_PerceptualSpeed,.~.+overall_gainscore)
og3_PerceptualSpeed=update(og2_PerceptualSpeed,.~overall_gainscore*group)

og1_EpisodicMemory=lm(GainEpisMem~group,data=gaingroupcast)
og2_EpisodicMemory=update(og1_EpisodicMemory,.~.+overall_gainscore)
og3_EpisodicMemory=update(og2_EpisodicMemory,.~overall_gainscore*group)

multreg_inter<-rbind(hregfxn(list(og1_WMn,og2_WMn,og3_WMn),1),
                     hregfxn(list(og1_WMspan,og2_WMspan,og3_WMspan),1),
                     hregfxn(list(og1_ReasACC,og2_ReasACC,og3_ReasACC),1),
                     hregfxn(list(og1_ReasRT,og2_ReasRT,og3_ReasRT),1),
                     hregfxn(list(og1_SelectiveAttention,og2_SelectiveAttention,og3_SelectiveAttention),1),
                     hregfxn(list(og1_DividedAttention,og2_DividedAttention,og3_DividedAttention),1),
                     hregfxn(list(og1_PerceptualSpeed,og2_PerceptualSpeed,og3_PerceptualSpeed),1),
                     hregfxn(list(og1_EpisodicMemory,og2_EpisodicMemory,og3_EpisodicMemory),1))


