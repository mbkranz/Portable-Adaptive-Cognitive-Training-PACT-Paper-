install.packages(c("reshape","lme4","plyr","pbkrtest","xlsx"))
library(lme4); library(reshape); library(lme4); library(pbkrtest);library(xlsx)
npwd<-'/Users/michaelkranz/Documents/Analyses/PACT_behavpaper_Pauline'
setwd(npwd)
np<-read.xlsx(file="USED_PACTBehTaskGainScores_081315.xlsx",sheetName="GainScores_forRscript")
np.melt<-melt(np,id=c("Subject","Group"))
np.melt<-rename(np.melt,c("value"="Score","variable"="Task"))
np.melt$Score<-as.numeric(np.melt$Score) #non numeric values were included (may want to just delete these in excel sheet)
#dataframe of each construct (change the  variable task names for tasks included in construct)
names.WMnback<-c("NBackThreeBackDPrime","NBackTwoBackDPrime","Dual2back","Dual3back") 
names.WMspan<-c("OSpanScore","RunSpanScore","SymSpanScore")
names.PSpeed<-c("SalthouseLetterCompCorrect","SalthousePatternCompCorrect","SalthouseDSSTCorrect")
names.AttControl=c("AntiSaccadeAcc","FlankerEffect","VigCounterBottomQuantileRT")
names.Reasoning=c("FormBoardsTotalCorrect","LetterSetsTotalCorrect",
"MatrixTotalCorrect","SpatialRelTotalCorrect","PaperFoldTotalCorrect","SalthouseShipleyTotalCorrect")
names.ReasoningRT=c("LetterSetsRT","MatrixRT","SpatialRelRT","PaperFoldRT")
names.EpisMem= c("SalthouseLogMemTotal","SalthousePairedAssociates","iPositionSwap")
names.DivAtt<-c("SalthouseTrailsBminusA","DodgeMaxLevelPassed","AttBlink")
#checks all Tasks are correctly typed
unique(np.melt[np.melt$Task %in% names.WMnback,"Task"])
unique(np.melt[np.melt$Task %in% names.WMspan,"Task"])
unique(np.melt[np.melt$Task %in% names.PSpeed,"Task"])
unique(np.melt[np.melt$Task %in% names.AttControl,"Task"])
unique(np.melt[np.melt$Task %in% names.Reasoning,"Task"])
unique(np.melt[np.melt$Task %in% names.ReasoningRT,"Task"])
unique(np.melt[np.melt$Task %in% names.EpisMem,"Task"])
unique(np.melt[np.melt$Task %in% names.DivAtt,"Task"])

lmm.fxn<-function(names.tasks){
  #takes task names to filter out np.melt (melted data frame of all np gain scores) and runs linear mixed models
  model.prelim<-lmer(Score~1+(1|Task)+(1|Subject),
                         data=np.melt[np.melt$Task %in% names.tasks,],na.action=na.omit)
  #intercepts only for random effects
  model.intercepts<-lmer(Score~Group+(1|Task)+(1|Subject),
                       data=np.melt[np.melt$Task %in% names.tasks,],na.action=na.omit,control=lmerControl(optimizer = "bobyqa"))
  #model.interceptsandslope<-lmer(Score~Group+(Group+1|Task)+(1|Subject),
   #                              data=np.melt[np.melt$Task %in% names.tasks,],na.action=na.omit,control=lmerControl(optimizer = "bobyqa"),REML=FALSE)
 # model.interceptsandslope2<-lmer(Score~Group+(Group-1|Task)+(1|Subject),
  #                              data=np.melt[np.melt$Task %in% names.tasks,],na.action=na.omit,control=lmerControl(optimizer = "bobyqa"),REML=FALSE)
  return(model.intercepts)
}

#gain score models
np.melt$Group<-as.factor(ifelse(np.melt$Group==1,"MF","AC")) #rename group to corresponding initials and make a factor
np.melt$Group<-relevel(np.melt$Group,ref="AC") #make active control reference variable of factor
#performance intercept only (see fxn above) on each construct.
lmm.WMnback<-lmm.fxn(names.WMnback)
lmm.WMspan<-lmm.fxn(names.WMspan)
lmm.PSpeed<-lmm.fxn(names.PSpeed)
lmm.AttControl<-lmm.fxn(names.AttControl)
lmm.Reasoning<-lmm.fxn(names.Reasoning)
lmm.ReasoningRT<-lmm.fxn(names.ReasoningRT)
lmm.DivAtt<-lmm.fxn(names.DivAtt)
lmm.EpisMem<-lmm.fxn(names.EpisMem)
# code taken from: http://mindingthebrain.blogspot.com/2014/02/three-ways-to-get-parameter-specific-p.html and made into functions


# Kenward-Roger approximation pvalues with pbkrtest package (Halekoh & Højsgaard,2014)-->method used in von Bastian and Eschen 2015
# get the KR-approximated degrees of freedom
#note: p.z= normal dist pval,p.kr=K-R approx pvals
pvals.fxn<-function(lmm){
# get p-values from the t-distribution using the t-values from model
coefs <- data.frame(coef(summary(lmm))) 
#parametric p values (assuming normality)
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
#K-R p values
df.KR <- get_ddf_Lb(lmm, fixef(lmm)) # degrees of freedom
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
#round p values
coefs$p.KR<-round(coefs$p.KR,digits=3)
coefs$p.z<-round(coefs$p.z,digits=3)
return(coefs)}

#Satterthwaite approximation (based on SAS PROC mixed )-->needs package lmerTest and to rerun lmer models (see website for more info)
# get Satterthwaite-approximated degrees of freedom
#pvals.satt.fxn<-function(lmm){
#coefs$df.Satt <- coef(summary(lmm))[, 3]
# get approximate p-values
#coefs$p.Satt <- coef(summary(lmm))[, 5])
#return(coefs)}

pvals.WMnback<-pvals.fxn(lmm.WMnback)
pvals.WMnback$construct<-"WMnback"
pvals.WMspan<-pvals.fxn(lmm.WMspan)
pvals.WMspan$construct<-"WMspan"
pvals.PSpeed<-pvals.fxn(lmm.PSpeed)
pvals.PSpeed$construct<-"PSpeed"
pvals.AttControl<-pvals.fxn(lmm.AttControl)
pvals.AttControl$construct<-"AttControl"
pvals.Reasoning<-pvals.fxn(lmm.Reasoning)
pvals.Reasoning$construct<-"Reasoning"
pvals.ReasoningRT<-pvals.fxn(lmm.ReasoningRT)
pvals.ReasoningRT$construct<-"ReasoningRT"
pvals.EpisMem<-pvals.fxn(lmm.EpisMem)
pvals.EpisMem$construct<-"EpisMem"
pvals.DivAtt<-pvals.fxn(lmm.DivAtt)
pvals.DivAtt$construct<-"DivAtt"

#rep(NA,each=4) to match with random effects columns
fixedeff<-rbind(pvals.WMnback,rep(NA,each=7),
                pvals.WMspan,rep(NA,each=7),
                pvals.Reasoning,rep(NA,each=7),
                pvals.ReasoningRT,rep(NA,each=7),
                pvals.AttControl,rep(NA,each=7),
                pvals.DivAtt,rep(NA,each=7),
                pvals.PSpeed,rep(NA,each=7),
                pvals.EpisMem)

write.csv(fixedeff,file="LinearMixedModelFixedEffects.csv")
raneff.fxn<-function(lmm){
  raneff<-data.frame(VarCorr(lmm))
  return(raneff)}
raneff.WMnback<-raneff.fxn(lmm.WMnback)
raneff.WMnback$construct<-"WMnback"
raneff.WMspan<-raneff.fxn(lmm.WMspan)
raneff.WMspan$construct<-"WMspan"
raneff.PSpeed<-raneff.fxn(lmm.PSpeed)
raneff.PSpeed$construct<-"PSpeed"
raneff.AttControl<-raneff.fxn(lmm.AttControl)
raneff.AttControl$construct<-"AttControl"
raneff.Reasoning<-raneff.fxn(lmm.Reasoning)
raneff.Reasoning$construct<-"Reasoning"
raneff.ReasoningRT<-raneff.fxn(lmm.ReasoningRT)
raneff.ReasoningRT$construct<-"ReasoningRT"
raneff.EpisMem<-raneff.fxn(lmm.EpisMem)
raneff.EpisMem$construct<-"EpisMem"
raneff.DivAtt<-raneff.fxn(lmm.DivAtt)
raneff.DivAtt$construct<-"DivAtt"

randomeff<-rbind(raneff.WMnback,
                 raneff.WMspan,
                 raneff.Reasoning,
                 raneff.ReasoningRT,
                 raneff.AttControl,
                 raneff.DivAtt,
                 raneff.PSpeed,
                 raneff.EpisMem)
write.csv(randomeff,file="LinearMixedModelRandomEffects.csv")


