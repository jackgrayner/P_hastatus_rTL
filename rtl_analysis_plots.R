####
# rTL analysis
####

library(ggplot2)
library(car)
library(ggbeeswarm)
library(patchwork)
library(ggeffects)
library(ggExtra)

#plot efficiency
eff<-read.csv('./tel_efficiency.csv',h=T)
eff<-eff[eff$primers %in% c("Tel","BDNF"),]
effplot<-ggplot(eff,aes(x=log10(dilution),y=Cp))+
  geom_point()+geom_smooth(method='lm',colour='black')+
  facet_wrap(.~primers,scales='free')+
  theme_bw()+
  labs(subtitle="Tel slope = -3.237, amplification factor = 2.04, efficiency = 103.66%\nBDNF slope = -3.431, amplification factor = 1.96, efficiency = 95.63%")
summary(lm(Cp~log10(dilution),dat=eff[eff$primers=="Tel",]))
summary(lm(Cp~log10(dilution),dat=eff[eff$primers=="BDNF",]))
ggsave('efficiency_plot.png',plot=effplot,dpi=600,height=4,width=6)
#

#read in filtered data
dat<-read.csv('/Users/Jack/Desktop/umd/manuscripts/iisage_telomere/all_tel_data_filtered_dryad.csv',h=T)
dat$sex<-factor(dat$SEX)
dat<-dat[!is.na(dat$CV_rtl),]#get rid of samples without replicate rTL
dat<-dat[dat$Season_sampled %in% c("Jan_23","Jan_24"),]#remove May '23 samples -- too few


#estimate ages of unaged samples (all female) based on toothwear scores 
summary(lm(Est.Age~TTH,data=dat[dat$sex=="F" & !dat$AgeSource=="ToothWear",]))
#plot relationship
ggplot(dat,aes(x=TTH,y=Est.Age,colour=sex))+scale_colour_manual(values=c("#fa8490","#9fd4fc"))+
 theme_minimal()+geom_point()+geom_smooth(method='lm')

#add estimate ages based on toothwear
dat$tthage<-NA
dat[dat$AgeSource=="ToothWear" & dat$sex=="F",]$Est.Age<-
 dat[dat$AgeSource=="ToothWear" & dat$sex=="F",]$TTH*3.1-0.6029

#plot rTL CoV histogram
CV_hist<-ggplot(dat,aes(x=CV_rtl))+geom_histogram(fill='#dddddd',colour='#555555')+theme_bw()+
  geom_vline(xintercept=0.5,colour='darkred',linetype='dashed')+
  geom_vline(xintercept=median(dat$CV_rtl),colour='blue',linetype='dashed')+
  xlim(c(0,1))+labs(tag="A")
# ggsave('rTL_CoV_hist.png',plot=CV_hist,dpi=300,height=3,width=4)

#plot primer CoV histograms
CV_tel<-ggplot(dat,aes(x=CV_tel))+
  geom_histogram(fill='#dddddd',colour='#555555',bins=100)+theme_bw()+
  geom_vline(xintercept=median(dat$CV_tel,na.rm = TRUE),colour='darkblue',linetype='dashed')+
  xlim(c(0,1))+labs(tag="B")
CV_BDNF<-ggplot(dat,aes(x=CV_BDNF))+geom_histogram(fill='#dddddd',colour='#555555',bins=100)+theme_bw()+
  geom_vline(xintercept=median(dat$CV_BDNF,na.rm = TRUE),colour='darkblue',linetype='dashed')+
  xlim(c(0,1))+labs(tag="C")
#ggsave('rTL_CoV_rTL_Tel_BDNF.png',plot=
#         CV_hist /
#         (CV_tel|CV_BDNF),dpi=300,height=5,width=5.5)


#remove samples with replicate rTL CoV > 0.5
dat<-dat[!dat$CV_rtl>0.5 & !is.na(dat$tlmean),]

#replace missing rtl1 / rtl2 values with rtl3 where applicable
dat[is.na(dat$rtl1),]$rtl1<-dat[is.na(dat$rtl1),]$rtl3
dat[is.na(dat$rtl2),]$rtl2<-dat[is.na(dat$rtl2),]$rtl3

#interplate correlation
cor.test(log(dat$rtl1),log(dat$rtl2))

#plot correlation between rtl1, rtl2
g.rtlcor<-ggplot(dat,aes(x=log2(rtl1),y=log2(rtl2),colour=plate))+geom_point(alpha=0.75,size=0.75)+
  geom_smooth(method='lm',se=FALSE,linewidth=0.5)+
  geom_smooth(data=dat,aes(x=log2(rtl1),y=log2(rtl2)),colour='black',method='lm',se=TRUE,linewidth=1.5)+
  xlab("Log2 rTL_1")+ylab("Log2 rTL_2")+
  theme_bw()+scale_colour_viridis(discrete=TRUE)+
  theme(legend.position='none')+
  xlim(c(-1.5,1.5))+ylim(c(-1.5,1.5))
#ggsave('rTL_replicate_correlation.png',plot=g.rtlcor,dpi=600,height=4,width=4)

#remove jan_23 vals for duplicate bands
dupes<-dat[duplicated(dat$band),]$band
dupes.df<-dat[dat$band %in% dupes,]
dupes.df$band<-factor(dupes.df$band)
dat<-dat[!(dat$band %in% dupes.df$band & dat$Season_sampled=="Jan_23"),]

#does it make a difference if we instead remove the jan_24 duplicate band values? - No
#dat<-dat[!duplicated(dat$band),]

#run forearm/Weight regressions
lm.fa<-lm(scale(FA)~sex+population,data=dat)
lm.wt<-lm(scale(WT)~sex+population+sex*poly(scale(Est.Age,scale=FALSE),2),data=dat)
Anova(lm.fa,type="II")
Anova(lm.wt,type="III")
summary(lm.fa)
summary(lm.wt)
#plot(lm.fa)#good
#plot(lm.wt)#good

#plot size and age across sexes
g.fa<-ggplot(dat,aes(x=sex,y=scale(FA),colour=sex))+
  scale_colour_manual(values=c("#fa8490","#9fd4fc"))+
  scale_fill_manual(values=c("#fa8490","#9fd4fc"))+
  geom_quasirandom(size=0.57)+
  stat_summary(fun.data=mean_sdl,fun.args=list(mult=1),
               colour='black',geom="errorbar",width=0.1)+
  stat_summary(fun="mean",shape=21,colour='black',aes(fill=sex))+
  theme_minimal()+theme(legend.position='none',axis.title.x=element_blank())+
  ylab("Forearm length")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='#fcfbfa',colour='#aaaaaa'))+
  labs(tag="A",caption=" ")
g.wt<-ggplot(dat,aes(x=sex,y=scale(WT),colour=sex))+
  scale_colour_manual(values=c("#fa8490","#9fd4fc"))+
  scale_fill_manual(values=c("#fa8490","#9fd4fc"))+
  geom_quasirandom(size=0.57)+
  stat_summary(fun.data=mean_sdl,fun.args=list(mult=1),
               colour='black',geom="errorbar",width=0.1)+
  stat_summary(fun="mean",shape=21,colour='black',aes(fill=sex))+  
  theme_minimal()+theme(legend.position='none',axis.title.x=element_blank())+
  ylab("Weight")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='#fcfbfa',colour='#aaaaaa'))+
  labs(tag="B",caption="mean ± SD")

#exploratory plot showing non-linear assocation between age and weight
g.wt.nonlinear<-ggplot(dat,aes(x=Est.Age,y=scale(WT),colour=sex))+
 scale_colour_manual(values=c("#fa8490","#9fd4fc"))+
 scale_fill_manual(values=c("#fa8490","#9fd4fc"))+
 geom_point(size=0.57)+
 geom_smooth(method='loess',span=1.5)+
 theme_minimal()+theme(legend.position='none',axis.title.x=element_blank())+
 ylab("Weight")+xlab("Est. Age")+
 theme(panel.grid=element_blank(),panel.background=element_rect(fill='#fcfbfa',colour='#aaaaaa'))
#ggsave('age_weight_loss.png',plot=g.wt.nonlinear,dpi=600,height=4,width=5.5)

#run rTL model
dat$log2tlmean<-log2(dat$tlmean)
dat$population<-relevel(factor(dat$population),ref="Tamana")
rtl.lm<-lm(scale(log2tlmean)~scale(FA,scale=FALSE)+scale(Est.Age,scale=FALSE)+
             sex+population+Season_sampled+
            scale(days_collection_to_assay,scale=FALSE)+
            scale(days_since_calibrator_extracted,scale=FALSE),
           #data=dat[dat$Est.Age<10,])
           data=dat)
#plot(rtl.lm)#looks fine
summary(rtl.lm)
Anova(rtl.lm,type="II")

# plot adjusted predicted values of rTL across sexes and ages, based on rtl.lm model
# we use the empirical marginilisation method. the other methods produce strange/uninterpretable plots,
# on scales very different from the observed data.
# the plot produced by the empirical marginilisation method is much more consistent with model results
# and removal of highly collinear variables in the model causes the other methods to produce
# plots much more similar to the empirical marginalisation plot.
# empirical marginilisation also appears to be recommended when you have skewed non-focal variables,
# which we do (e.g., sig difference in number of Jan_23 and Jan_24 samples)
g <- predict_response(rtl.lm, terms=c("Est.Age","sex"),back_transform = FALSE,margin='empirical',ci_level = 0.95) 
# g <- predict_response(rtl.lm, terms=c("Est.Age","sex"),back_transform = FALSE,margin='marginalmeans',ci_level = 0.95) 
# g <- predict_response(rtl.lm, terms=c("Est.Age","sex"),back_transform = FALSE,margin='mean_reference') 
# ggPredict(rtl.lm, terms=c("Est.Age","sex"),interactive=TRUE) 
g1<-data.frame(g)

# this could be used to adjust raw data points in plot by mean difference between season_sampled
# but we didn't do this for fig. 1
# diffyear<-mean(dat[dat$Season_sampled=="Jan_24" & dat$sex=="F",]$log2tlmean)-
#   mean(dat[dat$Season_sampled=="Jan_23" & dat$sex=="F",]$log2tlmean)
# dat1<-dat
# dat1[dat1$Season_sampled=="Jan_24",]$log2tlmean<-dat1[dat1$Season_sampled=="Jan_24",]$log2tlmean-diffyear
# summary(dat[dat$Season_sampled=="Jan_24",]$log2tlmean)
# summary(dat1[dat1$Season_sampled=="Jan_24",]$log2tlmean)

g.rtl<-ggplot() +
  geom_line(data=g1[g1$group=="M",],aes(x=x,y=predicted),colour="#5090bf",size=1) +
  geom_ribbon(data=g1[g1$group=="M",],linetype='dotted',
            aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high),alpha=0.2,colour="#5090bf",fill="#5090bf")+
  geom_rect(fill='#fcfbfa',aes(xmin=max(dat[dat$sex=="M",]$Est.Age),xmax=20,ymin=-2,ymax=2))+
  geom_line(data=g1[g1$group=="F",],aes(x=x,y=predicted),colour="#fa8490",size=1) +
  geom_ribbon(data=g1[g1$group=="F",],linetype='dotted',
             aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high),alpha=0.2,colour="#fa8490",fill="#fa8490")+
  geom_point(data=dat,aes(x=Est.Age,y=log2tlmean,colour=sex),alpha=0.6)+
 # geom_smooth(data=dat,aes(x=Est.Age,y=log2tlmean,colour=sex),alpha=0.6,method='lm')+
  scale_x_continuous(breaks=seq(1,20,2))+
  labs(x = "Age",y = "log2 rTL",color = "Sex") +
  theme_minimal()+scale_colour_manual(values=c("#fa8490","#5090bf"))+
  theme(panel.grid=element_blank(),legend.position='none',panel.background=element_rect(fill='#fcfbfa',colour='#aaaaaa'))
g.rtl

ggsave('fig1_size_rtl.png',dpi=600,height=7.5,width=6,
       plot=(g.fa | g.wt ) /
         ggMarginal(g.rtl,groupFill = TRUE,margins = 'x')+ labs(tag="C") + 
         plot_layout(heights = c(1, 1.75)))
