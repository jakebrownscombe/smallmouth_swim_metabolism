#Smallmouth bass swimming metabolism analysis 
#Journal of Fish Biology "Accelerometer-based swimming metabolism of smallmouth bass (Micropterus dolomieu)" 
#Jacob W. Brownscombe, Kurtis Smith, Paul Bzonek 


#load packages, R2 function, data ----
library(dplyr)
library(ggplot2)
library(shiny)
library(tidyverse)
library(plotly)
library(patchwork)
library(nlme)
library(AICcmodavg)
source("~/github/smallmouth_swim_metabolism/worksheets/functions.R") #R2 function

#load data 
swim.data.acc <- readRDS("~/github/smallmouth_swim_metabolism/data/swim.data.acc.rds")
swim.data.all <- readRDS("~/github/smallmouth_swim_metabolism/data/swim.data.all.rds")


#acc vs MO2 ----
head(swim.data.acc)
ggplot(swim.data.acc, aes(gs, log(O2.h.kg), col=temp.cat))+geom_point()+geom_smooth(method="lm")

#model
swim.data.acc$logMO2 <- log(swim.data.acc$O2.h.kg)
swim.data.acc$logmass.kg <- log(swim.data.acc$mass.kg)
swim.data.acc$logtemp <- log(swim.data.acc$temp)
swim.data.acc$loggs <- log(swim.data.acc$gs)


M1<-lme(logMO2 ~ gs + logmass.kg + logtemp + Sex +
          gs*logmass.kg + gs*logtemp + logmass.kg*logtemp,
        random = ~1 | Transmitter.ID, method = "REML", 
        data = swim.data.acc)# random intercept
M2<-lme(logMO2 ~ gs + logmass.kg + logtemp + Sex +
          gs*logmass.kg + gs*logtemp + logmass.kg*logtemp,
        random = ~gs -1 | Transmitter.ID, method = "REML", 
        data = swim.data.acc)# random slope

#this doesn't converge
#M3<-lme(logMO2 ~ gs + logmass.kg + logtemp + Sex +
#          gs*logmass.kg + gs*logtemp + logmass.kg*logtemp,
#        random = ~1 + gs | Transmitter.ID, method = "REML", 
#        data = swim.data.acc)# both
anova(M1,M2)
#intercept slightly better


M1.ML <-lme(logMO2 ~ gs + logmass.kg + logtemp + Sex +
              gs*logmass.kg + gs*logtemp + logmass.kg*logtemp,
            random = ~1 | Transmitter.ID, method = "ML", 
            data = swim.data.acc)# random intercept
drop1(M1.ML, test="Chisq")
#no interactions are helpful, neither is sex.

#all log
M1.final <- lme(logMO2 ~ gs + logmass.kg + logtemp,
                random = ~1 | Transmitter.ID, method = "REML", 
                data = swim.data.acc)

summary(M1.final)
rsquared.lme(list(M1.final)) 

#model val
E <- resid(M1.final)
F <- fitted(M1.final)
plot(E~F)
plot(E ~ swim.data.acc$gs)
plot(E ~ swim.data.acc$logtemp)
plot(E ~ swim.data.acc$mass.kg)


#show some predictions
mass <- expand.grid(mass.kg= seq(0.5, 3, 0.05))
temps <- expand.grid(temp= seq(1, 35, 1))
acc <- expand.grid(gs= seq(0.0, 3, 0.1))
preds.acc <- merge(mass, temps)
preds.acc <- merge(preds.acc, acc)
preds.acc$logmass.kg <- log(preds.acc$mass.kg)
preds.acc$logtemp <- log(preds.acc$temp)

#predict
preds.acc$MO2.pred <- exp(predict(M1.final, preds.acc, level = 0))

#select data to plot
preds.acc2 <- preds.acc %>% filter(temp==6 & mass.kg==0.5 | temp==6 & mass.kg==2 |
                                     temp==12 & mass.kg==0.5 | temp==12 & mass.kg==2 |
                                     temp==18 & mass.kg==0.5 | temp==18 & mass.kg==2 |
                                     temp==24 & mass.kg==0.5 | temp==24 & mass.kg==2 )
#plots
ggplot(swim.data.acc, aes(gs, O2.h.kg, col=as.factor(temp.cat)))+geom_point(pch=3)+
  geom_smooth(data=preds.acc2, aes(gs, MO2.pred, col=as.factor(temp), linetype=as.factor(mass.kg)), fill="NA")+
  scale_y_continuous(limits=c(0,400))+scale_x_continuous(limits=c(0,2.5))+
  scale_linetype(name="Mass (kg)")+
  scale_color_viridis_d(option="plasma", begin=0, end=0.9, name="ºC")+theme_bw()+
  labs(x="Acceleration (g)", y=bquote(~MO[2]* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), 
       title="Acceleration - Metabolism", subtitle="A) Acceleration + Temperature")+
  
  ggplot(preds.acc %>% filter(temp==15), aes(gs, mass.kg,  fill=MO2.pred))+geom_tile()+
  scale_fill_viridis_c(option="magma", name=bquote(~MO[2]))+theme_classic()+
  labs(x="Acceleration (g)", y="Mass (kg)", subtitle="B) Acceleration + Mass (15ºC)")+
  
  ggplot(preds.acc %>% filter(mass.kg==1), aes(gs, temp, fill=MO2.pred))+geom_tile()+
  scale_fill_viridis_c(option="magma", name=bquote(~MO[2]))+theme_classic()+
  labs(x="Acceleration (g)", y="Temperature (ºC)", subtitle="C) Acceleration + Temperature (1 kg fish)")+
  
  ggplot(preds.acc %>% filter(gs==1), aes(mass.kg, temp, fill=MO2.pred))+geom_tile()+
  scale_fill_viridis_c(option="magma", name=bquote(~MO[2]))+theme_classic()+
  labs(x="Mass (kg)", y="Temperature (ºC)", subtitle="D) Mass + Temperature (1 g acceleration)")










#acc vs speed ----
ggplot(swim.data.acc, aes(speed.corr, gs))+geom_point()+geom_smooth(method="lm", formula= y~poly(x,2))+
  scale_y_continuous(limits=c(0,3))+scale_x_continuous(limits=c(0,120))


#models
SM1<-lme(loggs ~ speed.corr + mass.kg + temp + Sex +
           speed.corr*mass.kg + speed.corr*temp + speed.corr*Sex,
         random = ~1 | Transmitter.ID, method = "REML", 
         data = swim.data.acc)# random intercept
SM2<-lme(loggs ~ speed.corr + mass.kg + temp + Sex +
           speed.corr*mass.kg + speed.corr*temp + speed.corr*Sex,
         random = ~speed.corr -1 | Transmitter.ID, method = "REML", 
         data = swim.data.acc)# random slope
SM3<-lme(loggs ~ speed.corr + mass.kg + temp + Sex +
           speed.corr*mass.kg + speed.corr*temp + speed.corr*Sex,
         random = ~1 + speed.corr | Transmitter.ID, method = "REML", 
         data = swim.data.acc)# both
anova(SM1,SM2, SM3)
#intercept 


SM1.ML <-lme(loggs ~ speed.corr + mass.kg + temp + Sex +
               speed.corr*mass.kg + speed.corr*temp + speed.corr*Sex,
             random = ~1 | Transmitter.ID, method = "ML", 
             data = swim.data.acc)
drop1(SM1.ML, test="Chisq")


SM2.ML <-lme(loggs ~ speed.corr + mass.kg + temp + Sex,
             random = ~1 | Transmitter.ID, method = "ML",
             data = swim.data.acc)# random intercept
drop1(SM2.ML, test="Chisq")


SM.final <-lme(loggs ~ speed.corr + mass.kg + temp + Sex,
               random = ~1 | Transmitter.ID, method = "REML",
               data = swim.data.acc)
summary(SM.final)
rsquared.lme(list(SM.final)) 

#intercept
exp(-1.69) #0.185 is resting 

#model validation 
E <- resid(SM.final)
F <- fitted(SM.final)
plot(E~F)
plot(E ~ swim.data.acc$speed.corr)
plot(E ~ swim.data.acc$temp)
plot(E ~ swim.data.acc$mass.kg)
plot(E ~ swim.data.acc$Sex)
#


#predictions
mass <- expand.grid(mass.kg= seq(0.5, 3, 0.05))
temps <- expand.grid(temp= seq(1, 35, 1))
speed <- expand.grid(speed.corr=seq(0,105,5))
preds.swim <- merge(mass, temps)
preds.swim <- merge(preds.swim, speed)
preds.swim <- merge(preds.swim, data.frame(Sex=as.factor(c("Male","Female"))))

preds.swim$gs.pred <- exp(predict(SM.final, preds.swim, level = 0))



#plots
swim.data.acc$mass.cat <- ifelse(swim.data.acc$mass.kg<=0.7, 0.5, 
                                 ifelse(swim.data.acc$mass.kg<=1.5, 1,
                                        ifelse(swim.data.acc$mass.kg<=2, 1.5,
                                               ifelse(swim.data.acc$mass.kg<=2.5, 2, 2.5))))

ggplot(swim.data.acc, aes(speed.corr, gs, col=temp.cat))+geom_point()+
  geom_smooth(data=preds.swim %>% filter(temp==6|temp==12|temp==18|temp==24), 
              aes(speed.corr, gs.pred, col=as.factor(temp)), fill="NA")+
  labs(x="Swimming Speed (cm/s)", y="Acceleration (g)", title="Swimming Speed - Acceleration", subtitle="A) Speed + Temperature")+
  scale_color_viridis_d(option="plasma", begin=0, end=0.9, name="Temperature")+theme_bw()+
  
  ggplot(swim.data.acc, aes(speed.corr, gs, col=mass.kg))+geom_point()+
  geom_smooth(data=preds.swim %>% filter(mass.kg==0.5 | mass.kg==1 | mass.kg==1.5 | mass.kg==2 | mass.kg==2.5), 
              aes(speed.corr, gs.pred, group=as.factor(mass.kg)), fill="NA")+
  labs(x="Swimming Speed (cm/s)", y="Acceleration (g)", subtitle="B) Speed + Mass")+theme_bw()+
  scale_color_viridis_c(begin=0, end=0.9, name="Mass (kg)")+
  
  ggplot(swim.data.acc, aes(speed.corr, gs, col=Sex))+geom_point()+
  geom_smooth(data=preds.swim, aes(speed.corr, gs.pred, col=Sex), fill="NA")+
  labs(x="Swimming Speed (cm/s)", y="Acceleration (g)", subtitle="C) Speed + Sex")+theme_bw()+
  
  ggplot(preds.swim %>% filter(speed.corr==50), aes(mass.kg, temp, fill=gs.pred))+geom_tile()+
  scale_fill_viridis_c(option="mako", name="Acc (g)")+
  labs(x="Mass (kg)", y="Temperature (ºC)", subtitle="D) Mass + Temperature (50 cm/s)")+theme_classic()
#











#speed vs metabolism ----
head(swim.data.all)
ggplot(swim.data.all, aes(speed.corr, log(O2.h.kg), col=temp.cat))+geom_point()+geom_smooth(method="lm")

#models
SMO1<-lme(logMO2 ~ speed.corr + logmass + temp + Sex +
            speed.corr*logmass + speed.corr*temp + speed.corr*Sex,
          random = ~1 | Transmitter.ID, method = "REML", 
          data = swim.data.all) # random intercept
SMO2<-lme(logMO2 ~ speed.corr + logmass + temp + Sex +
            speed.corr*logmass + speed.corr*temp + speed.corr*Sex,
          random = ~speed.corr -1 | Transmitter.ID, method = "REML", 
          data = swim.data.all)# random slope
SMO3<-lme(logMO2 ~ speed.corr + logmass + temp + Sex +
            speed.corr*logmass + speed.corr*temp + speed.corr*Sex,
          random = ~1 + speed.corr | Transmitter.ID, method = "REML", 
          data = swim.data.all)# both
anova(SMO1,SMO2,SMO3)
#intercept 

SMO1.1 <-lme(logMO2 ~ speed.corr + logmass + temp + Sex +
               speed.corr*logmass + speed.corr*temp + speed.corr*Sex,
             random = ~1 | Transmitter.ID, method = "ML", 
             data = swim.data.all)
drop1(SMO1.1, test="Chisq")

SMO1.2 <-lme(logMO2 ~ speed.corr + logmass + temp + Sex +
               speed.corr*logmass + speed.corr*temp,
             random = ~1 | Transmitter.ID, method = "ML", 
             data = swim.data.all)
drop1(SMO1.2, test="Chisq")
rsquared.lme(list(SMO1.2)) 

#retain sex for ease of use with other models.. have to drop temp interactions, crazy patterns. 
SMO.final <-lme(logMO2 ~ speed.corr + logmass + temp + logmass*speed.corr + Sex,
                random = ~1 | Transmitter.ID, method = "REML", 
                data = swim.data.all)
summary(SMO.final)
rsquared.lme(list(SMO.final)) 



#from non interacting model 
summary(lme(logMO2 ~ speed.corr + logmass + temp,
            random = ~1 | Transmitter.ID, method = "REML", 
            data = swim.data.all))

exp(4.08) #intercept
exp(0.038) #temp
exp(-0.27) #mass


#model val
E <- resid(SMO.final)
F <- fitted(SMO.final)
plot(E~F)
plot(E ~ swim.data.all$speed.corr)
plot(E ~ swim.data.all$temp)
plot(E ~ swim.data.all$mass.kg)


#preds
preds.swimMO2 <- merge(mass, temps)
preds.swimMO2 <- merge(preds.swimMO2, speed)
preds.swimMO2 <- merge (preds.swimMO2, data.frame(Sex=c("Male","Female")))
preds.swimMO2$logmass <- log(preds.swimMO2$mass.kg)
preds.swimMO2$logtemp <- log(preds.swimMO2$temp)
preds.swimMO2$MO2.pred <- exp(predict(SMO.final, preds.swimMO2, level = 0))
head(preds.swimMO2)


preds.swimMO2.2 <- preds.swimMO2 %>% filter(temp==6 & mass.kg==0.5 | temp==6 & mass.kg==2 |
                                              temp==12 & mass.kg==0.5 | temp==12 & mass.kg==2 |
                                              temp==18 & mass.kg==0.5 | temp==18 & mass.kg==2 |
                                              temp==24 & mass.kg==0.5 | temp==24 & mass.kg==2 )

#plots 
ggplot(swim.data.all, aes(speed.corr, O2.h.kg, col=as.factor(temp.cat)))+geom_point(pch=3)+
  geom_smooth(data=preds.swimMO2.2, aes(speed.corr, MO2.pred, col=as.factor(temp), linetype=as.factor(mass.kg)), fill="NA")+
  scale_color_viridis_d(option="plasma", begin=0, end=0.9, name="ºC")+
  scale_y_continuous(limits=c(0,480))+
  scale_linetype(name="Mass (kg)")+
  labs(x="Swimming Speed (cm/s)", y=bquote(~MO[2]* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), 
       title="Swimming Speed - Metabolism", subtitle="A) Speed + Temperature")+theme_bw()+
  
  ggplot(preds.swimMO2 %>% filter(temp==15), aes(speed.corr, mass.kg, fill=MO2.pred))+geom_tile()+scale_fill_viridis_c(option="magma", name=bquote(~MO[2]))+
  theme_classic()+labs(y="Fish Mass (kg)", x="Swimming Speed (cm/s)", subtitle="B) Speed + Mass (15ºC)")+
  
  ggplot(preds.swimMO2 %>% filter(mass.kg==1), aes(speed.corr, temp, fill=MO2.pred))+geom_tile()+scale_fill_viridis_c(option="magma", name=bquote(~MO[2]))+
  theme_classic()+labs(y="Temperature (ºC)", x="Swimming Speed (cm/s)", subtitle="B) Speed + Temperature (1 kg fish)")+
  
  ggplot(preds.swimMO2 %>% filter(speed.corr==50), aes(mass.kg, temp, fill=MO2.pred))+geom_tile()+scale_fill_viridis_c(option="magma", name=bquote(~MO[2]))+
  theme_classic()+labs(x="Fish Mass (kg)", y="Temperature (ºC)", subtitle="B) Temperature + Mass (50 cm/s)")

#









#Ucrit ----
head(swim.data.all)
ggplot(swim.data.all, aes(temp.cat, Ucrit))+geom_boxplot()+scale_y_continuous(limits=c(0, max(swim.data.all$Ucrit)+15))+
  theme_bw()+xlab("Temperature (ºC)")+ylab(bquote("Ucrit (cm"~s^'–1'*')'))
ggplot(swim.data.all, aes(temp, Ucrit))+geom_point()+geom_smooth(method='lm', formula=y ~ poly(x, 2))+
  theme_bw()+xlab("Temperature (ºC)")+ylab(bquote("Ucrit (cm"~s^'–1'*')'))
ggplot(swim.data.all, aes(log(temp), Ucrit))+geom_point()+geom_smooth(method='lm')+
  theme_bw()+xlab("Temperature (ºC)")+ylab(bquote("Ucrit (cm"~s^'–1'*')'))


#models
#models
UC1<-lme(Ucrit ~ poly(temp,2) + mass.kg +  Sex + 
           poly(temp,2)*mass.kg + poly(temp,2)*Sex,
         random = ~1 | Transmitter.ID, method = "REML", 
         data = swim.data.all) # random intercept
UC2<-lme(Ucrit ~ poly(temp,2) + mass.kg +  Sex + 
           poly(temp,2)*mass.kg + poly(temp,2)*Sex,
         random = ~speed.corr -1 | Transmitter.ID, method = "REML", 
         data = swim.data.all)# random slope
UC3<-lme(Ucrit ~ poly(temp,2) + mass.kg +  Sex + 
           poly(temp,2)*mass.kg + poly(temp,2)*Sex,
         random = ~1 + speed.corr | Transmitter.ID, method = "REML", 
         data = swim.data.all)# both
anova(UC1,UC2,UC3)

#intercept
UC1.1<-lme(Ucrit ~ poly(temp,2) + mass.kg +  Sex + 
             poly(temp,2)*mass.kg + poly(temp,2)*Sex,
           random = ~1 | Transmitter.ID, method = "ML", 
           data = swim.data.all)
drop1(UC1.1, test="Chisq")


UC.final<-lme(Ucrit ~ poly(temp,2) + mass.kg +  Sex + 
                poly(temp,2)*Sex,
              random = ~1 | Transmitter.ID, method = "REML", 
              data = swim.data.all)
summary(UC.final)
rsquared.lme(list(UC.final)) 

#model val
E <- resid(UC.final)
F <- fitted(UC.final)
plot(E~F)
plot(E ~ swim.data.all$speed.corr)
plot(E ~ swim.data.all$temp)
plot(E ~ swim.data.all$mass.kg)


#preds
preds.UC <- merge(mass, temps)
preds.UC <- merge(preds.UC, data.frame(Sex=c("Male","Female")))
preds.UC$logmass <- log(preds.UC$mass.kg)
preds.UC$logtemp <- log(preds.UC$temp)
preds.UC$Ucrit.pred <- predict(UC.final, preds.UC, level = 0)
head(preds.UC)


ggplot(swim.data.all, aes(temp, Ucrit, col=as.factor(mass.kg)))+geom_point()+
  geom_smooth(data=preds.UC %>% filter(mass.kg==0.5 | mass.kg==1 | mass.kg==1.5 | mass.kg==2 | mass.kg==2.5),
              aes(temp, Ucrit.pred, col=as.factor(mass.kg)), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote("Ucrit (cm"~s^'–1'*')'), title="Ucrit", subtitle="A) Temperature + Mass")+
  scale_y_continuous(limits=c(0,110))+scale_color_viridis_d(begin=0, end=0.9, name="Mass (kg)")+
  
  ggplot(swim.data.all, aes(temp, Ucrit, col=Sex))+geom_point()+
  geom_smooth(data=preds.UC %>% filter(mass.kg==1),
              aes(temp, Ucrit.pred, col=Sex), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote("Ucrit (cm"~s^'–1'*')'), subtitle="B) Temperature + Sex (1 kg fish)")+
  scale_y_continuous(limits=c(0,110))
#















#MMR ----

fish.metrics <- swim.data.all %>% group_by(Transmitter.ID, temp.cat) %>% 
  summarise( temp=mean(temp), mass.kg=mean(mass.kg), Ucrit=mean(Ucrit), MMR=max(O2.h.kg)) %>% as.data.frame()
fish.metrics$Sex <- swim.data.all$Sex[match(fish.metrics$Transmitter.ID, swim.data.all$Transmitter.ID)]
head(fish.metrics)

ggplot(fish.metrics, aes(temp, MMR, col=mass.kg))+geom_point()+geom_smooth(method='lm', formula=y ~ poly(x, 2))+
  theme_bw()+ylab("MMR O2.mg.kg.h")+xlab("Temperature (ºC)")+
  scale_y_continuous(limits=c(0,300))+scale_x_continuous(limits=c(0,25))

ggplot(fish.metrics, aes(log(temp), MMR, col=mass.kg))+geom_point()+geom_smooth(method='lm', formula=y ~ poly(x, 2))+
  theme_bw()+ylab("MMR O2.mg.kg.h")+xlab("Temperature (ºC)")


fish.metrics$logMMR <- log(fish.metrics$MMR)
fish.metrics$logtemp <- log(fish.metrics$temp)
fish.metrics$logmass <- log(fish.metrics$mass.kg)
fish.metrics$Sex <- as.factor(fish.metrics$Sex)
fish.metrics <- fish.metrics %>% droplevels()

mmr.lme <-lme(logMMR ~ poly(logtemp,2) + logmass + Sex +
                poly(logtemp,2)*logmass + poly(logtemp,2)*Sex,
              random = ~1 | Transmitter.ID, method = "ML", data = fish.metrics)# random intercept
drop1(mmr.lme, test="Chisq")


mmr.lme2 <-lme(logMMR ~ poly(logtemp,2) + logmass + Sex,
               random = ~1 | Transmitter.ID, method = "ML", data = fish.metrics)# random intercept
drop1(mmr.lme2, test="Chisq")

mmr.final <-lme(logMMR ~ poly(temp,2) + logmass + Sex,
                random = ~1 | Transmitter.ID, method = "REML", 
                data = fish.metrics)
summary(mmr.final)
rsquared.lme(list(mmr.final)) 

#preds
preds.swimMO2$MMR.pred <- exp(predict(mmr.final, preds.swimMO2, level = 0))
head(preds.swimMO2)

ggplot(fish.metrics, aes(temp, MMR, col=as.factor(mass.kg)))+geom_point()+
  geom_smooth(data=preds.swimMO2 %>% filter(mass.kg==0.7 | mass.kg==1 | mass.kg==1.5 | mass.kg==2),
              aes(temp, MMR.pred, col=as.factor(mass.kg)), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), title="MMR", subtitle="A) Temperature + Mass")+
  scale_y_continuous(limits=c(0,400))+scale_color_viridis_d(begin=0, end=0.9, name="Mass (kg)")+ 
  
  ggplot(fish.metrics, aes(temp, MMR, col=Sex))+geom_point()+
  geom_smooth(data=preds.swimMO2 %>% filter(mass.kg==1),
              aes(temp, MMR.pred, col=Sex), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), subtitle="B) Temperature + Sex")+
  scale_y_continuous(limits=c(0,400))




#plot Ucrit and MMR together 
swim.data.all$mass.r <- round(swim.data.all$mass.kg, 1)
fish.metrics$mass.r <- round(fish.metrics$mass.kg, 1)
#MMR
ggplot(fish.metrics, aes(temp, MMR, col=as.factor(mass.r)))+geom_point()+
  geom_smooth(data=preds.swimMO2 %>% filter(mass.kg==0.5 | mass.kg==1 | mass.kg==1.5 | mass.kg==2),
              aes(temp, MMR.pred, col=as.factor(mass.kg)), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), title="MMR", subtitle="A) Temperature + Mass")+
  scale_y_continuous(limits=c(0,400))+scale_color_viridis_d(begin=0, end=0.9, name="Mass (kg)", breaks = c("0.5","1","1.5","2"))+ 
  
  ggplot(fish.metrics, aes(temp, MMR, col=Sex))+geom_point()+
  geom_smooth(data=preds.swimMO2 %>% filter(mass.kg==1),
              aes(temp, MMR.pred, col=Sex), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), subtitle="B) Temperature + Sex")+
  scale_y_continuous(limits=c(0,400))+
  
  #Ucrit
  
  ggplot(swim.data.all, aes(temp, Ucrit, col=as.factor(mass.r)))+geom_point()+
  geom_smooth(data=preds.UC %>% filter(mass.kg==0.5 | mass.kg==1 | mass.kg==1.5 | mass.kg==2 ),
              aes(temp, Ucrit.pred, col=as.factor(mass.kg)), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote("Ucrit (cm"~s^'–1'*')'), title="Ucrit", subtitle="A) Temperature + Mass")+
  scale_y_continuous(limits=c(0,110))+
  scale_color_viridis_d(begin=0, end=0.9, name="Mass (kg)", breaks = c("0.5","1","1.5","2"))+
  
  ggplot(swim.data.all, aes(temp, Ucrit, col=Sex))+geom_point()+
  geom_smooth(data=preds.UC %>% filter(mass.kg==1),
              aes(temp, Ucrit.pred, col=Sex), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote("Ucrit (cm"~s^'–1'*')'), subtitle="B) Temperature + Sex (1 kg fish)")+
  scale_y_continuous(limits=c(0,110))









#metabolic metrics - MMR, RMR, AS ----

head(preds.swimMO2)
preds.RMR <- preds.swimMO2 %>% filter(speed.corr==0) %>% rename(RMR.pred = MO2.pred)
head(preds.RMR)

#RMR 
#pull from swim speed model 
ggplot()+
  geom_smooth(data=preds.RMR %>% filter(mass.kg==0.5 | mass.kg==1 | mass.kg==1.5 | mass.kg==2 | mass.kg==2.5),
              aes(temp, MMR.pred, col=as.factor(mass.kg)), fill="NA", linetype=2)+
  geom_smooth(data=preds.RMR %>% filter(mass.kg==0.5 | mass.kg==1 | mass.kg==1.5 | mass.kg==2 | mass.kg==2.5),
              aes(temp, RMR.pred, col=as.factor(mass.kg)), fill="NA")+
  theme_bw()+labs(x="Temperature (ºC)", y=bquote(~MO[2]* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), title="SMR + MMR")+
  scale_y_continuous(limits=c(0,400))+scale_color_viridis_d(option="magma", begin=0, end=0.9, name="Mass (kg)")


#AS
preds.RMR$AS <- preds.RMR$MMR.pred - preds.RMR$RMR.pred
preds.RMR$AS <- ifelse(preds.RMR$AS<0, 0, preds.RMR$AS)

ggplot(preds.RMR %>% filter(mass.kg==1))+
  geom_line(aes(temp, RMR.pred, linetype=Sex), col="black")+
  geom_line(aes(temp, MMR.pred, linetype=Sex), col="red")+
  geom_line(aes(temp, AS, linetype=Sex), col="blue")+
  annotate(geom = "label", 9, 200, label="MMR", col="red")+
  annotate(geom = "label", 26, 200, label="RMR", col="black")+
  annotate(geom = "label", 12, 50, label="AS", col="blue")+
  theme_bw()+
  labs(x="Temperature (ºC)", y=bquote(~MO[2]* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), title="Metabolic Metrics", subtitle="A) Metabolism (1 kg fish)")+
  
  ggplot(preds.RMR, aes(temp, mass.kg, fill=RMR.pred))+geom_tile()+
  theme_bw()+labs(x="Temperature (ºC)", y="Fish Mass (kg)", subtitle="B) Resting Metabolic Rate")+
  scale_fill_viridis_c(option="magma", name=bquote(~RMR*' '*O[2]))+scale_x_continuous(limits=c(0,35))+
  
  ggplot(preds.RMR, aes(temp, mass.kg, fill=MMR.pred))+geom_tile()+
  theme_bw()+labs(x="Temperature (ºC)", y="Fish Mass (kg)", subtitle="C) Maximum Metabolic Rate")+
  scale_fill_viridis_c(option="magma", name=bquote(~MMR*' '*O[2]))+scale_x_continuous(limits=c(0,35))+
  
  ggplot(preds.RMR, aes(temp, mass.kg, fill=AS))+geom_tile()+
  theme_bw()+labs(x="Temperature (ºC)", y="Fish Mass (kg)", subtitle="D) Aerobic scope")+
  scale_fill_viridis_c(option="magma", name=bquote(~AS*' '*O[2]))+scale_x_continuous(limits=c(0,35))

