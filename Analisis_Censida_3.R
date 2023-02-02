#Analisis Clusters SENSIDA
#Febrero 2022
#Analisis: Dr. Neftali Antonio Villa


library(tidyverse)
library(readxl)
library(dplyr)
library(mice)
library(haven)
library(epiR)
library(ggpubr)
library(factoextra)
library(FactoMineR)

setwd("/Users/nefoantonio/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/PROYECTOS/JAAF/CENSIDA - JAAF")
base<-read_sav("Base VIH_Fragilidad_TODOS_JUNTOS 05.01.18 MADRE.sav")

#####Recoding#####

base$id<-paste0(base$NUM_REG,"-",base$EdadCalc,"-",base$SEXO)
base<-base%>%filter(!duplicated(id))
base<-base%>%
  left_join(base%>%dplyr::select(id,TALLA)%>%filter(TALLA>=2)%>%
              mutate(TALLA=(as.numeric(TALLA)/100)),by="id")

base$TALLA.x[base$TALLA.x>=2.5]<-NA
base<-base%>%mutate(TALLA.x=coalesce(TALLA.x,TALLA.y))%>%
  dplyr::select(-c(TALLA.y))%>%
  rename("TALLA"="TALLA.x")
base$TALLA[base$TALLA<1]<-NA

base1<-base%>%filter(ZONA!=5)

#Edad 65 Años
base1$EDAD_65<-NULL
base1$EDAD_65[base1$EdadCalc>=65]<-1
base1$EDAD_65[base1$EdadCalc<65]<-0

#Fragilidad
base1$FRAGILIDAD_CRITERIOS[is.na(base1$FRAGILIDAD_CRITERIOS)]<-0

#Barthel
base1$Barthelcat_menor_90<-NULL
base1$Barthelcat_menor_90[base1$BARTHEL_TOTAL<=90]<-1
base1$Barthelcat_menor_90[base1$BARTHEL_TOTAL>90]<-0

#Cualquier Sindrome Discapacidad
base1$SX_DISCAPACIDAD_ANY<-NULL
base1$SX_DISCAPACIDAD_ANY[base1$RBcat==1 | base1$LAWTONCAT==1 | base1$Barthelcat_menor_90==1]<-1
base1$SX_DISCAPACIDAD_ANY[base1$RBcat!=1 & base1$LAWTONCAT!=1 & base1$Barthelcat_menor_90!=1]<-0
base1$SX_DISCAPACIDAD_ANY[is.na(base1$SX_DISCAPACIDAD_ANY)]<-0


#Deficit Cualquier Tipo
base1$DEFICIT_ANY<-NULL
base1$DEFICIT_ANY[base1$DEFICIT_VISUAL==1 | base1$DEFICIT_AUDITIVO==1]<-1
base1$DEFICIT_ANY[base1$DEFICIT_VISUAL!=1 & base1$DEFICIT_AUDITIVO!=1]<-0
base1$DEFICIT_ANY[is.na(base1$DEFICIT_ANY)]<-0

#Depression
base1$DEPRESSION_ANY<-NULL
base1$DEPRESSION_ANY[base1$GDScat==1]<-1
base1$DEPRESSION_ANY[base1$GDScat!=1]<-0
base1$DEPRESSION_ANY[is.na(base1$DEPRESSION_ANY)]<-0

#Falls
base1$CAIDAS_ANY<-NULL
base1$CAIDAS_ANY[base1$caidacat==1]<-1
base1$CAIDAS_ANY[base1$caidacat!=1]<-0
base1$CAIDAS_ANY[is.na(base1$CAIDAS_ANY)]<-0

#Cognitive Decline
base1$COGNITIVE_ANY<-NULL
base1$COGNITIVE_ANY[base1$MMSE<=23 | base1$EIDVIH_IHDScat==1]<-1
base1$COGNITIVE_ANY[base1$MMSE>23 & base1$EIDVIH_IHDScat!=1]<-0
base1$COGNITIVE_ANY[is.na(base1$COGNITIVE_ANY)]<-0

#POLIFARMACIA
base1$POLIFARMACIA_ANY<-NULL
base1$POLIFARMACIA_ANY[base1$NUM_FARMACOS>=3]<-1
base1$POLIFARMACIA_ANY[base1$NUM_FARMACOS<3]<-0
base1$POLIFARMACIA_ANY[is.na(base1$POLIFARMACIA_ANY)]<-0

#COMORBILIDADES
base1$MULTICOMORB_ANY<-NULL
base1$MULTICOMORB_ANY[base1$Comorbilidad>=3]<-1
base1$MULTICOMORB_ANY[base1$Comorbilidad<3]<-0
base1$MULTICOMORB_ANY[is.na(base1$MULTICOMORB_ANY)]<-0

#Incontinencia
base1$INCONTINENCIA<-NULL
base1$INCONTINENCIA[base1$ORINA<5 | base1$HECES<5]<-1
base1$INCONTINENCIA[base1$ORINA>=5 & base1$HECES>=5]<-0
base1$INCONTINENCIA[is.na(base1$INCONTINENCIA)]<-0

#Cummulative Syndromes

base1$SINDROMES_NUM<-as.numeric(base1$SX_DISCAPACIDAD_ANY)+
  as.numeric(base1$DEFICIT_ANY)+
  as.numeric(base1$DEPRESSION_ANY)+
  as.numeric(base1$CAIDAS_ANY)+
  as.numeric(base1$COGNITIVE_ANY)+
  as.numeric(base1$POLIFARMACIA_ANY)+
  as.numeric(base1$FRAGILIDAD_CRITERIOS>=1)+
  as.numeric(base1$INCONTINENCIA)

###----Cummulative and Any type of Geriatric Syndromes######

#Movilidad, Act. Instrumentadas, Act. Basicas
ncas <- table(base1$SX_DISCAPACIDAD_ANY)[2]; npop <- sum(!is.na(base1$SX_DISCAPACIDAD_ANY))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Sensorial Deficit
ncas <- table(base1$DEFICIT_ANY)[2]; npop <- sum(!is.na(base1$DEFICIT_ANY))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Depression
ncas <- table(base1$DEPRESSION_ANY)[2]; npop <- sum(!is.na(base1$DEPRESSION_ANY))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Falls
ncas <- table(base1$CAIDAS_ANY)[2]; npop <- sum(!is.na(base1$CAIDAS_ANY))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Cognitive Decline
ncas <- table(base1$COGNITIVE_ANY)[2]; npop <- sum(!is.na(base1$COGNITIVE_ANY))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev5<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Polifarmacia
ncas <- table(base1$POLIFARMACIA_ANY)[2]; npop <- sum(!is.na(base1$POLIFARMACIA_ANY))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev6<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Multicomorbidities
ncas <- table(base1$MULTICOMORB_ANY)[2]; npop <- sum(!is.na(base1$MULTICOMORB_ANY))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev7<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Fraility
ncas <- table(base1$FRAGILIDAD_CRITERIOS>=1)[2]; npop <- sum(!is.na(base1$FRAGILIDAD_CRITERIOS>=1))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev8<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Incontinencia
ncas <- table(base1$INCONTINENCIA)[2]; npop <- sum(!is.na(base1$INCONTINENCIA))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev9<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

sindromes.ger.df<-rbind(prev1,prev2,prev3,prev4,prev5,prev6,prev8,prev9)
sindromes.ger.df$group<-c("Physical Disability","Sensorial Deficit","Depression","Falls","Cognitive Decline","Polymedication","PreFraility or Fraility","Incontinency")
sindromes.ger.df$group<-factor(sindromes.ger.df$group,levels = c("Polymedication","Sensorial Deficit","Cognitive Decline","Physical Disability","PreFraility or Fraility","Falls","Depression","Incontinency"))

###Cummulative number of Syndromes
#SINDROMES
#>5
ncas <- table(base1$SINDROMES_NUM>=5)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM>=5))
tmp <- as.matrix(cbind(ncas, npop))
prev_MAS5<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                          conf.level = 0.95) * 100,2)


#0
ncas <- table(base1$SINDROMES_NUM==0)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==0))
tmp <- as.matrix(cbind(ncas, npop))
prev0<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#1
ncas <- table(base1$SINDROMES_NUM==1)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==1))
tmp <- as.matrix(cbind(ncas, npop))
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#2
ncas <- table(base1$SINDROMES_NUM==2)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==2))
tmp <- as.matrix(cbind(ncas, npop))
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#3
ncas <- table(base1$SINDROMES_NUM==3)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==3))
tmp <- as.matrix(cbind(ncas, npop))
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#4
ncas <- table(base1$SINDROMES_NUM==4)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==4))
tmp <- as.matrix(cbind(ncas, npop))
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#5
ncas <- table(base1$SINDROMES_NUM==5)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==5))
tmp <- as.matrix(cbind(ncas, npop))
prev5<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#6
ncas <- table(base1$SINDROMES_NUM==6)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==6))
tmp <- as.matrix(cbind(ncas, npop))
prev6<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#7
ncas <- table(base1$SINDROMES_NUM==7)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==7))
tmp <- as.matrix(cbind(ncas, npop))
prev7<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#8
ncas <- table(base1$SINDROMES_NUM==8)[2]; npop <- sum(!is.na(base1$SINDROMES_NUM==8))
tmp <- as.matrix(cbind(ncas, npop))
prev8<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

sindromes.cum.df<-rbind(prev0,prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)
sindromes.cum.df$group<-c(0,"1","2","3","4","5","6","7","8")

###----Description of Geriatric Syndromes------

#Discapacidad en Movilidad
ncas <- table(base1$RBcat)[2]; npop <- sum(!is.na(base1$RBcat))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Discapacidad Act. Instrumentadas

ncas <- table(base1$LAWTONCAT)[2]; npop <- sum(!is.na(base1$LAWTONCAT))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Discapacidad Act. Basicas

ncas <- table(base1$Barthelcat_menor_85)[2]; npop <- sum(!is.na(base1$Barthelcat_menor_85))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

discapacidad.df<-rbind(prev2,prev3,prev4)
discapacidad.df$group<-c("Mobility","Instrumental Activities","Basic Activities")
discapacidad.df$group<-factor(discapacidad.df$group,levels = c("Mobility","Instrumental Activities","Basic Activities"))
#Discapacidad Sensorial
#Visual

ncas <- table(base1$DEFICIT_VISUAL)[2]; npop <- sum(!is.na(base1$DEFICIT_VISUAL))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Audititva

ncas <- table(base1$DEFICIT_AUDITIVO)[2]; npop <- sum(!is.na(base1$DEFICIT_AUDITIVO))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prop.single", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)[,c(1,3,4)]

deficit.sens.df<-rbind(prev1,prev2)
deficit.sens.df$group<-c("Visual","Hearing")
deficit.sens.df$group<-factor(deficit.sens.df$group,levels = c("Visual","Hearing"))
#Depresivo
#GDS
ncas <- table(base1$GDScat)[2]; npop <- sum(!is.na(base1$GDScat))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Self-Reported
ncas <- table(base1$TDM)[2]; npop <- sum(!is.na(base1$TDM))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prop.single", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)[,c(1,3,4)]

depression.df<-rbind(prev1,prev2)
depression.df$group<-c("GDS ≥6 pts","Self-Reported")


#Caidas
ncas <- table(base1$caidacat)[2]; npop <- sum(!is.na(base1$caidacat))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Caidas ≥3
ncas <- table(base1$CAIDAS>=3)[2]; npop <- sum(!is.na(base1$CAIDAS>=3))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

caidad.df<-rbind(prev1,prev2)
caidad.df$group<-c("Any Reported Fall","≥3 Falls")
caidad.df$group<-factor(caidad.df$group,levels = c("Any Reported Fall","≥3 Falls"))

#Fx Cognitiva
ncas <- table(base1$MMSE<=23)[2]; npop <- sum(!is.na(base1$MMSE<=23))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Escala Internacional Demencia
ncas <- table(base1$EIDVIH_IHDScat)[2]; npop <- sum(!is.na(base1$EIDVIH_IHDScat))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

det.fx.cognitica.df<-rbind(prev1,prev2)
det.fx.cognitica.df$group<-c("MMSE ≤ 23","IHDS ≤ 10")

#Polifarmacia
#>3 farmacos
ncas <- table(base1$NUM_FARMACOS>=3)[2]; npop <- sum(!is.na(base1$NUM_FARMACOS>=3))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#>5 farmacos
ncas <- table(base1$NUM_FARMACOS>=5)[2]; npop <- sum(!is.na(base1$NUM_FARMACOS>=5))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

polifarmacia.df<-rbind(prev1,prev2)
polifarmacia.df$group<-c("≥3 Medications", "≥5 Medications")

#Comorbilidades
#≥1 comorbilidad
ncas <- table(base1$Comorbilidad>=1)[2]; npop <- sum(!is.na(base1$Comorbilidad>=1))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#≥3 comorbilidad
ncas <- table(base1$Comorbilidad>=3)[2]; npop <- sum(!is.na(base1$Comorbilidad>=3))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#≥5 comorbilidad
ncas <- table(base1$Comorbilidad>=5)[2]; npop <- sum(!is.na(base1$Comorbilidad>=5))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

comorbidities.df<-rbind(prev1,prev2,prev3)
comorbidities.df$group<-c("≥1","≥3","≥5")

#Fragilidad
#Componentes
#Perdida de Peso
ncas <- table(base1$DISMINUCION_PESO_KG)[2]; npop <- sum(!is.na(base1$DISMINUCION_PESO_KG))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Extenuacion
ncas <- table(base1$EXTENUACION)[2]; npop <- sum(!is.na(base1$EXTENUACION))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Prension
ncas <- table(base1$PRENSION)[2]; npop <- sum(!is.na(base1$PRENSION))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Marcha
ncas <- table(base1$MARCHA)[2]; npop <- sum(!is.na(base1$MARCHA))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Act. Fisica
ncas <- table(base1$ACTIVIDAD_FISICA)[2]; npop <- sum(!is.na(base1$ACTIVIDAD_FISICA))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev5<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

fragilidad.df<-rbind(prev1,prev2,prev3,prev4,prev5)
fragilidad.df$group<-c("Weight Loss","Exhaustion","Weakness","Slowness","Low Physical Activity")
fragilidad.df$group<-factor(fragilidad.df$group,levels = c("Exhaustion","Low Physical Activity","Weakness","Slowness","Weight Loss"))

#Incontinencia
#Urinaria
ncas <- table(base1$ORINA<5)[2]; npop <- sum(!is.na(base1$ORINA<5))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

#Heces 
ncas <- table(base1$HECES<5)[2]; npop <- sum(!is.na(base1$ACTIVIDAD_FISICA<5))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "exact", N = 500, design = 1, 
                      conf.level = 0.95) * 100,2)

incontinencia.df<-rbind(prev1,prev2)
incontinencia.df$group<-c("Urine","Feces")
incontinencia.df$group<-factor(incontinencia.df$group,levels = c("Urine","Feces"))
###----Poisson Regression Models-----

#Modelo Poisson
mod0<-glm(SINDROMES_NUM~scale(CD4_NADIR)+scale(CD4_ACTUAL)+scale(RNA_VIH_ACTUAL)+scale(EdadCalc)+SEXO+ESCOL.AÑOS,data=base1,family = "poisson")
summary(mod0)

#Escala sin Edad
base1$SINDROMES_NUM_SIN_EDAD_2<-glm(SINDROMES_NUM~EdadCalc,base1,family = "poisson")$residuals

#Normalizacion de Variables
base1$CD4_NADIR_NORM<-bestNormalize::orderNorm(base1$CD4_NADIR)$x.t
base1$CD4_ACTUAL_NORM<-bestNormalize::orderNorm(base1$CD4_ACTUAL)$x.t
base1$RNA_VIRAL_ACTUAL_NORM<-bestNormalize::orderNorm(base1$RNA_VIH_ACTUAL)$x.t

#Modelo de Regresion Lineal
mod1<-glm(SINDROMES_NUM_SIN_EDAD_2~CD4_NADIR_NORM+CD4_ACTUAL_NORM+RNA_VIH_ACTUAL+SEXO+ESCOL.AÑOS,data=base1)
options("scipen"=100, "digits"=4)
summary(mod1)
confint(mod1)

###----Cluster Analysis#####
#Depuration of Dataset
base3 <- base1 %>% dplyr::select(id,CD4_NADIR_NORM,CD4_ACTUAL_NORM,SINDROMES_NUM_SIN_EDAD_2)
base3.1<-base3%>%dplyr::select(-id)
base3.1 <- scale(base3.1)
base3.1 <- as.data.frame(base3.1)
base3.1_imp<-mice::mice(base3.1, m=5, maxit=5,seed = 123)
base3.1<-complete(base3.1_imp,1)

#Cluster Vizualization
set.seed(123)
k1<-kmeans(base3.1, 3, nstart = 50,iter.max = 15)
print(k1)
base3<-cbind(base3, cluster = k1$cluster)
base1_CLUSTER<-NULL
base1_CLUSTER<-base1%>%left_join(base3%>%dplyr::select(id,cluster),by = "id")
base1_CLUSTER$cluster<-factor(base1_CLUSTER$cluster,levels = c(2,1,3))
#fpc::clusterboot(base3.1,B=1000,bootmethod = "boot",clustermethod = fpc::kmeansCBI,k=3)

###----Figure 1----

Figure1A<-ggplot(sindromes.ger.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ 
  geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.0),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(-.1,80))+
  labs(fill="Type")+
  ggtitle("Geriatric Syndromes")

Figure1B<-ggplot(sindromes.cum.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.0),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(-.1,40))+
  labs(fill="")+
  ggtitle("Cumulative Syndromes")

Figure1<-ggarrange(Figure1A,Figure1B,ncol = 2,nrow = 1,labels = LETTERS[1:2])
ggsave(Figure1,
       filename = "Figure1.png", 
       width = 55, 
       height = 23,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)

###----Figure 2-----
cor.test(base1$SINDROMES_NUM_SIN_EDAD_2,base1$CD4_NADIR_NORM)
cor.test(base1$SINDROMES_NUM_SIN_EDAD_2,base1$CD4_ACTUAL_NORM)

Figure2A<-ggplot(base1,aes(CD4_NADIR_NORM,SINDROMES_NUM_SIN_EDAD_2))+
  geom_point(size=2)+
  stat_smooth(method = lm, formula = y ~ x, raw = TRUE, color="red", size=1)+
  ylab("Age-Independent \nCumulative Geriatric Syndromes")+
  xlab("CD4-Nadir (cells/μl) \nNormalized")+
  theme_pubclean()+
  annotate("text", x = 2, y = 2.4, label = "r = -0.126")+
  annotate("text", x = 2, y = 2.2, label = "95% CI: -0.223 to -0.026")+
  annotate("text", x = 2, y = 2.0, label = "p<0.05")+
  theme(
    axis.title.x = element_text(color="black", size=13, face="plain"),
    axis.title.y = element_text(color="black", size=13, face="plain"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"))

Figure2B<-ggplot(base1,aes(CD4_ACTUAL_NORM,SINDROMES_NUM_SIN_EDAD_2))+
  geom_point(size=2)+
  stat_smooth(method = lm, formula = y ~ x, raw = TRUE, color="red", size=1)+
  ylab("Age-Independent \nCumulative Geriatric Syndromes")+
  xlab("Last-CD4 Count (cells/μl) \nNormalized")+
  theme_pubclean()+
  annotate("text", x = 2, y = 2.4, label = "r = 0.8")+
  annotate("text", x = 2, y = 2.2, label = "95% CI: -0.080 to 0.098")+
  annotate("text", x = 2, y = 2.0, label = "p=0.800")+
  theme(
    axis.title.x = element_text(color="black", size=13, face="plain"),
    axis.title.y = element_text(color="black", size=13, face="plain"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"))

Figure2<-ggarrange(Figure2A,Figure2B,ncol = 2,nrow = 1,labels = LETTERS[1:2])

ggsave(Figure2,
       filename = "Figure_2.png", 
       bg = "white",
       width = 45, 
       height = 15,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)


###----Figure 3-----
my_comparisons <- list( c("1", "2"), c("1", "3"),
                        c("2", "3"))
#Characterization of Clustering
Figure3A<-fviz_cluster(k1, data=base3.1, ellipse.type="convex", 
                       palette="jama", ggtheme=theme_pubclean())+ theme(legend.position = "none")

#Boxplots


Figure3B.1<-ggplot(base1_CLUSTER,aes(factor(cluster),CD4_NADIR_NORM,fill=factor(cluster,labels=c("1","2","3"))))+
  geom_boxplot()+
  stat_compare_means(method="kruskal.test",label.y = 5.3)+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  ggsci::scale_fill_jama()+
  labs(fill="Cluster")+
  scale_x_discrete(labels = c('','',""))+ theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")

Figure3B.2<-ggplot(base1_CLUSTER,aes(factor(cluster),CD4_ACTUAL_NORM,fill=factor(cluster,labels=c("1","2","3"))))+
  geom_boxplot()+
  stat_compare_means(method="kruskal.test",label.y = 5.3)+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  ggsci::scale_fill_jama()+
  labs(fill="Cluster")+
  scale_x_discrete(labels = c('','',""))+ theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")

Figure3B.3<-ggplot(base1_CLUSTER,aes(factor(cluster),SINDROMES_NUM_SIN_EDAD_2,fill=factor(cluster,labels=c("1","2","3"))))+
  geom_boxplot()+
  stat_compare_means(method="kruskal.test",label.y = 5.3)+
  theme_pubclean()+
  xlab("")+
  ylab("Age-Independent \nCumulative Geriatric Syndromes")+
  ggsci::scale_fill_jama()+
  labs(fill="Cluster")+
  scale_x_discrete(labels = c('','',""))+ theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")

Figure3B<-ggarrange(Figure3B.1,Figure3B.2,Figure3B.3,ncol = 3,nrow = 1,labels = c("B","",""),common.legend = F)
Figure3<-ggarrange(Figure3A,Figure3B,ncol = 2,nrow = 1,labels = c("A",""),common.legend = T)
ggsave(Figure3,
       filename = "Figure_3.png", 
       bg = "white",
       width = 45, 
       height = 15,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)

###----Supplementary Figure 1------

Figure2A<-ggplot(discapacidad.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.0),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(-.1,50))+
  labs(fill="Type")+
  ggtitle("Disability")

Figure2B<-ggplot(deficit.sens.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.0),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(-.1,100))+
  labs(fill="Deficit")+
  ggtitle("Sensorial Deficit")

Figure2C<-ggplot(depression.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.0),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(0,20))+
  labs(fill="Method")+
  ggtitle("Depression")


Figure2D<-ggplot(caidad.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.0),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(0,50))+
  labs(fill="Method")+
  ggtitle("Falls")

Figure2E<-ggplot(det.fx.cognitica.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.5),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(0,70))+
  labs(fill="Method")+
  ggtitle("Cognitive Decline")

Figure2F<-ggplot(polifarmacia.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.5),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(0,100))+
  labs(fill="Method")+
  ggtitle("Polypharmacy")

Figure2G<-ggplot(comorbidities.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.5),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(0,75))+
  labs(fill="Cummulative")+
  ggtitle("Self-Reported Comorbidities")

Figure2H<-ggplot(fragilidad.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 5.5),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(0,40))+
  labs(fill="Component")+
  ggtitle("Fraility")+
  guides(fill = guide_legend(nrow = 2,byrow = TRUE))


Figure2I<-ggplot(incontinencia.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Paired")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 3.5),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(0,10))+
  labs(fill="Type")+
  ggtitle("Incontinency")


Figure2_UP<-ggarrange(Figure2A,Figure2B,Figure2C,ncol = 3,nrow = 1,labels = LETTERS[1:3])
Figure2_MID<-ggarrange(Figure2D,Figure2E,Figure2F,ncol = 3,nrow = 1,labels = LETTERS[4:6])
Figure2_DOWN<-ggarrange(Figure2G,Figure2H,Figure2I,ncol = 3,nrow = 1,labels = LETTERS[7:9])
Figure2<-ggarrange(Figure2_UP,Figure2_MID,Figure2_DOWN,ncol = 1,nrow = 3)

ggsave(Figure2,
       filename = "Supplementary_Figure_1.png", 
       width = 45, 
       height = 40,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)

###----Supplementary Figure 2-----
base2<-base1%>%dplyr::select(POLIFARMACIA_ANY,DEFICIT_ANY,COGNITIVE_ANY,
                             SX_DISCAPACIDAD_ANY,FRAGILIDAD_CRITERIOS,
                             CAIDAS_ANY,DEPRESSION_ANY,MULTICOMORB_ANY,INCONTINENCIA)
names(base2)<-c("Polypharmacy","Sensorial Deficit","Cognitive Decline",
                "Physical Disability","PreFraility or Fraility","Falls",
                "Depression","Multi-Comorbidities","Incontinency")

Sup_Figure2A.1<-ggplot(base1,aes(factor(POLIFARMACIA_ANY),CD4_NADIR_NORM,fill=factor(POLIFARMACIA_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Polypharmacy")

Sup_Figure2A.2<-ggplot(base1,aes(factor(POLIFARMACIA_ANY),CD4_ACTUAL_NORM,fill=factor(POLIFARMACIA_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Polypharmacy")

Sup_Figure2B.1<-ggplot(base1,aes(factor(DEFICIT_ANY),CD4_NADIR_NORM,fill=factor(DEFICIT_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Sensorial Deficit")

Sup_Figure2B.2<-ggplot(base1,aes(factor(DEFICIT_ANY),CD4_ACTUAL_NORM,fill=factor(DEFICIT_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Sensorial Deficit")

Sup_Figure2C.1<-ggplot(base1,aes(factor(COGNITIVE_ANY),CD4_NADIR_NORM,fill=factor(COGNITIVE_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Cognitive Decline")

Sup_Figure2C.2<-ggplot(base1,aes(factor(COGNITIVE_ANY),CD4_ACTUAL_NORM,fill=factor(COGNITIVE_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Cognitive Decline")

Sup_Figure2D.1<-ggplot(base1,aes(factor(SX_DISCAPACIDAD_ANY),CD4_NADIR_NORM,fill=factor(SX_DISCAPACIDAD_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Physical Disability")

Sup_Figure2D.2<-ggplot(base1,aes(factor(SX_DISCAPACIDAD_ANY),CD4_ACTUAL_NORM,fill=factor(SX_DISCAPACIDAD_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Physical Disability")

Sup_Figure2E.1<-ggplot(base1,aes(factor(FRAGILIDAD_CRITERIOS>=1),CD4_NADIR_NORM,fill=factor(FRAGILIDAD_CRITERIOS>=1,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("PreFraility or Fraility")

Sup_Figure2E.2<-ggplot(base1,aes(factor(FRAGILIDAD_CRITERIOS>=1),CD4_ACTUAL_NORM,fill=factor(FRAGILIDAD_CRITERIOS>=1,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("PreFraility or Fraility")

Sup_Figure2F.1<-ggplot(base1,aes(factor(CAIDAS_ANY),CD4_NADIR_NORM,fill=factor(CAIDAS_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Falls")

Sup_Figure2F.2<-ggplot(base1,aes(factor(CAIDAS_ANY),CD4_ACTUAL_NORM,fill=factor(CAIDAS_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Falls")

Sup_Figure2G.1<-ggplot(base1,aes(factor(DEPRESSION_ANY),CD4_NADIR_NORM,fill=factor(DEPRESSION_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Depression")

Sup_Figure2G.2<-ggplot(base1,aes(factor(DEPRESSION_ANY),CD4_ACTUAL_NORM,fill=factor(DEPRESSION_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Depression")

Sup_Figure2H.1<-ggplot(base1,aes(factor(MULTICOMORB_ANY),CD4_NADIR_NORM,fill=factor(MULTICOMORB_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Multi-Comorbidities")

Sup_Figure2H.2<-ggplot(base1,aes(factor(MULTICOMORB_ANY),CD4_ACTUAL_NORM,fill=factor(MULTICOMORB_ANY,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Multi-Comorbidities")

Sup_Figure2I.1<-ggplot(base1,aes(factor(INCONTINENCIA),CD4_NADIR_NORM,fill=factor(INCONTINENCIA,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Nadir (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Incontinency")

Sup_Figure2I.2<-ggplot(base1,aes(factor(INCONTINENCIA),CD4_ACTUAL_NORM,fill=factor(INCONTINENCIA,labels=c("No","Yes"))))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_pubclean()+
  xlab("")+
  ylab("CD4-Actual (cells/μl) \nNormalized")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill="")+
  scale_x_discrete(labels = c('',''))+
  ggtitle("Incontinency")

Sup_Figure2A<-ggarrange(Sup_Figure2A.1,Sup_Figure2A.2,ncol = 2,nrow = 1,labels = c("A",""),common.legend = T)
Sup_Figure2B<-ggarrange(Sup_Figure2B.1,Sup_Figure2B.2,ncol = 2,nrow = 1,labels = c("B",""),common.legend = T)
Sup_Figure2C<-ggarrange(Sup_Figure2C.1,Sup_Figure2C.2,ncol = 2,nrow = 1,labels = c("C",""),common.legend = T)
Sup_Figure2D<-ggarrange(Sup_Figure2D.1,Sup_Figure2D.2,ncol = 2,nrow = 1,labels = c("D",""),common.legend = T)
Sup_Figure2E<-ggarrange(Sup_Figure2E.1,Sup_Figure2E.2,ncol = 2,nrow = 1,labels = c("E",""),common.legend = T)
Sup_Figure2F<-ggarrange(Sup_Figure2F.1,Sup_Figure2F.2,ncol = 2,nrow = 1,labels = c("F",""),common.legend = T)
Sup_Figure2G<-ggarrange(Sup_Figure2G.1,Sup_Figure2G.2,ncol = 2,nrow = 1,labels = c("G",""),common.legend = T)
Sup_Figure2H<-ggarrange(Sup_Figure2H.1,Sup_Figure2H.2,ncol = 2,nrow = 1,labels = c("H",""),common.legend = T)
Sup_Figure2I<-ggarrange(Sup_Figure2I.1,Sup_Figure2I.2,ncol = 2,nrow = 1,labels = c("I",""),common.legend = T)
Sup_Figure2<-ggarrange(Sup_Figure2A,Sup_Figure2B,Sup_Figure2C,
                       Sup_Figure2D,Sup_Figure2E,Sup_Figure2F,
                       Sup_Figure2G,Sup_Figure2H,Sup_Figure2I,
                       ncol = 3,nrow = 3,common.legend = T)

ggsave(Sup_Figure2,
       filename = "Supplementary_Figure2.png", 
       bg = "white",
       width = 45, 
       height = 30,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)

###----Supplementary Figure 3----
#Evaluar Asumsiones
mod0<-glm(SINDROMES_NUM~EdadCalc,base1,family = "poisson")
summary(mod0)
sum(residuals(mod0, type = "pearson")^2) / mod0$df.residual
1-pchisq(mod0$deviance,499)


mod0.1<-MASS::glm.nb(SINDROMES_NUM~EdadCalc,base1)
require(pscl)
odTest(mod0.1)

plot(mod0)
lmtest::bptest(mod0)
hist(mod0$residuals)
nortest::ad.test(mod0$residuals)

par(mfrow = c(2, 3))
hist(mod0$residuals,xlab = "Residuals",main="Residuals from \nPoisson Regression Model")
plot(mod0)
plot(mod0,4)


###----Supplementary Figure 4-----
#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
data <- scale(base3.1)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})

par(mfrow=c(1, 3))
#Elbow Plot
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

#Calinski-Harabasz Criterion
library(fpc)
k2.1 <- fpc::kmeansruns(data,krange = 1:10, criterion = "ch")
k2.2 <- fpc::kmeansruns(data,krange = 1:10, criterion = "asw")
plot(1:10,k2.1$crit,type = "b",xlab="Number of Clusters", ylab="Calinski-Harabasz Criterion")
plot(1:10,k2.2$crit,type = "b",xlab="Number of Clusters", ylab="Average Silhouette Wight")



###----Descriptive characteristics (Table 1)####

columns <- c('Parameter',"All-Population (n=501)")

## Sex Male
sexo <- table(base1$SEXO,useNA = "always")[c(1)]
sexoprop <- round(prop.table(table(base1$SEXO,useNA = "always")),4)[c(1)]*100

Sexo<-`names<-`(data.frame("Men (%)",
                           matrix(c(paste(sexo,paste0('(',sexoprop,')'))),ncol = 1)),
                columns)

#Age
nortest::ad.test(base1$EdadCalc)
num1<-c(paste(round(median(base1$EdadCalc,na.rm = T ),4),
              paste0('(',round(quantile(base1$EdadCalc,na.rm = T,probs = c(0.25)),3),"-",
                     round(quantile(base1$EdadCalc,na.rm = T,probs = c(0.75)),3),')')))

Edad<-`names<-`(data.frame(matrix(c("Age (Years)",num1),ncol = 2)),columns)

## Education

nortest::ad.test(base1$ESCOL.AÑOS)
num1<-c(paste(round(median(base1$ESCOL.AÑOS,na.rm = T ),4),
              paste0('(',round(quantile(base1$ESCOL.AÑOS,na.rm = T,probs = c(0.25)),3),"-",
                     round(quantile(base1$ESCOL.AÑOS,na.rm = T,probs = c(0.75)),3),')')))

Educacion<-`names<-`(data.frame(matrix(c("Education (Years)",num1),ncol = 2)),columns)

#Comorbilidades
#Acumuladas
nortest::ad.test(base1$Comorbilidad)
num1<-c(paste(round(median(base1$Comorbilidad,na.rm = T ),4),
              paste0('(',round(quantile(base1$Comorbilidad,na.rm = T,probs = c(0.25)),3),"-",
                     round(quantile(base1$Comorbilidad,na.rm = T,probs = c(0.75)),3),')')))

Comorbilidades<-`names<-`(data.frame(matrix(c("Cummulative Comorbidities (Count)",num1),ncol = 2)),columns)

#Comorbilidades Categorias 
comor <- table(base1$comorbilidad2,useNA = "always")[c(1,2,3)]
comorprop <- round(prop.table(table(base1$comorbilidad2,useNA = "always")),4)[c(1,2,3)]*100

Comorbilidades_2<-`names<-`(data.frame(c("No-Comorbities (%)","1 Comorbidity (%)","≥2 Comorbities (%)"),
                                       matrix(c(paste(comor,paste0('(',comorprop,')'))),ncol = 1)),
                            columns)

#Diabetes 
diab <- table(base1$DM2,useNA = "always")[c(2)]
diabprop <- round(prop.table(table(base1$DM2,useNA = "always")),4)[c(2)]*100

Diabetes<-`names<-`(data.frame("Diabetes (%)",
                               matrix(c(paste(diab,paste0('(',diabprop,')'))),ncol = 1)),
                    columns)

#Hipertension 
has <- table(base1$HAS,useNA = "always")[c(2)]
hasprop <- round(prop.table(table(base1$HAS,useNA = "always")),4)[c(2)]*100

Hipertension<-`names<-`(data.frame("Arterial Hypertension (%)",
                                   matrix(c(paste(has,paste0('(',hasprop,')'))),ncol = 1)),
                        columns)

#Dislipidemia 
disli <- table(base1$DISLIPIDEMIA,useNA = "always")[c(2)]
disliprop <- round(prop.table(table(base1$DISLIPIDEMIA,useNA = "always")),4)[c(2)]*100

Dislipidemia<-`names<-`(data.frame("Dyslipidemia (%)",
                                   matrix(c(paste(disli,paste0('(',disliprop,')'))),ncol = 1)),
                        columns)

#Cancer 
cancer <- table(base1$CANCER_NOSIDA,useNA = "always")[c(2)]
cancerprop <- round(prop.table(table(base1$CANCER_NOSIDA,useNA = "always")),4)[c(2)]*100

Cancer_NOSIDA<-`names<-`(data.frame("Cancer (%)",
                                    matrix(c(paste(cancer,paste0('(',cancerprop,')'))),ncol = 1)),
                         columns)

#Isquemia 
EVC <- table(base1$EVC,useNA = "always")[c(2)]
EVCprop <- round(prop.table(table(base1$EVC,useNA = "always")),4)[c(2)]*100

EVC<-`names<-`(data.frame("Cardiovascular Disease (%)",
                          matrix(c(paste(EVC,paste0('(',EVCprop,')'))),ncol = 1)),
               columns)

#EPOC 
EPOC <- table(base1$EPOC,useNA = "always")[c(2)]
EPOCprop <- round(prop.table(table(base1$EPOC,useNA = "always")),4)[c(2)]*100

EPOC<-`names<-`(data.frame("COPD (%)",
                           matrix(c(paste(EPOC,paste0('(',EPOCprop,')'))),ncol = 1)),
                columns)

#Hepatopatia Crónica 
hepa <- table(base1$HEPAT_CRON,useNA = "always")[c(2)]
hepa.prop <- round(prop.table(table(base1$HEPAT_CRON,useNA = "always")),4)[c(2)]*100

Hepatopatia.Cronica<-`names<-`(data.frame("Chronic Hepatopathy (%)",
                                          matrix(c(paste(hepa,paste0('(',hepa.prop,')'))),ncol = 1)),
                               columns)

#OAD
oad <- table(base1$OAD,useNA = "always")[c(2)]
oad.prop <- round(prop.table(table(base1$HEPAT_CRON,useNA = "always")),4)[c(2)]*100

OAD<-`names<-`(data.frame("OAD (%)",
                          matrix(c(paste(oad,paste0('(',oad.prop,')'))),ncol = 1)),
               columns)

#Osteoartritis
osteo <- table(base1$OSTEOP,useNA = "always")[c(2)]
osteo.prop <- round(prop.table(table(base1$OSTEOP,useNA = "always")),4)[c(2)]*100

Osteoartritis<-`names<-`(data.frame("Osteoarthritis (%)",
                                    matrix(c(paste(osteo,paste0('(',osteo.prop,')'))),ncol = 1)),
                         columns)

#CKD
ckd <- table(base1$RENAL_CRON,useNA = "always")[c(2)]
ckd.prop <- round(prop.table(table(base1$RENAL_CRON,useNA = "always")),4)[c(2)]*100

CKD<-`names<-`(data.frame("Chronic Kidney Disease (%)",
                          matrix(c(paste(ckd,paste0('(',ckd.prop,')'))),ncol = 1)),
               columns)

#Depresion

depre <- table(base1$TDM,useNA = "always")[c(2)]
depre.prop <- round(prop.table(table(base1$TDM,useNA = "always")),4)[c(2)]*100

Depresion<-`names<-`(data.frame("Depression (%)",
                                matrix(c(paste(depre,paste0('(',depre.prop,')'))),ncol = 1)),
                     columns)

#Ansiedad
ansie <- table(base1$TAD,useNA = "always")[c(2)]
ansie.prop <- round(prop.table(table(base1$TAD,useNA = "always")),4)[c(2)]*100

Ansiedad<-`names<-`(data.frame("Anxiety (%)",
                               matrix(c(paste(ansie,paste0('(',ansie.prop,')'))),ncol = 1)),
                    columns)

#Peso
Peso.media<-round(mean(base1$PESO,na.rm = T),2)
Peso.sd<-round(sd(base1$PESO,na.rm = T),2)
PESO<-`names<-`(data.frame(matrix(c("Weight (Kg)",paste(Peso.media,"(","±",Peso.sd,")")),ncol = 2,byrow = T)),columns)

#Altura

talla.media<-round(mean(base1$TALLA,na.rm = T),2)
talla.sd<-round(sd(base1$TALLA,na.rm = T),2)
TALLA<-`names<-`(data.frame(matrix(c("Height (mts)",paste(talla.media,"(","±",talla.sd,")")),ncol = 2,byrow = T)),columns)

#IMC

nortest::ad.test(base1$IMC)
num1<-c(paste(round(median(base1$IMC,na.rm = T ),2),
              paste0('(',round(quantile(base1$IMC,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$IMC,na.rm = T,probs = c(0.75)),2),')')))

IMC<-`names<-`(data.frame(matrix(c("BMI (Kg/mts^2)",num1),ncol = 2)),columns)

#Hemoglobina
nortest::ad.test(base1$HEMOGLOB)
num1<-c(paste(round(median(base1$HEMOGLOB,na.rm = T ),2),
              paste0('(',round(quantile(base1$HEMOGLOB,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$HEMOGLOB,na.rm = T,probs = c(0.75)),2),')')))

Hemoglob<-`names<-`(data.frame(matrix(c("Hemoglobin (grs/ml)",num1),ncol = 2)),columns)

#Glucosa
nortest::ad.test(base1$GLUCOSA)
num1<-c(paste(round(median(base1$GLUCOSA,na.rm = T ),2),
              paste0('(',round(quantile(base1$GLUCOSA,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$GLUCOSA,na.rm = T,probs = c(0.75)),2),')')))

Gluco<-`names<-`(data.frame(matrix(c("Glucose (mg/ml)",num1),ncol = 2)),columns)

#Albumin

num1<-c(paste(round(median(base1$ALBUMINA,na.rm = T ),2),
              paste0('(',round(quantile(base1$ALBUMINA,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$ALBUMINA,na.rm = T,probs = c(0.75)),2),')')))

Albumin<-`names<-`(data.frame(matrix(c("Albumin (mg/ml)",num1),ncol = 2)),columns)

#Barthel

bartel <- table(base1$Barthelcat_menor_85,useNA = "always")[c(2)]
bartelprop <- round(prop.table(table(base1$Barthelcat_menor_85,useNA = "always")),4)[c(2)]*100

Barthel<-`names<-`(data.frame("Barthel-Scale <100pts (%)",
                              matrix(c(paste(bartel,paste0('(',bartelprop,')'))),ncol = 1)),
                   columns)

#Lawton

laton <- table(base1$LAWTONCAT,useNA = "always")[c(2)]
latonprop <- round(prop.table(table(base1$LAWTONCAT,useNA = "always")),4)[c(2)]*100

Lawton<-`names<-`(data.frame("Lawton-Scale Abnormal(%)",
                             matrix(c(paste(laton,paste0('(',latonprop,')'))),ncol = 1)),
                  columns)


#Caidas

caida <- table(base1$caidacat,useNA = "always")[c(2)]
caidaprop <- round(prop.table(table(base1$caidacat,useNA = "always")),4)[c(2)]*100

Caidas<-`names<-`(data.frame("History of ≥1 Fall (%)",
                             matrix(c(paste(caida,paste0('(',caidaprop,')'))),ncol = 1)),
                  columns)

#Rosow-Brews

rosow <- table(base1$ROSOW_BRESLAU_TOTAL,useNA = "always")[c(1:4)]
rosowprop <- round(prop.table(table(base1$ROSOW_BRESLAU_TOTAL,useNA = "always")),4)[c(1:4)]*100

Rosow<-`names<-`(data.frame(c("RB-0 pts (%)","RB-1 pts (%)","RB-2 pts (%)","RB-3 pts (%)"),
                            matrix(c(paste(rosow,paste0('(',rosowprop,')'))),ncol = 1)),
                 columns)

#Visual Deficit

def.visual <- table(base1$DEFICIT_VISUAL,useNA = "always")[c(2)]
def.visualprop <- round(prop.table(table(base1$DEFICIT_VISUAL,useNA = "always")),4)[c(2)]*100

Def.Visual<-`names<-`(data.frame(c("Visual Deficit (%)"),
                                 matrix(c(paste(def.visual,paste0('(',def.visualprop,')'))),ncol = 1)),
                      columns)


#Auditivo Deficit

def.auditivo <- table(base1$DEFICIT_AUDITIVO,useNA = "always")[c(2)]
def.auditivoprop <- round(prop.table(table(base1$DEFICIT_AUDITIVO,useNA = "always")),4)[c(2)]*100

Def.Auditivo<-`names<-`(data.frame(c("Hearing Loss (%)"),
                                   matrix(c(paste(def.auditivo,paste0('(',def.auditivoprop,')'))),ncol = 1)),
                        columns)

#Escala de Depresión

num1<-c(paste(round(median(base1$GDS,na.rm = T ),2),
              paste0('(',round(quantile(base1$GDS,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$GDS,na.rm = T,probs = c(0.75)),2),')')))

Escala_Depresion<-`names<-`(data.frame(matrix(c("Geriatric Depression Scale (pts)",num1),ncol = 2)),columns)

#MMSE

num1<-c(paste(round(median(base1$MMSE,na.rm = T ),2),
              paste0('(',round(quantile(base1$MMSE,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$MMSE,na.rm = T,probs = c(0.75)),2),')')))

MMSE<-`names<-`(data.frame(matrix(c("Mini-Mental Test (pts)",num1),ncol = 2)),columns)

#Reloj

num1<-c(paste(round(median(base1$RELOJ,na.rm = T ),2),
              paste0('(',round(quantile(base1$RELOJ,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$RELOJ,na.rm = T,probs = c(0.75)),2),')')))

Reloj<-`names<-`(data.frame(matrix(c("Clock-Test (pts)",num1),ncol = 2)),columns)

#FRED

num1<-c(paste(round(median(base1$PASE_TOTAL,na.rm = T ),2),
              paste0('(',round(quantile(base1$PASE_TOTAL,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$PASE_TOTAL,na.rm = T,probs = c(0.75)),2),')')))

PASE<-`names<-`(data.frame(matrix(c("PASE-test (pts)",num1),ncol = 2)),columns)

#Perdida de Peso
Perd.Peso <- table(base1$PERDIDA_PESO,useNA = "always")[c(2)]
Perd.Pesoprop <- round(prop.table(table(base1$PERDIDA_PESO,useNA = "always")),4)[c(2)]*100

Perdida_Peso_F<-`names<-`(data.frame(c("Weight-Loss (%)"),
                                     matrix(c(paste(Perd.Peso,paste0('(',Perd.Pesoprop,')'))),ncol = 1)),
                          columns)

#Extenuacion
Extenuacion <- table(base1$EXTENUACION,useNA = "always")[c(2)]
Extenuacionprop <- round(prop.table(table(base1$EXTENUACION,useNA = "always")),4)[c(2)]*100

Extenuacion_F<-`names<-`(data.frame(c("Exhaustion (%)"),
                                    matrix(c(paste(Extenuacion,paste0('(',Extenuacionprop,')'))),ncol = 1)),
                         columns)

#Prension
Prension <- table(base1$PRENSION,useNA = "always")[c(2)]
Prension.prop <- round(prop.table(table(base1$PRENSION,useNA = "always")),4)[c(2)]*100

Prension_Fuerza_F<-`names<-`(data.frame(c("Low Hand Grip Strength (%)"),
                                        matrix(c(paste(Prension,paste0('(',Prension.prop,')'))),ncol = 1)),
                             columns)

#Marcha
Marcha <- table(base1$MARCHA,useNA = "always")[c(2)]
Marcha.prop <- round(prop.table(table(base1$MARCHA,useNA = "always")),4)[c(2)]*100

Marcha_F<-`names<-`(data.frame(c("Low Gait-Speed (%)"),
                               matrix(c(paste(Marcha,paste0('(',Marcha.prop,')'))),ncol = 1)),
                    columns)


#Physical Activity
Act.Fisica <- table(base1$ACTIVIDAD_FISICA,useNA = "always")[c(2)]
Act.Fisica.prop <- round(prop.table(table(base1$ACTIVIDAD_FISICA,useNA = "always")),4)[c(2)]*100

Act_Fisica_F<-`names<-`(data.frame(c("Low Physical Activity (%)"),
                                   matrix(c(paste(Act.Fisica,paste0('(',Act.Fisica.prop,')'))),ncol = 1)),
                        columns)

#Fraility Categories
Fragilidad <- table(base1$FRAGILIDAD_CRITERIOS,useNA = "always")[c(1:3)]
Fragilidad.prop <- round(prop.table(table(base1$FRAGILIDAD_CRITERIOS,useNA = "always")),4)[c(1:3)]*100

Fragilicad_Comp<-`names<-`(data.frame(c("Normal (%)","Pre-Fragile (%)","Fragile (%)"),
                                      matrix(c(paste(Fragilidad,paste0('(',Fragilidad.prop,')'))),ncol = 1)),
                           columns)

#Escala de Demencia VIH

num1<-c(paste(round(median(base1$EIDVHI_IHDS,na.rm = T ),2),
              paste0('(',round(quantile(base1$EIDVHI_IHDS,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$EIDVHI_IHDS,na.rm = T,probs = c(0.75)),2),')')))

Demencia_VIH<-`names<-`(data.frame(matrix(c("VIH Dementia Scale (pts)",num1),ncol = 2)),columns)

#Table 1 
Table_1<-rbind(Sexo,Edad,Educacion,Comorbilidades_2,Diabetes,Hipertension,Dislipidemia,
               Cancer_NOSIDA,EVC,EPOC,Hepatopatia.Cronica,OAD,Osteoartritis,CKD,Depresion,Ansiedad,
               PESO,TALLA,IMC,Hemoglob,Gluco,Albumin)
Table1_Flex_FINAL<-flextable::align(flextable::flextable(Table_1,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table1_Flex_FINAL,path="Table_1_FINAL.docx")

#Table 2
Table_2<-rbind(Barthel,Lawton,Caidas,Rosow,Def.Visual,Def.Auditivo,Escala_Depresion,MMSE,
               Reloj,PASE,Perdida_Peso_F,Extenuacion_F,Prension_Fuerza_F,Marcha_F,Act_Fisica_F,
               Fragilicad_Comp,Demencia_VIH)
Table2_Flex_FINAL<-flextable::align(flextable::flextable(Table_2,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table2_Flex_FINAL,path="Table_2_FINAL.docx")

###----Characterization of clusters (Table 3)----

columns <- c('Parameter',"Cluster 1 (n=120)","Cluster 2 (n=183)","Cluster 3 (n=198)")

base1_CLUSTER$cluster_REC<-NULL
base1_CLUSTER$cluster_REC[base1_CLUSTER$cluster==2]<-1
base1_CLUSTER$cluster_REC[base1_CLUSTER$cluster==1]<-2

#Sexo
sexo.cl <- table(base1_CLUSTER$SEXO,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
sexo.cl.prop <- round(prop.table(table(base1_CLUSTER$SEXO,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Sexo_clusters<-`names<-`(data.frame("Women (%)",
                                    matrix(c(paste(sexo.cl,paste0('(',sexo.cl.prop,')'))),ncol = 3)),
                         columns)

#Age
num1<-c(paste(round(mean(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$EdadCalc,na.rm = T),2),
              paste0('(',"±",round(sd(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$EdadCalc,na.rm = T),2),")")))

num2<-c(paste(round(mean(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$EdadCalc,na.rm = T),2),
              paste0('(',"±",round(sd(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$EdadCalc,na.rm = T),2),")")))

num3<-c(paste(round(mean(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$EdadCalc,na.rm = T),2),
              paste0('(',"±",round(sd(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$EdadCalc,na.rm = T),2),")")))

Edad_cluster<-`names<-`(data.frame(matrix(c("Age (Years)",num1,num2,num3),ncol = 4)),columns)

#Educacion
num1<-c(paste(round(mean(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$ESCOL.AÑOS,na.rm = T),2),
              paste0('(',"±",round(sd(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$ESCOL.AÑOS,na.rm = T),2),")")))

num2<-c(paste(round(mean(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$ESCOL.AÑOS,na.rm = T),2),
              paste0('(',"±",round(sd(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$ESCOL.AÑOS,na.rm = T),2),")")))

num3<-c(paste(round(mean(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$ESCOL.AÑOS,na.rm = T),2),
              paste0('(',"±",round(sd(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$ESCOL.AÑOS,na.rm = T),2),")")))

Educacion_cluster<-`names<-`(data.frame(matrix(c("Education (Years)",num1,num2,num3),ncol = 4)),columns)


#Comorbidities (0 Comorb)
Comorbidities_0 <- table(base1_CLUSTER$Comorbilidad==0,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Comorbidities_0.prop <- round(prop.table(table(base1_CLUSTER$Comorbilidad==0,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Comorbidities_cero<-`names<-`(data.frame("No-Comorbidities (%)",
                                         matrix(c(paste(Comorbidities_0,paste0('(',Comorbidities_0.prop,')'))),ncol = 3)),
                              columns)


#Comorbidities (1 Comorb)
Comorbidities_1 <- table(base1_CLUSTER$Comorbilidad==1,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Comorbidities_1.prop <- round(prop.table(table(base1_CLUSTER$Comorbilidad==1,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Comorbidities_one<-`names<-`(data.frame("1 Comorbidity (%)",
                                        matrix(c(paste(Comorbidities_1,paste0('(',Comorbidities_1.prop,')'))),ncol = 3)),
                             columns)

#Comorbidities (≥2 Comorb)
Comorbidities_2 <- table(base1_CLUSTER$Comorbilidad>=2,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Comorbidities_2.prop <- round(prop.table(table(base1_CLUSTER$Comorbilidad>=2,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Comorbidities_two<-`names<-`(data.frame("≥2 Comorbidities (%)",
                                        matrix(c(paste(Comorbidities_2,paste0('(',Comorbidities_2.prop,')'))),ncol = 3)),
                             columns)

#Diabetes
diabetes.cl <- table(base1_CLUSTER$DM2,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
diabetes.cl.prop <- round(prop.table(table(base1_CLUSTER$DM2,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Diab_clusters<-`names<-`(data.frame("Diabetes (%)",
                                    matrix(c(paste(diabetes.cl,paste0('(',diabetes.cl.prop,')'))),ncol = 3)),
                         columns)

#Arterial Hypertension
hipertension.cl <- table(base1_CLUSTER$HAS,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
hipertension.cl.prop <- round(prop.table(table(base1_CLUSTER$HAS,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
HAS_clusters<-`names<-`(data.frame("Arterial Hypertension (%)",
                                   matrix(c(paste(diabetes.cl,paste0('(',diabetes.cl.prop,')'))),ncol = 3)),
                        columns)

#Dyslipidemia
Dislipidemia.cl <- table(base1_CLUSTER$DISLIPIDEMIA,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Dislipidemia.cl.prop <- round(prop.table(table(base1_CLUSTER$DISLIPIDEMIA,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Dislipidemia_clusters<-`names<-`(data.frame("Dyslipidemia (%)",
                                            matrix(c(paste(Dislipidemia.cl,paste0('(',Dislipidemia.cl.prop,')'))),ncol = 3)),
                                 columns)

#Cancer
Cancer.cl <- table(base1_CLUSTER$CANCER_NOSIDA,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Cancer.cl.prop <- round(prop.table(table(base1_CLUSTER$CANCER_NOSIDA,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Cancer_clusters<-`names<-`(data.frame("Cancer (%)",
                                      matrix(c(paste(Cancer.cl,paste0('(',Cancer.cl.prop,')'))),ncol = 3)),
                           columns)

#Cardiovascular
CVD.cl <- table(base1_CLUSTER$EVC,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
CVD.cl.prop <- round(prop.table(table(base1_CLUSTER$EVC,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
CVD_clusters<-`names<-`(data.frame("Cardiovascular Disease (%)",
                                   matrix(c(paste(CVD.cl,paste0('(',CVD.cl.prop,')'))),ncol = 3)),
                        columns)

#COPD
COPD.cl <- table(base1_CLUSTER$EPOC,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
COPD.cl.prop <- round(prop.table(table(base1_CLUSTER$EPOC,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
COPD_clusters<-`names<-`(data.frame("COPD (%)",
                                    matrix(c(paste(COPD.cl,paste0('(',COPD.cl.prop,')'))),ncol = 3)),
                         columns)

#Chronic Hepatopathy 
Hepat.cl <- table(base1_CLUSTER$HEPAT_CRON,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Hepat.cl.prop <- round(prop.table(table(base1_CLUSTER$HEPAT_CRON,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Hepat_clusters<-`names<-`(data.frame("Chronic Hepatopathy (%)",
                                     matrix(c(paste(Hepat.cl,paste0('(',Hepat.cl.prop,')'))),ncol = 3)),
                          columns)

#OAD
OAD.cl <- table(base1_CLUSTER$OAD,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
OAD.cl.prop <- round(prop.table(table(base1_CLUSTER$OAD,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
OAD_clusters<-`names<-`(data.frame("OAD (%)",
                                   matrix(c(paste(OAD.cl,paste0('(',OAD.cl.prop,')'))),ncol = 3)),
                        columns)

#Osteoarthritis
Osteoart.cl <- table(base1_CLUSTER$OSTEOP,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Osteoart.cl.prop <- round(prop.table(table(base1_CLUSTER$OSTEOP,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Osteoart_clusters<-`names<-`(data.frame("OAD (%)",
                                        matrix(c(paste(Osteoart.cl,paste0('(',Osteoart.cl.prop,')'))),ncol = 3)),
                             columns)


#CKD
CKD.cl <- table(base1_CLUSTER$RENAL_CRON,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
CKD.cl.prop <- round(prop.table(table(base1_CLUSTER$RENAL_CRON,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
CKD_clusters<-`names<-`(data.frame("Chronic Kidney Disease (%)",
                                   matrix(c(paste(CKD.cl,paste0('(',CKD.cl.prop,')'))),ncol = 3)),
                        columns)

#TAD
TDM.cl <- table(base1_CLUSTER$TDM,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
TDM.cl.prop <- round(prop.table(table(base1_CLUSTER$TDM,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
TDM_clusters<-`names<-`(data.frame("Depression (%)",
                                   matrix(c(paste(TDM.cl,paste0('(',TDM.cl.prop,')'))),ncol = 3)),
                        columns)

#TAD
TAD.cl <- table(base1_CLUSTER$TAD,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
TAD.cl.prop <- round(prop.table(table(base1_CLUSTER$TAD,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
TAD_clusters<-`names<-`(data.frame("Anxiety (%)",
                                   matrix(c(paste(TAD.cl,paste0('(',TAD.cl.prop,')'))),ncol = 3)),
                        columns)


#BMI
num1<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$IMC,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$IMC,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$IMC,na.rm = T,prob=c(0.75)),")"))

num2<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$IMC,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$IMC,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$IMC,na.rm = T,prob=c(0.75)),")"))

num3<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$IMC,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$IMC,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$IMC,na.rm = T,prob=c(0.75)),")"))

IMC_cluster<-`names<-`(data.frame(matrix(c("BMI (Kg)",num1,num2,num3),ncol = 4)),columns)

#Hemoglobin
num1<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$HEMOGLOB,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$HEMOGLOB,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$HEMOGLOB,na.rm = T,prob=c(0.75)),")"))

num2<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$HEMOGLOB,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$HEMOGLOB,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$HEMOGLOB,na.rm = T,prob=c(0.75)),")"))

num3<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$HEMOGLOB,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$HEMOGLOB,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$HEMOGLOB,na.rm = T,prob=c(0.75)),")"))

Hemoglobin_cluster<-`names<-`(data.frame(matrix(c("Hemoglobin (grs/ml)",num1,num2,num3),ncol = 4)),columns)


#Glucose
num1<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$GLUCOSA,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$GLUCOSA,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$GLUCOSA,na.rm = T,prob=c(0.75)),")"))

num2<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$GLUCOSA,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$GLUCOSA,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$GLUCOSA,na.rm = T,prob=c(0.75)),")"))

num3<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$GLUCOSA,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$GLUCOSA,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$GLUCOSA,na.rm = T,prob=c(0.75)),")"))

Glucose_cluster<-`names<-`(data.frame(matrix(c("Glucose (mg/dl)",num1,num2,num3),ncol = 4)),columns)


#Albumin
num1<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$ALBUMINA,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$ALBUMINA,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$ALBUMINA,na.rm = T,prob=c(0.75)),")"))

num2<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$ALBUMINA,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$ALBUMINA,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$ALBUMINA,na.rm = T,prob=c(0.75)),")"))

num3<-c(paste(median(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$ALBUMINA,na.rm = T),"(",
              quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$ALBUMINA,na.rm = T,prob=c(0.25)),
              "-",quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$ALBUMINA,na.rm = T,prob=c(0.75)),")"))

Albumin_cluster<-`names<-`(data.frame(matrix(c("Albumin (mg/dl)",num1,num2,num3),ncol = 4)),columns)

#Polyfarmacia
Polyfarm <- table(base1_CLUSTER$POLIFARMACIA_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
Polyfarm.prop <- round(prop.table(table(base1_CLUSTER$POLIFARMACIA_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Polyfarmacia<-`names<-`(data.frame("Polypharmacy (%)",
                                   matrix(c(paste(Polyfarm,paste0('(',Polyfarm.prop,')'))),ncol = 3)),
                        columns)


#Deficit Sensorial
def.sen <- table(base1_CLUSTER$DEFICIT_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
def.sen.prop <- round(prop.table(table(base1_CLUSTER$DEFICIT_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Deficit_Sensorial<-`names<-`(data.frame("Sensorial Deficit (%)",
                                        matrix(c(paste(def.sen,paste0('(',def.sen.prop,')'))),ncol = 3)),
                             columns)

#Deficit Cognitivo
def.cog <- table(base1_CLUSTER$COGNITIVE_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
def.cog.prop <- round(prop.table(table(base1_CLUSTER$COGNITIVE_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Deficit_Cognitivo<-`names<-`(data.frame("Cognitive Decline (%)",
                                        matrix(c(paste(def.cog,paste0('(',def.cog.prop,')'))),ncol = 3)),
                             columns)

#Deficit Fisico
def.fisico <- table(base1_CLUSTER$SX_DISCAPACIDAD_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
def.fisico.prop <- round(prop.table(table(base1_CLUSTER$SX_DISCAPACIDAD_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Deficit_Fisico<-`names<-`(data.frame("Physical Disability (%)",
                                     matrix(c(paste(def.fisico,paste0('(',def.fisico.prop,')'))),ncol = 3)),
                          columns)

#Fragilidad
fragilidad <- table(base1_CLUSTER$FRAGILIDAD_CRITERIOS,base1_CLUSTER$cluster_REC,useNA = "always")[c(1,2,3,5,6,7,9,10,11)]
fragilidad.prop <- round(prop.table(table(base1_CLUSTER$FRAGILIDAD_CRITERIOS,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(1,2,3,5,6,7,9,10,11)]*100
Fragilidad<-`names<-`(data.frame(c("Normal (%)","PreFraility (%)","Fraility (%)"),
                                 matrix(c(paste(fragilidad,paste0('(',fragilidad.prop,')'))),ncol = 3)),
                      columns)

#Caidas
caida <- table(base1_CLUSTER$CAIDAS_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
caida.prop <- round(prop.table(table(base1_CLUSTER$CAIDAS_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Caidas<-`names<-`(data.frame("Falls (%)",
                             matrix(c(paste(caida,paste0('(',caida.prop,')'))),ncol = 3)),
                  columns)

#Caidas
depre <- table(base1_CLUSTER$DEPRESSION_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
depre.prop <- round(prop.table(table(base1_CLUSTER$DEPRESSION_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Depre<-`names<-`(data.frame("Depression (%)",
                            matrix(c(paste(depre,paste0('(',depre.prop,')'))),ncol = 3)),
                 columns)

#Multicomorbilidades
multicom <- table(base1_CLUSTER$MULTICOMORB_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
multicom.prop <- round(prop.table(table(base1_CLUSTER$MULTICOMORB_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Multicomorb<-`names<-`(data.frame("Multicomorbidities (%)",
                                  matrix(c(paste(multicom,paste0('(',multicom.prop,')'))),ncol = 3)),
                       columns)

#Incontinencia
inconti <- table(base1_CLUSTER$MULTICOMORB_ANY,base1_CLUSTER$cluster_REC,useNA = "always")[c(2,5,8)]
inconti.prop <- round(prop.table(table(base1_CLUSTER$MULTICOMORB_ANY,base1_CLUSTER$cluster_REC,useNA = "always"),2),4)[c(2,5,8)]*100
Incontinencia<-`names<-`(data.frame("Incontinency (%)",
                                    matrix(c(paste(inconti,paste0('(',inconti.prop,')'))),ncol = 3)),
                         columns)

Table_4<-rbind(Edad_cluster,Sexo_clusters,Edad_cluster,Educacion_cluster,Comorbidities_cero,Comorbidities_one,Comorbidities_two,Diab_clusters,
               HAS_clusters,Dislipidemia_clusters,Cancer_clusters,CVD_clusters,COPD_clusters,Hepat_clusters,OAD_clusters,Osteoart_clusters,
               CKD_clusters,TDM_clusters,TAD_clusters,IMC_cluster,Hemoglobin_cluster,Glucose_cluster,Albumin_cluster)
Table4_Flex_FINAL<-flextable::align(flextable::flextable(Table_4,cwidth=4),align="center",part="all")%>%flextable::autofit()

flextable::save_as_docx(Table4_Flex_FINAL,path="Table_4_FINAL.docx")

table(base1_CLUSTER$ZONA)
prop.table(table(base1_CLUSTER$ZONA))*100
kruskal.test(base1_CLUSTER$ALBUMINA~base1_CLUSTER$cluster_REC)
stats::chisq.test(table(base1_CLUSTER$TAD,base1_CLUSTER$cluster_REC),correct = T,simulate.p.value = T)
base1$RNA_VIH_ACTUAL
median(base1$RNA_VIH_ACTUAL,na.rm = T)
quantile(base1$RNA_VIH_ACTUAL,na.rm = T,prob=c(0.25,0.75))







table(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$RNAcargaCAT)
table(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$RNAcargaCAT)
table(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$RNAcargaCAT)

prop.table(table(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$RNAcargaCAT))*100
prop.table(table(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$RNAcargaCAT))*100
prop.table(table(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$RNAcargaCAT))*100


summary(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$CD4_ACTUAL)
summary(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$CD4_ACTUAL)
summary(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$CD4_ACTUAL)

quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==1,]$CD4_ACTUAL,probs = c(0.25,0.75),na.rm = T)
quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==2,]$CD4_ACTUAL,probs = c(0.25,0.75),na.rm = T)
quantile(base1_CLUSTER[base1_CLUSTER$cluster_REC==3,]$CD4_ACTUAL,probs = c(0.25,0.75),na.rm = T)
