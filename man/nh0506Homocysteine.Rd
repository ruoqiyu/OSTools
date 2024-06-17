\name{nh0506Homocysteine}
\alias{nh0506Homocysteine}
\docType{data}
\title{
Homocysteine and Smoking
}
\description{
NHANES 2005-2006 data on smoking and homocysteine levels in adults.
}
\usage{data("nh0506Homocysteine")}
\format{
  A data frame with 2759 observations on the following 11 variables.
  \describe{
    \item{\code{X}}{Row number, 1 to 2759}
    \item{\code{SEQN}}{NHANES identification number}
    \item{\code{z}}{Smoking status, 1 = daily smoker, 0 = never smoker}
    \item{\code{female}}{1 = female, 0 = male}
    \item{\code{age}}{Age in years, >=20, capped at 85}
    \item{\code{black}}{1=black race, 0=other}
    \item{\code{education}}{Level of education}
    \item{\code{povertyr}}{Ratio of family income to the poverty level, capped at 5 times poverty}
    \item{\code{bmi}}{BMI or body-mass-index}
    \item{\code{cigsperday30}}{Cigarettes smoked per day, 0 for never smokers}
    \item{\code{cotinine}}{Blood cotinine level, a biomarker of recent exposure to tobacco}
    \item{\code{homocysteine}}{Level of homocysteine}
  }
}
\details{
The following code constructed the data as used here.
Attention is confined to adults, excluding children.
Also, people who have smoked in the past, but do not
now smoke daily, are excluded.
A moderate number of individuals with missing
povertyr, cotinine or homocysteine were excluded.

library(foreign)
DEMO<-read.xport("DEMO_D.XPT")
HCY<-read.xport("HCY_D.XPT")
SMQ<-read.xport("SMQ_D.XPT")
BMX<-read.xport("BMX_D.XPT")
COT<-read.xport("COT_D.XPT")
d<-merge(DEMO,HCY,by="SEQN",all.x=TRUE)
d<-merge(d,SMQ,by="SEQN",all.x=TRUE)
d<-merge(d,COT,by="SEQN",all.x=TRUE)
d<-merge(d,BMX,by="SEQN",all.x=TRUE)
rm(DEMO,HCY,SMQ,COT,BMX)
SEQN<-d$SEQN
age<-d$RIDAGEYR
race<-d$RIDRETH1
black<-1*(race==4)
hispanic<-1*((race==1)|(race==2))
female<-1*(d$RIAGENDR==2)
education<-d$DMDEDUC2
education[education>6]<-NA
povertyr<-d$INDFMPIR
homocysteine<-d$LBXHCY
bmi<-d$BMXBMI
cotinine<-d$LBXCOT
cigs100life<-d$SMQ020
cigs100life[cigs100life>3]<-NA
cigs100life<-(cigs100life==1)*1
smokenow<-1*(d$SMQ040<2.5)
smokenow[cigs100life==0]<-0
cigsdays30<-d$SMD641
cigsdays30[cigsdays30>32]<-NA
cigsdays30[smokenow==0]<-0
cigsperday30<-d$SMD650
cigsperday30[cigsperday30>100]<-NA
cigsperday30[smokenow==0]<-0
dailysmoker<-1*((cigs100life==1)&(cigsdays30==30)&(smokenow==1)&(cigsperday30>=10))
neversmoker<-1*((cigs100life==0)&(smokenow==0))
z<-dailysmoker
z[(neversmoker==0)&(dailysmoker==0)]<-(-999)
ds<-data.frame(SEQN,female,age,black,education,povertyr,bmi,
               homocysteine,cotinine,
               cigs100life,smokenow,cigsdays30,cigsperday30,
               dailysmoker,neversmoker,z)
use<-age>=20
ds1<-ds[use,]
use<-complete.cases(ds1)
ds1$z[ds1$z==-999]<-NA
ds2<-ds1[use&!is.na(ds1$z),]
rm(SEQN,female,age,black,hispanic,education,povertyr,
   homocysteine,cotinine,
   cigs100life,smokenow,cigsdays30,cigsperday30,
   dailysmoker,neversmoker,z,race,use,bmi)
ds2<-ds2[order(1-ds2$z),]
attach(ds2)
nh0505Homocysteine<-data.frame(SEQN,z,female,age,black,education,
   povertyr,bmi,cigsperday30,cotinine,homocysteine)
write.csv(nh0506Homocysteine,"nh0506Homocysteine.csv")

}
\source{
From the NHANES web page, for NHANES 2005-2006
}
\references{
US National Health and Nutrition Examination Survey, 2005-2006.  From the US Center for Health Statistics.
}
\examples{
data(nh0506Homocysteine)
summary(nh0506Homocysteine)
}
\keyword{datasets}
