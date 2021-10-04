install.packages("dlm")
install.packages("R0")
install.packages("EpiEstim")
install.packages("TSdist")
library(dlm)
require(dplyr)
library(TSdist)

coronavirus <- read.csv("https://raw.githubusercontent.com/RamiKrispin/coronavirus/master/csv/coronavirus.csv", header = T)
vaccination <- read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv", header = T)
population <- read.csv("https://raw.githubusercontent.com/datasets/population/master/data/population.csv", header = T)

df <- coronavirus %>%
  dplyr::filter(country == "Poland") %>%
  dplyr::group_by(date, type) %>%
  dplyr::summarise(total = sum(cases, na.rm = TRUE)) %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = total
  ) %>%
  dplyr::arrange(date) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(active = confirmed - death - recovered) %>%
  dplyr::mutate(
    confirmed_cum = cumsum(confirmed),
    death_cum = cumsum(death),
    recovered_cum = cumsum(recovered),
    active_cum = cumsum(active)
  )

print(which(df$confirmed!=0)[1])




#przygotowanie danych o szczepieniach i wybranie kolumn z data i wyszczepionymi
vaccinationDF <- vaccination %>% 
  dplyr::filter(location == "Poland") %>%
  dplyr::select(3,6)

df <- left_join(df, vaccinationDF, by=c('date'))

#Zastapienie pierwszej wartosci NA zerem 
df$people_fully_vaccinated[1] <- 0

#zastapienie pozostaluch NA poprzednia wartoscia
library(zoo)
df$people_fully_vaccinated <- na.locf(df$people_fully_vaccinated)


population <- population %>%
  dplyr::group_by(population$Country.Name) %>%
  slice(n())%>%
  ungroup()

print(max(population$Year))
####################################################################################
################################# Program glowny ###################################
####################################################################################


library(zoo)
library(reshape2)
library(ggplot2)

R0Country <- function(countryName){
  df <- coronavirus %>%
    dplyr::filter(country == countryName) %>%
    dplyr::group_by(date, type) %>%
    dplyr::summarise(total = sum(cases, na.rm = TRUE)) %>%
    tidyr::pivot_wider(
      names_from = type,
      values_from = total
    ) %>%
    dplyr::arrange(date) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(active = confirmed - death - recovered) %>%
    dplyr::mutate(
      confirmed_cum = cumsum(confirmed),
      death_cum = cumsum(death),
      recovered_cum = cumsum(recovered),
      active_cum = cumsum(active)
    )
  
  vaccinationDF <- vaccination %>% 
    dplyr::filter(location == countryName) %>%
    dplyr::select(3,6)
  
  df <- left_join(df, vaccinationDF, by=c('date'))
  
  df$people_fully_vaccinated[1] <- 0

  df$people_fully_vaccinated <- na.locf(df$people_fully_vaccinated)

  countryPop <- population %>%
    dplyr::filter(population$Country.Name == countryName) %>%
    slice(n())%>%
    ungroup()
  
  print(countryPop$Value)
  
  N <- countryPop$Value
  xrang <- which(df$confirmed!=0)[1]:nrow(df)
  Infected <- df$active_cum[xrang]
  Recover <- df$recovered_cum[xrang]+df$death_cum[xrang]
  Vaccined <- df$people_fully_vaccinated[xrang]
  Suspect <- N -Infected - Recover - Vaccined
  init<-c(S = N-max(Infected[1], 1), I = max(Infected[1], 1), R = 0)
  
  Opt <- optim(c(0.5, 0.5),
               RSS,
               method = "L-BFGS-B",
               lower = c(0, 0),
               upper = c(1, 1)
  )
  
  Opt$message
  Opt_par <- setNames(Opt$par, c("beta", "gamma"))
  Opt_par
  beta <- Opt$par[[1]]
  gamma <- Opt$par[[2]]
  
  
  ww <- SIRfunction(N, Infected, Recover, Vaccined, df$death_cum, Opt_par[[1]], Opt_par[[2]])
  
  #makePlot(ww)
  
  ww$R0
}


wx1 = matrix(c(0), 10,10)
rownames(wx1) <- colnames(wx1) <- c("Poland", "China", "Germany", "France", "Ukraine", "Italy", "Canada", "Brazil", "Israel", "Australia")
wx2 = list()
for(i in 1:10)  wx2[i] <- R0Country(rownames(wx1)[[i]])

for(i in 1:9){
  for(j in (i+1):10){
    wx1[i,j] <- TSDistances(wx2[i], wx2[j], distance ="dtw", sigma = 10)
    wx1[j,i] = wx1[i,j]
    print((wx1[j, i]))
  }
}
wx3 = hclust(as.dist(wx1))

plot(wx3)

wx3den <- as.dendrogram(wx3)

plot(wx3den,type = "triangle", ylab = "Height")

plot(wx3den, xlab = "Height", horiz = TRUE)


library("ape")

plot(as.phylo(wx3), type = "unrooted", cex = 0.6, no.margin = TRUE)

plot(as.phylo(wx3), type = "fan")


colors = c("red", "blue", "green", "black")
clus4 = cutree(wx3, 4)
plot(as.phylo(wx3), type = "fan", tip.color = colors[clus4], label.offset = 0.2, cex = 0.8)
plot(as.phylo(wx3), type = "unrooted",  tip.color = colors[clus4], cex = 0.6, no.margin = TRUE)
################dla 15 krajow

wx12 = matrix(c(0), 15,15)
rownames(wx12) <- colnames(wx12) <- c("Poland", "China", "Germany", "France", "Ukraine", "Italy", "Canada", "Brazil", "Israel", "Australia", "Japan", "Morocco", "Sweden", "Hungary", "Latvia")
wx22 = list()
for(i in 1:15)  wx22[i] <- R0Country(rownames(wx12)[[i]])

for(i in 1:14){
  for(j in (i+1):15){
    wx12[i,j] <- TSDistances(wx22[i], wx22[j], distance ="dtw", sigma = 10)
    wx12[j,i] = wx12[i,j]
    print((wx12[j, i]))
  }
}

wx32 = hclust(as.dist(wx12))

plot(wx32)

wx32den <- as.dendrogram(wx32)

plot(wx32den,type = "triangle", ylab = "Height")

plot(wx32den, xlab = "Height", horiz = TRUE)


library("ape")

plot(as.phylo(wx32), type = "unrooted", cex = 0.6, no.margin = TRUE)

plot(as.phylo(wx32), type = "fan")

colors = c("red", "blue", "green", "black")
clus4 = cutree(wx32, 4)
plot(as.phylo(wx32), type = "fan", tip.color = colors[clus4], label.offset = 0.2, cex = 0.8)
plot(as.phylo(wx32), type = "unrooted",  tip.color = colors[clus4], cex = 0.6, no.margin = TRUE)

#library("Rpart") ladniejszy wykres <- nie da sie zastosowac
#tsdist korelacja
#porównac beta gamma
#porownac po epiestim

#poprobować rozne odleglosci








setwd("E:/studia/seminarium/R0")
saveRDS(df, file = "sirDF.tab")

df<-as.data.frame(readRDS("sirDF.tab"))


SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = xrang, func = SIR, parms = parameters)
  fit <- out[, 3]
  sum((Infected - fit)^2)
}


##############     SIR W FUNKCJI     ##################
library(deSolve)

N <- 37846605
xrang <- 43:nrow(df)
Infected <- df$active_cum[xrang]
Recover <- df$recovered_cum[xrang]+df$death_cum[xrang]
Vaccined <- df$people_fully_vaccinated[xrang]
Suspect <- N -Infected - Recover - Vaccined
init<-c(S = N-max(Infected[1], 1), I = max(Infected[1], 1), R = 0)

Opt <- optim(c(0.5, 0.5),
             RSS,
             method = "L-BFGS-B",
             lower = c(0, 0),
             upper = c(1, 1)
)

Opt$message
Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par
beta <- Opt$par[[1]]
gamma <- Opt$par[[2]]


SIRfunction <- function(N, I, R, S,Vacc,D, Opt1, Opt2){
  nt <- length(I)
  
  grI <- I[2:nt]/I[1:nt-1]-1
  epidemia <-which(!is.na(grI) && grI != Inf)
  if(min(epidemia==1)){
    epidemia <- setdiff(epidemia, Opt1)
  }
  I <- I[epidemia]
  R <- R[epidemia]+D[epidemia]
  Vacc <- Vacc[epidemia]
  S <- N - I - R - Vacc
  
  
  
  init<-c(S = N-max(Infected[1], 1), I = max(Infected[1], 1), R = 0)
  
  
  RS1 <- function(v,w,gamma, grJ) {
    m1<-dlm(FF=gamma, V=v, GG=1, W=w, m0=0, C0=1)
    wyn<-dlmFilter(grJ,m1)
    sd(wyn$f-grJ) 
  }
  
  blad = 1.0e10
  v0 = 1
  w0 = 1
  
  grI1<-grI[epidemia]
  zakres1<-c(min(epidemia)-1,epidemia)
  Susp <- N-I-R-df$death_cum - Vacc
  Susp<-Susp[zakres1]
  dat1<-df$dat[zakres1]
  dan<-c(v=0,w=0,fw=Inf)
  for (v in df$seq1) for (w in df$seq2) {
    w1<-RS1(v,w,gamma,grI1)  
    if (w1<dan['fw']) 
      dan<-c(V=v,W=w,fw=w1)
  }
  
  model <- dlm(FF = gamma, V=v0, GG=1, W=w0, m0 = 4, C0=1e2) #petla dla roznych wartosci
  wyn <- dlmSmooth(grI, model)
 
  #View(wyn)
  
  mse.list=dlmSvd2var(wyn$U.S, wyn$D.S)
  se = 1.96*t(sapply(mse.list, FUN = function(x) sqrt(diag(x)))*sqrt(N/S))
  
  
  wyn01 = data.frame(t=df$date[zakres1], R0=(wyn$s+1)*N/S, R0_UCL=(wyn$s+1)*N/S+se[1,],
                     R0_LCL = (wyn$s+1)*N/S-se[1,], grI=(c(0, grI)+1)*N/S)
  wyn01
}
  ww<-SIRfunction(N, Infected, Recover, Suspect, Vaccined, df$death_cum, Opt_par[[1]], Opt_par[[2]])
  
makePlot <- function(ww){

  wwTo1 = 0
  ww1plus = 0
  for(row in 1:nrow(ww)){
    r0val <- ww[row, "R0"]
    if(r0val <=1)
      wwTo1<-wwTo1+1
  }
  ww1plus <- nrow(ww)-wwTo1
  ww$t = as.numeric(1:nrow(ww))
  melted<- melt(ww[,c(1,2,5)], id.vars="t")
  ggplot(melted, aes(x=t, y=value, color = variable)) + theme_bw() + geom_line()
}




##############     R0     ##################

library(R0)
check.incid(df$confirmed)


W1<-data.frame(dates = as.Date(df$date), I = as.numeric(df$confirmed), num = seq(1, nrow(df)))

W1 <- W1[W1$I!=0, ]

t<-W1$dates

GT.flu <- generation.time("gamma", c(4.7, 1))
rs <- estimate.R(W1$I, GT = GT.flu, methods = c("EG"))
WR0 <- c(EGR = rs$estimates$EG$R, EGR_LCL = rs$estimates$EG$conf.int[[1]],
         EGR_UCL = rs$estimates$EG$conf.int[[2]])



rs <- estimate.R(W1$I, GT = GT.flu, methods = c("ML"))
WR0 <- c(WR0, MLR = rs$estimates$ML$R, MLR_LCL = rs$estimates$ML$conf.int[[1]],
         MLR_UCL = rs$estimates$ML$conf.int[[2]])



rs <- estimate.R(W1$I, GT = GT.flu, methods = c("SB"))
tw1<-floor((nrow(W1)-length(rs$estimates$SB$R))/2)
tw2<-tw1+length(rs$estimates$SB$R)-1
wrt = data.frame(dat = t, SBR = c(NA), SBR_LCL = c(NA), SBR_UCL = c(NA))
wrt[tw1:tw2, 2] <- rs$estimates$SB$R
wrt[tw1:tw2, 3] <- rs$estimates$SB$conf.int[[1]]
wrt[tw1:tw2, 4] <- rs$estimates$SB$conf.int[[2]]



WR0 <- c(WR0, EGR = rs$estimates$EG$R, EGR_LCL = rs$estimates$EG$conf.int[[1]],
         EGR_UCL = rs$estimates$EG$conf.int[[2]])


GT_obj = generation.time("empirical", val = df$confirmed)



est.R0.TD(df$confirmed, GT = GT_obj,begin=1, end=end, nsim=1000)



##############  EpiEstim  ##################
library(EpiEstim)
library(incidence)

data("Flu2009")

W3<-estimate_R(incid = df$confirmed, method = c("non_parametric_si"), 
               config = make_config(list(t_start = 49, t_end = nrow(df), si_distr = Flu2009$si_distr)))

wynEpiEstim = c(ER = w3$R$'Mean(R)', ER_LCL = W3$R$'Mean(R)'-1.96*W3$R$'Std(R)',
                ER_UCL =W3$R$'Mean(R)'+1.96*W3$R$'Std(R)')

View(df$confirmed)


tmp <- data.frame(date = df$date[xrang],I = df$confirmed[xrang])

res_parametric_si <- estimate_R(tmp$I,
                                method = "parametric_si",
                                config = make_config(list(
                                  mean_si = 2.6,
                                  std_si = 1.5)))

plot(as.incidence(tmp$I, dates = tmp$date ))
head(res_parametric_si$R)

plot(res_parametric_si,"R")











##############    SIR     ##################
#########   TYLKO DO WGLADU   ##############

N <- 37846605
xrang <- 43:nrow(df)
Infected <- df$active_cum[xrang]
Recover <- df$recovered_cum[xrang]+df$death_cum[xrang]
Suspect <- N -Infected - Recover
init<-c(S = N-max(Infected[1], 1), I = max(Infected[1], 1), R = 0)

library(deSolve)
Opt <- optim(c(0.5, 0.5),
             RSS,
             method = "L-BFGS-B",
             lower = c(0, 0),
             upper = c(1, 1)
)

Opt$message
Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par
beta <- Opt$par[[1]]
gamma <- Opt$par[[2]]


library(dlm)

nt <- nrow(df)

grI <- df$confirmed_cum[44:nt]/df$confirmed_cum[43:(nt-1)]-1

blad = 1.0e10
v0 = 1
w0 = 1
for(V in sek1)for(W in sek2){
  dlm(FF = gamma, V = v, GG = 1, W = w, m0 = 4, C0 = 1e2)
  wyn <- dlmFilter(grI, m1)
  if(sd(wyn$f - grI)<blad){
    w0 = w
    v0 = v
    blad = sd(wyn$f - grI)
  }
  
}
model <- dlm(FF = gamma, V=v0, GG=1, W=w0, m0 = 4, C0=1e2) #petla dla roznych wartosci
wyn <- dlmSmooth(grI, model)

mse.list=dlmSvd2var(wyn$U.S, wyn$D.S)
se = 1.96*t(sapply(mse.list, FUN = function(x) sqrt(diag(x)))*sqrt(N/Suspect))
wyn01 = data.frame(t=df$date[zakres], R0=(wyn$s+1)*N/Suspect, R0_UCL=(wyn1$s+1)*N/Suspect+se[1,],
                   R0_LCL = (wyn1$s+1)*N/Suspect-se[1,], grI=(c(0, grI)+1)*N/Suspect)


wyn1<- (wyn$s+1)*N/Suspect[43:nt]
library(ggplot2)
betadf<- data.frame(t = df$date[43:nt], x=wyn1, typ = c("beta"))
ggplot(betadf, aes(x = as.Date(t), y = x, color = typ)) + geom_line()
plot(wyn1, type = "l")

#wynik funkcji blad, w, v, beta, gamma

makePlot<-function(dataToProcess){
  N <- 37846605
  xrang <- 43:nrow(dataToProcess)
  Infected <- dataToProcess$active_cum[xrang]
  Recover <- dataToProcess$recovered_cum[xrang]+dataToProcess$death_cum[xrang]
  Suspect <- N -Infected - Recover
  init<-c(S = N-max(Infected[1], 1), I = max(Infected[1], 1), R = 0)
  
  library(deSolve)
  Opt <- optim(c(0.5, 0.5),
               RSS,
               method = "L-BFGS-B",
               lower = c(0, 0),
               upper = c(1, 1)
  )
  # check for convergence
  Opt$message
  Opt_par <- setNames(Opt$par, c("beta", "gamma"))
  Opt_par
  beta <- Opt$par[[1]]
  gamma <- Opt$par[[2]]
  
  nt <- nrow(dataToProcess)
  
  grI <- dataToProcess$confirmed_cum[44:nt]/dataToProcess$confirmed_cum[43:(nt-1)]-1
  model <- dlm(FF = gamma, V=1, GG=1, W=1,m0 = 4, C0=1e2)
  wyn <- dlmSmooth(grI, model)
  
  wyn1<- (wyn$s+1)*N/Suspect[43:nt]
  library(ggplot2)
  betadf<- data.frame(t = df$date[43:nt], x=wyn1, typ = c("beta"))
  ggplot(betadf, aes(x = as.Date(t), y = x, color = typ)) + geom_line()
  plot(wyn1, type = "l")

}
makePlot(df)