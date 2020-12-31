#######################################################################
# OTTO DAHLIN - Stock and Watson 3*3 Dimensional VAR Replication
######################################################################


library(writexl)
library(ggplot2)
library(dplyr)
library(AER)
library(lmtest)
library(tseries)
library(urca)
library(dynlm)
library(sandwich)
library(readxl)
library(forecast)
library(xts)
library(vars)
library(zoo)
library(timeSeries)
library(quantmod)
library(mFilter)
library(seasonal)
library(lubridate)
library(CARS)
library(car)


# Reading in the TXT file/data.
data <- read.delim("sw2001.txt", header= TRUE)
data

#########################################################################
# CONVERTING TO TIME SERIES
#########################################################################

# INFLATION
Inflation <- ts(data[, 2], frequency = 4)
Inflation
# Convert to ts INflation
Inflation.ts <- ts(Inflation, frequency=4, start=c(1960,1))
Inflation.ts
is.ts(Inflation.ts)
ts.plot(Inflation.ts, main="Inflation", ylab="Percentage %")

# UNEMPLOYMENT
Unemployment <- ts(data[, 3], frequency = 4)
Unemployment
# Convert to ts Unemployment
Unemployment.ts <- ts(Unemployment, frequency=4, start=c(1960,1))
Unemployment.ts
is.ts(Unemployment.ts)
ts.plot(Unemployment.ts, main="Unemployment", ylab="Percentage %")

# Federal Funds Rate(FFR)
Fedfunds <- ts(data[, 4], frequency = 4)
Fedfunds
# Convert to ts Unemployment
Fedfunds.ts <- ts(Fedfunds, frequency=4, start=c(1960,1))
Fedfunds.ts
is.ts(Fedfunds.ts)
ts.plot(Fedfunds.ts, main="Federal Funds Rate (FFR)", ylab="Percentage %")

#########################################################################

# Data Transformation.

#########################################################################


# KPSS-test Stationarity Test <INFLATION> in levels

ur.kpss(Inflation.ts, type = "tau")@teststat
ur.kpss(Inflation.ts, type ="tau")@cval
# Not stationary. Need diff

# First Difference of INFLATION levels series
d.Inflation.ts <- diff(Inflation.ts)
ur.kpss(d.Inflation.ts, type = "tau")@teststat
ur.kpss(d.Inflation.ts, type ="tau")@cval
# Stationary!

# Diff log (growth rates interpretation.)
d.log.Inflation.ts <- diff(log(Inflation.ts))
ur.kpss(d.log.Inflation.ts, type = "tau")@teststat
ur.kpss(d.log.Inflation.ts, type ="tau")@cval


# KPSS-test Stationarity Test <UNEMPLOYMENT> in levels

ur.kpss(Unemployment.ts, type = "tau")@teststat
ur.kpss(Unemployment.ts, type ="tau")@cval
# Not stationary

# First Difference of UNEMPLOYMENT levels series
d.Unemployment.ts <- diff(Unemployment.ts)
ur.kpss(d.Unemployment.ts, type = "tau")@teststat
ur.kpss(d.Unemployment.ts, type ="tau")@cval

# Diff log (growth rates interpretation.)
d.log.Unemployment.ts <- diff(log(Unemployment.ts))
ur.kpss(d.log.Unemployment.ts, type = "tau")@teststat
ur.kpss(d.log.Unemployment.ts, type ="tau")@cval


# FED FUNDS RATE
ur.kpss(Fedfunds.ts, type = "tau")@teststat
ur.kpss(Fedfunds.ts, type ="tau")@cval

#########################################################

# VAR-engineering

#########################################################

# Recursive ordering with FFR last while Inflation first 
# followed by Unemployment.
var.SW <- cbind(window(Inflation.ts, start = c(1960,1)),
                window(Unemployment.ts, start = c(1960,1)),
                window(Fedfunds.ts, start = c(1960,1))) 
var.SW
# Changing column names
dimnames(var.SW)[[2]] <- c("Inflation","Unemployment","Fedfunds")
var.SW

# VAR Select with a constant term just as SW(2001)
VARselect(var.SW, type ="const") # AIC: 10 lags while BIC/SIC: 2 lags

# SC/BIC says 2 lags
# AIC says 10 lags
var.SW.levels <- VAR(var.SW, p = 4, type ="const") 
summary(var.SW.levels)
serial.test(var.SW.levels) 


roots(var.SW.levels) # roots stable. all roots within the unit disk.


#################################################################################
# IRF plots.
#
# SW(2001) Replication. All Variables in Levels because one can not mix
# transformed series and series in levels in one single VAR system.
#
# IRF Plots. In total Stock Watson has generated
# a 3*3 IRF plots.
# 
# Generate 66% Confidence intervals just as SW(2001) for each impulse response
#################################################################################
# Plot 1 :  Response of Inflation from shock to Inflation
irf(var.SW.levels, impulse ="Inflation", response = "Inflation", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Inflation", response = "Inflation", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: Inflation Shock to Inflation")

# Plot 2 : Inflation shock to Unemployment
irf(var.SW.levels, impulse ="Inflation", response = "Unemployment", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Inflation", response = "Unemployment", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: Inflation Shock to Unemployment")


# Plot 3 : Inflation shock to FED FUnds Rate
irf(var.SW.levels, impulse ="Inflation", response = "Fedfunds", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Inflation", response = "Fedfunds", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: Inflation Shock to FFR")


#######
# Plot 4 : Unemployment shock to Inflation
irf(var.SW.levels, impulse ="Unemployment", response = "Inflation", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Unemployment", response = "Inflation", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: Unemployment Shock to Inflation")


# Plot 5 : Unemployment shock to Unemployment
irf(var.SW.levels, impulse ="Unemployment", response = "Unemployment", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Unemployment", response = "Unemployment", n.ahead = 24, seed = 4654,ci=0.66 ),
     ylab="Percent",main ="Orthogonal IRF: Unemployment Shock to Unemployment")


# Plot 6 : Unemployment shock to Fed Funds Rate
irf(var.SW.levels, impulse ="Unemployment", response = "Fedfunds", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Unemployment", response = "Fedfunds", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: Unemployment Shock to FFR")

#######

# Plot 7 : Fed Funds rate shock to Inflation
irf(var.SW.levels, impulse ="Fedfunds", response = "Inflation", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Fedfunds", response = "Inflation", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: FFR Shock to Inflation")


# Plot 8 : Fed Funds rate shock to Unemployment
irf(var.SW.levels, impulse ="Fedfunds", response = "Unemployment", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Fedfunds", response = "Unemployment", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: FFR Shock to Unemployment")


# Plot 9 : Fed Funds rate shock to FEd Funds Rate
irf(var.SW.levels, impulse ="Fedfunds", response = "Fedfunds", n.ahead = 24, seed = 4654)
plot(irf(var.SW.levels, impulse ="Fedfunds", response = "Fedfunds", n.ahead = 24, seed = 4654, ci=0.66),
     ylab="Percent",main ="Orthogonal IRF: FFR shock to FFR")


#############################################################
# FEVD of model "var.SW.levels"

# With ordering: Inflation, Unemployment, Fed Funds Rate

# Similar to SW(2001) we utilize forecast horizon 12

#############################################################

# The below computed Variance decompositions are broadly 
# very similar to the computed Variance decompositions in SW(2001)
# in Tables 1B.i - 1B.iii from SW paper.

# Variance decompositions
fevd(var.SW.levels, n.ahead=12)$Inflation
fevd(var.SW.levels, n.ahead=12)$Unemployment
fevd(var.SW.levels, n.ahead=12)$Fedfunds

#############################################################
##############################################################
# END
#############################################################
#############################################################




