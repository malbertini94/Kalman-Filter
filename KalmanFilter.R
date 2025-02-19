library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(knitr)
library(readxl)
library(ivreg)

## Set up function

KalmanFilter <- function(Y, W_00, P_00, Q, R, f, H, t) {
  # initialize Vector
  W =c(W_00)
  P =c(P_00)
  # Compute W_tt and P_tt
  for (i in 1:t) {
    mu1  <- f%*%W[i]
    mu2  <- t(H) %*% mu1
    S_11 <- f^2 %*% P[i] + Q
    S_22 <- t(H)%*%S_11%*%H + R
    S_12 <- S_11%*%H
    S_21 <- t(S_12)
    K    <- S_12 %*% solve(S_22)
    eta  <- matrix(c(Y[i]), nrow=1, ncol=1) - mu2
    W_tt <- mu1 + K %*% eta
    P_tt <- S_11 - K %*% t(H) %*% S_11
    # Update Vectors
    W = append(W, W_tt)
    P = append(P, P_tt)
    }
  # Prepare Results
  results <-as.data.frame(cbind(W[2:t], P[2:t]))
  KF <<- results
}

# 1) Y is the data;
# 2) W_00 and P_00 are initial conditions for $W_t$;
# 3) H is the coefficient of the space equation, i.e. 1;
# 4) Q, R are variance covariance matrices respectively of the state and of the space
# 5) f is $\phi_1$
#6) t is the number of recursions.

The parameters H, Q, R, f can be fully understood from Time Serie Analysis, James D. Hamilton (chp 13, p. 372).

## Example

# Create a time series where outcome W_t is an arima
dates <- seq(as.Date("1916-1-1"), as.Date("2015-1-1"), by = "years")

Wt <- arima.sim(list(order=c(1,0,0), ar=.5), n=100)

# Yt=W_t+noise 
Yt  <- 2*Wt+rnorm(100)

Yt <- data.frame(dates,Yt)
colnames(Yt) <- c("dates","Y")

## Visualize the filter

# The unfiltered process
ggplot() + 
  geom_line(data = Yt, aes(x = dates, y = Y), color = "red") 

# apply the filter
KalmanFilter(Y=Yt$Y, W_00=0, P_00=1/(1-0.5^2), Q=1, R=1, f=0.5, H=2, t=100)

# Visualize the filtered process
Yt <- Yt[-100,]

mydata <- data.frame(Yt,KF)

ggplot() + 
  geom_line(data = mydata, aes(x = dates, y = Y), color = "red") +
  geom_line(data = mydata, aes(x = dates, y = V1), color = "blue")

# Where the blue process is our estimated signal. As by definition, indeed Y was driven by the process of W plus noise.

