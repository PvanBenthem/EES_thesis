################################# TANG #########################################

# first run bedum05_correct_mfr_&_P_sl
library(viridisLite)
library(viridis)
library(tidyverse)
library(grid)

##############################UNITS ###########################################
#R <-               Radius well [m]
#z <-               depth [m]
#r <-               inner radius [m]
#t <-               insulation thickness [m]
#V <-               Velocity [m/s]
#P <-               Pressure [Pa]
#dh_down  <-        hydraulic diameter [m]
#dh_up    <-        hydraulic diameter [m]
#mu <-              dynamic viscosity [Pa s]
#ro <-              density [kg/m^3]
#k <-               thermal conductivity [W/m K]
#Cp <-              specific heat [J/kg K]
#alpha <-           formation diffusivity [m^2 s^-1]
########################### Variables #########################################
#V <- 14.06588 #Velocity [m/s]
#V <- seq(from = 0.13, to = 20.01, by = 0.00399277)\
T_in <- 293         #Temperature ingoing fluid in Kelvin [K]
T_eio <- 283.1
t <- 300
g_G <- 100/4979
L <- 4979

Vf <- (27.5-0.13)/L
V <- seq(from = 0.13, to = 27.5-Vf, by = Vf)

#other variables |
alpha <- 10^-6

#####################CREATING WELL DATAFRAME ###################################
# dat1 <- data.frame(1:4979)
# dat1 <- dat1 %>%
#   rename(z=X1.4979)
# 
# dat1$rock_type <- c("a")
# 
# dat1 <- dat1 %>% 
#   rows_update(tibble(z = 1:700, rock_type = "sandstone")) %>%
#   rows_update(tibble(z = 701:1000, rock_type = "claystone")) %>%
#   rows_update(tibble(z = 1001:1900, rock_type = "limestone")) %>%
#   rows_update(tibble(z = 1901:2450, rock_type = "claystone")) %>%
#   rows_update(tibble(z = 2451:4600, rock_type = "salt")) %>%
#   rows_update(tibble(z = 4601:4979, rock_type = "sandstone"))
# 
# dat1 <- dat1 %>%
#   mutate(cp = case_when(
#     startsWith(rock_type, "san") ~ 775,
#     startsWith(rock_type, "cl") ~ 860,
#     startsWith(rock_type, "li") ~ 680,
#     startsWith(rock_type, "sal") ~ 880
#   ))
# 
# dat1 <- dat1 %>%
#   mutate(lambda = case_when(
#     startsWith(rock_type, "san") ~ 2.75,
#     startsWith(rock_type, "cl") ~ 1.10,
#     startsWith(rock_type, "li") ~ 2.45,
#     startsWith(rock_type, "sal") ~ 6
#   ))
# 
# dat1 <- dat1 %>%
#   mutate(ro = case_when(
#     startsWith(rock_type, "san") ~ 2640,
#     startsWith(rock_type, "cl") ~ 2680,
#     startsWith(rock_type, "li") ~ 2760,
#     startsWith(rock_type, "sal") ~ 2160
#  ))
dat1 <- data.frame(1:4979)
dat1 <- dat1 %>%
  rename(z=X1.4979)
dat1$z <- 1

dat1$r_ti <- 0.425 * 0.0254     #r  
dat1$r_to <- 0.475 * 0.0254
dat1$r_in <- 0.675 * 0.0254
dat1$r_ci <- 1.525 * 0.0254              #R 
dat1$r_co <- 1.75 * 0.0254
dat1$r_wb <- 5.624046 * 0.0254 #average of r_wb below
# dat1 <- dat1 %>%
#   rows_update(tibble(z = 1:91, r_wb = 15 * 0.0254 )) %>%
#   rows_update(tibble(z = 92:1030, r_wb = 8 * 0.0254 )) %>%
#   rows_update(tibble(z = 1031:2545, r_wb = 6.125 * 0.0254 )) %>%
#   rows_update(tibble(z = 2546:4580, r_wb = 4.25 * 0.0254 )) %>%
#   rows_update(tibble(z = 4581:4979, r_wb = 3 * 0.0254 ))
dat1$t    <- 0.25 * 0.0254
dat1$Tw <- seq(from = 283.1, to= 383, by = 0.02006427 )
dat1$T_ei <- T_eio + g_G * dat1$z


dat1 <- dat1 %>%
  mutate(A_i = pi*((r_ci^2) - ((r_in)^2))) %>%
  mutate(A_p = pi*r_ti^2) %>%
  mutate(de_p = 2 * r_ti) %>%
  mutate(de_i = 2* sqrt((r_ci^2)-(r_in^2)))

dat1$avg_cp <- 817.6793
dat1$avg_lambda <- 3.817483
dat1$avg_ro <- 2461.249

datP <- data.frame(1:4979)
datP <- datP %>%
  rename(z=X1.4979)

datP$z <- 1
datP$V <- V

################## CREATING FLUID DATAFRAME ####################################
dat2 <- data.frame(1:4979)
dat2 <- dat2 %>%
  rename(z = X1.4979)  
dat2$z <- 1

dat2$V <- V
dat2$delta_z <- c(1)
dat2$P_i <- c(2)
dat2$dh_i <- 2 * sqrt((dat1$r_ci^2)-((dat1$r_ti+dat1$t)^2))
dat2$dh_p <- 2 * dat1$r_ti

dat2 <- dat2 %>%
  mutate(P_i = P_i *100000)

# water K 293 (liquid)

dat2$ro <- 762
dat2$mu <- 0.0033
dat2$Cp <- 2.5 * 10^3
dat2$k <- 0.13


#for second loop ###############333333

dat3 <- data.frame(1:4979)
dat3 <- dat3 %>%
  rename(z = X1.4979)  

dat3$V <- V
dat3$delta_z <- c(1)
dat3$P_i <- c(2)
dat3$dh_i <- 2 * sqrt((dat1$r_ci^2)-((dat1$r_ti+dat1$t)^2))
dat3$dh_p <- 2 * dat1$r_ti

dat3 <- dat2 %>%
  mutate(P_i = P_i *100000)

# water K 343 (liquid)

dat3$ro <- 762
dat3$mu <- 0.0033
dat3$Cp <- 2.5 * 10^3
dat3$k <- 0.13



########################  FUNCTIONS     ########################################

get_Re_i <- function(dat2, datP) {
  
  Re <- 0
  Re <- (dat2$ro*datP$V*dat2$dh_i)/dat2$mu
  dat2$Re_i <- Re
}

get_Re_p <- function(dat2, datP ) {
  
  Re <- 0
  Re <- (dat2$ro*datP$V*dat2$dh_p)/dat2$mu
  dat2$Re_i <- Re
}

get_Pr <- function(dat2) {
  
  Pr <- 0 
  Pr <- (dat2$Cp*dat2$mu)/dat2$k
  dat2$Pr <- Pr
}

get_hf_i <- function(dat1,dat2){
  hf_i <- 0
  hf_i <- (0.023*dat2$k*(dat2$Re_i^0.8)*(dat2$Pr^0.4))/dat1$de_i
  
  dat2$hf_i <- hf_i
}

get_hf_p <- function(dat1,dat2){
  hf_p <- 0
  hf_p <- (0.023*dat2$k*(dat2$Re_p^0.8)*(dat2$Pr^0.4))/dat1$de_p
  
  dat2$hf_p <- hf_p
}

get_U_to <- function(dat1, dat2, kt = 57, kin = 0.023){
  U_to <- 0
  U_to <- 1/( (dat1$r_to/(dat1$r_ti*dat2$hf_p)) + ((dat1$r_to*log(dat1$r_to/dat1$r_ci))/kt)+
                ((dat1$r_to*log(dat1$r_in/dat1$r_to))/kin) + (dat1$r_to/(dat1$r_in*dat2$hf_p)) )
  dat2$U_to <- U_to
}

get_U_co <- function(dat1, dat2, kc = 57, kcem = 2.3 ){
  U_co <- 0
  U_co <- 1/( (dat1$r_co/(dat1$r_ci * dat2$hf_i)) + ((dat1$r_co*log(dat1$r_co/dat1$r_ci))/kc) +
                (dat1$r_co*log(dat1$r_wb/dat1$r_co)/kcem)  )
  dat2$U_co <- U_co
}


get_M <- function(dat1,dat2, datP){
  M <- 0 
  M <- (dat1$A_p*dat2$ro*dat2$Cp*datP$V)/(2*pi*dat1$r_ti*dat2$U_to)
  dat2$M <- M
}

get_tau_D <- function(dat1){
  tau_D <- 0 
  tau_D <- (dat1$avg_lambda * t)/(dat1$avg_ro * dat1$avg_cp * (dat1$r_wb^2))
  dat1$tau_D <- tau_D
}

get_T_D <- function(dat1, dat2) {
  T_D <- 0
  T_D <- log((exp(-0.2*dat1$tau_D))+(1.5-0.3719*(exp(-dat1$tau_D))*sqrt(dat1$tau_D)))
  dat2$T_D <- T_D
}

get_N <- function(dat1, dat2, datP, ke = 1.8){ #Not yet completed waiting for TD
  N <- 0 
  N <- (dat1$A_p*dat2$ro*dat2$Cp*datP$V*(ke + dat1$r_co*dat2$U_co*dat2$T_D))/
    (2*pi*dat1$r_co*dat2$U_co*ke)
  dat2$N <- N
}

get_l_1 <- function(dat2){
  l_1 <- 0 
  l_1 <- -(1/(2*dat2$N)) + (1/(2*dat2$N))*sqrt(1+((4*dat2$N)/dat2$M))
  dat2$l_1 <- l_1
}

get_l_2 <- function(dat2){
  l_2 <- 0 
  l_2 <- -(1/(2*dat2$N)) - (1/(2*dat2$N))*sqrt(1+((4*dat2$N)/dat2$M))
  dat2$l_2 <- l_2
}

get_a <- function(dat2, T_in) {
  a <- 0 
  a <- -(((T_in - T_eio)*dat2$l_2*exp(dat2$l_2*L) + g_G * (1 - dat2$l_2*dat2$M))/
           (dat2$l_1*exp(dat2$l_1*L)*(1-dat2$l_2*dat2$M)-dat2$l_2*exp(dat2$l_2*L) *(1-dat2$l_1*dat2$M))   )
  dat2$a <- a
}

get_b <- function(dat2,  T_in){
  b <- 0 
  b <- ((T_in-T_eio)*dat2$l_1*exp(dat2$l_1*L)+g_G*(1-dat2$l_1*dat2$M))/
    (dat2$l_1*exp(dat2$l_1*L)*(1-dat2$l_2*dat2$M)-dat2$l_2*exp(dat2$l_2*L) *(1-dat2$l_1*dat2$M)) 
  dat2$b <- b
}

get_Tp <- function(dat1, dat2){
  Tp <- 0 
  Tp <- dat2$a * exp(dat2$l_1*dat2$z) + dat2$b * exp(dat2$l_2 * dat2$z) +
    dat2$M * g_G + dat1$T_ei
  dat2 <- Tp
}

get_Ti <- function(dat1, dat2) {
  Ti <- 0 
  Ti <- (1 - dat2$l_1*dat2$M)*dat2$a*exp(dat2$l_1*dat2$z) + (1 - dat2$l_2*dat2$M) *
    dat2$b * exp(dat2$l_2 * dat2$z) + dat1$T_ei
  dat1$Ti <- Ti
}

get_a2 <- function(dat2, datP) {
  a <- 0 
  a <- -(((datP$Tp - T_eio)*dat2$l_2*exp(dat2$l_2*L) + g_G * (1 - dat2$l_2*dat2$M))/
           (dat2$l_1*exp(dat2$l_1*L)*(1-dat2$l_2*dat2$M)-dat2$l_2*exp(dat2$l_2*L) *(1-dat2$l_1*dat2$M))   )
  dat2$a <- a
}

get_b2 <- function(dat2,  datP){
  b <- 0 
  b <- ((datP$Tp-T_eio)*dat2$l_1*exp(dat2$l_1*L)+g_G*(1-dat2$l_1*dat2$M))/
    (dat2$l_1*exp(dat2$l_1*L)*(1-dat2$l_2*dat2$M)-dat2$l_2*exp(dat2$l_2*L) *(1-dat2$l_1*dat2$M)) 
  dat2$b <- b
}

get_mfr <- function (dat1, dat2, datP){
  mfr <- 0 
  mfr <- dat2$ro * dat1$A_p * datP$V
  datP$mfr <- mfr
}



################ Calculating well properties###################################

#fluid coefficients 
dat2$Re_i <- get_Re_i(dat2, datP)
dat2$Re_p <- get_Re_p(dat2, datP)
dat2$Pr <- get_Pr(dat2)
dat2$hf_i <- get_hf_i(dat1,dat2)
dat2$hf_p <- get_hf_p(dat1,dat2)
dat2$U_to <- get_U_to(dat1,dat2)
dat2$U_co <- get_U_co(dat1,dat2)


#time calculations (dimensionless)
dat1$tau_D <- get_tau_D(dat1)
dat2$T_D <- get_T_D(dat1, dat2)

#Calculating constants 
dat2$M <- get_M(dat1,dat2, datP)
dat2$N <- get_N(dat1, dat2, datP)
dat2$l_1 <- get_l_1(dat2)
dat2$l_2 <- get_l_2(dat2)
dat2$a <- get_a(dat2, T_in)
dat2$b <- get_b(dat2, T_in)

#Get Temperatures of wells
dat1$Tp <- get_Tp(dat1, dat2)
dat1$Ti <- get_Ti(dat1, dat2)

#################Second loop 
datP$Tp <- dat1$Tp
#fluid coefficients 
dat3$Re_i <- get_Re_i(dat3, datP)
dat3$Re_p <- get_Re_p(dat3, datP)
dat3$Pr <- get_Pr(dat3)
dat3$hf_i <- get_hf_i(dat1,dat3)
dat3$hf_p <- get_hf_p(dat1,dat3)
dat3$U_to <- get_U_to(dat1,dat3)
dat3$U_co <- get_U_co(dat1,dat3)


#time calculations (dimensionless)
dat1$tau_D <- get_tau_D(dat1)
dat3$T_D <- get_T_D(dat1, dat3)

#Calculating constants 
dat3$M <- get_M(dat1,dat3, datP)
dat3$N <- get_N(dat1, dat3, datP)
dat3$l_1 <- get_l_1(dat3)
dat3$l_2 <- get_l_2(dat3)
dat3$a <- get_a2(dat3, datP)
dat3$b <- get_b2(dat3, datP)

#Get Temperatures of wells
dat1$Tp2 <- get_Tp(dat1, dat3)
datP$Tp2 <- dat1$Tp2
dat1$Ti2 <- get_Ti(dat1, dat3) 

#get mfr
datP$mfr <- get_mfr(dat1, dat2, datP)



#########################TEST AREA ############################################
#Mass flow rate at outflow of well: density * velocity* Area ()
mf_out <- (dat2$ro[1]/1000) * V * dat2$dh_p[1] 
#Energy within the fluid (outflow temperature - input temperature): cp*(T_out-T-in)
E_1 <- dat3$Cp[1] * (dat1$Tp[1] - T_in)
E_2 <- dat3$Cp[1] * (dat1$Tp2[1] - T_in)
#Power output: mass flow rate * energy
P_1 <- mf_out * E_1
P_2 <- mf_out* E_2
## P in this case (V = 2.8) is 15.7 kW. with a V = 10, P=41.5 kW, with V = 20 P = 54.6 kW, V=100 P=69 kW -> the increase is non linear
# It also makes now a large difference if one or two wells are used (see case of V=20)
# There will probably be a physical limit to how fast the fluid can flow through the well (pumping power?, turbulent flow?)


#Mass flow rate at outflow of well: density * velocity* Area ()
#datP$mf_out <- dat2$ro * datP$V * dat2$dh_p
#Energy within the fluid (outflow temperature - input temperature): cp*(T_out-T-in)
datP$E_1 <- dat3$Cp[1] * (datP$Tp - T_in)
datP$E_2 <- dat3$Cp[1] * (datP$Tp2 - T_in)
#Power output: mass flow rate * energy
#datP$P_1 <- (datP$mf_out * datP$E_1)/1000
#datP$P_2 <- (datP$mf_out * datP$E_2)/1000


tmp <- datP

datP$E_1 <- dat3$Cp * (datP$Tp - T_in)
datP$E_2 <- dat3$Cp * (datP$Tp2 - T_in)

datP$P_1 <- (datP$mfr * datP$E_1)
datP$P_2 <- (datP$mfr * datP$E_2)


############################# PLOTS ###########################################

datp4 <- datP$V
datp4 <- as.data.frame(datp4)
datp4 <- datp4 %>%
  rename(V = datp4) 
datp4$Tp2 <- datP$Tp2 -273
datp4$Tp <- datP$Tp - 273
datp4$mfr <- datP$mfr



# datp4 <- datp4 %>%
#   rename(V = V, "Fluid production pipe" = Tp, "Mass flow rate (kg/s)" = mfr) %>%
#   gather(key = 'variable', value = "value", -z)

p4 <- ggplot(datp4, aes(x=mfr, y = Tp)) +
  geom_line(color = "steelblue", size=1.3) + #aes(color = variable)
  #scale_color_manual(values = c("steelblue","darkred","darkgreen")) +
  # scale_y_reverse() + 
  #ggtitle("Outlet temperature for different working fluid inlet velocities and two wells") +
  ylab("Temperature (°C)") + xlab("Mass Flow rate (kg/s)") +
  geom_vline(xintercept = 1.191936, linetype="dotted",  color = "darkred", size=1) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10), position = "top")  +
  theme_bw(base_size = 16) +
  theme(text=element_text(size=30)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r =0 , b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r =20 , b = 0, l = 0)))

p4

#### ploting Production . 
datp5 <- datP$mfr
datp5 <- as.data.frame(datp5)
datp5 <- datp5 %>%
  rename(mfr = datp5) 

datp5$P_1 <- datP$P_1 / 1000
datp5$P_2 <- datP$P_2 / 1000


datp5 <- datp5 %>%
  rename(mfr = mfr, "First well" = P_1, "Second well" = P_2) %>%
  gather(key = 'variable', value = "value", -mfr)

p5 <- ggplot(datp5, aes(x=mfr, y=value)) +
  geom_line(aes(color = variable),size = 1.3)+ #)aes(color = variable)) + 
  scale_color_manual(values = c("steelblue","darkred")) +
  # geom_line(fill = NA) +
  #scale_y_reverse() + 
  #ggtitle("Production 1 and 2 loops") +
  xlab("Mass Flow rate (kg/s)") + ylab("Thermal power (kW)") +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10), position = "top")  +
  theme_bw(base_size = 16) +
  theme(text=element_text(size=30)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r =0 , b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r =20 , b = 0, l = 0)))

p5

################## plot for oil vs water 

datp6 <- datPwater$mfr
datp6 <- as.data.frame(datp6)
datp6 <- datp6 %>%
  rename(mfr = datp6) 

datp6$P_1_oil <- datP$P_1 / 1000
#datp6$P_2 <- datP$P_2 / 1000
datp6$P_1_water <- datPwater$P_1 /1000


datp6 <- datp6 %>%
  rename(mfr = mfr, "Water" = P_1_water, "Oil" = P_1_oil) %>%
  gather(key = 'Fluid', value = "value", -mfr)

p6 <- ggplot(datp6, aes(x=mfr, y=value)) +
  geom_line(aes(color = Fluid),size = 1.3)+ #)aes(color = variable)) + 
  scale_color_manual(values = c("steelblue","darkred")) +
  # geom_line(fill = NA) +
  #scale_y_reverse() + 
  #ggtitle("Production 1 and 2 loops") +
  xlab("Mass Flow rate (kg/s)") + ylab("Thermal power (kW)") +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10), position = "top")  +
  theme_bw(base_size = 16) +
  theme(text=element_text(size=30)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r =0 , b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r =20 , b = 0, l = 0)))

p6





