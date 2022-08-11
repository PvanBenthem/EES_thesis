################################# TANG #########################################

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

V <- 	3.859789 #Velocity [m/s]
T_in <- 293         #Temperature ingoing fluid in Kelvin [K]
T_eio <- 283.1
Tb <- 473
t <- 300
L <- 6010
g_G <- (Tb-T_eio)/L



#other variables 
alpha <- 10^-6

#####################CREATING WELL DATAFRAME ###################################
dat1 <- data.frame(1:6010)
dat1 <- dat1 %>%
  rename(z=X1.6010)

dat1$rock_type <- c("a")

dat1 <- dat1 %>% 
  rows_update(tibble(z = 1:800, rock_type = "limestone")) %>%
  rows_update(tibble(z = 801:1700, rock_type = "sandstone")) %>%
  rows_update(tibble(z = 1701:2014, rock_type = "claystone")) %>%
  rows_update(tibble(z = 2015:2691, rock_type = "salt")) %>%
  rows_update(tibble(z = 2692:2741, rock_type = "claystone")) %>%
  rows_update(tibble(z = 2742:3117, rock_type = "sandstone")) %>%
  rows_update(tibble(z = 3118:3975, rock_type = "claystone")) %>%
  rows_update(tibble(z = 3975:6010, rock_type = "sandstone"))

dat1 <- dat1 %>%
  mutate(cp = case_when(
    startsWith(rock_type, "san") ~ 775,
    startsWith(rock_type, "cl") ~ 860,
    startsWith(rock_type, "li") ~ 680,
    startsWith(rock_type, "sal") ~ 880
  ))

dat1 <- dat1 %>%
  mutate(lambda = case_when(
    startsWith(rock_type, "san") ~ 2.75,
    startsWith(rock_type, "cl") ~ 1.10,
    startsWith(rock_type, "li") ~ 2.45,
    startsWith(rock_type, "sal") ~ 6
  ))

dat1 <- dat1 %>%
  mutate(ro = case_when(
    startsWith(rock_type, "san") ~ 2640,
    startsWith(rock_type, "cl") ~ 2680,
    startsWith(rock_type, "li") ~ 2760,
    startsWith(rock_type, "sal") ~ 2160
  ))


dat1$cp <- mean(dat1$cp)
dat1$lambda <- mean(dat1$lambda)
dat1$ro <- mean(dat1$ro)


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
#   rows_update(tibble(z = 4581:6010, r_wb = 3 * 0.0254 )) 
dat1$t    <- 0.25 * 0.0254
dat1$Tw <- seq(from = T_eio, to= (Tb-g_G), by = g_G )
dat1$T_ei <- T_eio + g_G * dat1$z


dat1 <- dat1 %>%
  mutate(A_i = pi*((r_ci^2) - ((r_in)^2))) %>%
  mutate(A_p = pi*r_ti^2) %>%
  mutate(de_p = 2 * r_ti) %>%
  mutate(de_i = 2* sqrt((r_ci^2)-(r_in^2)))



################## CREATING FLUID DATAFRAME ####################################
dat2 <- data.frame(1:6010)
dat2 <- dat2 %>%
  rename(z = X1.6010)  

dat2$V <- V
dat2$delta_z <- c(1)
dat2$P_i <- c(2)
dat2$dh_i <- 2 * sqrt((dat1$r_ci^2)-((dat1$r_ti+dat1$t)^2))
dat2$dh_p <- 2 * dat1$r_ti

dat2 <- dat2 %>%
  mutate(P_i = P_i *100000)

# for water K 293 (liquid)

dat2$ro <- 998.61
dat2$mu <- 0.0010012
dat2$Cp <- 4.18200 * 10^3
dat2$k <- 0.59887

########################  FUNCTIONS     ########################################

get_Re_i <- function(dat2, V = dat2$V) {
  
  Re <- 0
  Re <- (dat2$ro*V*dat2$dh_i)/dat2$mu
  dat2$Re_i <- Re
}

get_Re_p <- function(dat2, V = dat2$V ) {
  
  Re <- 0
  Re <- (dat2$ro*V*dat2$dh_p)/dat2$mu
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


get_M <- function(dat1,dat2, V = dat2$V){
  M <- 0 
  M <- (dat1$A_p*dat2$ro*dat2$Cp*V)/(2*pi*dat1$r_ti*dat2$U_to)
  dat2$M <- M
}

get_tau_D <- function(dat1){
  tau_D <- 0 
  tau_D <- (dat1$lambda * t)/(dat1$ro * dat1$cp * (dat1$r_wb^2))
  dat1$tau_D <- tau_D
}

get_T_D <- function(dat1, dat2) {
  T_D <- 0
  T_D <- log((exp(-0.2*dat1$tau_D))+(1.5-0.3719*(exp(-dat1$tau_D))*sqrt(dat1$tau_D)))
  dat2$T_D <- T_D
}

get_N <- function(dat1, dat2, V = dat2$V, ke = 1.8){ #Not yet completed waiting for TD
  N <- 0 
  N <- (dat1$A_p*dat2$ro*dat2$Cp*V*(ke + dat1$r_co*dat2$U_co*dat2$T_D))/
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

get_a <- function(dat2) {
  a <- 0 
  a <- -(((T_in - T_eio)*dat2$l_2*exp(dat2$l_2*L) + g_G * (1 - dat2$l_2*dat2$M))/
           (dat2$l_1*exp(dat2$l_1*L)*(1-dat2$l_2*dat2$M)-dat2$l_2*exp(dat2$l_2*L) *(1-dat2$l_1*dat2$M))   )
  dat2$a <- a
}

get_b <- function(dat2){
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



################ Calculating well properties###################################

#fluid coefficients 
dat2$Re_i <- get_Re_i(dat2)
dat2$Re_p <- get_Re_p(dat2)
dat2$Pr <- get_Pr(dat2)
dat2$hf_i <- get_hf_i(dat1,dat2)
dat2$hf_p <- get_hf_p(dat1,dat2)
dat2$U_to <- get_U_to(dat1,dat2)
dat2$U_co <- get_U_co(dat1,dat2)


#time calculations (dimensionless)
dat1$tau_D <- get_tau_D(dat1)
dat2$T_D <- get_T_D(dat1, dat2)

#Calculating constants 
dat2$M <- get_M(dat1,dat2)
dat2$N <- get_N(dat1, dat2)
dat2$l_1 <- get_l_1(dat2)
dat2$l_2 <- get_l_2(dat2)
dat2$a <- get_a(dat2)
dat2$b <- get_b(dat2)

#Get Temperatures of wells
dat1$Tp <- get_Tp(dat1, dat2)
dat1$Ti <- get_Ti(dat1, dat2)

############################# PLOTS ###########################################


datp1 <- dat1$z 
datp1 <- as.data.frame(datp1)
datp1 <- datp1 %>%
  rename(z = datp1) 
datp1$Rock <- dat1$Tw -273
datp1$T_down <- dat1$Ti -273
datp1$T_up <- dat1$Tp - 273

datp1 <- datp1 %>%
  rename(z = z, Rock = Rock, "Fluid downward" = T_down, "Fluid upward" = T_up) %>%
  gather(key = 'variable', value = "value", -z)

p1 <- ggplot(datp1, aes(x=value, y = z)) +
  geom_point(aes(color = variable),size=0.01) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(values = c("#0072B2","#D55E00","#999999")) +
  scale_y_reverse() +# ggtitle("Water and rock temperature versus depth") +
  xlab("Temperature (°C)") + ylab("Depth (m)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), position = "top")  +
  theme_bw(base_size = 16) +
  theme(text=element_text(size=30)) +
  theme(axis.title.x = element_text(margin = margin(t = 0, r =0 , b = 20, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r =20 , b = 0, l = 0)))

p1


dat1$T_i_C <- dat1$Ti -273

#########################TEST AREA ############################################

