###########################Thesis Geothermal ###################################

library(tidyverse)
library(reshape2)

#Constants
#kt <-              total heat transfer coefficient [W/m^2 K]
#lambda <-          thermal conductivity s from rock [w/m K]
#ro <-              density [kg/m^3]
#cp <-              specific heat capacity [J/kg K]
#as <-              thermal diffustivity [m^2/s]
#rw <-              external radius of  [m]
#t'  <-              time [s]

#h  <-              convective heat transfer coefficient
#mu <-              dynamic viscosity [Pa s]
#lambda_f<-         thermal conductivity fluid [W/m K]
#v  <-              kinamtic viscosity [m^2/s]
#D  <-              Diameter [m]
#rc <-              casing radius [m]


####################making a dataset for BEDUM-05 ##############################


dat1 <- data.frame(1:4979)
dat1 <- dat1 %>%
  rename(z=X1.4979)

dat1$rock_type <- c("a")

dat1 <- dat1 %>% 
  rows_update(tibble(z = 1:700, rock_type = "sandstone")) %>%
  rows_update(tibble(z = 701:1000, rock_type = "claystone")) %>%
  rows_update(tibble(z = 1001:1900, rock_type = "limestone")) %>%
  rows_update(tibble(z = 1901:2450, rock_type = "claystone")) %>%
  rows_update(tibble(z = 2451:4600, rock_type = "salt")) %>%
  rows_update(tibble(z = 4601:4979, rock_type = "sandstone"))

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

dat1$Tw <- seq(from = 283.1, to= 383, by = 0.02006427 )
dat1$delta_z <- c(1)

dat1$D <- c(0)
dat1 <- dat1 %>% 
  rows_update(tibble(z = 1:2000, D = 7)) %>%
  rows_update(tibble(z = 2001:4000, D = 6.8)) %>%
  rows_update(tibble(z = 4001:4979, D = 6.7)) 

dat1 <- dat1 %>%
  mutate(D = D * 0.0254) %>%
  mutate(rw = D * 0.5)

dat1$Tf_down <-c(0)
dat1$Tf_down[1] <- 293
dat1$Q_down <-c(0)

###################### fluid data frame ########################################

fluid <- c("water", "diathermic oil")
ro <- c(1000, 762)
lambda_f <- c(0.67, 0.13)
cp <- c(4186, 2500)
rc <- c( 7, 7)
v <- c(1, 3.3 )
mu <- c(8.90*10^-4, 4*10^-4)

dat2 <- data.frame(fluid, ro, lambda_f, cp, rc, v, mu)
dat2 <- dat2 %>%
  mutate(rc = rc * 0.0254 *0.5) 
dat2 <- dat2 %>%
  slice(1)

##################### formulas ################################################

#as thermal diffusitivity usded for Rs. 
get_as <- function(dat){
  as <- 0 
  as <- dat$lambda/(dat$ro*dat$cp)
  return(as)
}

get_RS <- function(dat, t=1){
  rs <- 0
  
  as <- get_as(dat)
  rs <- (1/(2*dat$lambda))*log((2*sqrt(as*t))/dat$rw)
  dat$rs <- rs
  
  return(dat)
}

#rs function for fig 3. 
rs2 <- function(dat, t){
  rs <- 0
  
  as <- get_as(dat)
  rs <- (1/(2*dat$lambda))*log((2*sqrt(as*t))/dat$rw)
  
  return(rs)
  
}

#Calculating convective heat transfer coefficient
get_h <- function(dat){
  h  <- 0
  Pr <- 0
  Re <- 0 
  
  Pr <- dat$cp*dat$mu/dat$lambda_f
  Re <- (dat$ro*dat$v*(dat$rc*2))/dat$mu
  h <- ((0.023*dat$lambda_f*(Re^0.8))*(Pr^0.4))/(2*dat$rc)
  return(h)
} 


#Ra thermal resistance due to the heat tranfer by convection into the pipe
get_Ra <- function(dat){
  Ra <- 0
  h  <- 0
  
  h <- get_h(dat)
  Ra <- 1/(2*dat$rc*h)
  
  #df1$Ra <- Ra 
  
  return(Ra)
}

#kt total heat transfer coefficient. 
get_kt <- function(dat1, dat2, t = 0.0435){
  kt <- 0 
  
  as <- get_as(dat1)
  h  <- get_h(dat2) 
  
  kt <- 1/((dat1$D/(2*dat1$lambda))*log((4*sqrt(as*t))/dat1$D)+(1/h))
  return(kt)
  
}


#Calculating Q_down and Tf_down 
get_Q_and_Tf <- function(dat1){

  end <- nrow(dat1)-1
for (i in 1:end) {
  newQ <- 2*pi*dat1$rw[i]*dat1$kt[i]*(dat1$Tw[i]-dat1$Tf_down[i])*dat1$delta_z[i]
  print(newQ)
  
  newT <- dat1$Tw[i+1]-(newQ/(2*pi*dat1$rw[i+1]*dat1$kt[i+1]*dat1$delta_z[i+1]))
  dat1$Tf_down[i+1] <- newT
  print(newT)
  
  # dat1$Q_down <- newQ
  # dat1$Tf_down <- newT
  
  
  return(dat1$Q_down)
  return(dat1$Tf_down)
  }
}

############################## Recreating fig. 3 ###############################

time <- seq(0,378683112, by=15778463)

for_t <- dat1 %>%
  slice(200, 1200, 2200, 3200, 4200)

inf1 <- rs2(for_t[1,], time)
inf2 <- rs2(for_t[2,], time)
inf3 <- rs2(for_t[3,], time)
inf4 <- rs2(for_t[4,], time)


inf_radius <- data.frame(inf1,inf2,inf3,inf4)
inf_radius[1,] <- 0
inf_radius <- inf_radius %>%
 mutate(years = seq(0,12, by=0.5)) 

inf_radius_trans <- melt(inf_radius, id.vars = c("years"))

inf_radius_trans$variable <- recode_factor(inf_radius_trans$variable, inf1 = "Sandstone at depth 200 m", 
                                inf2 = "Limestone at depth 1200 m", inf3 = "Claystone at depth 2200 m", 
                                inf4 = "Salt at depth 3200 m")


p1 <- ggplot(inf_radius_trans, aes(x=years, y = value, fill=variable, color = variable)) +
  geom_line(size=1) +# ggtitle("Influence radius over time") +
  xlab("Years") + ylab("influence radius (m)")  + theme_bw(base_size = 16) +
  theme(text=element_text(size=30))

p1

######################### Calculations for Well ################################

dat1$kt <- get_kt(dat1,dat2)

dat3 <- dat1
dat1 <- dat3

#dat1 <- get_Q_and_Tf(dat1)
# 
# dat1$Tf_down <- 0
# dat1$Tf_down[1] <- 293.0
# dat1$Q_down <- 0
# 

#end <- nrow(dat1)-1 
#end <- 3
# #dat1 <- 
# 


dat1$Tf_down <- 0
dat1$Tf_down[1] <- 293.0
dat1$Q_down <- 0
dat1$Q_down[1] <- 2*pi*dat1$rw[1]*dat1$kt[1]*(dat1$Tw[1]-dat1$Tf_down[1])*dat1$delta_z[1]



end <- 4
end <- nrow(dat1)-1
for (i in 2:end) {
  #Tw <- dat1$Tw[i]
  
  newT <- dat1$Tw[i]-(dat1$Q_down[i-1]/(2*pi*dat1$rw[i]*dat1$kt[i]*dat1$delta_z[i]))
  dat1$Tf_down[i-1] <- newT
  print(newT)
  
  newQ <- 2*pi*dat1$rw[i]*dat1$kt[i]*(dat1$Tw[i]-newT)*dat1$delta_z[i]
  print(newQ)
  dat1$Q_down[i] <- newQ
  
  #dat1$Q_down[i] <- 2*pi*dat1$rw*dat1$kt*(dat1$Tw-dat1$Tf_down[i])*dat1$delta_z
  #newQ <- dat1$Q_down[i]
  #j <- i+1
  #dat1$Tf_down[j] <- dat1$Tw-(newQ/(2*pi*dat1$rw*dat1$kt*dat1$delta_z))
  # newT <-dat1$Tf_down[j]
  
}



ggplot(dat1, aes(x = z, y = Q_down)) + geom_line()
ggplot(dat1, aes(x = z, y = Tf_down)) + geom_line()








#################Calculations where T changes#################
##########THIS IS WORKING NOW FOR, HOwever it should be inverted################
dat1 <- dat3

#end <- 3
end <- nrow(dat1)-1
for (i in 2:end) {
  #Tw <- dat1$Tw[i]
  
  
  newQ <- 2*pi*dat1$rw[i]*dat1$kt[i]*(dat1$Tw[i-1]-dat1$Tf_down[i-1])*dat1$delta_z[i]
  print(newQ)
  dat1$Q_down[i] <- newQ
  
  newT <- dat1$Tw[i+1]-(dat1$Q_down[i]/(2*pi*dat1$rw[i+1]*dat1$kt[i+1]*dat1$delta_z[i+1]))
  dat1$Tf_down[i] <- newT
  print(newT)
  
  #dat1$Q_down[i] <- 2*pi*dat1$rw*dat1$kt*(dat1$Tw-dat1$Tf_down[i])*dat1$delta_z
  #newQ <- dat1$Q_down[i]
  #j <- i+1
  #dat1$Tf_down[j] <- dat1$Tw-(newQ/(2*pi*dat1$rw*dat1$kt*dat1$delta_z))
  # newT <-dat1$Tf_down[j]
  
}

ggplot(dat1, aes(x = z, y = Q_down)) + geom_line()
ggplot(dat1, aes(x = z, y = Tf_down)) + geom_line()

##################Johannes' code #####################################
for (i in 2:nrow(dat1)) {
  dat1$Q_down[i]<- 2*pi*dat1$rw[i]*dat1$kt[i]*(dat1$Tw[i-1]-dat1$Tf_down[i-1])
  
  dat1$Tf_down[i] <- dat1$Tw[i]-(dat1$Q_down[i]/(2*pi*dat1$rw[i]*dat1$kt[i]))
}

ggplot(dat1, aes(x = z, y = Q_down)) + geom_line()
ggplot(dat1, aes(x = z, y = Tf_down)) + geom_line()

#####################################################################

dat1 <- dat3
dat1$Q_down[1] <- 2*pi*dat1$rw[1]*dat1$kt[1]*(dat1$Tw[1]-dat1$Tf_down[1])*dat1$delta_z[1]
dat1$Tf_down[1] <- 0

end <- nrow(dat1)-1
for (i in 2:end) {
  
  newT <- dat1$Tw[i-1]-(dat1$Q_down[i-1]/(2*pi*dat1$rw[i]*dat1$kt[i]*dat1$delta_z[i]))
  dat1$Tf_down[i] <- newT
  print(newT)
  
  newQ <- 2*pi*dat1$rw[i+1]*dat1$kt[i+1]*(dat1$Tw[i+1]-dat1$Tf_down[i])*dat1$delta_z[i+1]
  print(newQ)
  dat1$Q_down[i] <- newQ
  
  #original
   # newQ <- 2*pi*dat1$rw[i]*dat1$kt[i]*(dat1$Tw[i-1]-dat1$Tf_down[i-1])*dat1$delta_z[i]
   # print(newQ)
   # dat1$Q_down[i] <- newQ
   # 
   # newT <- dat1$Tw[i+1]-(dat1$Q_down[i]/(2*pi*dat1$rw[i+1]*dat1$kt[i+1]*dat1$delta_z[i+1]))
   # dat1$Tf_down[i] <- newT
   # print(newT)
   # 

  
}

ggplot(dat1, aes(x = z, y = Q_down)) + geom_line()
ggplot(dat1, aes(x = z, y = Tf_down)) + geom_line()





