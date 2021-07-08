library(tidyverse)

# boundary conditions
wtr_profile <- c(25,25,25,25,25,24,23,20,18,17,16,15,15,15,15,15,15)
depth_profile <- seq(0,length(wtr_profile)-1,1)
area_profile <- seq(from = 150000, to = 0, length.out = length(wtr_profile))
volume_profile <- area_profile * mean(diff(depth_profile))

calc_dens <- function(wtemp){
  dens = 999.842594 + (6.793952 * 10^-2 * wtemp) - (9.095290 * 10^-3 *wtemp^2) + (1.001685 * 10^-4 * wtemp^3) - (1.120083 * 10^-6* wtemp^4) + (6.536336 * 10^-9 * wtemp^5)
  return(dens)
}

# inflow parameters
cd <- 0.016
slope <- 1.5
angle <- 65
Q_in_0 <- 3
wtr_in_0 <- 15

# initial entrainment
Ri_in <- (cd * (1 + 0.21 * sqrt(cd) * sin(angle))) / (sin(angle) * tan(slope))
E_in <- 1.6 * (cd^(3/2)) / Ri_in

# entrainment algorithm according to Antenucci (2005), Ayla (2014), Fischer (1979), Hipsey (2019)
g_in = 9.81 * (calc_dens(wtr_in_0) - calc_dens(wtr_profile[1])) / calc_dens(wtr_profile[1])
delta_z_in <- (2 * (Ri_in / g_in) * (Q_in_0 / tan(angle))^2)^(1/5)

delta_z = rep(NA, 100)
delta_z[1] <- delta_z_in
delta_x = rep(NA, 100)
delta_Q = rep(NA, 100)
Q_in <- rep(NA, 100)
Q_in[1] <- Q_in_0
h <- rep(NA, 100)
h[1] <- 0
h <- max(depth_profile) - h
z <- rep(NA, 100)
z[1] <- 0
z = max(depth_profile)  - z
partial_z <- rep(NA, 100)
wtr_in <- rep(NA, 100)
wtr_in[1] <- wtr_in_0
hs <- max(depth_profile) -min(depth_profile)
for (j in 2:100){
  partial_z[j] = (hs - z[j-1]) - depth_profile[ which.min(((max(depth_profile) - depth_profile) - z[j-1])^2) + 1]
  delta_x[j] = partial_z[j] / sin(slope)
  delta_z[j] = 1.2 * E_in * delta_x[j] + delta_z[j-1]
  z[j] <- z[j-1] + delta_x[j] * sin(angle)
  
  delta_Q[j] <- Q_in[j-1] * (abs(delta_z[j] / delta_z[j-1])^(5/3) - 1)
  Q_in[j] <- Q_in[j-1] + delta_Q[j]

  Mtotal <- Q_in[j] * calc_dens(wtr_in[j-1]) + 
    delta_Q[j] * calc_dens((wtr_profile[which.min(((max(depth_profile) - depth_profile )- 
                                                            (abs(z[j])))^2)])
                             )
  
  wtr_in[j] <- (wtr_in[j-1] * Q_in[j] * calc_dens(wtr_in[j-1]) +
                  delta_Q[j] * calc_dens((wtr_profile[which.min(((max(depth_profile) - depth_profile )-  (abs(z[j])))^2)])) * 
                  (wtr_profile[which.min(((max(depth_profile) - depth_profile )- (abs(z[j])))^2)])) /
    Mtotal
  
  if (calc_dens(wtr_in[j]) <= calc_dens(wtr_profile[which.min(((max(depth_profile) - depth_profile )- 
                                                            (abs(z[j])))^2)])){
    break
  }
}

print(depth_profile[which.min(((max(depth_profile) - depth_profile )- 
                                 (abs(z[j])))^2)])

# theoretical and experimental entrainment depth according to Wells & Nadarajah (2009)
B = 9.81 * (calc_dens(wtr_in_0) - calc_dens(wtr_profile[1])) / calc_dens(wtr_profile[1]) * Q_in_0
wtemp <- data.frame('Depth' = depth_profile, 'Temp' = wtr_profile, 'Density' =
                      calc_dens(wtr_profile))
wtempLong = wtemp %>% mutate(dp = Density - lag(Density)) %>% 
  mutate(dp = ifelse(Depth == 0, NA, dp))
wtempLong = wtempLong %>% mutate(finite.dp = lead(Density) - Density) %>% 
  mutate(finite.dp = ifelse(Depth == max(Depth), Density - lag(Density), finite.dp)) 
wtempLong = wtempLong %>% mutate(buoyFreq = sqrt(9.81/998.2 * (dp/1))) %>% 
  mutate(finite.buoyFreq = sqrt(9.81/998.2 * (finite.dp/1)))

# Find the max depth
maxDepths = wtempLong %>% 
  filter(finite.buoyFreq == max(finite.buoyFreq,na.rm = T)) 
N = maxDepths$finite.buoyFreq
  
Z_min =  (3 - 1) * B^(1/3) / N
Z_max = (3 + 1) * B^(1/3) / N

ggplot(wtemp, aes(x = Density, y = Depth)) +
  geom_line() +
  scale_y_reverse()  +
  geom_hline(aes(yintercept=
               depth_profile[which.min((calc_dens(wtr_profile) - calc_dens(wtr_in_0))^2)],
               linetype = 'Neutral buoyancy'),
             color = "green") + 
  geom_hline(aes(yintercept=
                   depth_profile[which.min(((max(depth_profile) - depth_profile )- 
                                              (abs(z[j])))^2)], linetype = 'GLM algorithm'), 
             color = "red") + 
  geom_hline(aes(yintercept=
               Z_min,
               linetype = 'Wells and Nadarajah (2009)'),
             color = "blue") + 
  geom_hline(aes(yintercept=
               Z_max,
             linetype = 'Wells and Nadarajah (2009)'),
             color = "blue") + 
  scale_linetype_manual(name = "", values = c(2, 2, 2, 2), 
                        guide = guide_legend(override.aes = list(color = c("red", "green",
                                                                           'blue')))) +
  ggtitle('Entrainment depth under stratified conditions') +
  theme_minimal() +
  theme(legend.position="bottom") 
