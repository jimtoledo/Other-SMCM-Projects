# Primitive 2-dimensional COVID-19 PDE Model as reaction-diffusion system in an L x L km space
# Jim Toledo

rm(list=ls())

library(deSolve)
library(ReacTran)

#model parameters:
birth_rate <- 0 # number births per person per day
symp_rate <- 1/7 # inverse of incubation period (1/days)
asymp_rec_rate <- 1/6 # asymptomatic recovery rate (1/days)
inf_rec_rate <- 1/24 # infected recovery rate (1/days)
asymp_cont_rate <- 0.00005 # asymptomatic contact rate (1/person/day)
sympt_cont_rate <- 0.00005 # symptomatic contact rate (1/person/day)
inf_mort_rate <- 1/160 # infected mortality rate (1/days)
gen_mort_rate <- 0 # non-covid death rate (1/days)
A <- 1000 # Allee effect (number of people)

#diffusivity km^2/day

s_diffusivity <- 0.02

e_diffusivity <- 0.01

i_diffusivity <- 0.00001

r_diffusivity <- 0.01


#boundary conditions (Neumann)
left_flux = 0 # number of people/km^2/day
right_flux = 0
bottom_flux = 0
top_flux = 0


#grid definition (discretizing space domain)
x.N <- 20 # number of cells in x-direction
y.N <- 20 # number of cells in y-direction
x.L <- 1000 # domain size x-direction (km)
y.L <- 1000 # domain size y-direction (km)
dx <- x.L/x.N   # cell size x-direction (km)
dy <- y.L/y.N   # cell size y-direction (km)


#initial conditions

avg_pop_density <- 50 # number people per km^2 (in this model, we assume uniform population distribution)

#in a real model, the data would be scanned from some file and outputted to the matrices
susceptible <- matrix(data = dx*dy*avg_pop_density, nrow = x.N, ncol = y.N, byrow = FALSE)
exposed <- matrix(data = 0, nrow = x.N, ncol = y.N, byrow = FALSE)
infected <- matrix(data = 0, nrow = x.N, ncol = y.N, byrow = FALSE)
recovered <- matrix(data = 0, nrow = x.N, ncol = y.N, byrow = FALSE)
dead <- matrix(data = 0, nrow = x.N, ncol = y.N, byrow = FALSE)

#for initial condition, assume we start out with 0.1 percent infected in the center region
susceptible[c(x.N/2, x.N/2+1), c(y.N/2, y.N/2+1)] <- 0.999*dx*dy*avg_pop_density
infected[c(x.N/2, x.N/2+1), c(y.N/2, y.N/2+1)] <- 0.001*dx*dy*avg_pop_density

#initial condition if virus started from corner of the region
#susceptible[1, 1] <- 0.995*dx*dy*avg_pop_density
#infected[1, 1] <- 0.005*dx*dy*avg_pop_density

#random initial condition test
#for(x in 1:x.N) {
#  for(y in 1:y.N) {
#    random_dens <- runif(1, 20, 200) #random density between 20 and 200 people per square km
#    random_inf_rate <- runif(1, 0, 0.01) #random initial infected percentage between 0 and 1 percent
#    susceptible[x, y] <- (1-random_inf_rate)*dx*dy*random_dens
#    infected[x, y] <- random_inf_rate*dx*dy*random_dens
#  }
#}

#model PDE system:
seird <- function(t, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    #get state variable matrices from the state vector
    s <- matrix(nrow = x.N, ncol = y.N, byrow = FALSE, data = state[1:(x.N*y.N)])
    e <- matrix(nrow = x.N, ncol = y.N, byrow = FALSE, data = state[(x.N*y.N+1):(2*x.N*y.N)])
    i <- matrix(nrow = x.N, ncol = y.N, byrow = FALSE, data = state[(2*x.N*y.N+1):(3*x.N*y.N)])
    r <- matrix(nrow = x.N, ncol = y.N, byrow = FALSE, data = state[(3*x.N*y.N+1):(4*x.N*y.N)])
    d <- matrix(nrow = x.N, ncol = y.N, byrow = FALSE, data = state[(4*x.N*y.N+1):(5*x.N*y.N)])
    
    n <- s+e+i+r
    
    #diffusivity terms
    s_diff <- tran.2D(C = s, flux.x.up = left_flux, flux.x.down = right_flux, flux.y.up = bottom_flux, flux.y.down = top_flux, 
                           D.x = s_diffusivity, D.y = s_diffusivity, dx = dx, dy = dy)
    e_diff <- tran.2D(C = e, flux.x.up = left_flux, flux.x.down = right_flux, flux.y.up = bottom_flux, flux.y.down = top_flux,
                           D.x = e_diffusivity, D.y = e_diffusivity, dx = dx, dy = dy)
    i_diff <- tran.2D(C = i, flux.x.up = 0, flux.x.down = 0, flux.y.up = 0, flux.y.down = 0,
                           D.x = i_diffusivity, D.y = i_diffusivity, dx = dx, dy = dy) #assume infected persons never enter
    r_diff <- tran.2D(C = r, flux.x.up = left_flux, flux.x.down = right_flux, flux.y.up = bottom_flux, flux.y.down = top_flux,
                           D.x = r_diffusivity, D.y = r_diffusivity, dx = dx, dy = dy)
    
    #SEIRD equations
    dS <- birth_rate*n - (1-(A/n))*sympt_cont_rate*s*i - (1-(A/n))*asymp_cont_rate*s*e - gen_mort_rate*s + s_diff$dC
    dE <- (1-(A/n))*sympt_cont_rate*s*i + (1-(A/n))*asymp_cont_rate*s*e - symp_rate*e - asymp_rec_rate*e - gen_mort_rate*e + e_diff$dC
    dI <- symp_rate*e - inf_mort_rate*i - inf_rec_rate*i - gen_mort_rate*i + i_diff$dC
    dR <- inf_rec_rate*i + asymp_rec_rate*e - gen_mort_rate*r + r_diff$dC
    dD <- inf_mort_rate*i
    
    return(list(c(dS, dE, dI, dR, dD)))
  })
  
}

#set parameter vector
parameters <- c(birth_rate = birth_rate, symp_rate = symp_rate, asymp_rec_rate = asymp_rec_rate, inf_rec_rate = inf_rec_rate, 
                asymp_cont_rate = asymp_cont_rate, sympt_cont_rate = sympt_cont_rate, inf_mort_rate = inf_mort_rate, gen_mort_rate = gen_mort_rate, A = A,
                s_diffusivity = s_diffusivity, e_diffusivity = e_diffusivity, i_diffusivity = i_diffusivity, r_diffusivity = r_diffusivity,
                left_flux = left_flux, right_flux = right_flux, bottom_flux = bottom_flux, top_flux = top_flux)

#set state vector
state <- c(s = susceptible, e = exposed, i = infected, r = recovered, d = dead)

#set times for 2-D ode solver
times <- 0:120

out <- ode.2D(y = state, func = seird, times = times, parms = parameters, nspec = 5, dimens = c(x.N, y.N), lrw = 10000000)

#get model data at specific times
data_at_0 <- matrix(nrow = x.N, ncol = y.N*5, byrow = FALSE, subset(out, time == 0))
data_at_5 <- matrix(nrow = x.N, ncol = y.N*5, byrow = FALSE, subset(out, time == 5))
data_at_30 <- matrix(nrow = x.N, ncol = y.N*5, byrow = FALSE, subset(out, time == 30))
data_at_60 <- matrix(nrow = x.N, ncol = y.N*5, byrow = FALSE, subset(out, time == 60))
data_at_120 <- matrix(nrow = x.N, ncol = y.N*5, byrow = FALSE, subset(out, time == 120))

#get infected data
infected_data_at_0 <- data_at_0[,(2*y.N+1):(3*y.N)]
infected_data_at_5 <- data_at_5[,(2*y.N+1):(3*y.N)]
infected_data_at_30 <- data_at_30[,(2*y.N+1):(3*y.N)]
infected_data_at_60 <- data_at_60[,(2*y.N+1):(3*y.N)]
infected_data_at_120 <- data_at_120[,(2*y.N+1):(3*y.N)]

#get death data
death_data_at_0 <- data_at_0[,(4*y.N+1):(5*y.N)]
death_data_at_5 <- data_at_5[,(4*y.N+1):(5*y.N)]
death_data_at_30 <- data_at_30[,(4*y.N+1):(5*y.N)]
death_data_at_60 <- data_at_60[,(4*y.N+1):(5*y.N)]
death_data_at_120 <- data_at_120[,(4*y.N+1):(5*y.N)]

#get susceptible data
susc_data_at_0 <- data_at_0[,1:y.N]
susc_data_at_5 <- data_at_5[,1:y.N]
susc_data_at_30 <- data_at_30[,1:y.N]
susc_data_at_60 <- data_at_60[,1:y.N]
susc_data_at_120 <- data_at_120[,1:y.N]

#get exposed data
exp_data_at_0 <- data_at_0[,(y.N+1):(2*y.N)]
exp_data_at_5 <- data_at_5[,(y.N+1):(2*y.N)]
exp_data_at_30 <- data_at_30[,(y.N+1):(2*y.N)]
exp_data_at_60 <- data_at_60[,(y.N+1):(2*y.N)]
exp_data_at_120 <- data_at_120[,(y.N+1):(2*y.N)]

#get recovered data
rec_data_at_0 <- data_at_0[,(3*y.N+1):(4*y.N)]
rec_data_at_5 <- data_at_5[,(3*y.N+1):(4*y.N)]
rec_data_at_30 <- data_at_30[,(3*y.N+1):(4*y.N)]
rec_data_at_60 <- data_at_60[,(3*y.N+1):(4*y.N)]
rec_data_at_120 <- data_at_120[,(3*y.N+1):(4*y.N)]

#active cases = exposed + infected
active_cases_at_0 <- infected_data_at_0 + exp_data_at_0
active_cases_at_5 <- infected_data_at_5 + exp_data_at_5
active_cases_at_30 <- infected_data_at_30 + exp_data_at_30
active_cases_at_60 <- infected_data_at_60 + exp_data_at_60
active_cases_at_120 <- infected_data_at_120 + exp_data_at_120

#total covid cases = exposed + infected + recovered + dead
total_cases_at_0 <- infected_data_at_0 + rec_data_at_0 + death_data_at_0 + exp_data_at_0
total_cases_at_5 <- infected_data_at_5 + rec_data_at_5 + death_data_at_5 + exp_data_at_5
total_cases_at_30 <- infected_data_at_30 + rec_data_at_30 + death_data_at_30 + exp_data_at_30
total_cases_at_60 <- infected_data_at_60 + rec_data_at_60 + death_data_at_60 + exp_data_at_60
total_cases_at_120 <- infected_data_at_120 + rec_data_at_120 + death_data_at_120 + exp_data_at_120

#plotting total covid cases
windows() #this command opens a new window for the graphs (only on Windows machines). you could also skip this and just have the plot to show up on the plot pane
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((total_cases_at_0)/(dx*dy)), plot.title = title(main = "Initial COVID Cases")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((total_cases_at_5)/(dx*dy)), plot.title = title(main = "Total COVID Cases after 5 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((total_cases_at_30)/(dx*dy)), plot.title = title(main = "Total COVID Cases after 30 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((total_cases_at_60)/(dx*dy)), plot.title = title(main = "Total COVID Cases after 60 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((total_cases_at_120)/(dx*dy)), plot.title = title(main = "Total COVID Cases after 120 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))


#plotting total covid deaths
windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((death_data_at_5)/(dx*dy)), plot.title = title(main = "COVID Deaths after 5 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((death_data_at_30)/(dx*dy)), plot.title = title(main = "COVID Deaths after 30 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((death_data_at_60)/(dx*dy)), plot.title = title(main = "COVID Deaths after 60 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((death_data_at_120)/(dx*dy)), plot.title = title(main = "COVID Deaths after 120 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

#plotting active covid cases
windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((active_cases_at_5)/(dx*dy)), plot.title = title(main = "Active COVID Cases after 5 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((active_cases_at_30)/(dx*dy)), plot.title = title(main = "Active COVID Cases after 30 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((active_cases_at_60)/(dx*dy)), plot.title = title(main = "Active COVID Cases after 60 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((active_cases_at_120)/(dx*dy)), plot.title = title(main = "Active COVID Cases after 120 Days")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))


#steady state solution

steady_state <- steady.2D(y = state, func = seird, time = 0, parms = parameters, nspec = 5, dimens = c(x.N, y.N), lrw = 10000000, method = "runsteady")

steady_matrix <- matrix(nrow = x.N, ncol = y.N*5, byrow = FALSE, data = unlist(steady_state))

steady_infected <- steady_matrix[,(2*y.N+1):(3*y.N)]

steady_death <- steady_matrix[,(4*y.N+1):(5*y.N)]

steady_susc <- steady_matrix[,1:y.N]

steady_exp <- steady_matrix[,(y.N+1):(2*y.N)]

steady_rec <- steady_matrix[,(3*y.N+1):(4*y.N)]

steady_total <- steady_infected + steady_rec + steady_death + steady_exp

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((steady_total)/(dx*dy)), plot.title = title(main = "COVID Cases")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))

windows()
filled.contour(x = seq(dx, x.L, by = dx), y = seq(dy, y.L, by = dy), z = log((steady_death)/(dx*dy)), plot.title = title(main = "COVID Deaths")
               , key.title = title(main = "log\n(#/km^2)"), zlim = c(-5, 5))
