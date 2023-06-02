
##----------------------- PARENTAL LOSS FROM DRUGS/FIREARMS --------------------
##
##  Code projecting the age distribution of dead/alive parents over the
##  age of their child, and years using the matrix kinship model
##  developed by Hal Caswell. The projection is performed 
##  independently for each ethnic group.
##
## 
##  Author: Benjamin Schlüter
##  Date: May 2023
##------------------------------------------------------------------------------
##
##  Notes
## ------
##
## 
##------------------------------------------------------------------------------

rm(list = ls())



## Load packages ---------------------------------------------------------------

## Install/load packages
packages <- c("tidyverse", "here")
for(p in packages){
        if(!require(p,character.only = TRUE)) install.packages(p)
        library(p,character.only = TRUE)
}



## Functions -------------------------------------------------------------------

## Script containing the functions computing
## the main elements required to run the kin dynamics
source(here("code", "fct_matrix_kinship_model.R"))



## Load data -------------------------------------------------------------------

lt.US <- readRDS(here("data_private", "lt_US.rda")) %>% 
        ## Compute mx by causes
        mutate(mx_other = (n_deaths - (n_drug + n_firearm))/pop,
               mx_drug = n_drug/pop,
               mx_firearm = n_firearm/pop)

fx.US <- readRDS(here("data", "fx_US.rda"))



## Key dimensions for code ---------------------------------------------------------

## Consider three different "causes" of death
alpha <- 3
## Ages 
ages <- unique(lt.US$age)
omega <- length(ages)
## Years 
years <- unique(lt.US$year)
n.years <- length(years)
## Races considered: only look at most populous 
## because no smoothing involved in px &
## zero death counts (zero mx) create issue in computing
## the kin dynamics
races <- c("black", "white", "hispanic", "total")
n.races <- length(races)



## Container for dynamics of parents of Focal ----------------------------------

d_x_t <- array(NA,
               dim = c(2*(omega + alpha*omega), 
                       omega, 
                       n.years, 
                       n.races),
               dimnames = list("Parents" = 1:(2*(omega + alpha*omega)),
                               "Focal age" = ages,
                               "Year" = years,
                               "Race" = races))



## Run the dynamics ------------------------------------------------------------

## These lists allow to select part of Parent vector 
## alive/dead while looping
rows <- list("alive" = 1:(2*omega),
             "dead" = ((2*omega)+1):nrow(d_x_t)
)

## Boundary conditions at t=0 (year=2000)
## -> Assume vital rate as if stable population
## because Focal is a randomly drawn individual
## in the population but he/she could be any age.
## d(0,0) consists of the dist. of the ages 
## of the parents of offspring  in the year 2000
## at Focal's birth.

## START LOOP ON RACE
for (r in races) {
        
        ## Boundary conditions
        ## Distribution of parents' age of offspring
        d_x_t[rows[["alive"]], 1, 1, r] <- get_pi(ages, fx.US, lt.US, min(years), r)
        ## Assumption of no dead parents at birth of Focal
        d_x_t[rows[["dead"]], 1, 1, r] <- rep(0, (2*(alpha*omega)))
        
        ## Get d(x,0) for all x assuming stable pop in 
        ## the year 2000 (Focal could be any age in the
        ## year 2000)
        U_tilde <- get_U_tilde(lt.US, min(years), r, cum = F)
        
        for (x in 2:omega) {
                
                d_x_t[, x, 1, r] <- U_tilde %*% d_x_t[, (x-1), 1, r]
        }
        
        ## Create the dynamics by looping over the years
        
        ## START LOOP ON YEARS
        
        for (y in 2:n.years) {
                
                y.n <- years[y]
                
                ## Boundary conditions
                ## Focal birth:
                ## Distribution of parents' age at offspring
                d_x_t[rows[["alive"]], 1, y, r] <- get_pi(ages, fx.US, lt.US, y.n, r)
                ## Assumption of no dead parents at birth of Focal
                d_x_t[rows[["dead"]], 1, y, r] <- rep(0, (2*(alpha*omega)))
                
                ## Dynamics of parents
                U_tilde <- get_U_tilde(lt.US, y.n, r, cum = F)
                
                d_x_t[, 2:omega, y, r] <- U_tilde %*% d_x_t[, 1:(omega-1), (y-1), r]
        }
}



## Save dynamics of parents ----------------------------------------------------

saveRDS(d_x_t,
        here("data_private", "d_x_t.rda"))




