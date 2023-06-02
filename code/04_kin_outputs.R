
##----------------------- PARENTAL LOSS FROM DRUGS/FIREARMS --------------------
##
##  Code generating the outputs from the kin dynamics of parents.
##
## 
##  Author: Benjamin Schlüter
##  Date: May 2023
##------------------------------------------------------------------------------
##
##  Notes
## ------
##  - Equations (Eq.) in comments refer to the ones found in PAA 2023 paper of 
##    Caswell, Margolis & Verdery. They are there to help readers
##    grasp the code.
## 
##------------------------------------------------------------------------------

rm(list = ls())



## Load packages ---------------------------------------------------------------

## Install/load packages
packages <- c("tidyverse", "ggplot2", "here", "MASS", "viridis")
for(p in packages){
        if(!require(p,character.only = TRUE)) install.packages(p)
        library(p,character.only = TRUE)
}



## Load data -------------------------------------------------------------------

lt.US <- readRDS(here("data_private", "lt_US.rda"))
d_x_t <- readRDS(here("data_private", "d_x_t.rda"))



## Key dimensions for code ---------------------------------------------------------

## Consider three different "causes" of death
alpha <- 3

## Ages 
ages <- unique(lt.US$age)
omega <- length(ages)

## Years 
years <- unique(lt.US$year)
n.years <- length(years)

## Sex
sex <- unique(lt.US$sex)
n.sex <- length(sex)

## Races considered: only look at most populous 
## because no smoothing involved in px &
## zero death counts (zero mx) create issue in computing
## the kin dynamics
races <- c("black", "white", "hispanic", "total")
n.races <- length(races)

## rows alive/dead in d_x_t
rows.stage <- list(
        "alive" = 1:(2*omega),
        "dead" = ((2*omega)+1):nrow(d_x_t)
        )

## rows female/male dead in d_x_t
rows.dead.sex <- list(
        "female" = ((2*omega)+1):((((2*omega)+1) + (alpha*omega))-1),
        "male" = ((((2*omega)+1) + (alpha*omega))):nrow(d_x_t)
        )



## Nber of deaths by cause and age ---------------------------------------------

## Focus on three age class : 1 yo, <5 yo, and <18 yo
age.class <- c("1 yo", "<5 yo", "<18 yo")
n.age.class <- length(age.class)
upper.col <- c(2, 5, 18)


## Container
nber.dth <- array(NA, 
                  dim = c(alpha, 
                          n.sex,
                          n.age.class,
                          n.years,
                          n.races),
                  dimnames = list("cause" = c("other", "drugs", "firearms"),
                                  "sex" = c("mother", "father"),
                                  "focal age" = age.class,
                                  "year" = years,
                                  "race" = races))

## Sum over ages of parents death by cause (Eq. 16)
for (r in races) {
        for (y in 1:n.years) {
                for (s in 1:n.sex) {
                        for (x in 1:n.age.class) {
                                
                                ## Perform summation of age distribution of parent's death  
                                ## over age classes of Focal (1 yo, <5 yo, <18 yo), by sex
                                d_a_t <- apply(d_x_t[rows.dead.sex[[s]], 1:upper.col[x], ,], c(1,3,4), sum)
                                
                                ## Vector allowing to sum over all parent's ages
                                ## by cause
                                sum.vec <- (t(rep(1, omega)) %x% diag(alpha))
                                
                                ## Sum over all parent's age by cause
                                nber.dth[, s, x, y, r] <- sum.vec %*% d_a_t[, y, r]
                        }
                }
        }
}



## Proportion deaths by cause --------------------------------------------------

## Container
prop.dth <- array(NA, 
                  dim = c(alpha, 
                          n.age.class,
                          n.years,
                          n.races),
                  dimnames = list("cause" = c("other", "drugs", "firearms"),
                                  "focal age" = age.class,
                                  "year" = years,
                                  "race" = races))

## Sum mother and father death to compute proportions
nber.dth.parents <- apply(nber.dth, c(1,3,4,5), sum)

## Eq.20
for (r in races) {
        for (y in 1:n.years) {
                
                ## Compute Moore-Penrose pseudo-inverse
                sumDth <- t(rep(1, alpha)) %*% nber.dth.parents[, , y, r] 
                D <- diag(c(sumDth))
                MPpi <- ginv(D)
                ## Get prop by cause
                prop.dth[, , y, r] <- nber.dth.parents[, , y, r] %*% MPpi
        }
}

## Proportional death of parents by drugs & firearms over the years
df.prop.dth <- as.data.frame.table(prop.dth) %>%
        rename("prop" = Freq) %>%
        mutate(year = as.character(year) %>% as.numeric)

## Save outputs
saveRDS(df.prop.dth,
        here("data_private", "df_prop_dth.rda"))



## Prob. of bereavement --------------------------------------------------------

## Container both sex
prob.ber.parents <- array(NA, 
                  dim = c(alpha, 
                          n.age.class,
                          n.years,
                          n.races),
                  dimnames = list("cause" = c("other", "drugs", "firearms"),
                                  "focal age" = age.class,
                                  "year" = years,
                                  "race" = races))
## Eq. 22
for (r in races) {
        for (y in 1:n.years) {
                
                ## Get prop by cause
                prob.ber.parents[, , y, r] <- 1 - (1 - (nber.dth.parents[, , y, r]/2))^2
        }
}

## Probability of parental loss by drugs or firearms over the year
df.prob.ber.parents <- as.data.frame.table(prob.ber.parents) %>%
        rename("prob" = Freq) %>%
        mutate(year = as.character(year) %>% as.numeric)

## Save outputs
saveRDS(df.prob.ber.parents,
        here("data_private", "df_prob_ber_parents.rda"))


## Check: Probability by cause approximately additive?

# ## Sum deaths over drugs and firearms
# nber.dth.parents.sum <- apply(nber.dth.parents[2:3,,,], c(2,3,4), sum)
# 
# ## Container both sex over drugs and firearms combined
# prob.ber.parents <- array(NA, 
#                       dim = c(n.age.class,
#                               n.years,
#                               n.races),
#                       dimnames = list("focal age" = age.class,
#                                       "year" = years,
#                                       "race" = races))
# ## Eq. 22
# for (r in races) {
#         for (y in 1:n.years) {
#                 
#                 ## Get prop with drugs and firearms death summed
#                 prob.ber.parents[, y, r] <- 1 - (1 - (nber.dth.parents.sum[ , y, r]/2))^2
#         }
# }
# ## Probability of parental loss summing over drugs and firearms over the year
# df.prob.ber.parents.2 <- as.data.frame.table(prob.ber.parents) %>%
#         rename("prob" = Freq) %>%
#         mutate(year = as.character(year) %>% as.numeric)
# 
# 
# ## Plot diagonal ?
# df.prob.ber.parents.1 %>% 
#         left_join(df.prob.ber.parents.2,
#                   by = c("focal.age", "year", "race")) %>% 
#         mutate(same = (prob.x == prob.y)) %>% 
#         ggplot(aes(x = prob.x, y = prob.y)) +
#         facet_wrap( ~ race,
#                    scales = "free") +
#         geom_point() +
#         geom_abline(slope = 1) +
#         theme_bw() +
#         labs(x = "Addition of cause-specific prob. to get prob firearms+drugs",
#              y = "Addition of cause-specific deaths to get prob firearms+drugs")


## Container by sex
prob.ber <- array(NA, 
                      dim = c(alpha, 
                              n.sex,
                              n.age.class,
                              n.years,
                              n.races),
                      dimnames = list("cause" = c("other", "drugs", "firearms"),
                                      "sex" = c("mother", "father"),
                                      "focal age" = age.class,
                                      "year" = years,
                                      "race" = races))

for (r in races) {
        for (y in 1:n.years) {
                for (s in 1:n.sex) {
                        
                        ## Get prob by cause and sex
                        prob.ber[, s, , y, r] <- 1 - (1 - (nber.dth[, s, , y, r]))
                }
        }
}

## Probability of parental loss by drugs or firearms over the year
df.prob.ber <- as.data.frame.table(prob.ber) %>%
        rename("prob" = Freq) %>%
        mutate(year = as.character(year) %>% as.numeric)

## Save outputs
saveRDS(df.prob.ber,
        here("data_private", "df_prob_ber.rda"))


## Check: sum of cause-specific prob over all focal's ages ?= 1
## -> Difference always < 0.06
## -> More exact to sum death then compute prob of summed causes
##    In that later case, difference always < 0.025

# ## Container
# nber.dth.check <- array(NA, 
#                         dim = c(alpha, 
#                                 omega,
#                                 n.sex,
#                                 n.years,
#                                 n.races),
#                         dimnames = list("cause" = c("other", "drugs", "firearms"),
#                                         "focal age" = ages,
#                                         "sex" = c("mother", "father"),
#                                         "year" = years,
#                                         "race" = races))
# 
# ## Sum over ages of parents death by cause and sex
# for (r in races) {
#         for (y in 1:n.years) {
#                 for (s in 1:2) {
#                         for (x in 1:omega) {
#                                 
#                                 nber.dth.check[, x, s, y, r] <- (t(rep(1, omega)) %x% diag(alpha)) %*% d_x_t[rows.dead.sex[[s]], x, y, r]
#                         }
#                 }
#         }
# }
# 
# ## Sum mother and father, and all causes death to compute proportions
# nber.dth.parents.check <- apply(nber.dth.check, c(1,4,5), sum)
# 
# ## Exactly 2 only for the year 2000 as it is the only real cohort
# ## due to the assumption of stable population in that year.
# ## High mortality due to Covid19 in the year 2020 leads to > 2 !
# apply(nber.dth.check, c(4,5), sum)
# 
# ## Container both sex
# prob.ber.check <- array(NA, 
#                         dim = c(alpha,
#                                 n.years,
#                                 n.races),
#                         dimnames = list("cause" = c("other", "drugs", "firearms"),
#                                         "year" = years,
#                                         "race" = races))
# 
# ## Eq. 22
# for (r in races) {
#         for (y in 1:n.years) {
#                 
#                 ## Get prop by cause
#                 prob.ber.check[, y, r] <- 1 - (1 - (nber.dth.parents.check[, y, r]/2))^2
#         }
# }
# 
# ## Probability of parental loss by drugs or firearms over the year
# as.data.frame.table(prob.ber.check) %>%
#         rename("prob" = Freq) %>%
#         mutate(year = as.character(year) %>% as.numeric) %>% 
#         group_by(year, race) %>% 
#         summarise(prob = sum(prob)) %>% 
#         ggplot(aes(x = year, y = prob)) +
#         geom_point() +
#         theme_bw()


