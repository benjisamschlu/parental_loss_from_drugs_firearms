
##----------------------- PARENTAL LOSS FROM DRUGS/FIREARMS --------------------
##
##  Code performing the explanatory data analysis.
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
packages <- packages <- c("tidyverse", "ggplot2", "here", "openxlsx", "viridis")
for(p in packages){
        if(!require(p,character.only = TRUE)) install.packages(p)
        library(p,character.only = TRUE)
}




## Functions -------------------------------------------------------------------

source(here("code", "fct_matrix_kinship_model.R"))



## Load data ------------------------------------------------------------------------

lt.US <- readRDS(here("data_private", "lt_US.rda"))
fx.US <- readRDS(here("data", "fx_US.rda"))
mac.US <- readRDS(here("data", "mac_US.rda"))


## Mortality data from Matt
df.dth <- readRDS(here("data_private",  
                       "state_year_age-sex-race_drug-opioid-mortality.RDS")) %>% 
        rename("age" = age_years) %>% 
        dplyr::select(!c(st_fips, division))



## Key dimensions for code ---------------------------------------------------------

## Ages 
ages <- unique(lt.US$age)
omega <- length(ages)
## Years (do not account for Covid19 year)
years <- unique(lt.US$year)
n.years <- length(years)
## Races considered
races <- unique(lt.US$race_eth)
n.races <- length(races)



## Data checks related to the call with Math -----------------------------------

## Create only one data set is possible?
deathsU5.by.drugs.or.firearms <- lt.US %>% 
        filter(age<5) %>% 
        group_by(year, race_eth) %>% 
        summarise_at(vars(n_deaths:n_firearm), sum)

deaths80plus.by.drugs.or.firearms <- lt.US %>% 
        filter(age>=80) %>% 
        group_by(year, race_eth) %>% 
        summarise_at(vars(n_deaths:n_firearm), sum)

deathsTOT.by.drugs.or.firearms <- lt.US %>% 
        group_by(year, race_eth) %>% 
        summarise_at(vars(n_drug:n_firearm), sum)

## Save all into an Excel files
list.datasets <- list("Under 5 yo" = deathsU5.by.drugs.or.firearms, 
                      "80+ yo" = deaths80plus.by.drugs.or.firearms,
                      "All ages" = deathsTOT.by.drugs.or.firearms)
write.xlsx(list.datasets, file = here("rmds", "table_deaths.xlsx"))



## EDA mortality ---------------------------------------------------------------

## MX all-causes 
lt.US %>% 
        ggplot(aes(x = age, y = mx, group = year)) +
        facet_grid(sex ~ race_eth) +
        geom_line(aes(col = year)) +
        theme_bw() +
        scale_y_log10() +
        labs(y = "mx (log-scale)")

## Mx all-causes other than drugs & firearms
lt.US %>% 
        mutate(mx_other = (n_deaths - (n_drug + n_firearm))/pop) %>% 
        ggplot(aes(x = age, y = mx_other, group = year)) +
        facet_grid(sex ~ race_eth) +
        geom_line(aes(col = year)) +
        theme_bw() +
        scale_y_log10() +
        labs(y = "mx (log-scale)")

## Mx drugs
lt.US %>% 
        mutate(mx_drug = n_drug/pop) %>% 
        ggplot(aes(x = age, y = mx_drug, group = year)) +
        facet_grid(sex ~ race_eth) +
        geom_line(aes(col = year)) +
        theme_bw() +
        labs(y = "mx drugs")

## Mx firearms
lt.US %>% 
        mutate(mx_firearm = n_firearm/pop) %>% 
        ggplot(aes(x = age, y = mx_firearm, group = year)) +
        facet_grid(sex ~ race_eth) +
        geom_line(aes(col = year)) +
        theme_bw() +
        labs(y = "mx drugs")

## Life expectancy at birth
lt.US %>% 
        filter(age == 0) %>% 
        ggplot(aes(x = year, y = ex, group = race_eth, col = race_eth)) +
        facet_wrap(~ sex) +
        geom_line(linewidth = 1) +
        geom_line(aes(y = ifelse(race_eth == "total", ex, NA)), linewidth = 2) +
        theme_bw() +
        scale_color_manual(values = c("aian" = "#440154FF",
                                      "api" = "#3B528BFF",
                                      "black" = "#21908CFF",
                                      "hispanic" = "#5DC863FF",
                                      "white" = "#FDE725FF",
                                      "total" = "#000000")) +
        labs(y = expression(e^0),
             x = "Year",
             col = "Race/ethnic group")



## EDA Fertility -------------------------------------------------------------------

## Modelled difference in MAC
## vs observed MAC
ggplot(mac.US) +
        geom_point(aes(x = year, y = diff.mac)) +
        geom_line(aes(x = year, y = fit, group = 1)) +
        theme_minimal()

## Interpolated fertility schedules across races
fx.US %>% 
        ggplot(aes(x = age, group = interaction(sex, year), col = sex)) +
        facet_wrap(~ race_eth) +
        geom_line(aes(y = fx), linewidth = 1) +
        theme_bw() +
        labs(y = "fx")

## Age of the parents of offspring
pi <- array(NA,
            dim = c(omega, 2, n.years, n.races),
            dimnames = list("age" = ages,
                            "sex" = c("female", "male"),
                            "year" = years,
                            "race" = races))

## Compute age distribution of parents of offspring
for (y in 1:n.years) {
        for (r in races) {
                
                y.n <- years[y]
                pi.temp <- get_pi(ages, fx.US, lt.US, y.n, r)
                pi[, 1, y, r] <- pi.temp[1:omega]
                pi[, 2, y, r] <- pi.temp[(omega+1):(2*omega)]
        }
}
## Convert to df
df.pi <- as.data.frame.table(pi) %>% 
        rename("d_x" = Freq) %>% 
        mutate(age = as.character(age) %>% as.numeric,
               year = as.character(year) %>% as.numeric)
## Plot
df.pi %>% 
        filter(year %in% c(2000, 2005, 2010, 2015, 2020),
               age %in% 10:65,
               race != "total") %>% 
        ggplot(aes(x = age, y = d_x, group = race)) +
        facet_grid(sex ~ year) +
        geom_line(aes(col = race), linewidth = 1.2) +
        theme_bw() +
        scale_color_manual(values = c("aian" = "#440154FF",
                                      "api" = "#3B528BFF",
                                      "black" = "#21908CFF",
                                      "hispanic" = "#5DC863FF",
                                      "white" = "#FDE725FF",
                                      "total" = "#000000")) +
        labs(y = "Age of parents at offspring")
