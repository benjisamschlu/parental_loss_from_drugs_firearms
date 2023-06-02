
##----------------------- PARENTAL LOSS FROM DRUGS/FIREARMS --------------------
##
##  Code visualizing the generated outputs of the parents' kin dynamics.
##
## 
##  Author: Benjamin Schlüter
##  Date: June 2023
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
packages <- c("tidyverse", "ggplot2", "here", "viridis")
for(p in packages){
        if(!require(p,character.only = TRUE)) install.packages(p)
        library(p,character.only = TRUE)
}



## Load data -------------------------------------------------------------------

df.prop.dth <- readRDS(here("data_private", "df_prop_dth.rda"))
df.prob.ber.parents <- readRDS(here("data_private", "df_prob_ber_parents.rda"))
df.prob.ber <- readRDS(here("data_private", "df_prob_ber.rda"))



## Figures ---------------------------------------------------------------------

## Proportional parents deaths by cause
df.prop.dth %>% 
        filter(cause != "other") %>% 
        ggplot(aes(x = year, y = prop, group = race)) +
        facet_grid(cause ~ focal.age) +
        geom_line(aes(col = race), linewidth = 1) +
        theme_bw() +
        scale_color_manual(values = c("aian" = "#440154FF",
                                      "api" = "#3B528BFF",
                                      "black" = "#21908CFF",
                                      "hispanic" = "#5DC863FF",
                                      "white" = "#FDE725FF",
                                      "total" = "#000000")) +
        labs(y = "Proportion of parents' death",
             x = "Year",
             col = "Race/ethnic group")

## Prob of bereavement
## Sum cause-specific probability of losing a parent
df.prob.ber.parents.1 <- df.prob.ber.parents %>% 
        filter(cause != "other") %>% 
        group_by(focal.age, year, race) %>% 
        summarise(prob = sum(prob)) %>% 
        ungroup() %>% 
        mutate(cause = "drugs & firearms")

## Cause-specific probability of losing a parent       
bind_rows(df.prob.ber.parents %>% 
                  filter(cause != "other"),
          ## Add drugs & firearms together
          df.prob.ber.parents.1) %>% 
        mutate(cause = factor(cause,
                              levels = c("drugs", "firearms", "drugs & firearms"))) %>% 
        ggplot(aes(x = year, y = prob, group = race)) +
        facet_grid(cause ~ focal.age) +
        geom_line(aes(col = race), linewidth = 1) +
        theme_bw() +
        scale_color_manual(values = c("aian" = "#440154FF",
                                      "api" = "#3B528BFF",
                                      "black" = "#21908CFF",
                                      "hispanic" = "#5DC863FF",
                                      "white" = "#FDE725FF",
                                      "total" = "#000000")) +
        scale_y_log10() +
        labs(y = "Probability of losing one or more parent (log-scale)",
             x = "Year",
             col = "Race/ethnic group")

## Prob of bereavement by parent
## Sum cause-specific probability of losing a parent
df.prob.ber.sum <- df.prob.ber %>% 
        filter(cause != "other") %>% 
        group_by(sex, focal.age, year, race) %>% 
        summarise(prob = sum(prob)) %>% 
        ungroup() %>% 
        mutate(cause = "drugs & firearms")

bind_rows(df.prob.ber %>% 
                  filter(cause != "other",
                         race != "total"),
          ## Add drugs & firearms together
          df.prob.ber.sum) %>% 
        mutate(cause = factor(cause,
                              levels = c("drugs", "firearms", "drugs & firearms"))) %>% 
        ggplot(aes(x = year, y = prob, group = interaction(race, sex))) +
        facet_grid(cause ~ focal.age) +
        geom_line(aes(col = race, linetype = sex), linewidth = 1) +
        theme_bw() +
        scale_color_manual(values = c("aian" = "#440154FF",
                                      "api" = "#3B528BFF",
                                      "black" = "#21908CFF",
                                      "hispanic" = "#5DC863FF",
                                      "white" = "#FDE725FF",
                                      "total" = "#000000")) +
        labs(y = "Probability of losing one or more parent",
             x = "Year",
             col = "Race/ethnic group")

## Same as above but on log-scale
bind_rows(df.prob.ber %>% 
                  filter(cause != "other",
                         race != "total"),
          ## Add drugs & firearms together
          df.prob.ber.sum) %>% 
        mutate(cause = factor(cause,
                              levels = c("drugs", "firearms", "drugs & firearms"))) %>% 
        ggplot(aes(x = year, y = prob, group = interaction(race, sex))) +
        facet_grid(cause ~ focal.age) +
        geom_line(aes(col = race, linetype = sex), linewidth = 1) +
        theme_bw() +
        scale_color_manual(values = c("aian" = "#440154FF",
                                      "api" = "#3B528BFF",
                                      "black" = "#21908CFF",
                                      "hispanic" = "#5DC863FF",
                                      "white" = "#FDE725FF",
                                      "total" = "#000000")) +
        scale_y_log10() +
        labs(y = "Probability of losing one or more parent (log-scale)",
             x = "Year",
             col = "Race/ethnic group")


