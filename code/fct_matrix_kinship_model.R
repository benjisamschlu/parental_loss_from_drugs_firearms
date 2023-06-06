
##------------- FUNCTIONS FOR MATRIX KINSHIP MODEL -----------------------------
## 
## 
##  Author: Benjamin Schlüter
##  Date: April 2023
##
##
## -----------------------------------------------------------------------------



## Load packages ---------------------------------------------------------------

packages <- c("tidyverse", "Matrix")
invisible( lapply(packages, library, character.only = TRUE))



## Functions -------------------------------------------------------------------

## Construct U matrix from a lifetable
get_U <- function(s, y, r, lifeTable) {
        
        ## Extract survival probs 
        px <- lifeTable %>% 
                filter(sex == s,
                       year == y,
                       race_eth == r) %>% 
                pull(px)
        ## Dim of A
        omega <- length(px)
        ## Creation of U
        U <- matrix(0, 
                    nrow = omega, 
                    ncol = omega)
        ## Survival prob on subdiagonal
        U[row(U)-col(U)==1] <- head(px,-1) # ASSUMPTION: Everybody die at last age
        
        return(U)
}

## Construct \hat{M} matrix
get_M_hat <- function(s, y, r, lifeTable) {
        
        ## Extract death probs 
        qx <- lifeTable %>% 
                filter(sex == s,
                       year == y,
                       race_eth == r) %>% 
                pull(qx)
        
        ## Create hazard matrix, dim(H)=(alpha*omega)
        H <- lifeTable %>% 
                filter(sex == s,
                       year == y,
                       race_eth == r) %>% 
                dplyr::select(mx_other, mx_drug, mx_firearm) %>% 
                as.matrix %>% 
                t()
        
        ## Create summed hazard
        h_tilde <- t(rep(1, nrow(H))) %*% H 
        
        ## Compute M (Caswell et al (2023))
        M <- H %*% solve(diag(c(h_tilde))) %*% diag(qx)
        
        ## Store columns of M as a list of vectors
        M.cols <- lapply(1:ncol(H), function(j) return(M[,j]))
        
        ## Create M_hat using the vectors as elements of the block diagonal
        M_hat <- bdiag(M.cols)
        
        return(M_hat)
}

## Construct U_tilde from U and M
get_U_tilde <- function(lifeTable, y, r, cum = F) {
        
        
        ## U for both sexes
        U <- lapply(c("female", "male"), get_U, y, r, lifeTable)
        
        ## M_hat for both sexes
        M_hat <- lapply(c("female", "male"), get_M_hat, y, r, lifeTable)
        
        ## Dimensions
        omega <- dim(M_hat[[1]])[2]
        alphaomega <- dim(M_hat[[1]])[1]
        
        ## Construct block matrix, block by block
        ## Upper-left block: survival probs.
        block_UL <- bdiag(U)
        ## Lower-left block: death probs. by cause
        block_LL <- bdiag(M_hat)
        ## Upper-right block: death can't become alive 
        zeros <- matrix(0,nrow = omega, ncol = alphaomega)
        block_UR <- bdiag(list(zeros, zeros))
        ## Lower-right block: 
        if (cum) {
                I <- diag(alphaomega)
                block_LR <- bdiag(list(I, I))
        } else {
                zeros <- matrix(0, nrow = alphaomega, ncol = alphaomega)
                block_LR <- bdiag(list(zeros, zeros))
        }
        ## Combine
        block_U <- cbind(block_UL, block_UR)
        block_L <- cbind(block_LL, block_LR)
        
        U_tilde <- rbind(block_U, block_L)
        
        return(as.matrix(U_tilde))
}

## Construct F matrix from fx
get_F <- function(s, y, r, ages, asfr) {
        
        omega <- length(ages)
        ## Extract asfr
        fx <- asfr %>% 
                filter(sex == s,
                       year == y,
                       race_eth == r) %>% 
                pull(fx)
        ## Reproduction ages present in asfr
        ages.repro <- asfr %>% 
                filter(sex == s,
                       year == y,
                       race_eth == r) %>%
                pull(age)
        ## Creation of F
        F. <- matrix(0, 
                    nrow = omega, 
                    ncol = omega)
        ## ASFR on 1st row
        F.[1, (ages.repro + 1)] <- fx 
        
        return(F.)
}

## Dist. of the ages of the parents of offspring 
## NOTE: Using norm() for 1-norm does not lead to 
##      take the sum of the absolute values of a vector.
##      -> We compute it manually
get_pi <- function(ages, asfr, lifeTable, y, r) {
        
        
        PI <- sapply(c("female", "male"), function(s) {
                
                ## Fertility matrix
                F. <- get_F(s, y, r, ages, asfr)
                ## Population vector
                z <- lifeTable %>% 
                        filter(sex == s,
                               year == y,
                               race_eth == r) %>% 
                        arrange(age) %>% 
                        pull(pop)
                ## Population age structure
                z <- z/sum(z)
                ## Compute distribution of ages at offspring
                num <- t(F.[1, ]) * z
                denom <- sum( abs(num) )
                pi <- num / denom
                
                return(pi)
        },
        simplify = T,
        USE.NAMES = T)

        return(c(PI))
}

## Parents death at a population level 
## -> weight the age-specific frequency of parents death
##    by the survivorship column of a cohort life table
##    with a starting population equal to the initial
##    size of each child birth cohort. Look at male
##    and female children together (life table both sexes
##    combined)
get_df_nber_children <- function(y, r, lifeTable, asfr) {
        
        ## Focus on age group <5 to be able to go until year 2017
        ## with a cohort perspective
        y.window <- y:(y+3)
        
        ## Get period mx both sexes combined
        mx.matrix <- lifeTable %>% 
                filter(year %in% y.window,
                       race_eth == r) %>% 
                ## Sum over sexes
                group_by(year, race_eth, age) %>% 
                summarise(n_deaths = sum(n_deaths),
                          pop = sum(pop)) %>% 
                ungroup() %>% 
                ## Mx two sexes
                mutate(mx = n_deaths/pop) %>%
                ## Convert into a matrix
                dplyr::select(year, age, mx) %>% 
                pivot_wider(names_from = year, values_from = mx) %>% 
                dplyr::select(!age) %>% 
                as.matrix()
        
        ## Take period rates as cohort rates using diagonal
        mx.cohort <- diag(mx.matrix)[1:4] 
        
        ## Build cohort life table for age <5 yo
        ax = c(0.07 + 1.7*mx.cohort[1], rep(0.5, 3))
        n = rep(1, 4)
        qx = (n * mx.cohort)/(1 + (n - ax) * mx.cohort)
        px = 1 - qx
        ## Fraction surviving
        lx = cumprod(c(1, px))
        
        ## Get total birth using female fx and pop
        birth <- asfr %>% 
                filter(sex == "female",
                       ## Birth cohort
                       year == y,
                       race_eth == r) %>% 
                left_join(lt.US %>% 
                                  dplyr::select(year, race_eth, sex, age, pop),
                          by = c("year", "race_eth", "sex", "age")) %>% 
                mutate(birth = fx * pop) %>% 
                group_by(year) %>% 
                summarise(birth = sum(birth)) %>% 
                pull(birth)
        
        ## Compute nber of surviving children
        surv.children <- lx * birth
        
        ## Container
        child.losing.parent <- array(NA,
                                     dim = c(alpha, 5, n.sex),
                                     dimnames = list("cause" = c("other", "drugs", "firearms"),
                                                     "focal age" = 0:4,
                                                     "sex" = c("mother", "father")))
        for (s in 1:n.sex) {
                
                ## Vector allowing to sum over all parent's ages
                ## by cause
                sum.vec <- (t(rep(1, omega)) %x% diag(alpha))
                
                ## Sum over all parent's age by cause
                dist.dth.cause <- sum.vec %*% d_x_t[rows.dead.sex[[s]], 1:5, which(y == years), r] 
                
                ## Loop over the ages: 0 -> 4 yo
                child.losing.parent[,,s] <- sapply(1:5, function(j) {
                        
                        surv.children[j] * dist.dth.cause[,j]
                })
                
        }
        out <- as.data.frame.table(child.losing.parent) %>% 
                rename("N" = Freq) %>% 
                mutate(race = r,
                       year = y)
        return(out)
}


