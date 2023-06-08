# Packages ----------------------------------------------------------------
library(DiceKriging)
library(SLHD)

library(ggpubr)
library(latex2exp)
library(metR)
library(tidyverse)

source("Code/Simulation.R")

# Models ------------------------------------------------------------------
Covariate_model <- function(n) {
  matrix(runif(n*2), ncol = 2)
}

Propensity_model <- function(x) {
  -2 + 2*x[1]*x[2]
}

Outcome_model <- function(x, a, sd = 0) {
  3/4*exp(-1/4*(9*x[1] - 2)^2 - 1/4*(9*x[2] - 2)^2
            ) + 3/4*exp(-1/49*(9*x[1] + 1)^2 - 1/10*(9*x[2] + 1)^2
                          ) + 1/2*a*exp(-1/4*(9*x[1] - 7)^2 - 1/4*(9*x[2] - 3)^2
                                      ) - 1/5*a*exp(-(9*x[1] - 4)^2 - (9*x[2] - 7)^2) + rnorm(1, 0, sd)
}



# Visualization -----------------------------------------------------------
x <- seq(0, 1, 0.01)
X <- cbind(rep(x, length(x)), rep(x, each = length(x)))

tibble(
  x1 = rep(x, length(x)),
  x2 = rep(x, each = length(x)),
  Treatment = mapply(Outcome_model, asplit(X, 1), 1, 0),
  Control = mapply(Outcome_model, asplit(X, 1), 0, 0)
) %>%
  pivot_longer(cols = Treatment:Control, names_to = "Group", values_to = "Y") %>%
  ggplot() +
  geom_contour_fill(aes(x1, x2, z = Y), bins = 25) +
  labs(x = TeX("$x_1$"), y = TeX("$x_2$"), title = "") +
  scale_fill_gradient2(midpoint = 0.4, low = "steelblue", mid = "white", high = "darkred") +
  facet_grid(. ~ Group) +
  theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    legend.title = element_blank()
  )

ggsave("Visualizations/Franke.png", width = 8, height = 4)

tibble(
  x1 = rep(x, length(x)),
  x2 = rep(x, each = length(x)),
  Treatment = mapply(Outcome_model, asplit(X, 1), 1, 0),
  Control = mapply(Outcome_model, asplit(X, 1), 0, 0),
  TE = Treatment - Control
) %>%
  ggplot() +
  geom_contour_fill(aes(x1, x2, z = TE), bins = 25) +
  labs(x = "", y = "", title = "Treatment effect") +
  scale_fill_gradient2(midpoint = 0.4, low = "steelblue", mid = "white", high = "darkred") +
  theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    legend.title = element_blank()
  )


tibble(
  x1 = rep(x, length(x)),
  x2 = rep(x, each = length(x)),
  Propensity = inv_logit(mapply(Propensity_model, asplit(X, 1)))
) %>%
  ggplot() +
  geom_contour_fill(aes(x1, x2, z = Propensity), bins = 15) +
  labs(x = TeX("$x_1$"), y = TeX("$x_2$"), title = "") +
  scale_fill_gradient2(midpoint = 0.3, low = "white", mid = "firebrick1", high = "firebrick4") +
  theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    legend.title = element_blank()
  )

ggsave("Visualizations/Franke_propensity.png", width = 5, height = 4)

# True values -------------------------------------------------------------

set.seed(0)
n_test <- 100
X_test <- maximinSLHD(1, n_test, 2)$StandDesign
ATE_true <- mean(mapply(Outcome_model, asplit(X_test, 1), 1, 0
                        ) - mapply(Outcome_model, asplit(X_test, 1), 0, 0))

Propensity_true <- inv_logit(mapply(Propensity_model, asplit(X_test, 1)))
ATTE_true <- sum(Propensity_true/sum(Propensity_true)*(
  mapply(Outcome_model, asplit(X_test, 1), 1, 0) - mapply(Outcome_model, asplit(X_test, 1), 0, 0)))
ATO_true <- sum(Propensity_true*(1 - Propensity_true)/sum(Propensity_true*(1 - Propensity_true))*(
  mapply(Outcome_model, asplit(X_test, 1), 1, 0) - mapply(Outcome_model, asplit(X_test, 1), 0, 0)))


# Scenario 1 --------------------------------------------------------------

set.seed(1)
ATE_random <- Simulation_random(B = 50, n_train = 100,
                                Covariate_model = Covariate_model, Outcome_model = Outcome_model,
                                X_test = X_test)
ATE_VR <- Simulation1_VR(B = 50, n_train = 100, n_init = 20,
                          Covariate_model = Covariate_model, Outcome_model = Outcome_model,
                          X_test = X_test)
ATE_ALC <- Simulation1_ALC(B = 50, n_train = 100, n_init = 20,
                           Covariate_model = Covariate_model, Outcome_model = Outcome_model,
                           X_test = X_test)

cbind(ATE_random, ATE_VR, ATE_ALC) %>% 
  apply(2, Eval, true = ATE_true, simplify = T) %>% 
  map(~ select(as.data.frame(.x), bias, RMSE)) %>%
  enframe %>%
  unnest(value, keep_empty = TRUE)

write_csv(tmp, "Outputs/ATE_1.csv")

# Scenario 2A -------------------------------------------------------------

set.seed(1)
ATE_random <- Simulation_random(B = 50, n_train = 100,
                                Covariate_model = Covariate_model, Outcome_model = Outcome_model,
                                X_test = X_test)
ATE_VR <- Simulation2_VR(B = 50, n_train = 100, n_init = 20, n_cand = 500,
                         Covariate_model = Covariate_model, Outcome_model = Outcome_model,
                         X_test = X_test)
ATE_ALC <- Simulation2_ALC(B = 50, n_train = 100, n_init = 20, n_cand = 500,
                          Covariate_model = Covariate_model, Outcome_model = Outcome_model,
                          X_test = X_test)

tmp <- cbind(ATE_random, ATE_VR, ATE_ALC) %>% 
  apply(2, Eval, true = ATE_true, simplify = T) %>% 
  map(~ select(as.data.frame(.x), bias, RMSE)) %>%
  enframe %>%
  unnest(value, keep_empty = TRUE)

write_csv(tmp, "Outputs/ATE_2A.csv")

# Scenario 2B -------------------------------------------------------------

set.seed(1)
ATE_random <- Simulation_random(B = 50, n_train = 100,
                                Covariate_model = Covariate_model, 
                                Propensity_model = Propensity_model,
                                Outcome_model = Outcome_model,
                                X_test = X_test)
ATE_VR <- Simulation2_VR(B = 50, n_train = 100, n_init = 20, n_cand = 500,
                         Covariate_model = Covariate_model,
                         Propensity_model = Propensity_model,
                         Outcome_model = Outcome_model,
                         X_test = X_test)
ATE_ALC <- Simulation2_ALC(B = 50, n_train = 100, n_init = 20, n_cand = 500,
                           Covariate_model = Covariate_model,
                           Propensity_model = Propensity_model,
                           Outcome_model = Outcome_model,
                           X_test = X_test)

tmp <- cbind(ATE_random, ATE_VR, ATE_ALC) %>% 
  apply(2, Eval, true = ATE_true, simplify = T) %>% 
  map(~ select(as.data.frame(.x), bias, RMSE)) %>%
  enframe %>%
  unnest(value, keep_empty = TRUE)

write_csv(tmp, "Outputs/ATE_2B.csv")

set.seed(1)
ATTE_random <- Simulation_random(B = 50, n_train = 200,
                                 Covariate_model = Covariate_model, 
                                 Propensity_model = Propensity_model,
                                 Outcome_model = Outcome_model,
                                 X_test = X_test, target = "ATTE")
ATTE_VR <- Simulation2_VR(B = 50, n_train = 200, n_init = 50, n_cand = 1000,
                          Covariate_model = Covariate_model,
                          Propensity_model = Propensity_model,
                          Outcome_model = Outcome_model,
                          X_test = X_test, target = "ATTE")
ATTE_ALC <- Simulation2_ALC(B = 50, n_train = 200, n_init = 50, n_cand = 1000,
                            Covariate_model = Covariate_model,
                            Propensity_model = Propensity_model,
                            Outcome_model = Outcome_model,
                            X_test = X_test, target = "ATTE")

tmp <- cbind(ATTE_random, ATTE_VR, ATTE_ALC) %>% 
  apply(2, Eval, true = ATTE_true, simplify = T) %>% 
  map(~ select(as.data.frame(.x), bias, RMSE)) %>%
  enframe %>%
  unnest(value, keep_empty = TRUE)

write_csv(tmp, "Outputs/ATTE_2B.csv")

set.seed(2)
ATO_random <- Simulation_random(B = 50, n_train = 200,
                                Covariate_model = Covariate_model, 
                                Propensity_model = Propensity_model,
                                Outcome_model = Outcome_model,
                                X_test = X_test, target = "ATO")
ATO_VR <- Simulation2_VR(B = 50, n_train = 200, n_init = 50, n_cand = 1000,
                         Covariate_model = Covariate_model,
                         Propensity_model = Propensity_model,
                         Outcome_model = Outcome_model,
                         X_test = X_test, target = "ATO")
ATO_ALC <- Simulation2_ALC(B = 50, n_train = 200, n_init = 50, n_cand = 1000,
                           Covariate_model = Covariate_model,
                           Propensity_model = Propensity_model,
                           Outcome_model = Outcome_model,
                           X_test = X_test, target = "ATO")

tmp <- cbind(ATO_random, ATO_VR, ATO_ALC) %>% 
  apply(2, Eval, true = ATO_true, simplify = T) %>% 
  map(~ select(as.data.frame(.x), bias, RMSE)) %>%
  enframe %>%
  unnest(value, keep_empty = TRUE)

write_csv(tmp, "Outputs/ATO_2B.csv")

# Scenario 3 --------------------------------------------------------------

set.seed(3)
Effect_random <- Simulation3_random(B = 50, n_train = 50,
                                    Covariate_model = Covariate_model, 
                                    Propensity_model = Propensity_model,
                                    Outcome_model = Outcome_model)
Effect_greedy <- Simulation3_greedy(B = 50, n_train = 50, n_init = 10, n_cand = 1000,
                                    Covariate_model = Covariate_model,
                                    Propensity_model = Propensity_model,
                                    Outcome_model = Outcome_model)
Effect_UCB <- Simulation3_UCB(B = 50, n_train = 50, n_init = 10, n_cand = 1000,
                              Covariate_model = Covariate_model,
                              Propensity_model = Propensity_model,
                              Outcome_model = Outcome_model, c = 0.0001)

tibble(
  Method = factor(rep(c("Random", "Greedy", "UCB"), each = 50), levels = c("Random", "Greedy", "UCB")),
  Effect = c(Effect_random, Effect_greedy, Effect_UCB)
) %>% 
ggplot() +
  geom_boxplot(aes(x = Method, y = Effect, color = Method)) +
  labs(x = "", y = "") +
  scale_color_manual(values = c("black", "steelblue", "firebrick")) +
  theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    legend.position = "None"
  )
  
ggsave("Visualizations/Franke_Sim3.png", width = 5, height = 4)

  