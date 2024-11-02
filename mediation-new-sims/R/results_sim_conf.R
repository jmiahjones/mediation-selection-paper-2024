library(tidyverse)
library(qs)
conf_results <- qread("results/conf_results.qs")

conf_results %>% unnest(data) %>%
  select(-ests,-kappa_summary) %>%
  unnest(sel_result) %>%
  mutate(kappa_str = case_when(is.na(kappa) & cv_type == "mse" ~ "CV MSE",
                               is.na(kappa) & cv_type == "bic" ~ "CV BIC",
                               TRUE ~ sprintf("\u03ba=%i", kappa)),
         kappa = factor(kappa_str),
         lam_type = factor(lam_type, levels=c("fixed","local"),
                           labels=c("Fixed-asymptotic \u03bb",
                                    "Local-asymptotic \u03bb")),
         sel_type = factor(sel_idx, levels=c(1:7),
                           labels=c(
                             "Mediator 1",
                             "Mediator 2",
                             rep("Outcome Predictor",3),
                             "Exposure Effect",
                             "Noise"
                           )),
         n = factor(n, labels=sprintf("n=%i", c(200,400)))) %>%

  ggplot(aes(x=sel_type, y=Shat)) +
  geom_col(aes(fill=lam_type), position=position_dodge()) +
  theme_bw() +
  scale_y_continuous(limits=c(0,1), labels = ~sprintf("%i%%", .x*100)) +
  facet_grid(cols=vars(kappa), rows=vars(n)) +
  labs(y="Selection Probability",
       fill="Lambda Rate",
       x="Variable Type") +
  theme(axis.text.x = element_text(hjust = 1, angle=15),
        legend.position = "bottom")

library(cowplot)


CV_kappa_stat <- conf_results %>% unnest(data) %>% select(-sel_result, -ests) %>%
  unnest(kappa_summary) %>% filter(is.na(kappa)); CV_kappa_stat

bias_table <- conf_results %>% unnest(data) %>%
  select(-sel_result, -kappa_summary) %>% unnest(ests) %>%
  mutate(NIE_star = 2, bias=abs(NIE_hat_mean - NIE_star)) %>%
  select(n:cv_type, NIE_star, NIE_hat_mean, bias)

bias_table %>%
  select(n, lam_type, cv_type, bias)

p1 <-
conf_results %>% unnest(data) %>%
  select(-ests,-kappa_summary) %>%
  unnest(sel_result) %>%
  mutate(kappa_str = case_when(is.na(kappa) & cv_type == "mse" ~ "CV MSE",
                               is.na(kappa) & cv_type == "bic" ~ "CV BIC",
                               TRUE ~ sprintf("\u03ba=%i", kappa)),
         kappa = factor(kappa_str),
         lam_type = factor(lam_type, levels=c("fixed","local"),
                           labels=c("Fixed-asymptotic \u03bb",
                                    "Local-asymptotic \u03bb")),
         sel_type = factor(sel_idx, levels=c(1:7),
                           labels=c(
                             "Mediator 1",
                             "Mediator 2",
                             rep("Outcome Predictor",3),
                             "Exposure Effect",
                             "Noise"
                           )),
         n = factor(n, labels=sprintf("n=%i", c(200,400)))) %>%
  filter(n=="n=200") %>%
  ggplot(aes(y=Shat)) +
  # geom_point(aes(color=lam_type)) +
  geom_col(aes(x=sel_type, fill=lam_type), position=position_dodge()) +
  theme_bw() +
  scale_y_continuous(limits=c(0,1), labels = ~sprintf("%i%%", .x*100)) +
  facet_grid(cols=vars(kappa)) +
  labs(y="Selection Probability",
       fill="Lambda Rate",
       x="Variable Type") +
  theme(axis.text.x = element_text(hjust = 1, angle=15),
        legend.position = "none")

p2 <- bias_table %>%
  mutate(kappa_str = case_when(is.na(kappa) & cv_type == "mse" ~ "CV MSE",
                               is.na(kappa) & cv_type == "bic" ~ "CV BIC",
                               TRUE ~ sprintf("\u03ba=%i", kappa)),
         kappa = factor(kappa_str),
         lam_type = factor(lam_type, levels=c("fixed","local"),
                           labels=c("Fixed-asymptotic \u03bb",
                                    "Local-asymptotic \u03bb")),
         n = factor(n, labels=sprintf("n=%i", c(200,400)))) %>%
  filter(n=="n=200") %>%
  ggplot(aes(x=lam_type, y=bias)) +
  geom_col(aes(fill=lam_type), position=position_dodge()) +
  theme_bw() +
  facet_grid(cols=vars(kappa)) +
  labs(y="Bias",
       fill="Lambda Rate") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")

p2
plot_grid(p1,p2, nrow=2, labels=LETTERS[1:2])
ggsave("assumption-violation.png", width=1750, height=750, units="px",
       dpi = 700)
