library(tidyverse)
library(qs)
conf_results <- qread("results/too_small_results.qs")

kappa_colors <- c("magenta3", "steelblue1")

p1 <- conf_results %>% unnest(data) %>%
  select(-ests,-kappa_summary) %>%
  unnest(sel_result) %>%
  filter(is.na(cv_type), lam_type=="local") %>%
  mutate(
    kappa_str = paste0("\u03ba=",kappa),
         kappa = factor(kappa_str),
         sel_type = factor(sel_idx, levels=c(1:5),
                           labels=c(
                             "Mediator 1",
                             "Mediator 2",
                             "Outcome Predictor",
                             "Exposure Effect",
                             "Noise"
                           )),
         n = factor(n, labels=sprintf("n=%i", c(2,4)*100))) %>%

  ggplot(aes(x=sel_type, y=Shat)) +
  geom_col(aes(fill=kappa), position=position_dodge()) +
  theme_bw() +
  scale_y_continuous(limits=c(0,1), labels = ~sprintf("%i%%", .x*100)) +
  scale_fill_manual(values=kappa_colors) +
  facet_grid(cols=vars(n)) +
  labs(y="Selection Probability",
       fill="\u03ba",
       x="Variable Type") +
  theme(axis.text.x = element_text(hjust = 1, angle=15),
        legend.position = "bottom",
        legend.title = element_blank())

CV_kappa_stat <- conf_results %>% unnest(data) %>% select(-sel_result, -ests) %>%
  unnest(kappa_summary) %>% filter(is.na(kappa)); CV_kappa_stat

bias_table <- conf_results %>% unnest(data) %>%
  select(-sel_result, -kappa_summary) %>% unnest(ests) %>%
  mutate(NIE_star = 1 + n**(-3/2),
         bias=NIE_hat_mean - NIE_star) %>%
  select(n:cv_type, NIE_star, NIE_hat_mean, bias)

p1

p2 <- bias_table %>%
  filter(is.na(cv_type), lam_type=="local") %>%
  mutate(rtnbias=sqrt(n)*abs(bias),
         kappa_str = paste0("\u03ba=",kappa),
         kappa=factor(kappa_str)) %>%
  ggplot(aes(x=n, y=rtnbias)) +
  geom_line(aes(color=kappa, group=kappa)) +
  theme_bw() +
  scale_color_manual(values=kappa_colors) +
  labs(y="\u221an x bias", color="\u03ba") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

p1
p2




p1 <-
  conf_results %>% unnest(data) %>%
  select(-ests,-kappa_summary) %>%
  unnest(sel_result) %>%
  mutate(kappa_str = paste0("\u03ba=",kappa),
         kappa = factor(kappa_str),
         lam_type = factor(lam_type, levels=c("fixed","local"),
                           labels=c("Fixed-asymptotic \u03bb",
                                    "Local-asymptotic \u03bb")),
         sel_type = factor(sel_idx, levels=c(1:5),
                           labels=c(
                             "Mediator 1",
                             "Mediator 2",
                             "Outcome Predictor",
                             "Exposure Effect",
                             "Noise"
                           )),
         n = factor(n, labels=sprintf("n=%i", c(200,400)))) %>%
  ggplot(aes(y=Shat)) +
  geom_col(aes(x=sel_type, fill=kappa), position=position_dodge()) +
  theme_bw() +
  scale_fill_manual(values=kappa_colors) +
  scale_y_continuous(limits=c(0,1), labels = ~sprintf("%i%%", .x*100)) +
  facet_grid(cols=vars(n)) +
  labs(y="Selection Probability",
       fill="Lambda Rate",
       x="Variable Type") +
  theme(axis.text.x = element_text(hjust = 1, angle=15),
        legend.position = "none")

p2 <- bias_table %>%
  mutate(kappa_str = paste0("\u03ba=",kappa),
         kappa = factor(kappa_str),
         lam_type = factor(lam_type, levels=c("fixed","local"),
                           labels=c("Fixed-asymptotic \u03bb",
                                    "Local-asymptotic \u03bb")),
         n = factor(n, labels=sprintf("n=%i", c(200,400)))) %>%
  ggplot(aes(x=lam_type, y=bias)) +
  geom_col(aes(fill=kappa), position=position_dodge()) +
  scale_fill_manual(values=kappa_colors) +
  theme_bw() +
  facet_grid(cols=vars(n)) +
  labs(y="Bias",
       fill="Lambda Rate") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")
plot_grid(p1,p2,nrow=2)
