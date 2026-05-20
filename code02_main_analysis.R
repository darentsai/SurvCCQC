source("code01_data_preprocess.R")

library(tidyverse)
library(survival)  # Cox model
library(coxme)     # Cox frailty model
library(survminer) # survival curve
library(patchwork)

#############################
#--- Cox (Frailty) Model ---#
#############################

haz_ratio <- function(x, var) {
  coef <- coef(summary(x))
  ci <- confint(x)
  coef <- coef[grep(var, rownames(coef)), , drop = FALSE]
  ci <- ci[grep(var, rownames(ci)), , drop = FALSE]
  data.frame(HR = exp(coef[, "coef"]),
             HR_L = exp(ci[, 1]),
             HR_U = exp(ci[, 2]),
             pval = coef[, colnames(coef) %in% c("p", "Pr(>|z|)")])
}

hr1 <- ca_dat %>%
  group_by(pri_site) %>%
  group_modify(~ {
    form <- Surv(surv_m, is_dead) ~ pass
    mod <- coxph(form, data = .x)
    haz_ratio(mod, "pass")
  }) %>%
  ungroup()

hr2 <- ca_dat %>%
  group_by(pri_site) %>%
  group_modify(~ {
    form <- Surv(surv_m, is_dead) ~ pass + hosplv + sex + agegp + fnstage + region + do_surg + do_chem + do_rt
    if(.y$pri_site %in% c("2_breast", "5_cervix", "6_uterine", "7_ovary")) {
      form <- update(form, . ~ . - sex)
    }
    mod <- coxph(form, data = .x)
    haz_ratio(mod, "pass")
  }) %>%
  ungroup()

hr3 <- ca_dat %>%
  group_by(pri_site) %>%
  group_modify(~ {
    form <- Surv(surv_m, is_dead) ~ pass + (1 | hospcode)
    mod <- coxme(form, data = .x)
    haz_ratio(mod, "pass")
  }) %>%
  ungroup()

hr4 <- ca_dat %>%
  group_by(pri_site) %>%
  group_modify(~ {
    form <- Surv(surv_m, is_dead) ~ pass + hosplv + sex + agegp + fnstage + region + do_surg + do_chem + do_rt + (1 | hospcode)
    if(.y$pri_site %in% c("2_breast", "5_cervix", "6_uterine", "7_ovary")) {
      form <- update(form, . ~ . - sex)
    }
    mod <- coxme(form, data = .x)
    haz_ratio(mod, "pass")
  }) %>%
  ungroup()

hr_dat <- bind_rows(hr1, hr2, hr3, hr4, .id = "model")

##################
#--- Table 1. ---#
##################

#--- Hazard Ratio ---#

hr_tab <- hr_dat %>%
  mutate(pval_star = symnum(pval, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', '')),
         hr_txt = sprintf("%.2f (%.2f-%.2f)%s", HR, HR_L, HR_U, pval_star)) %>% 
  pivot_wider(id_cols = pri_site, names_from = model, values_from = hr_txt)

#--- 5-year Cumulative Probability of Death ---#

risk5y_tab <- ca_dat %>%
  group_by(pri_site) %>%
  group_modify(~ {
    km_fit <- survfit(Surv(surv_m, is_dead) ~ pass, data = .x) %>%
      summary(times = 60)
    data.frame(strata = km_fit$strata,
               n = km_fit$n,
               cumrisk = (1 - km_fit$surv) * 100)
  }) %>%
  ungroup()

#--- Merge ---#

tab1 <- hr_tab %>%
  mutate(strata = "pass=1") %>% 
  left_join(risk5y_tab, ., by = join_by(pri_site, strata))

# openxlsx::write.xlsx(tab1, "Table_1.xlsx")

###################
#--- Figure 2. ---#
###################

ca_lookup <- c("1_colorectal" = "Colorectal Cancer",
               "2_breast" = "Female Breast Cancer",
               "3_lung" = "Lung Cancer",
               "4_liver" = "Liver Cancer",
               "5_cervix" = "Cervical Cancer",
               "6_uterine" = "Uterine Cancer",
               "7_ovary" = "Ovarian Cancer",
               "8_bladder" = "Bladder Cancer")

p <- hr_dat %>%
  mutate(model_txt1 = factor(model, labels = c("bold(Crude)",
                                               "bold('Adjusted for ' * confounders^a)",
                                               "bold('Adjusted for ' * frailty^b)",
                                               "bold(atop('Adjusted for ' * confounders^a, 'and ' * frailty^b))")),
         model_txt2 = factor(model, labels = c("bold('Crude')",
                                               "bold('Confounder-adjusted')",
                                               "bold('Frailty-adjusted')",
                                               "bold('Fully adjusted')"))) %>% 
  ggplot(aes(HR, pri_site)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "gray50") +
  geom_errorbar(aes(xmin = HR_L, xmax = HR_U), width = 0.2) +
  geom_point(aes(colour = model, shape = model), size = 1.1, fill = "white", stroke = 0.3) +
  facet_grid(~ model_txt2, labeller = label_parsed) +
  scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 2), limits = c(0.5, 2),
                labels = ~ round(.x, 2)) +
  scale_y_discrete(limits = rev, labels = ca_lookup) +
  scale_colour_manual(values = c("darkblue", "darkblue", "red4", "red4")) +
  scale_shape_manual(values = c(21, 16, 22, 15)) +
  labs(x = "\nHazard Ratio (HR)", y = NULL) +
  theme_minimal(base_size = 8, base_line_size = 0.25) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 7),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        legend.position = "none")

# ggsave("figure/Figure_2.pdf", p, width = 7, height = 2.8, dpi = 600)

#################################
#--- Supplementary Materials ---#
#################################

#--- Survival Curve by Stage ---#

survplot_by_stage <- function(data, stage) {
  subdat <- data %>%
    filter(fnstage == .env$stage)
  survplot <- subdat %>%
    group_by(pri_site) %>%
    group_map(~ {
      site <- .y$pri_site
      label <- letters[match(site, names(ca_lookup))]
      hr <- exp(coef(coxph(Surv(surv_m, is_dead) ~ pass, data = .x)))
      km_fit <- survfit(Surv(surv_m, is_dead) ~ pass, data = .x)
      ggsurvplot(
        km_fit, data = .x,
        censor = FALSE, conf.int = TRUE, pval = TRUE, pval.method = TRUE,
        legend.labs = c("Non-certified hospitals", "Certified hospitals"),
        legend.title = "", xlab = "Year", ylab = "Survival probability",
        title = sprintf("(%s) %s", label, ca_lookup[site]),
        break.x.by = 12, break.y.by = 0.2, xscale = 12, surv.scale = "percent",
        xlim = c(0, 60), ylim = c(0, 1),
        size = 0.5, font.main = 12, font.tickslab = 9,
        pval.method.coord = c(0, 0.05), pval.coord = c(12, 0.05),
        pval.size = 3, pval.method.size = 3
      )$plot +
        geom_text(x = 0, y = 0.15, label = sprintf("HR = %.2f", hr),
                  hjust = 0, size = 3)
    })
  
  guide_area() /
    wrap_plots(survplot, ncol = 2, axis_titles = "collect") +
    plot_layout(guides = "collect", heights = c(1, 15)) &
    theme(legend.position = "top")
}

p1 <- survplot_by_stage(ca_dat, stage = "1")
p2 <- survplot_by_stage(ca_dat, stage = "2")
p3 <- survplot_by_stage(ca_dat, stage = "3")
p4 <- survplot_by_stage(ca_dat, stage = "4")

# ggsave("figure/Figure_S1.jpeg", p1, width = 7, height = 10, dpi = 600)
# ggsave("figure/Figure_S2.jpeg", p2, width = 7, height = 10, dpi = 600)
# ggsave("figure/Figure_S3.jpeg", p3, width = 7, height = 10, dpi = 600)
# ggsave("figure/Figure_S4.jpeg", p4, width = 7, height = 10, dpi = 600)

#--- Baseline Characteristics ---#

baseline_char_table <- function(data, var) {
  
  if(is.numeric(pull(data, {{ var }}))) {
    tab <- data %>%
      mutate(level = NA) %>% 
      summarise(x = {
        q <- round(quantile({{ var }}, c(0.25, 0.50, 0.75), na.rm = TRUE))
        sprintf("%d\n(%d-%d)", q[2], q[1], q[3])
      }, .by = c(pri_site, pass, level)) %>% 
      pivot_wider(id_cols = level, names_from = c(pri_site, pass), values_from = x,
                  names_sort = TRUE, names_expand = TRUE)
  } else {
    tab <- data %>%
      mutate(level = {{ var }}) %>% 
      count(pri_site, pass, level) %>% 
      pivot_wider(id_cols = level, names_from = c(pri_site, pass), values_from = n, values_fill = 0,
                  id_expand = TRUE, names_sort = TRUE, names_expand = TRUE) %>%
      mutate(across(-1, ~ sprintf("%s\n(%.1f%%)",
                                  format(.x, big.mark = ",", trim = TRUE),
                                  .x / sum(.x, na.rm = TRUE) * 100)))
  }
  
  cbind(var = deparse(substitute(var)), tab)
}

tabS2 <- bind_rows(baseline_char_table(mutate(ca_dat, all = 'all'), all),
                   baseline_char_table(ca_dat, age),
                   baseline_char_table(ca_dat, sex),
                   baseline_char_table(ca_dat, fnstage),
                   baseline_char_table(ca_dat, hosplv),
                   baseline_char_table(ca_dat, region),
                   baseline_char_table(mutate(ca_dat, do_surg = factor(do_surg)), do_surg),
                   baseline_char_table(mutate(ca_dat, do_chem = factor(do_chem)), do_chem),
                   baseline_char_table(mutate(ca_dat, do_rt = factor(do_rt)), do_rt),
                   baseline_char_table(mutate(ca_dat, do_target = factor(do_target)), do_target),
                   baseline_char_table(mutate(ca_dat, do_palli = factor(do_palli)), do_palli))

# openxlsx::write.xlsx(tabS2, "Table_S2.xlsx")
