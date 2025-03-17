# Make plots for manuscript
library(tidyverse)
library(splines2)
library(foreach)
library(ggplot2)
library(latex2exp)
library(ggnewscale)
rm(list = ls())
gc(verbose = FALSE)
load("biocard_result_group20nonzeros.RData")
load("AD_Diagnose_Age.rda")
load("DECAGE.rda")

bionames <- c(
  "MMSE", "LogMem", "DSST", "ENT-THICK", "HIPPO", "ENT-VOL",
  "MTL", "SPARE-AD", "t-tau", "p-tau181", "AB42/AB40 Ratio"
)

# Age effect from age 50 to 90 years old
age_59 <- c(50, 90)
basis_59 <- ibs(age_59,
  knots = knot.list[[1]], Boundary.knots = boundary.knot,
  degree = 2, intercept = TRUE
)
basis_59 <- basis_59[, 3:(ncol(basis_59) - 2)]
incre_59 <-
  foreach(i = 1:K, .combine = rbind) %do%
  {
    curves <- basis_59 %*% coefs[-(1:4), i, (R / 2 + 1):R]
    incre <- apply(curves, 2, function(x) x[2] - x[1])
    incre_est <- c(mean = mean(incre), HDInterval::hdi(incre))
    incre_est
  } %>%
  as.data.frame() %>%
  mutate_all(round, digits = 2) %>%
  mutate(biomarker = bionames) %>%
  remove_rownames() %>%
  column_to_rownames("biomarker")

# Make Point-CI plot for covariate effect ----
indice <- (R / 2 + 1):R
summary_fixed <- apply(
  coefs[2:4, , indice], c(1, 2),
  function(x) c(mean = mean(x), HDInterval::hdi(x))
)
dimnames(summary_fixed)[[2]] <- c("APOE4 vs Non-APOE4", "Female vs Male", "Education")
dimnames(summary_fixed)[[1]] <- c("mean", "lower", "upper")
dimnames(summary_fixed)[[3]] <- bionames
# Education effect from 12 to 16 years
summary_fixed[, 3, ] <- (
  summary_fixed[, 3, ] * (max(df$education) - min(df$education)) / 2
)
# 3D-array containing credible intervals
summary_fixed <- abind::abind(summary_fixed,
  Age = array(t(incre_59), dim = c(3, 1, 11)),
  along = 2
)

sf <- t(Reduce(cbind, list(
  summary_fixed[, , 1],
  summary_fixed[, , 2],
  summary_fixed[, , 3],
  summary_fixed[, , 4],
  summary_fixed[, , 5],
  summary_fixed[, , 6],
  summary_fixed[, , 7],
  summary_fixed[, , 8],
  summary_fixed[, , 9],
  summary_fixed[, , 10],
  summary_fixed[, , 11]
)))
# Print LaTeX-style table of fixed effects
sf <- data.frame(sf,
  variable = rownames(sf), biomarker = rep(bionames, each = 4),
  x = rep(c(-0.15, -0.05, 0.05, 0.15), 11) + rep(1:11, each = 4),
  row.names = NULL
)
sink("table.txt")
sf %>%
  mutate(text = sprintf("% 5.2f(% 5.2f,% 5.2f)", mean, lower, upper)) %>%
  select(biomarker, text, variable) %>%
  pivot_wider(names_from = biomarker, values_from = text) %>%
  t() %>%
  janitor::row_to_names(1) %>%
  xtable::xtable() %>%
  print()
sink()

# Make Goodness-of-fit Plot for select variables
fitted_values <- array(0, dim(Y))
colnames(fitted_values) <- bionames
colnames(Y) <- bionames
for (i in 1:K) {
  fit <- covar.list[[i]] %*% coefs[, i, indice] + REs[df$ID, i, indice]
  fitted_values[, i] <- apply(fit, 1, mean)
  fitted_values[is.na(Y[, i]), i] <- NA
}
fitted_values <- as.data.frame(fitted_values)
Y$type <- "Observed"
Y$id <- df$ID
Y$age <- df$ageori
fitted_values$type <- "Fitted"
fitted_values$id <- df$ID
fitted_values$age <- df$ageori
select_biom <- c(-1, -2, -4, -5, -7, -8, -10, -11)
gof_data <- rbind(Y, fitted_values)[, select_biom]
pdf("Goodness_of_fit_Grid3.pdf")
print(
  gof_data %>% gather("Biomarker", "Value", DSST, "ENT-VOL", "t-tau") %>%
    mutate(Biomarker = ifelse(Biomarker == 't-tau', 'Total tau', Biomarker)) %>% 
    ggplot() +
    geom_line(aes(x = age, y = Value, group = interaction(id, type), color = type),
      alpha = 0.3
    ) +
    facet_wrap(vars(Biomarker), ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_line(colour = "grey75"),
      panel.grid.minor = element_line(colour = "grey90")
    ) +
    scale_color_manual(values = c(Observed = "black", Fitted = "darkred"), name = NULL) +
    scale_x_continuous(limits = c(50, 100)) +
    scale_y_continuous(name = 'Standardized Abnormality Scores')
)
dev.off()
# Make CI plot for Standardized Biomarker Trajectory
agex <- seq(from = min(boundary.knot), to = max(boundary.knot), by = 0.1)
basis <- ibs(agex,
  knots = knot.list[[1]], Boundary.knots = boundary.knot,
  degree = 2, intercept = TRUE
)
basis <- basis[, 3:(ncol(basis) - 2)]
spline_std <-
  foreach(i = 1:K, .combine = rbind) %do%
  {
    curves <- basis %*% coefs[-(1:4), i, indice]
    curves <- apply(curves, 2, function(x) x / max(x))
    frame <- apply(curves, 1, function(x) c(value = mean(x), HDInterval::hdi(x))) %>%
      t() %>%
      as.data.frame() %>%
      mutate(age = agex, Biomarker = bionames[i])
    frame
  } %>% mutate(Biomarker = factor(Biomarker, levels = bionames), type = "Curve")

spline_std <- rbind(spline_std,
  foreach(i = 1:K, .combine = rbind) %do%
    {
      curves <- basis %*% coefs[-(1:4), i, indice]
      curves <- apply(curves, 2, function(x) x / max(x))
      curve_mean <- apply(curves, 1, mean)
      frame <- apply(curves, 2, function(x) agex[max(which(diff(x, differences = 2) >= 0)) + 1]) %>%
        mean() %>%
        as.data.frame() %>%
        mutate(Biomarker = bionames[i]) %>%
        rename(age = 1)
      frame$value <- approx(agex, curve_mean, frame$age)$y
      frame$lower <- NA
      frame$upper <- NA
      frame
    } %>% mutate(Biomarker = factor(Biomarker, levels = bionames), type = "Mark"),
  fill = TRUE
) %>% filter(!is.na(Biomarker))

colors0 <- c(
  "darkred", "red", "orange", "yellow", "green", "darkgreen", "turquoise4",
  "royalblue", "blue", "purple", "black"
)
names(colors0) <- bionames
print_labels <- unname(TeX(c(
  "A$\\beta$ ratio", "p-tau181", "Total tau", "ENT-THICK", "ENT-VOL",
  "MTL", "SPARE-AD", "HIPPO", "DSST", "MMSE", "LogMem"
)))
names(print_labels) <- bionames

inflect_order <- spline_std %>%
  filter(type == "Mark") %>%
  select(age) %>%
  as.matrix() %>%
  as.vector() %>%
  order(decreasing = FALSE)
bionames <- bionames[inflect_order]
colors0 <- colors0[inflect_order]
spline_std <- mutate(spline_std, Biomarker = factor(Biomarker, levels = bionames))
# spline_std <- mutate(
#   spline_std, Biomarker = case_when(
#     Biomarker == 'AB42/AB40 Ratio' ~ 'Abeta ratio',
#     Biomarker == 't-tau' ~ 'Total tau',
#     .default = Biomarker
#   )
# )
p1 <-
  ggplot() +
  geom_ribbon(aes(x = age, ymin = lower, ymax = upper, group = Biomarker, fill = Biomarker),
              alpha = 0.1,
              data = subset(spline_std, type == "Curve"&Biomarker%in%bionames[1:3])
  ) +
  geom_line(aes(x = age, y = value, group = Biomarker, color = Biomarker), 
            data = subset(spline_std, type == "Curve"&Biomarker%in%bionames[1:3])) +
  geom_point(aes(x = age, y = value, color = Biomarker), 
             data = subset(spline_std, type == "Mark"&Biomarker%in%bionames[1:3]), 
             shape = "I", size = 4) +
  scale_color_manual(values = colors0[1:3], labels = print_labels[1:3], 
                     breaks = bionames[1:3], name = 'Biomarker',
                     guide = guide_legend(order=1)) +
  scale_fill_manual(values = colors0[1:3], labels = print_labels[1:3], 
                    breaks = bionames[1:3], name = 'Biomarker',
                    guide = guide_legend(order=1)) +
  
  new_scale_fill() + new_scale_color() +
  geom_ribbon(aes(x = age, ymin = lower, ymax = upper, group = Biomarker, fill = Biomarker),
              alpha = 0.1,
              data = subset(spline_std, type == "Curve"&Biomarker%in%bionames[4:8])
  ) +
  geom_line(aes(x = age, y = value, group = Biomarker, color = Biomarker), 
            data = subset(spline_std, type == "Curve"&Biomarker%in%bionames[4:8])) +
  geom_point(aes(x = age, y = value, color = Biomarker), 
             data = subset(spline_std, type == "Mark"&Biomarker%in%bionames[4:8]), 
             shape = "I", size = 4) +
  scale_color_manual(values = colors0[4:8], labels = print_labels[4:8], 
                     breaks = bionames[4:8], name = ' ',
                     guide = guide_legend(order=2)) +
  scale_fill_manual(values = colors0[4:8], labels = print_labels[4:8], 
                    breaks = bionames[4:8], name = ' ',
                    guide = guide_legend(order=2)) +
  
  new_scale_fill() + new_scale_color() +
  geom_ribbon(aes(x = age, ymin = lower, ymax = upper, group = Biomarker, fill = Biomarker),
              alpha = 0.1,
              data = subset(spline_std, type == "Curve"&Biomarker%in%bionames[9:11])
  ) +
  geom_line(aes(x = age, y = value, group = Biomarker, color = Biomarker), 
            data = subset(spline_std, type == "Curve"&Biomarker%in%bionames[9:11])) +
  geom_point(aes(x = age, y = value, color = Biomarker), 
             data = subset(spline_std, type == "Mark"&Biomarker%in%bionames[9:11]), 
             shape = "I", size = 4) +
  scale_color_manual(values = colors0[9:11], labels = print_labels[9:11], 
                     breaks = bionames[9:11], name = '  ',
                     guide = guide_legend(order=3)) +
  scale_fill_manual(values = colors0[9:11], labels = print_labels[9:11], 
                    breaks = bionames[9:11], name = '  ',
                    guide = guide_legend(order=3)) +
  
  
  scale_x_continuous(limits = c(50, 100)) +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Abnormality") +
  theme(
    legend.justification = c(0, 0.5),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 16),
    title = element_text(size = 16)
  )

p2 <-
  rbind(cbind(age = AD_Diagnose_Age, status = "Dementia"), cbind(age = decs, status = "Onset")) %>%
  as.data.frame() %>%
  filter(!is.na(age)) %>%
  mutate(age = as.double(age)) %>%
  ggplot() +
  geom_density(aes(x = age, group = status, fill = status), bounds = c(50, 100), alpha = 0.3) +
  scale_x_continuous(limits = c(50, 100)) +
  scale_fill_manual(
    values = c("Dementia" = "black", "Onset" = "grey"),
    labels = c(
      sprintf("Dementia(n=%d)", sum(!is.na(AD_Diagnose_Age))),
      sprintf("Symptom Onset(n=%d)", sum(!is.na(decs)))
    )
  ) +
  theme_bw() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 16),
    title = element_text(size = 16),
    legend.justification = c(0, 0.5)
  )
pdf("SplineStd.pdf", width = 14)
print(
  cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(7, 3))
)
dev.off()
