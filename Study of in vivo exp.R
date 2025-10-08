#Make analysis of in vivo study
library(readxl)
library(dplyr)
library(tidyverse)

#Upload data
data <- read_excel("/Users/racheleniccolai/Desktop/Statistics_in_vivo/LSS_invivo_2.xlsx", sheet = "tumorsize")
data_surv <- read_excel("/Users/racheleniccolai/Desktop/Statistics_in_vivo/LSS_invivo_2.xlsx", sheet = "Sheet1")

data2 <-data %>% separate(
  Group,
  into = c("cell_line", "treatment"),
  sep = "\\s*[._-]?\\s*(?=(?:vehicle|vehicle|GSK)$)",
  remove = FALSE,
  extra = "merge",
  fill = "right"
) %>%
  
  mutate(treatment = ifelse(treatment == "vehicle", "vehicle", treatment))
data2$treatment[is.na(data2$treatment)] <- "vehicle"

data2$cell_line <- as.factor(data2$cell_line)
data2$ID <- as.factor(data2$ID)
data2$treatment <- as.factor(data2$treatment)

data2 <- data2 %>%  group_by(ID) %>%                        # Group by mouse ID
  arrange(Day) %>%                           # Ensure data is ordered by date for each mouse
  mutate(first_volume = first(Size),          # Get the first tumor volume for each mouse
         normalized_volume = Size / first_volume)    # Normalize the tumor volume

data2$cell_line <- sub("^\\s*\\d+\\.\\s*", "", data2$cell_line)

data2 <- droplevels(data2)
table(data2$cell_line, data2$treatment)

data2$treatment <- factor(data2$treatment, levels = c("vehicle", "GSK"))  # <- and here
data2$cell_line <- factor(data2$cell_line, levels = c("LSS KO", "Oci-Ly1"))  # <- and here

p_1 <- ggplot(data2, aes(x = Day, y = Size, group = ID, color = Group)) +
  geom_line() +
  labs(title="Tumor Volume per Animal by Group",
       x="Day", y="Tumor Volume") +
  theme_minimal() + theme(legend.position="right")

plotly::ggplotly(p_1)

p_2 <- ggplot(data2, aes(x = Day, y = log(normalized_volume), group = ID, color = Group)) +
  geom_line() +
  labs(title="Normalized Tumor Volume per Animal by Group",
       x="Day", y="Log(Normalized Volume)") +
  theme_minimal() + theme(legend.position="right")

plotly::ggplotly(p_2)


library(nlme) 

fit <- lme(
  fixed = log(normalized_volume) ~ 1 + Day + Day:treatment + Day:cell_line + Day:treatment:cell_line,
  random = ~ 1 | ID,
  data   = data2, 
  method = "REML"
)

table <- summary(fit)$tTable
DT::datatable(round(table, 4))

# -------------------------------------------------------------------------
# Build a combined "Group" = treatment × cell_line for plotting & coloring
# -------------------------------------------------------------------------
data2$treatment  <- droplevels(factor(data2$treatment))
data2$cell_line  <- droplevels(factor(data2$cell_line))
data2$Group      <- interaction(data2$treatment, data2$cell_line, sep = " • ", drop = TRUE)

# --- base spaghetti plot (individual animals) ---
p_2 <- ggplot(data2, aes(x = Day,
                         y = log(normalized_volume),
                         group = ID,                # random-effect unit
                         color = Group)) +
  geom_line(alpha = 0.35) +
  labs(title = "Normalized Tumor Volume per Animal",
       subtitle = "Colored by treatment × cell line",
       x = "Day", y = "Log(Normalized Volume)") +
  theme_minimal() +
  theme(legend.position = "right")

# -------------------------------------------------------
# Population (fixed-effects) lines from the mixed model
# -------------------------------------------------------
# prediction grid over Day × treatment × cell_line
Day_seq <- seq(min(data2$Day, na.rm = TRUE),
               max(data2$Day, na.rm = TRUE),
               length.out = 100)

tlev <- levels(data2$treatment)
clev <- levels(data2$cell_line)

pred_df <- expand.grid(Day = Day_seq,
                       treatment = tlev,
                       cell_line = clev,
                       KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)

# create Group in pred_df to match plot mapping
pred_df$treatment <- factor(pred_df$treatment, levels = tlev)
pred_df$cell_line <- factor(pred_df$cell_line, levels = clev)
pred_df$Group     <- interaction(pred_df$treatment, pred_df$cell_line, sep = " • ", drop = TRUE)

# population-level predictions (no random effects)
pred_df$fit <- as.numeric(predict(fit, newdata = pred_df, level = 0))  # nlme::lme

# --- overlay the lines ---
p_2_lines <- p_2 +
  geom_line(data = pred_df,
            aes(x = Day, y = fit, 
                color = Group, 
                group = Group),
            linewidth = 1.3,
            inherit.aes = FALSE) +
  scale_color_manual(values = c("#00BFC4","#C77CFF","#F8766D","#7CAE00")) +
  labs(subtitle = "Solid lines = mixed-model population fits per treatment × cell line")


# palette = c("#F8766D","#00BFC4","#C77CFF", "#7CAE00"),


# show it
# p_2_lines
# Optional interactive:
plotly::ggplotly(p_2_lines)

# one row per subject
surv_df <- data2 %>%
  group_by(ID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(
    treatment = droplevels(factor(treatment)),
    cell_line = droplevels(factor(cell_line)),
    Group = interaction(treatment, cell_line, sep = " • ", drop = TRUE),
    status = 1L                       # all events occurred
  ) %>%
  select(ID, treatment, cell_line, Day, status, Group)

library(survival)

# Kaplan–Meier (will drop to 0 eventually; no censor ticks)
surv_obj <- Surv(time = surv_df$Day, event = surv_df$status)

# Cox PH with interaction
cox_fit <- coxph(surv_obj ~ treatment * cell_line, data = surv_df)
summary(cox_fit)

exp(cbind(HR = coef(cox_fit), confint(cox_fit)))

cox.zph(cox_fit)

library(survminer)

km_fit <- survfit(surv_obj ~ Group, data = surv_df)

# KM (no built-in median lines)
km_plot <- ggsurvplot(
  km_fit, data = surv_df, risk.table = TRUE, conf.int = FALSE,
  pval = TRUE, ggtheme = theme_minimal(),
  xlab = "Time (days)", ylab = "Survival probability",
  surv.median.line = "none",
  xlim = c(0, 35),
  palette = c("#00BFC4", "#C77CFF","#F8766D","#7CAE00")  
)

surv_df <- surv_df[surv_df$Day <= 35, ]

# Per-group medians
med <- survminer::surv_median(km_fit)  # columns: strata, median, etc.

km_plot$plot <- km_plot$plot +
  geom_segment(data = med,
               aes(x = 0, xend = median, y = 0.5, yend = 0.5, color = strata)) +
  geom_segment(data = med,
               aes(x = median, xend = median, y = 0.5, yend = 0, color = strata),
               linetype = "dashed")

print(km_plot)


scale_x_continuous(limits = c(0, 35), breaks = seq(0, 35, 10))

# 1) Pairwise log-rank tests (BH-adjusted p-values)
pw_lr <- pairwise_survdiff(Surv(Day, status) ~ Group,
                           data = surv_df, p.adjust.method = "bonferroni")
pw_lr$p.value  # matrix of adjusted p-values


