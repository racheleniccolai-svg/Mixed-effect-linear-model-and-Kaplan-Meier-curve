# Mixed effect linear model and Kaplan-Meier curve
R script for mixed-effects linear modeling of tumor volume and Kaplan–Meier survival analysis in in vivo experiments.

# Loading and preprocessing 
data <- read_excel("data/LSS_invivo_2.xlsx", sheet = "tumorsize")
data_surv <- read_excel("data/LSS_invivo_2.xlsx", sheet = "Sheet1")

data2 <-data %>% separate(
Group,
into = c("cell_line", "treatment"),
sep = "\\s*[._-]?\\s*(?=(?:vehicule|vehicle|GSK)$)",
remove = FALSE,
extra = "merge",
fill = "right"
) %>%
mutate(treatment = ifelse(treatment== "vehicule", "vehicle", treatment))
data2$treatment[is.na(data2$treatment)] <- "vehicle"
data2$cell_line <- as.factor(data2$cell_line)
data2$ID <- as.factor(data2$ID)
data2$treatment <- as.factor(data2$treatment)

data2 <- data2 %>% group_by(ID) %>%        # Group by mouse ID
arrange(Day) %>%                           # Ensure data is ordered by date for each mouse
mutate(first_volume = first(Size),         # Get the first tumor volume for each mouse
normalized_volume = Size / first_volume)   # Normalize the tumor volume

data2$cell_line <- sub("ˆ\\s*\\d+\\.\\s*", "", data2$cell_line)
data2 <- droplevels(data2)
table(data2$cell_line, data2$treatment)

# EDA: Trend pairwise comparison
p_1 <- ggplot(data2, aes(x = Day, y = Size, group = ID, color = Group)) +
geom_line() +
labs(title="Tumor Volume per Animal by Group",
x="Day", y="Tumor Volume") +
theme_minimal() + theme(legend.position="right")

#plotly::ggplotly(p_1)
p_1

p_2 <- ggplot(data2, aes(x = Day, y = log(normalized_volume), group = ID, color = Group)) +
geom_line() +
labs(title="Normalized Tumor Volume per Animal by Group",
x="Day", y="Log(Normalized Volume)") +
theme_minimal() + theme(legend.position="right")

#plotly::ggplotly(p_2)
p_2

fit <- lme(
random =~ 1 | ID,
data = data2,
method = "REML"
fixed = log(normalized_volume)~ 1 + Day + Day:treatment + Day:cell_line + Day:treatment:cell_line,)

table <- summary(fit)$tTable
#DT::datatable(round(table, 4))
print(round(table, 4))







