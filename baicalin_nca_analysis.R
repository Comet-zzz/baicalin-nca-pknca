# Non-compartmental analysis (NCA) parameter estimation and graphical visualization

# Load required packages
library(PKNCA)
library(ggplot2)
library(dplyr)

# 1. Data entry
data_ug <- data.frame(
  subject = 1,
  time = c(0.17, 0.33, 0.5, 1, 2, 4, 6, 8, 10, 12, 16, 24),
  conc_ugL = c(320.5, 680.2, 1120.8, 1587.9, 1420.3, 980.6, 1150.4, 890.2, 720.5, 580.3, 310.7, 120.4)
)

# Unit conversion: μg/L → mg/L and select relevant columns
data <- data_ug %>%
  mutate(conc = conc_ugL / 1000) %>%1
  select(subject, time, conc)

# 2. Create concentration and dose objects
conc_obj <- PKNCAconc(data, conc ~ time | subject)
dose_data <- data.frame(subject = 1, time = 0, dose = 200)  # Oral dose: 200 mg/kg
dose_obj <- PKNCAdose(dose_data, dose ~ time | subject)

# 3. Combine and perform NCA calculations
nca_data <- PKNCAdata(conc_obj, dose_obj)
nca_result <- pk.nca(nca_data)

# 4. View complete results
summary(nca_result)

# 5. Extract selected parameters (without relying on PPUNITS)
params <- as.data.frame(nca_result)

# Uncomment to view column names
# colnames(params)

# Filter and retain only TESTCD and the numeric result
params_selected <- params %>%
  filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "aucinf.obs", 
                         "half.life", "lambda.z", "cl.obs", "vz.obs")) %>%
  select(PPTESTCD, PPORRES)

print(params_selected)

# 6. Generate concentration–time plots
p_linear <- ggplot(data, aes(x = time, y = conc)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, color = "red") +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       title = "Baicalin concentration-time profile (oral 200 mg/kg)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p_semilog <- ggplot(data, aes(x = time, y = conc)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, color = "red") +
  scale_y_log10() +
  labs(x = "Time (h)", y = "Concentration (mg/L) [log scale]",
       title = "Semi-log plot") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(p_linear)
print(p_semilog)

# Save plots to files
ggsave("baicalin_linear.png", p_linear, width = 6, height = 4)
ggsave("baicalin_semilog.png", p_semilog, width = 6, height = 4)
print(params_selected)


# Terminal slope diagnostic: assessing the rationale for λz selection

# Extract terminal data points (≥ 12 h) for regression
tail_data <- data %>% filter(time >= 12)

# Fit linear regression: ln(concentration) ~ time
lm_fit <- lm(log(conc) ~ time, data = tail_data)
r_squared <- summary(lm_fit)$r.squared

# Diagnostic plot: semi-log scatter with fitted regression line
p_diagnostic <- ggplot(tail_data, aes(x = time, y = log(conc))) +
  geom_point(size = 2, color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", alpha = 0.3) +
  labs(x = "Time (h)", y = "ln(Concentration) (mg/L)",
       title = paste0("Terminal phase regression (R² = ", round(r_squared, 4), ")")) +
  theme_minimal() +
  annotate("text", x = max(tail_data$time) * 0.7, y = max(log(tail_data$conc)) * 0.9,
           label = paste("λz =", round(-coef(lm_fit)[2], 3), "1/h"), size = 4)

print(p_diagnostic)
ggsave("terminal_slope_diagnostic.png", p_diagnostic, width = 5, height = 4)

