library(dplyr)
library(ggplot2)

# ── 1. Load & clean ───────────────────────────────────────────────────────────
df <- bind_rows(lapply(nep_first1000, function(x) mutate(x, across(everything(), as.character)))) %>%
  mutate(
    NEP    = trimws(NEP),
    sal_ppt = as.numeric(sal_ppt),
    temp_c  = as.numeric(temp_c),
    do_mgl  = as.numeric(do_mgl),
    ph     = as.numeric(ph),
    ph_T   = as.numeric(ph_T),
    pH     = ifelse(!is.na(ph_T), ph_T, ph)
  ) %>%
  mutate(across(c(do_mgl, pH, sal_ppt, temp_c), ~ ifelse(. == -999, NA, .)))

# ── 2. Order NEPs west → east ─────────────────────────────────────────────────
site_order <- df %>%
  mutate(lon = as.numeric(lon)) %>%
  group_by(NEP) %>%
  summarise(med_lon = median(lon, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_lon) %>%
  pull(NEP)

df$NEP <- factor(df$NEP, levels = site_order)

# ── 3. Colour palette (one per NEP) ──────────────────────────────────────────
nep_pal <- setNames(
  colorRampPalette(c("#08306b","#2171b5","#6baed6","#006d2c",
                     "#41ab5d","#a1d99b","#a50f15","#fb6a4a"))(nlevels(df$NEP)),
  levels(df$NEP)
)

# ── 4. Theme ──────────────────────────────────────────────────────────────────
ecdf_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.5, "cm")
  )

# ── 5. DO ECDF ────────────────────────────────────────────────────────────────
p_do <- df %>%
  filter(!is.na(do_mgl)) %>%
  ggplot(aes(x = do_mgl, colour = NEP)) +
  stat_ecdf(linewidth = 0.8) +
  scale_colour_manual(values = nep_pal, name = "NEP") +
  labs(
    title    = "Empirical Cumulative Distribution — Dissolved Oxygen",
    subtitle = "All sites within each NEP combined  ·  n = 1000 per NEP",
    x        = "DO (mg/L)",
    y        = "Cumulative probability"
  ) +
  ecdf_theme

ggsave("fig_ecdf_do.png", p_do, width = 10, height = 7, dpi = 180)
message("Saved fig_ecdf_do.png")

# ── 6. pH ECDF ────────────────────────────────────────────────────────────────
p_ph <- df %>%
  filter(!is.na(pH)) %>%
  ggplot(aes(x = pH, colour = NEP)) +
  stat_ecdf(linewidth = 0.8) +
  scale_colour_manual(values = nep_pal, name = "NEP") +
  labs(
    title    = "Empirical Cumulative Distribution — pH",
    subtitle = "ph_T used where available  ·  All sites within each NEP combined  ·  n = 1000 per NEP",
    x        = "pH",
    y        = "Cumulative probability"
  ) +
  ecdf_theme

ggsave("fig_ecdf_ph.png", p_ph, width = 10, height = 7, dpi = 180)
message("Saved fig_ecdf_ph.png")

# ── 7. Salinity ECDF ──────────────────────────────────────────────────────────
p_sal <- df %>%
  filter(!is.na(sal_ppt)) %>%
  ggplot(aes(x = sal_ppt, colour = NEP)) +
  stat_ecdf(linewidth = 0.8) +
  scale_colour_manual(values = nep_pal, name = "NEP") +
  labs(
    title    = "Empirical Cumulative Distribution — Salinity",
    subtitle = "All sites within each NEP combined  ·  n = 1000 per NEP",
    x        = "Salinity (ppt)",
    y        = "Cumulative probability"
  ) +
  ecdf_theme

ggsave("fig_ecdf_sal.png", p_sal, width = 10, height = 7, dpi = 180)
message("Saved fig_ecdf_sal.png")

# ── 8. Temperature ECDF ───────────────────────────────────────────────────────
p_temp <- df %>%
  filter(!is.na(temp_c)) %>%
  ggplot(aes(x = temp_c, colour = NEP)) +
  stat_ecdf(linewidth = 0.8) +
  scale_colour_manual(values = nep_pal, name = "NEP") +
  labs(
    title    = "Empirical Cumulative Distribution — Temperature",
    subtitle = "All sites within each NEP combined  ·  n = 1000 per NEP",
    x        = "Temperature (°C)",
    y        = "Cumulative probability"
  ) +
  ecdf_theme

ggsave("fig_ecdf_temp.png", p_temp, width = 10, height = 7, dpi = 180)
message("Saved fig_ecdf_temp.png")