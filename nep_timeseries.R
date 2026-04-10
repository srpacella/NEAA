library(dplyr)
library(ggplot2)
library(lubridate)

# ── 1. Load & clean ───────────────────────────────────────────────────────────
df <- bind_rows(lapply(nep_first1000, function(x) mutate(x, across(everything(), as.character)))) %>%
  mutate(
    NEP       = trimws(NEP),
    site_code = trimws(site_code),
    do_mgl    = as.numeric(do_mgl),
    ph        = as.numeric(ph),
    ph_T      = as.numeric(ph_T),
    pH        = ifelse(!is.na(ph_T), ph_T, ph)
  ) %>%
  mutate(across(c(do_mgl, pH), ~ ifelse(. == -999, NA, .)))

df$datetime <- parse_date_time(df$datetime_utc,
                               orders = c("ymd HMS", "ymd HM", "mdy HM", "mdy HMS"),
                               tz = "UTC", quiet = TRUE)

# ── 2. Assign a local site index within each NEP (resets per NEP) ─────────────
# site_idx = 1, 2, 3... based on alphabetical order of site_code within each NEP
df <- df %>%
  group_by(NEP) %>%
  mutate(site_idx = as.factor(as.integer(factor(site_code)))) %>%
  ungroup()

cat("Max sites in one NEP:", max(as.integer(df$site_idx), na.rm = TRUE), "\n")

# ── 3. Order NEPs west → east ─────────────────────────────────────────────────
site_order <- df %>%
  mutate(lon = as.numeric(lon)) %>%
  group_by(NEP) %>%
  summarise(med_lon = median(lon, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_lon) %>%
  pull(NEP)

df$NEP <- factor(df$NEP, levels = site_order)

# ── 4. Colour palette — fixed set cycling per site index ─────────────────────
n_sites <- max(as.integer(df$site_idx), na.rm = TRUE)
base_pal <- c("#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3","#A65628","#F781BF","#999999","#000000")
site_pal <- setNames(
  rep_len(base_pal, n_sites),
  as.character(1:n_sites)
)

# ── 5. Theme ──────────────────────────────────────────────────────────────────
ts_theme <- theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "grey92"),
    axis.text.x      = element_text(size = 7.5, angle = 30, hjust = 1),
    axis.text.y      = element_text(size = 8),
    axis.title.y     = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9.5, colour = "grey40"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7.5),
    legend.key.size  = unit(0.4, "cm")
  )

# ── 6. pH time series ─────────────────────────────────────────────────────────
p_ph <- df %>%
  filter(!is.na(pH), !is.na(datetime)) %>%
  ggplot(aes(x = datetime, y = pH, colour = site_idx, label = site_code)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_manual(values = site_pal, name = "Site (local rank)") +
  scale_x_datetime(date_labels = "%b '%y") +
  facet_wrap(~ NEP, ncol = 3, scales = "free_x") +
  labs(title    = "pH Time Series — 15 National Estuary Programs",
       subtitle = "ph_T used where available, otherwise ph  ·  colour restarts per panel (site 1, 2, ... by site_code alphabetical order)",
       x = NULL, y = "pH") +
  ts_theme +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

ggsave("fig_ph_timeseries.png", p_ph, width = 14, height = 17, dpi = 180)
message("Saved fig_ph_timeseries.png")

# ── 7. DO time series ─────────────────────────────────────────────────────────
p_do <- df %>%
  filter(!is.na(do_mgl), !is.na(datetime)) %>%
  ggplot(aes(x = datetime, y = do_mgl, colour = site_idx)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_manual(values = site_pal, name = "Site (local rank)") +
  scale_x_datetime(date_labels = "%b '%y") +
  facet_wrap(~ NEP, ncol = 3, scales = "free_x") +
  labs(title    = "Dissolved Oxygen Time Series — 15 National Estuary Programs",
       subtitle = "colour restarts per panel (site 1, 2, ... by site_code alphabetical order)",
       x = NULL, y = "DO (mg/L)") +
  ts_theme +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

ggsave("fig_do_timeseries.png", p_do, width = 14, height = 17, dpi = 180)
message("Saved fig_do_timeseries.png")
