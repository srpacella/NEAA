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
  mutate(across(c(do_mgl, pH), ~ ifelse(. == -999, NA, .))) %>%
  filter(!is.na(do_mgl), !is.na(pH))

# ── 2. Order NEPs west → east ─────────────────────────────────────────────────
site_order <- df %>%
  mutate(lon = as.numeric(lon)) %>%
  group_by(NEP) %>%
  summarise(med_lon = median(lon, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_lon) %>%
  pull(NEP)

df$NEP <- factor(df$NEP, levels = site_order)

# ── 3. Local site index (resets per NEP) ─────────────────────────────────────
df <- df %>%
  group_by(NEP) %>%
  mutate(site_idx = as.factor(as.integer(factor(site_code)))) %>%
  ungroup()

# ── 4. Colour palette ─────────────────────────────────────────────────────────
base_pal <- c("#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3","#A65628","#F781BF","#999999","#000000")
n_sites  <- max(as.integer(df$site_idx), na.rm = TRUE)
site_pal <- setNames(rep_len(base_pal, n_sites), as.character(1:n_sites))

# ── 5. Theme ──────────────────────────────────────────────────────────────────
pp_theme <- theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "grey92"),
    axis.text        = element_text(size = 8),
    axis.title       = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9.5, colour = "grey40"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7.5),
    legend.key.size  = unit(0.4, "cm")
  )

# ── 6. DO vs pH ───────────────────────────────────────────────────────────────
p <- df %>%
  ggplot(aes(x = pH, y = do_mgl, colour = site_idx)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = site_pal, name = "Site") +
  facet_wrap(~ NEP, ncol = 3, scales = "free") +
  labs(
    title    = "DO vs. pH — 15 National Estuary Programs",
    subtitle = "ph_T used where available  ·  sites ordered west → east  ·  colour = site_code (resets per panel)",
    x        = "pH",
    y        = "DO (mg/L)"
  ) +
  pp_theme +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

ggsave("fig_do_vs_ph.png", p, width = 14, height = 17, dpi = 180)
message("Saved fig_do_vs_ph.png")