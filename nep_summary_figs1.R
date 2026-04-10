library(dplyr)
library(ggplot2)
library(gridExtra)

# ── 1. Combine list — coerce all columns to character first to avoid type conflicts
df <- bind_rows(lapply(nep_first1000, function(x) mutate(x, across(everything(), as.character)))) %>%
  mutate(NEP     = trimws(NEP),
         temp_c  = as.numeric(temp_c),
         sal_ppt = as.numeric(sal_ppt),
         do_mgl  = as.numeric(do_mgl),
         ph      = as.numeric(ph),
         ph_T    = as.numeric(ph_T),
         lon     = as.numeric(lon),
         pH      = ifelse(!is.na(ph_T), ph_T, ph)) %>%
  # Remove -999 sentinel values
  mutate(across(c(temp_c, sal_ppt, do_mgl, ph, ph_T, pH), ~ ifelse(. == -999, NA, .)))

cat("Total rows:", nrow(df), "\n")
cat("NEPs:", paste(sort(unique(df$NEP)), collapse = ", "), "\n")
cat("Non-NA — temp:", sum(!is.na(df$temp_c)),
    "| sal:", sum(!is.na(df$sal_ppt)),
    "| do:",  sum(!is.na(df$do_mgl)),
    "| pH:",  sum(!is.na(df$pH)), "\n")

# ── 2. Order sites west → east by median longitude ───────────────────────────
site_order <- df %>%
  group_by(NEP) %>%
  summarise(med_lon = median(lon, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_lon) %>%
  pull(NEP)

df$NEP <- factor(df$NEP, levels = site_order)

# ── 3. Colour palette ─────────────────────────────────────────────────────────
pal <- setNames(
  colorRampPalette(c("#08306b","#2171b5","#6baed6","#006d2c",
                     "#41ab5d","#a1d99b","#a50f15","#fb6a4a"))(nlevels(df$NEP)),
  levels(df$NEP)
)

# ── 4. Theme + plot function ──────────────────────────────────────────────────
nep_theme <- theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 40, hjust = 1, size = 9),
    axis.title         = element_text(size = 11),
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = 12)
  )

make_violin <- function(var, ylab, title) {
  d <- df[!is.na(df[[var]]), ]
  ggplot(d, aes(x = NEP, y = .data[[var]], fill = NEP)) +
    geom_violin(trim = TRUE, alpha = 0.75, colour = "white", linewidth = 0.3) +
    geom_boxplot(width = 0.1, outlier.size = 0.6, outlier.alpha = 0.3,
                 colour = "grey20", fill = "white", linewidth = 0.35) +
    scale_fill_manual(values = pal) +
    labs(title = title, x = NULL, y = ylab) +
    nep_theme
}

# ── 5. Build & save individual plots ─────────────────────────────────────────
p_temp <- make_violin("temp_c",  "Temperature (°C)",        "Temperature")
p_sal  <- make_violin("sal_ppt", "Salinity (ppt)",          "Salinity")
p_do   <- make_violin("do_mgl",  "Dissolved Oxygen (mg/L)", "Dissolved Oxygen")
p_ph   <- make_violin("pH",      "pH",                      "pH  (ph_T preferred)")

ggsave("fig_temp.png", p_temp, width = 12, height = 5.5, dpi = 180)
ggsave("fig_sal.png",  p_sal,  width = 12, height = 5.5, dpi = 180)
ggsave("fig_do.png",   p_do,   width = 12, height = 5.5, dpi = 180)
ggsave("fig_ph.png",   p_ph,   width = 12, height = 5.5, dpi = 180)

# ── 6. Four-panel summary via gridExtra ──────────────────────────────────────
png("fig_nep_summary_4panel.png", width = 20, height = 11, units = "in", res = 180)
grid.arrange(p_temp, p_sal, p_do, p_ph, ncol = 2,
             top = "National Estuary Programs — Water Quality Summary (n = 1000 per site)")
dev.off()

message("Done.")