###############################################################################
##  Application
##  "Plots of Annual Unemployment Rates (AURs) for (DZ, CA, IT, JP, USA, EG)"
##  Source : World Bank WDI, indicator SL.UEM.TOTL.ZS, 1991-2024
##  Style : fig05.pdf  (p1 colored lines, p2 colored box plots + points,
##                      p3 single grouped-bar ACF)
###############################################################################

# ── 1.  download database ────────────────────────────────────────────────────
meta <- data.frame(
  iso2c = c("DZ", "CA", "IT", "JP", "US", "EG"),
  label = c("Algeria", "Canada", "Italy", "Japan", "USA", "Egypt"),
  abbr  = c("DZ", "CA", "IT", "JP", "USA", "EG"),
  stringsAsFactors = FALSE
)

dat <- WDI(
  indicator = "SL.UEM.TOTL.ZS",
  country   = meta$iso2c,
  start = 1991, end = 2024,
  extra = FALSE
) |>
  rename(unemp = SL.UEM.TOTL.ZS) |>
  filter(!is.na(unemp)) |>
  left_join(meta, by = "iso2c") |>
  mutate(
    label = factor(label, levels = meta$label),
    abbr  = factor(abbr,  levels = meta$abbr)
  ) |>
  arrange(label, year)

pal_col <- c(
  Algeria = "#0000CD",
  Canada  = "#CC0000",
  Italy   = "#006400",
  Japan   = "#5555FF",
  USA     = "#FF4500",
  Egypt   = "#8B0000"
)

# ── Thème commun (style des graphes de simulation SRE) ───────────────────────
theme_fig <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid         = element_blank(),
      panel.grid.major   = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor   = element_line(color = "grey90", linewidth = 0.25),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.75),
      panel.spacing      = unit(0.02, "lines"),
      legend.position    = "top",
      legend.title       = element_blank(),
      legend.text        = element_text(size = 11, face = "bold"),
      legend.key.width   = unit(1.1, "cm"),
      legend.box         = "horizontal",
      legend.box.margin  = margin(t = -8),
      legend.box.spacing = unit(0.2, "cm"),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black")
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# PANNEAU 1  –  time series  (lignes colorées, sans formes)
# ══════════════════════════════════════════════════════════════════════════════
pal_col_abbr <- pal_col
names(pal_col_abbr) <- c("DZ", "CA", "IT", "JP", "USA", "EG")

p1 <- ggplot(dat, aes(x = year, y = unemp, color = abbr)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = pal_col_abbr, name = NULL) +
  scale_x_continuous(
    breaks = c(1991, 1995, 2000, 2005, 2010, 2015, 2020, 2025),
    limits = c(1991, 2025)
  ) +
  scale_y_continuous(breaks = seq(0, 35, by = 5), limits = c(0, 35)) +
  labs(x = "Year", y = "Annual Unemployment Rate") +
  guides(color = guide_legend(nrow = 1)) +
  theme_fig() +
  theme(plot.margin = margin(1, 13, 1, 13))
p1

# ══════════════════════════════════════════════════════════════════════════════
# PANNEAU 2  –  box plot  (boîtes colorées + points)
# ══════════════════════════════════════════════════════════════════════════════
p2 <- ggplot(dat, aes(x = abbr, y = unemp)) +
  stat_boxplot(geom = "errorbar", width = 0.28, linewidth = 0.45) +
  geom_boxplot(aes(fill = abbr),
               color         = "black",
               outlier.shape = NA,
               linewidth     = 0.45,
               width         = 0.55,
               alpha         = 0.65) +
  geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.5, color = "grey25") +
  scale_fill_manual(values = pal_col_abbr, guide = "none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5), limits = c(0, 35)) +
  labs(x = NULL, y = "Value") +
  theme_fig() +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        plot.margin = margin(1, 13, 1, 13))
p2

# ══════════════════════════════════════════════════════════════════════════════
# PANNEAU 3  –  Grille ACF / PACF 
# ══════════════════════════════════════════════════════════════════════════════

names(pal_col) <- c("DZ", "CA", "IT", "JP", "USA", "EG")
ordre_pays <- c("DZ", "CA", "IT", "JP", "USA", "EG")

cor_data <- dat %>%
  mutate(abbr = factor(abbr, levels = ordre_pays)) %>% 
  group_by(abbr) %>%
  do({
    a <- acf(.$unemp, lag.max = 12, plot = FALSE)
    p <- pacf(.$unemp, lag.max = 12, plot = FALSE)
    data.frame(
      lag = 0:12,
      acf = c(1, as.numeric(a$acf[-1, 1, 1])),
      pacf = c(1, as.numeric(p$acf[, 1, 1]))
    )
  }) %>%
  ungroup() %>%
  pivot_longer(cols = c(acf, pacf), names_to = "type", values_to = "valeur") %>%
  mutate(type = toupper(type))

p3 <- ggplot(cor_data, aes(x = factor(lag), y = valeur, fill = abbr)) +
  geom_bar(stat = "identity", alpha = 0.65, width = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = c(ci, -ci), color = "grey40", linetype = "dashed", linewidth = 0.4) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_fill_manual(values = pal_col) + 
  facet_grid(type ~ abbr) + 
  labs(x = "Lag", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    panel.spacing.x = unit(0, "lines"),
    panel.spacing.y = unit(0.2, "lines"),
    plot.margin = margin(0, -7, 0, 30) 
  )

p3
