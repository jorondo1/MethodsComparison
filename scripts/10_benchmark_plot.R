# Resource benchmark: per-sample CPU cost of each method against input size, and peak memory of each step.
# Input: the output folder of scripts/benchmark/benchmark.sh (samples.tsv, bases.tsv, jobs.tsv, sacct.tsv).
# Output in the same folder:
#   resource_use.tsv         one row per sample x step: CPU-hours, wall time, peak memory
#   cost_model.tsv           per method: overhead per sample (a), CPU-hours per Gbp (b), 95 % CIs, peak memory
#   benchmark_resources.pdf  A: CPU-hours per sample against Gbp, with the fit; B: peak memory per step
# Usage, from the repo root: Rscript scripts/10_benchmark_plot.R

pacman::p_load(tidyverse, patchwork, magrittr)
source("scripts/0_config.R") # methods (colours)

dir <- "Out/benchmark"

# sacct labels (jobs.tsv, Kraken step names) to display names
tools <- c(
  sylph_sketch = "sylph sketch", SYref = "sylph RefSeq", SY220 = "sylph GTDB 220",
  sourmash_sketch = "sourmash sketch", SMref = "sourmash RefSeq", SM220 = "sourmash GTDB 220",
  KB10 = "Kraken-Bracken RefSeq 0.10", KB45 = "Kraken-Bracken RefSeq 0.45", KB90 = "Kraken-Bracken RefSeq 0.90",
  KG10 = "Kraken-Bracken GTDB 0.10", KG45 = "Kraken-Bracken GTDB 0.45", KG90 = "Kraken-Bracken GTDB 0.90",
  MOTU3 = "mOTUs 3", MOTU4 = "mOTUs 4",
  MPA22 = "MetaPhlAn 2022", MPA23 = "MetaPhlAn 2023", MPA25 = "MetaPhlAn 2025"
)

# Methods as a user runs them (see Cost model below), in figure order, with the profile object whose colour they take
# (methods table, 0_config.R)
benchmark_methods <- c(
  `sylph RefSeq` = "SYLPH_RefSeq_Seq", `sylph GTDB 220` = "SYLPH_gtdb_rs220_Seq",
  `sourmash RefSeq` = "SM_RefSeq_20250528", `sourmash GTDB 220` = "SM_gtdb-rs220-rep",
  `Kraken-Bracken RefSeq` = "KB45", `Kraken-Bracken GTDB` = "KB45_GTDB",
  `mOTUs 3` = "MOTUS3", `mOTUs 4` = "MOTUS4",
  `MetaPhlAn 2022` = "MPA_db2022", `MetaPhlAn 2023` = "MPA_db2023", `MetaPhlAn 2025` = "MPA_db2025")
method_colours <- setNames(methods$colour[match(benchmark_methods, methods$Database)], names(benchmark_methods))

# Kraken's index sits in /dev/shm, charged to the copy step, so a sample's step may not show it (GB)
index_gb <- c(KB = 85, KG = 502)

# Read inputs -----------------------------------------------------------------
samples <- read_tsv(
  file.path(dir, "samples.tsv"), show_col_types = FALSE,
  col_names = c("Dataset", "Sample", "p1", "p2", "u1", "u2", "Bytes")) %>%
  mutate(Row = row_number())

bases <- read_tsv(file.path(dir, "bases.tsv"), show_col_types = FALSE)

# Bases per read file, with its sample. Gbp per sample counts all 4 read files (sylph only reads the 2 paired ones).
sample_files <- samples %>%
  pivot_longer(c(p1, p2, u1, u2), values_to = "file") %>%
  left_join(bases, by = "file")

input <- summarise(
  sample_files, Gbp = sum(sum_len) / 1e9, 
  .by = c(Row, Dataset, Sample))

jobs <- read_tsv(
  file.path(dir, "jobs.tsv"), 
  show_col_types = FALSE, col_names = c("Label", "Job"),
  col_types = "cc") # every submission: a rerun of some rows (Kraken, OOM tasks) adds a job to its label

# sacct output (slurm logs job usage)
sacct <- read_delim(
  file.path(dir, "sacct.tsv"), delim = "|", show_col_types = FALSE,
  col_types = cols(.default = "c"))

# Unit conversions ------------------------------------------------------------

# [D-]HH:MM:SS or MM:SS.mmm to seconds
to_seconds <- function(x) {
  map_dbl(x, function(s) {
    days <- if (str_detect(s, "-")) as.numeric(str_extract(s, "^\\d+")) else 0
    parts <- as.numeric(str_split_1(str_remove(s, "^\\d+-"), ":"))
    days * 86400 + sum(parts * 60^rev(seq_along(parts) - 1))
  })
}

# 123K, 45M, 6G to GB
to_gb <- function(x) {
  value <- as.numeric(str_remove(x, "[KMGT]$"))
  unit <- str_extract(x, "[KMGT]$")
  value * c(K = 1 / 1024^2, M = 1 / 1024, G = 1, T = 1024)[unit]
}

# One row per sample x step ---------------------------------------------------

# Array tasks: 123_7.batch, 7 = row of samples.tsv
array_tasks <- sacct %>%
  filter(str_detect(JobID, "^\\d+_\\d+\\.batch$")) %>%
  tidyr::extract(JobID, c("Job", "Row"), "^(\\d+)_(\\d+)\\.", remove = FALSE) %>%
  inner_join(jobs, by = "Job")

# Kraken: srun steps named <label>_<row>, for example KG45_7
kraken_steps <- sacct %>%
  filter(str_detect(JobName, "^K[BG]\\d+_\\d+$"), str_remove(JobID, "\\..*") %in% jobs$Job) %>%
  tidyr::extract(JobName, c("Label", "Row"), "^(K[BG]\\d+)_(\\d+)$", remove = FALSE)

# A rerun row: its last job wins. That also drops the 2026-11-10 KB/KG steps that ran on 0 reads (COMPLETED anyway:
# nfs-ip34 paths not visible in the container), rerun on 2026-11-12.
use <- bind_rows(array_tasks, kraken_steps) %>%
  filter(State == "COMPLETED") %>%
  mutate(Row = as.integer(Row), Job_n = as.numeric(str_extract(JobID, "^\\d+"))) %>%
  slice_max(Job_n, n = 1, by = c(Label, Row)) %>%
  inner_join(input, by = "Row") %>%
  mutate(
    Tool = factor(tools[Label], levels = tools),
    Wall_h = as.numeric(ElapsedRaw) / 3600,
    CPU_h = to_seconds(TotalCPU) / 3600,
    Mem_GB = to_gb(MaxRSS),
    Index_GB = index_gb[str_sub(Label, 1, 2)],
    # add the index only if the step's memory can't already include it
    Mem_GB = if_else(!is.na(Index_GB) & Mem_GB < Index_GB, Mem_GB + Index_GB, Mem_GB)
  ) %>%
  select(Tool, Dataset, Sample, Gbp, CPUs = AllocCPUS, Wall_h, CPU_h, Mem_GB, Node = NodeList)

write_tsv(use, file.path(dir, "resource_use.tsv"))


# Cost model ------------------------------------------------------------------
# A run = what a user runs for one sample:
#   - sylph, sourmash: sketch + profile (the sketch is needed, and costs more than sylph's profiling);
#   - Kraken: Kraken2 + Bracken at one confidence. The 3 confidences cost the same, so the 3 runs per sample are pooled.
# Per run: CPU_h = a + b * Gbp. a = overhead per sample (e.g. loading the MetaPhlAn index), b = CPU-hours per Gbp.
#   Gamma GLM, identity link: the spread grows with the cost, so the deepest samples don't dominate the fit.
# Left out: Kraken's index copy to /dev/shm, done once per dataset run on a node (< 0.15 CPU-hours).

runs <- use %>%
  mutate(Method = str_remove(Tool, " 0\\.\\d+$"), Sketcher = word(Tool, 1), Key = paste(Dataset, Sample))

sketches <- runs %>%
  filter(str_detect(Tool, "sketch")) %>%
  select(Sketcher, Key, Sketch_h = CPU_h, Sketch_GB = Mem_GB)

runs <- runs %>%
  filter(!str_detect(Tool, "sketch")) %>%
  left_join(sketches, by = c("Sketcher", "Key")) %>%
  mutate(CPU_h = CPU_h + coalesce(Sketch_h, 0),
         Mem_GB = pmax(Mem_GB, Sketch_GB, na.rm = TRUE),
         Method = factor(Method, levels = names(method_colours))) %>%
  select(Method, Key, Dataset, Sample, Gbp, CPU_h, Mem_GB)

fit_cost <- function(d) {
  start <- pmax(coef(lm(CPU_h ~ Gbp, d)), c(1e-4, 1e-6))
  f <- tryCatch(suppressWarnings(glm(CPU_h ~ Gbp, family = Gamma("identity"), data = d, start = start)),
                error = \(e) NULL)
  if (!is.null(f) && coef(f)[1] >= 0) return(f)
  # a < 0 (Kraken-Bracken RefSeq: no measurable overhead): through the origin
  suppressWarnings(glm(CPU_h ~ 0 + Gbp, family = Gamma("identity"), data = d, start = start[2]))
}

by_method <- split(runs, runs$Method)
fits <- map(by_method, fit_cost)

# Wald 95 % intervals; a = 0 without interval for a fit through the origin
model <- imap_dfr(fits, \(f, m) {
  ci <- confint.default(f)
  has_a <- "(Intercept)" %in% rownames(ci)
  tibble(Method = m,
         a = if (has_a) coef(f)[["(Intercept)"]] else 0,
         a_lo = if (has_a) ci["(Intercept)", 1] else NA, a_hi = if (has_a) ci["(Intercept)", 2] else NA,
         b = coef(f)[["Gbp"]], b_lo = ci["Gbp", 1], b_hi = ci["Gbp", 2],
         n_samples = n_distinct(by_method[[m]]$Key),
         Mem_GB_median = median(by_method[[m]]$Mem_GB), Mem_GB_max = max(by_method[[m]]$Mem_GB))
}) %>%
  mutate(Method = factor(Method, levels = names(method_colours))) %>%
  arrange(Method)

model %>%
  mutate(across(where(is.numeric) & !n_samples, \(x) signif(x, 3))) %>%
  write_tsv(file.path(dir, "cost_model.tsv"))

# A: CPU-hours per sample, with the fit ------------------------------------------

number <- scales::label_number(drop0trailing = TRUE, big.mark = ",")

# Each method's line over the Gbp range it was measured on; ribbon = 95 % CI of the mean, built on the log scale
# (delta method: SE of log(fit) = SE / fit) so it stays positive on the log axes
fit_lines <- imap_dfr(fits, \(f, m) {
  g <- range(by_method[[m]]$Gbp)
  grid <- tibble(Gbp = 10^seq(log10(g[1]), log10(g[2]), length.out = 50))
  p <- predict(f, grid, se.fit = TRUE)
  mutate(grid, Method = m, CPU_h = p$fit,
         lo = p$fit * exp(-1.96 * p$se.fit / p$fit), hi = p$fit * exp(1.96 * p$se.fit / p$fit))
}) %>%
  mutate(Method = factor(Method, levels = names(method_colours)))

fit_plot <- ggplot(runs, aes(Gbp, CPU_h)) +
  geom_point(colour = "grey45", size = 0.7) +
  geom_ribbon(data = fit_lines, aes(ymin = lo, ymax = hi, fill = Method), alpha = 0.25) +
  geom_line(data = fit_lines, aes(colour = Method), linewidth = 0.6) +
  scale_x_log10(labels = number) +
  scale_y_log10(labels = number) +
  scale_colour_manual(values = method_colours, aesthetics = c("colour", "fill")) +
  labs(x = "Sample size (Gbp, post-QC)", y = "CPU-hours per sample") +
  theme_classic() +
  theme(strip.background = element_blank()); fit_plot

# B: peak memory per step ---------------------------------------------------------
# One point per sample and step. Kraken: the highest of the 3 confidences (they differ by < 1 GB), index included.
# Sketches get their own rows (in grey: each serves 2 databases).

steps <- c("sylph sketch", "sylph RefSeq", "sylph GTDB 220", "sourmash sketch", "sourmash RefSeq", "sourmash GTDB 220",
           "Kraken-Bracken RefSeq", "Kraken-Bracken GTDB", "mOTUs 3", "mOTUs 4", "MetaPhlAn 2022", "MetaPhlAn 2023", "MetaPhlAn 2025")
step_colours <- c(method_colours, `sylph sketch` = "grey50", `sourmash sketch` = "grey50")

memory <- use %>%
  mutate(Step = str_remove(Tool, " 0\\.\\d+$")) %>%
  summarise(Mem_GB = max(Mem_GB), .by = c(Step, Dataset, Sample)) %>%
  mutate(Step = factor(Step, levels = steps))

mem_plot <- ggplot(memory, aes(x = Mem_GB, y = Step, colour = Step)) +
  geom_point(position = position_jitter(width = 0, height = 0.2, seed = 1), size = 1, alpha = 0.5, shape = 21) +
  stat_summary(fun = median, geom = "point", shape = "|", size = 6, colour = "black") +
  scale_x_log10(labels = number, limits = c(NA, 1000), breaks = c(1, 10, 100, 1000)) +
  scale_y_discrete(limits = rev(steps)) +
  scale_colour_manual(values = step_colours, guide = "none") +
  labs(x = "Peak memory (GB)", y = NULL) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))

resources <- free(fit_plot, type = "label") / mem_plot +
  plot_layout(widths = c(3, 2), guides = "collect") +
  plot_annotation(tag_levels = "A"); resources

ggsave(file.path(dir, "benchmark_resources.pdf"), 
       plot = resources, bg = 'white', units = 'px',
       width = 2600, height = 2200)


memory %>% 
  group_by(Step) %>% 
  summarise(Median = median(Mem_GB),
            Max = max(Mem_GB))
  
  
  
  
  
  
  
  
  
  
  
  
  
