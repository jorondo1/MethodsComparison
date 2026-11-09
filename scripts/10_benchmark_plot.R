# Resource benchmark: CPU-hours and peak memory per tool x database x dataset.
# Input: the output folder of scripts/benchmark/benchmark.sh (samples.tsv, bases.tsv, jobs.tsv, sacct.tsv)
# Output in the same folder: resource_use.tsv, summary.tsv, benchmark_ridges.pdf, benchmark_scaling.pdf
# Usage: Rscript scripts/10_benchmark_plot.R [folder]

library(tidyverse)
library(ggridges)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "Out/benchmark"

# Top to bottom in the ridge plots
tools <- c(
  sylph_sketch = "sylph sketch", SYref = "sylph RefSeq", SY220 = "sylph GTDB 220",
  sourmash_sketch = "sourmash sketch", SMref = "sourmash RefSeq", SM220 = "sourmash GTDB 220",
  KB10 = "Kraken RefSeq 0.10", KB45 = "Kraken RefSeq 0.45", KB90 = "Kraken RefSeq 0.90",
  KG10 = "Kraken GTDB 0.10", KG45 = "Kraken GTDB 0.45", KG90 = "Kraken GTDB 0.90",
  MOTU3 = "mOTUs 3", MOTU4 = "mOTUs 4",
  MPA22 = "MetaPhlAn 2022", MPA23 = "MetaPhlAn 2023", MPA25 = "MetaPhlAn 2025"
)

# Kraken's index sits in /dev/shm, charged to the copy step, so a sample's step may not show it (GB)
index_gb <- c(KB = 85, KG = 502)

# Read inputs -----------------------------------------------------------------

samples <- read_tsv(file.path(dir, "samples.tsv"), show_col_types = FALSE,
                    col_names = c("Dataset", "Sample", "p1", "p2", "u1", "u2", "Bytes")) %>%
  mutate(Row = row_number())

bases <- read_tsv(file.path(dir, "bases.tsv"), show_col_types = FALSE)

# Gbp per sample, all 4 read files. Sylph only reads the 2 paired files.
input <- samples %>%
  pivot_longer(c(p1, p2, u1, u2), values_to = "file") %>%
  left_join(bases, by = "file") %>%
  summarise(Gbp = sum(sum_len) / 1e9, .by = c(Row, Dataset, Sample))

jobs <- read_tsv(file.path(dir, "jobs.tsv"), show_col_types = FALSE, col_names = c("Label", "Job"),
                 col_types = "cc") %>%
  slice_tail(n = 1, by = Label) # the last submission of each label

sacct <- read_delim(file.path(dir, "sacct.tsv"), delim = "|", show_col_types = FALSE,
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

# One row per sample x tool ---------------------------------------------------

# Array tasks: 123_7.batch, 7 = row of samples.tsv
array_tasks <- sacct %>%
  filter(str_detect(JobID, "^\\d+_\\d+\\.batch$")) %>%
  extract(JobID, c("Job", "Row"), "^(\\d+)_(\\d+)\\.", remove = FALSE) %>%
  inner_join(jobs, by = "Job")

# Kraken: srun steps named <label>_<row>, for example KG45_7
kraken_steps <- sacct %>%
  filter(str_detect(JobName, "^K[BG]\\d+_\\d+$"), str_remove(JobID, "\\..*") %in% jobs$Job) %>%
  extract(JobName, c("Label", "Row"), "^(K[BG]\\d+)_(\\d+)$", remove = FALSE)

use <- bind_rows(array_tasks, kraken_steps) %>%
  filter(State == "COMPLETED") %>%
  mutate(Row = as.integer(Row)) %>%
  inner_join(input, by = "Row") %>%
  mutate(
    Tool = factor(tools[Label], levels = rev(tools)),
    Wall_h = as.numeric(ElapsedRaw) / 3600,
    CPU_h = to_seconds(TotalCPU) / 3600,
    CPU_h_per_Gbp = CPU_h / Gbp,
    Mem_GB = to_gb(MaxRSS),
    Index_GB = index_gb[str_sub(Label, 1, 2)],
    # add the index only if the step's memory can't already include it
    Mem_GB = if_else(!is.na(Index_GB) & Mem_GB < Index_GB, Mem_GB + Index_GB, Mem_GB)
  ) %>%
  select(Tool, Dataset, Sample, Gbp, CPUs = AllocCPUS, Wall_h, CPU_h, CPU_h_per_Gbp, Mem_GB, Node = NodeList)

write_tsv(use, file.path(dir, "resource_use.tsv"))

use %>%
  summarise(across(c(Gbp, Wall_h, CPU_h, CPU_h_per_Gbp, Mem_GB), list(median = median, min = min, max = max)),
            n = n(), .by = c(Tool, Dataset)) %>%
  arrange(Tool, Dataset) %>%
  write_tsv(file.path(dir, "summary.tsv"))

# Plots -----------------------------------------------------------------------

ridge <- function(value, x_label) {
  ggplot(use, aes(x = {{ value }}, y = Tool, fill = Dataset, colour = Dataset)) +
    geom_density_ridges(alpha = 0.3, jittered_points = TRUE, point_size = 0.8) +
    scale_x_log10() +
    labs(x = x_label, y = NULL) +
    theme_classic()
}

# Raw cost, cost per Gbp, memory
ridges <- ridge(CPU_h, "CPU-hours per sample") +
  (ridge(CPU_h_per_Gbp, "CPU-hours per Gbp") + theme(axis.text.y = element_blank())) +
  (ridge(Mem_GB, "Peak memory (GB)") + theme(axis.text.y = element_blank())) +
  plot_layout(guides = "collect")
ggsave(file.path(dir, "benchmark_ridges.pdf"), ridges, width = 12, height = 7)

# Fixed cost against scaling with depth
scaling <- ggplot(use, aes(x = Gbp, y = CPU_h)) +
  geom_point(aes(colour = Dataset), size = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ fct_rev(Tool)) +
  labs(x = "Input (Gbp)", y = "CPU-hours") +
  theme_classic()
ggsave(file.path(dir, "benchmark_scaling.pdf"), scaling, width = 10, height = 8)
