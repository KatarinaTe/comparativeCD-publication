#!/usr/bin/env Rscript
# emission_state_plot.R — chromatin state emission-probability / genome-annotation-enrichment plot.
# Source: 06_cCRE/emission_state_UU.R — literal reproduction, only the setwd()/hardcoded
# filenames replaced with CLI args. Confirmed as the real historical generator of the deposited
# emission_states_plot.pdf/.png (ggsave's width=16, height=7, dpi=300 exactly matches their pixel
# dimensions, 4800x2100 — a duplicate Python implementation, emission_state_UU.ipynb, produces a
# visually similar but differently-sized figure and was not used here for that reason).

args <- commandArgs(trailingOnly = TRUE)
probs_file <- args[1]
annot_file <- args[2]
pdf_out    <- args[3]
png_out    <- args[4]

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# ── Load data ──────────────────────────────────────────────────────────────────
probs <- read.table(probs_file, header=TRUE, sep='\t')
annot <- read.table(annot_file, header=TRUE, sep='\t', check.names=FALSE)

# Keep presence=1 rows only
probs1 <- probs %>% filter(presence == 1)
marks_order <- c('ATAC', 'H3K27ac', 'H3K4me3', 'H3K27me3')
probs1$mark <- factor(probs1$mark, levels=marks_order)

# ── State annotations ──────────────────────────────────────────────────────────
state_labels <- data.frame(
  state = 1:9,
  label = c('Promoter', 'Enhancer', 'Open chromatin', 'Promoter',
            'Promoter', 'Promoter', 'No signal', 'Enhancer', 'Repressed')
)
state_colors <- c(
  'Promoter'       = '#F4A7B9',
  'Enhancer'       = '#4A90D9',
  'Open chromatin' = '#C8860A',
  'No signal'      = '#FFFFFF',
  'Repressed'      = '#5A8A3C'
)

# ── Panel 1: Emission probabilities ───────────────────────────────────────────
p1 <- ggplot(probs1, aes(x=mark, y=factor(state, levels=rev(1:9)), fill=prob)) +
  geom_tile(color='white', linewidth=1) +
  scale_fill_gradient(
    low='white', high='#1A3A8F',
    limits=c(0,1), name='Emission\nprobability',
    breaks=c(0, 0.5, 1),
    labels=c('0', '0.5', '1'),
    guide=guide_colorbar(
      barwidth       = unit(3, 'cm'),
      barheight      = unit(0.4, 'cm'),
      title.position = 'top',
      title.hjust    = 0.5,
      label.theme    = element_text(size=8),
      ticks          = FALSE
    )
  )+
  scale_x_discrete(position='bottom') +
  scale_y_discrete(name='State (Emission order)') +
  labs(x='Mark') +
  theme_minimal(base_size=11) +
  theme(
    axis.text.x      = element_text(angle=45, hjust=1, size=10),
    axis.text.y      = element_text(size=10),
    panel.grid       = element_blank(),
    legend.position  = 'bottom',
    legend.key.width = unit(0.5, 'cm'),
    plot.margin      = margin(5, 2, 5, 5)
  )

# ── Panel 2a: Genome % only ────────────────────────────────────────────────────
genome_df <- annot %>%
  rename(state = `State (Emission order)`) %>%
  select(state, `Genome %`)

p2a <- ggplot(genome_df, aes(x='Genome %', y=factor(state, levels=rev(1:9)), fill=`Genome %`)) +
  geom_tile(color='white', linewidth=1) +
  scale_fill_gradient(low='#FFF5EB', high='#8C2D04',
                      name='Genome %') +
  scale_x_discrete(position='bottom') +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(
    axis.text.x      = element_text(angle=45, hjust=1, size=9),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = 'bottom',
    legend.key.width = unit(0.3, 'cm'),
    plot.margin      = margin(5, 1, 5, 1)
  )

# ── Panel 2b: Annotation enrichment ───────────────────────────────────────────
annot_long <- annot %>%
  rename(state = `State (Emission order)`) %>%
  select(-`Genome %`) %>%
  pivot_longer(-state, names_to='annotation', values_to='enrichment')

annot_long$annotation <- gsub('\\.cf4\\.bed\\.gz', '', annot_long$annotation)
annot_order <- c('RefSeqTSS', 'RefSeqTSS2kb', 'RefSeqExon', 'RefSeqTES', 'RefSeqGene')
annot_long$annotation <- factor(annot_long$annotation, levels=annot_order)

p2b <- ggplot(annot_long, aes(x=annotation, y=factor(state, levels=rev(1:9)), fill=enrichment)) +
  geom_tile(color='white', linewidth=1) +
  scale_fill_gradientn(
    colours = c('#F7FBFF', '#C6DBEF', '#6BAED6', '#2171B5', '#08306B'),
    values  = scales::rescale(c(0, 5, 10, 20, 41)),
    limits  = c(0, 41),
    name    = 'Fold\nenrichment'
  )+
  labs(x='Genome annotation', y=NULL) +
  theme_minimal(base_size=11) +
  theme(
    axis.text.x      = element_text(angle=45, hjust=1, size=9),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = 'bottom',
    legend.key.width = unit(0.5, 'cm'),
    plot.margin      = margin(5, 2, 5, 1)
  )

# ── Panel 3: State label bars ──────────────────────────────────────────────────
state_label_df <- state_labels %>%
  mutate(state_rev = factor(state, levels=rev(1:9)))

p3 <- ggplot(state_label_df, aes(x=1, y=state_rev, fill=label)) +
  geom_tile(color='white', linewidth=1, width=0.9, height=0.9) +
  geom_text(aes(label=label,
                color=ifelse(label %in% c('Enhancer','Repressed','Open chromatin'),
                             'white', 'black')),
            size=3.2, fontface='bold') +
  scale_fill_manual(values=state_colors, guide='none') +
  scale_color_identity() +
  scale_x_continuous(expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 2)
  )

# ── Combine panels ─────────────────────────────────────────────────────────────
combined <- p1 + p2a + p2b + p3 +
  plot_layout(widths=c(2, 0.6, 3, 1)) +
  plot_annotation(
    title='Chromatin state emission probabilities and genomic enrichment',
    theme=theme(plot.title=element_text(size=12, face='bold', hjust=0.5))
  )

# ── Save ───────────────────────────────────────────────────────────────────────
ggsave(pdf_out, combined, width=16, height=7, device='pdf')
ggsave(png_out, combined, width=16, height=7, dpi=300)
