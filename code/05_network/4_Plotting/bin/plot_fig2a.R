#!/usr/bin/env Rscript
# plot_fig2a.R -- Fig. 2a: genes captured per network hierarchy.
#
# Source: originally a standalone script (remake_fig2A.R) that read a live Google Sheet
# (hardcoded to a collaborator's Google account) -- not reproducible outside that one
# account, and the sheet itself is just a hand-copy of the same numbers
# build_network_gene_counts.py now computes live. Converted 2026-09-27: reads that script's
# network_gene_counts.csv output directly, no live external auth. Plot logic unchanged.
#
# Replaces the hardcoded (captured, total) literals plot_network_overlap.R used to draw this
# panel with -- see that script's header for the real bug this closes (Depression/dogCD's
# total was coded assuming the wrong z_comb threshold).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]  # network_gene_counts.csv, from BUILD_NETWORK_GENE_COUNTS
plot_out    <- args[2]

counts <- read.csv(counts_file, check.names = FALSE)

df <- counts |>
  rename(
    net      = Network,
    type     = Type,
    n_tot = `genes in network`,
    n_hier   = `genes in hierarchy`,
    seed_tot = `seeds total`,
    dog_seed_tot   = `dog seeds total`,
    hum_seed_tot   = `human seeds total`,
    dog_seed_net  = `dog seeds in network`,
    hum_seed_net  = `human seeds in network`,
    seed_tot_net = `seeds in network`,
    dog_seed_hier = `dog seeds in hierarchy`,
    hum_seed_hier = `human seeds in hierarchy`
  ) |>
  # seeds in network but not carried into the hierarchy-filtered subnetwork -- these counts are
  # too small (0-4 genes) to show as a visually distinguishable bar segment at this plot's scale,
  # so a text annotation states the count directly instead, only when it's actually nonzero
  mutate(seed_hier_net = dog_seed_hier + hum_seed_hier,
         seed_not_hier = seed_tot_net - seed_hier_net,
         seed_label = ifelse(seed_not_hier > 0,
                              paste0(seed_tot_net,"/",seed_tot, "\nseeds (", seed_not_hier, " not in hier.)"),
                              paste0(seed_tot_net,"/",seed_tot, "\nseeds"))) |>
  # top-to-bottom order on the y axis; ggplot draws the first level at the bottom, hence rev()
  mutate(net = factor(net, levels = rev(c("dogCD","OCD","OCD/dogCD", "Depression","Depression/dogCD", "Schizophrenia","Schizophrenia/dogCD"))))

p <- ggplot(df,aes(y=net,fill=species))
p <- p + geom_col(aes(x=n_tot,color=species),alpha=0.3, linewidth=0.2, width=2/3) # genes in network, faint bar with solid outline -- drawn
                                                                                   # first/underneath so it doesn't wash out the seed colours below
p <- p + geom_col(aes(x=n_hier), width=2/3)                                  # genes in hierarchy, drawn on top
# seed genes in network, darkens the start of the bar -- human-seed width drawn first (always the
# wider or equal segment, since dog_seed_net + hum_seed_net = seed_tot_net), then dog-seed width
# drawn on top of it; since dog_seed_net is always <= seed_tot_net this leaves [0, dog_seed_net]
# showing the dog-seed colour and (dog_seed_net, seed_tot_net] showing the human-seed colour --
# i.e. within the cross-species rows' seed segment, which species each seed gene came from.
# Drawn last/on top so neither colour gets tinted by the faint n_tot overlay above. Darker shades
# of the species colours, not the same hues used for n_hier/n_tot -- single-species rows (all-dog
# or all-human seeds) would otherwise blend invisibly into their own same-coloured base bar.
p <- p + geom_col(aes(x=seed_tot_net), fill="#003E7A", alpha=0.85, width=2/3)
p <- p + geom_col(aes(x=dog_seed_net), fill="#7A1414", alpha=0.85, width=2/3)
p <- p + geom_text(aes(x=seed_tot_net, label=seed_label,lineheight=0.9),      # seed genes in network, just past the black bar --
                   hjust=0, nudge_x=5, colour="white", size=1.25)             # states how many of those aren't in the hierarchy, when nonzero
p <- p + geom_text(aes(x=n_tot, label=n_tot), hjust=-0.15, colour="grey30", size=1.5)
p <- p + geom_text(aes(x=n_hier, label=paste0(n_hier, " (", round(100*n_hier/n_tot), "%)")),
                   hjust=1.1, colour="white", size=1.5)
p <- p + scale_x_continuous(breaks = c(0, 200, 400, 600),
                            expand = expansion(mult = c(0, 0.07)))          # room for the outside labels
p <- p + scale_fill_manual(values = c(human = "#0074E6", dog = "#D54040", cross = "#D04900"))
p <- p + scale_colour_manual(values = c(human = "#0074E6", dog = "#D54040", cross = "#D04900"), guide = "none")
p <- p + labs(x = "number of genes", y = NULL)
p <- p + theme_cowplot(font_size=6, line_size=0.25)                        # axis lines and ticks half the default 0.5
p <- p + theme(axis.ticks.y = element_blank(),                               # ticks on x axis only
               axis.ticks.length.x = unit(3, "pt"),
               axis.line.y=element_blank(),
               legend.position = "none")

ggsave(plot=p,filename=plot_out,width=100,height=50,units="mm")
