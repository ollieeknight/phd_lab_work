library(ggtext)
library(ggsci)
library(ggrepel)
library(gt)
library(patchwork)
library(pals)
library(ggthemes)
library(ggpubr)
library(ggridges)
library(BuenColors)
library(ggstatsplot)
library(rcartocolor)
library(Nebulosa)
library(broom)
library(RColorBrewer)

#rcartocolor::display_carto_all(type = 'diverging', colorblind_friendly = TRUE)

# Define color palettes and scales
cols <- c()
cols$origin <- c('#C4140A', '#4DABAC', '#238B45')
cols$celltype <- c('#2C5D9D', '#BFA387', '#2C9D4C', '#9D2E2C')
cols$synovial <- c(pal_jco('default')(6), cols$celltype[2:4])
cols$peripheral <- c(tol()[1:4], '#BFA387', '#2C9D4C')
cols$total <- c(cols$peripheral[1:4], cols$synovial[1:7], cols$peripheral[5], cols$synovial[8], cols$peripheral[6])
cols$donor <- tol()
cols$phase <- pal_uchicago('default')(3)
cols$disease <- c('cornsilk4', 'cornflowerblue', 'cyan4', 'coral3')
cols$disease2 <- c('#DC571D', '#008080', '#F1D96E')
cols$celltype2 <- pal_npg(palette = 'nrc')(5)
cols$ILCs <- c('darkblue', 'deepskyblue4', 'deepskyblue2', '#D55882', '#492351', '#EE7744',
               "#BFA387", "#2C9D4C", "#9D2E2C")

scale <- c()
scale$rna <- colorRampPalette(c('#f3d8e1', '#8da0cf', '#1a69d0', '#353979', '#300b36'))(100)
scale$adt <- colorRampPalette(c("#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#125B49"))(100)
scale$adt2 <- colorRampPalette(c('#DEDDD3', '#ACB9AF', '#699B84', '#305949', '#012015'))(100)
scale$dotplot <- colorRampPalette(c('#FFF5E0', '#FFDE73', '#F9BC4F', '#EE7744', '#2C7C79', '#4C4D78'))(100)
scale$score1 <- colorRampPalette(c('#EBE2BC', '#FFBE5F', '#FB7D2D', '#DB242E', '#003147'))(100)
scale$score2 <- colorRampPalette(c('#F4EACB', '#EDA64E', '#E6403C', '#9B242E', '#55281E'))(100)
scale$score3 <- colorRampPalette(c('#FDFBF8', '#F1E6E0', '#FF8F97', '#E44458', '#2E2364', '#150640'))(100)
scale$score4 <- colorRampPalette(c('#F6E5D9', '#FFBF92', '#FF808A', '#D55882', '#7F376D', '#492351'))(100)
scale$heatmap <- ocean.balance(100)

themes <- c()

themes$UMAP_axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, 'npc'),
  trunc_upper = unit(0.15, 'npc')
)

themes$UMAP_theme <- list(
  guides(
    x = themes$UMAP_axis, 
    y = themes$UMAP_axis
  ),
  labs(x = 'UMAP_1', y = 'UMAP_2'),
  theme(
    axis.title.x = element_markdown(size = 8, hjust = 0),
    axis.title.y = element_markdown(size = 8, hjust = 0),
    axis.text.x = element_markdown(size = 0),
    axis.text.y = element_markdown(size = 0),
    axis.ticks = element_blank(),
    legend.position = 'right',
    axis.line = element_line(arrow = arrow(length = unit(3.0, 'npc'), type = 'closed')),
    plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc')
  ),
  coord_fixed())

themes$feature_plot_theme <- list(
  guides(
    x = themes$UMAP_axis, 
    y = themes$UMAP_axis
  ),
  labs(x = 'UMAP_1', y = 'UMAP_2'),
  theme(
    axis.title.x = element_markdown(size = 8, hjust = 0),
    axis.title.y = element_markdown(size = 8, hjust = 0),
    axis.text.x = element_markdown(size = 0),
    axis.text.y = element_markdown(size = 0),
    axis.ticks = element_blank(),
    legend.position = 'right',
    axis.line = element_line(arrow = arrow(length = unit(3.0, 'npc'), type = 'closed')),
    plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc')
  ),
  coord_fixed())

themes$violin_split_theme <- list(
  theme(strip.text.x = element_text(size = 16, face = 'italic'),
        plot.title = element_markdown(size = 16, colour = 'black', face = 'italic'),
        axis.text.x = element_markdown(size = 12, colour = 'black'),
        axis.title.y = element_markdown(size = 14, colour = 'black'),
        axis.text.y = element_markdown(size = 12, colour = 'black'),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc'),
        strip.background.x = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(size = 10)
        )
)

themes$ridge <- list(
  theme(strip.text.x = element_text(size = 16, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(size = 10, colour = 'black', angle = 45, vjust = 0.6),
        axis.text.y = element_text(size = 10, colour = 'black'),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc'),
        strip.background.x = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(size = 10)
  )
)

themes$bar_split_theme <- list(
  theme(strip.text.x = element_text(size = 16, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(size = 10, colour = 'black'),
        axis.text.y = element_text(size = 10, colour = 'black'),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc'),
        strip.background.x = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(size = 10)
  )
)

themes$bar_split_theme_gene <- list(
  theme(strip.text.x = element_text(size = 16, face = 'italic'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(size = 10, colour = 'black'),
        axis.text.y = element_text(size = 10, colour = 'black'),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc'),
        strip.background.x = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(size = 10)
  )
)

themes$dotplot <- list(
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = 'italic', colour = 'black'), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, linetype = 'dashed',
                                        colour = 'grey'),
        axis.line = element_line(colour = 'black'),
        plot.margin = margin(0, 0, 0, 0, 'cm')),
    guides(size = guide_legend(order = 1, title = 'Percent\nexpressed'),
           color = guide_colorbar(order = 2, title = 'Average\nexpression', 
                                  label.theme = element_text(size = 10))),
    theme(plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc'))
  
)

themes$featureplot <- list(
  guides(x = themes$UMAP_axis, y = themes$UMAP_axis),
    labs(title = NULL, x = 'UMAP_1', y = 'UMAP_2'),
    coord_fixed(), theme(axis.title.x = element_text(size = 9, hjust = 0),
                          axis.title.y = element_text(size = 9, hjust = 0),
                          axis.text.x = element_text(size = 0),
                          axis.text.y = element_text(size = 0),
                          axis.ticks = element_blank(),
                          axis.line = element_line(arrow = arrow(length = unit(2.5, 'npc'), type = 'closed')),
                          plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'npc'),
                          legend.title = element_text(size = 10))
  
)
