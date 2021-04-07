## Code for: 
## White TE, Umbers KDL (2021) Meta-analytic evidence for quantitative honesty 
## in aposematic signals. Proc. R. Soc. B.

## Reset
  rm(list = ls())

## Libraries
  library(tidyverse)
  library(metafor)
  library(PRISMAstatement)
  library(ape)
  library(rotl)
  library(gt)
  library(multcomp)
  library(cowplot)
  library(magick)
  library(ggpubr)
  library(metaAidR)
  
## Functions
  # MLMR builder with common parameters, for simplicity
  mlmr <- function(dat, variable){
    rma.mv(yi, 
           vi,
           mods = ~ variable - 1,
           random = list(~ 1 | obs, ~1 | study),
           method = "REML", 
           data = dat)
  }
  
## Data
  coldat <- read.csv('../data/data.csv')
  
## Processing
  # Need n > 4 to estimate sampling variance
  coldat <- filter(coldat, n > 4)
  
  # Take absolute values when 'hue' is the colour measure as it's directionless
  coldat2 <- coldat  # Save untransformed data for sensitivity analysis 
  coldat[which(coldat$col_var == 'hue'), ]['r'] <- abs(coldat[which(coldat$col_var == 'hue'), ]['r'])
  
  # Relevel factors for convenience
  coldat$col_var <- factor(coldat$col_var, levels = c('hue', 'saturation', 'brightness', 'area', 'internal_contrast', 'external_contrast'))
  coldat$study_type <- factor(coldat$study_type, levels = c('among_individuals', 'among_populations', 'among_species'))
  
## Check if phylo guff is necessary
  # Replate some synonymous species names
  coldat[which(coldat$species == 'Epipedobates_bilinguis'), ]$species <- "Ameerega_bilinguis"
  coldat[which(coldat$species == 'Epipedobates_parvulus'), ]$species <- "Ameerega_parvula"
  coldat[which(coldat$species == 'Epipedobates_hahneli'), ]$species <- "Ameerega_hahneli"
  coldat[which(coldat$species == 'Parasemia_plantaginis'), ]$species <- "Arctia_plantaginis"
  
  # Subset to include only species-level data
  coldat_spec <- filter(coldat, study_type == 'among_individuals')
  
  # Estimate tree, compute branch lengths and correlation matrix
  taxa <- tnrs_match_names(names = gsub("_", " ", unique(coldat_spec$species)), do_approximate_matching = FALSE)
  # tree <- tol_induced_subtree(ott_ids = ott_id(taxa), label_format = 'name')
  # tree <- compute.brlen(tree, power = 0.5)  # Estimate branch lengths
  tree <- read.tree("../data/tree.tre")
  phylo <- vcv(tree, corr = TRUE) # Correlation matrix
  
  # Check for name mis-matches
  setdiff(row.names(phylo), unique(coldat$species))

## Modelling
  
  # Calculate effects
  coldat <- escalc("ZCOR", ri = r, ni = n, data = coldat, append = TRUE)
  coldat_spec <- escalc("ZCOR", ri = r, ni = n, data = coldat_spec, append = TRUE)

  # Create a variance-covariance matrix to account for reapeated measurements from groups
  control_vcv <- make_VCV_matrix(coldat, V = 'vi', cluster = 'group', rho = 0.5)
  control_vcv_spec <- make_VCV_matrix(coldat_spec, V = 'vi', cluster = 'group', rho = 0.5)
  
  # Phylogenetic null model
  all_null_phylo <- rma.mv(yi, 
                           control_vcv_spec,
                           random = list(~ 1 | obs, ~1 | species),
                           R = list(species = phylo),
                           method = "REML",
                           slab = paste(author, year, sep=", "),
                           data = coldat_spec)
  summary(all_null_phylo)
  I2(all_null_phylo, coldat_spec$vi, obs = 'obs')   
  
  # Study null model
  all_null_study <- rma.mv(yi, 
                            control_vcv_spec,
                            random = list(~ 1 | obs, ~1 | study),
                            method = "REML",
                            slab = paste(author, year, sep=", "),
                            data = coldat_spec)
  summary(all_null_study)
  I2(all_null_study, coldat_spec$vi, obs = 'obs')

  # Is phylo necessary?
  AIC(all_null_study, all_null_phylo)
  abs(AIC(all_null_study) - AIC(all_null_phylo))
    
  # Full null model
  models <- list()  # A list for tidily storing all the models
  I_sq <- list()  # A list for tidily storing I2's
  
  models$all_null <- rma.mv(yi, 
                            control_vcv,
                            random = list(~ 1 | obs, ~1 | study),
                            method = "REML",
                            slab = paste(author, year, sep=", "),
                            data = coldat)
  summary(models$all_null)
  I_sq$all_null <- I2(models$all_null, coldat$vi, obs = 'obs')   
  
  # Check publication bias
  # Eggers test
  m_eggers <- regtest(rma(yi = residuals(models$all_null), sei = 1/sqrt(1/models$all_null$vi)), model = "lm")
  
  # Funnel plot
  png('../figs/fig_funnel.png', width = 21, height = 21, units = 'cm', res = 300)
    par(mar = c(6, 5, 4, 2))
    funnel(rma(yi = residuals(models$all_null), sei = 1/sqrt(1/models$all_null$vi)), 
           xlab = 'Effect size (Fisher\'s Z)',
           ylab = 'Standard error',
           xlim = c(-3, 3),
           digits = 2,
           back = 'whitesmoke')
  dev.off()  
  
## MLMR: Moderators
  
  # Study type (intraspecific, interpopulation, interspecific)
  models$m_type <- mlmr(coldat, coldat$study_type)
  summary(models$m_type)
  I_sq$m_type <- I2(models$m_type, coldat$vi, obs = 'obs')
  summary(glht(models$m_type, linfct = contrMat(table(coldat$study_type), type = "Tukey"), test = adjusted(type = 'bonferroni')))
  
  # Colour measure (hue, saturation, brightness, internal contrast, external contrast)
  models$m_col <- mlmr(coldat, coldat$col_var)
  summary(models$m_col)
  I_sq$m_col <- I2(models$m_col, coldat$vi, obs = 'obs')
  summary(glht(models$m_col, linfct = contrMat(table(coldat$col_var), type = "Tukey"), test = adjusted(type = 'bonferroni')))
  
  # Defence measure (chemical/physiological measure versus defence assay)
  models$m_tox <- mlmr(coldat, coldat$tox_measure)
  summary(models$m_tox)
  I_sq$m_tox <- I2(models$m_tox, coldat$vi, obs = 'obs')
  summary(glht(models$m_tox, linfct = contrMat(table(coldat$tox_measure), type = "Tukey"), test = adjusted(type = 'bonferroni')))
  
  # Taxa (intraspecific, interpopulation, interspecific)
  models$m_class <- mlmr(coldat, coldat$class)
  summary(models$m_class)
  I_sq$m_class <- I2(models$m_class, coldat$vi, obs = 'obs')
  summary(glht(models$m_class, linfct = contrMat(table(coldat$class), type = "Tukey"), test = adjusted(type = 'bonferroni')))
  
  # Tidy summary effects for plotting
  effects <- data.frame(name = str_replace_all(str_remove(do.call(c, lapply(models, function(x) row.names(x$b))), 'variable'), fixed("_"), " ") ,
                        b = do.call(c, lapply(models, function(x) x$b)),
                        lci = do.call(c, lapply(models, function(x) x$ci.lb)),
                        uci = do.call(c, lapply(models, function(x) x$ci.ub)))
  effects$name[1] <- 'overall'
  effects$n <- c(nrow(coldat),
                 length(which(coldat$study_type == 'among_individuals')),
                 length(which(coldat$study_type == 'among_populations')),
                 length(which(coldat$study_type == 'among_species')),
                 length(which(coldat$col_var == 'hue')),
                 length(which(coldat$col_var == 'saturation')),
                 length(which(coldat$col_var == 'brightness')),
                 length(which(coldat$col_var == 'area')),
                 length(which(coldat$col_var == 'internal_contrast')),
                 length(which(coldat$col_var == 'external_contrast')),
                 length(which(coldat$tox_measure == 'physiological')),
                 length(which(coldat$tox_measure == 'bioassay')),
                 length(which(coldat$class == 'amphibia')),
                 length(which(coldat$class == 'gastropoda')),
                 length(which(coldat$class == 'insecta')))
  
  # Back-transform to r
  effects$r <- transf.ztor(effects$b)
  effects$r.lci <- transf.ztor(effects$lci)
  effects$r.uci <- transf.ztor(effects$uci)

## Effects table
  effects_tab <- 
    data.frame(
    effects$name,
    effects$r,
    effects$r.lci,
    effects$r.uci,
    effects$n
  )
  effects_tab$I2 <- c(round(I_sq$all_null[2,1], 2),
                      round(I_sq$m_type[2,1], 2),
                      ' ',
                      ' ',
                      round(I_sq$m_col[2,1], 2),
                      ' ',
                      ' ',
                      ' ',
                      ' ',
                      ' ',
                      round(I_sq$m_tox[2,1], 2),
                      ' ',
                      round(I_sq$m_class[2,1], 2),
                      ' ',
                      ' ')
  effects_tab$effects.name[1] <- "intercept-only"
  e_table <-
    effects_tab %>%
    gt() %>%
    fmt_number(
      c(2:4),
      decimals = 2
    ) %>%
    cols_label(
      effects.name = " ",
      effects.r = "r",
      effects.r.lci = "lower CI",
      effects.r.uci = "upper CI",
      effects.n = "n",
      I2 = "I2"
    ) %>%
    tab_row_group(
      group = "Taxon",
      rows = 13:15
    ) %>% 
    tab_row_group(
      group = "Defence measure",
      rows = 11:12
    ) %>%
    tab_row_group(
      group = "Signal component",
      rows = 5:10
    ) %>%
    tab_row_group(
      group = "Scale",
      rows = 2:4
    ) %>%
    tab_row_group(
      group = "Overall",
      rows = 1
    )
  e_table

## Plots

# Main forest
   
  effects$name <- factor(effects$name, levels = c('insecta', 'amphibia', 'gastropoda',
                                                  'physiological', 'bioassay', 
                                                  'external contrast', 'internal contrast', 'area', 'brightness', 'saturation', 'hue',
                                                  'among species', 'among populations', 'among individuals', 
                                                  'overall'))
  effects$label <- c(' ',
                     'Scale', ' ', ' ', 
                     'Signal component', ' ', ' ', ' ', ' ', ' ',
                     'Defence measure', ' ', 
                     ' ', 'Taxon', ' ')
  
ggplot(data = effects, aes(x = name, y = r, ymin = r.lci, ymax = r.uci, label = label)) +
            geom_pointrange() +
            geom_hline(yintercept = 0, linetype = 2) +
            geom_errorbar(aes(ymin = r.lci, ymax = r.uci), width = 0.5, cex = 0.7) + 
            geom_text(y = -1.8,
                      vjust = -1.5,
                      hjust = 0,
                      size = 3,
                      aes(fontface = 'bold')) +
            xlab(' ') + 
            ylab("Pearson's r (95% CI)") +
            ylim(-1, 1) +
            theme_bw() +
            theme(plot.title = element_text(size = 16, face = "bold"),
                  axis.text.x = element_text(face = "bold"),
                  axis.text.y = element_text(size = 7.5),
                  axis.title = element_text(size = 12, face = "bold"),
                  strip.text.y = element_text(hjust = 0, vjust = 1, angle = 180, face = "bold"),
                  plot.margin = margin(1, 3, 1, 2, "cm"),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()) +
            coord_flip(clip = 'off')
  
  ggsave('../figs/fig_mod_forest.tiff', units = 'cm', width = 18,  height = 18)

  
# Effects summary
  
  # Guide images
  ind <- image_read('../figs/pics/ind.png')
  pop <- image_read('../figs/pics/pop.png')
  spec <- image_read('../figs/pics/spec.png')
  frog <- image_read('../figs/pics/frog.png')
  bee <- image_read('../figs/pics/bee.png')
  snail <- image_read('../figs/pics/snail.png')
  needle <- image_read('../figs/pics/needle.png')
  phys <- image_read('../figs/pics/phys.png')
  internal <- image_read('../figs/pics/internal.png')
  external <- image_read('../figs/pics/external.png')
  saturation <- image_read('../figs/pics/saturation.png')
  brightness <- image_read('../figs/pics/brightness.png')
  hue <- image_read('../figs/pics/hue.png')
  area <- image_read('../figs/pics/area.png')
  
  effects_scale <- effects[2:4, ]
  effects_scale$name <- factor(effects_scale$name, levels = c('among individuals', 'among populations', 'among species'))
  
  gg_scale <- ggplot(effects_scale, aes(x = name, y = n)) +
                geom_bar(stat = 'identity') +
                xlab('') +
                ylab('') +
                ylim(0, 100) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1.1),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank()) +
  draw_image(ind, x = 0.5, y = 89, scale = 29) +
  draw_image(pop, x = 1.5, y = 89, scale = 29) +
  draw_image(spec, x = 2.5, y = 89, scale = 29)
  
  effects_taxa <- effects[13:15, ]
  effects_taxa$name <- factor(effects_taxa$name, levels = c('insecta', 'amphibia', 'gastropoda'))
  
  gg_taxa <- ggplot(effects_taxa, aes(x = name, y = n)) +
              geom_bar(stat = 'identity') +
              xlab('') +
              ylab(' ') +
              ylim(0, 100) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1.1),
                    axis.text.y = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank()) +
              draw_image(bee, x = 0.5, y = 98, scale = 10) +
              draw_image(frog, x = 1.5, y = 98, scale = 12) + 
              draw_image(snail, x = 2.5, y = 98, scale = 10)
  
  effects_col <- effects[5:10, ]
  effects_col$name <- factor(effects_col$name, levels = c('internal contrast', 'external contrast', 'saturation', 'brightness', 'hue', 'area'))
  
  gg_col <- ggplot(effects_col, aes(x = name, y = n)) +
              geom_bar(stat = 'identity') +
              xlab('') +
              ylab('') +
              ylim(0, 100) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1.1),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank()) +
  draw_image(internal, x = 0.5, y = 95, scale = 16) +
  draw_image(external, x = 1.5, y = 95, scale = 16) +
  draw_image(saturation, x = 2.5, y = 95, scale = 16) +
  draw_image(brightness, x = 3.5, y = 95, scale = 16) +
  draw_image(hue, x = 4.5, y = 95, scale = 16) +
  draw_image(area, x = 5.5, y = 95, scale = 16)
  
  effects_tox <- effects[11:12, ]
  effects_tox$name <- factor(effects_tox$name, levels = c('bioassay', 'physiological'))
  
  gg_tox <- ggplot(effects_tox, aes(x = name, y = n)) +
              geom_bar(stat = 'identity') +
              xlab('') +
              ylab('') +
              ylim(0, 100) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1.1),
                    axis.text.y = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank()) +
            draw_image(needle, x = 0.5, y = 100, scale = 8) +
            draw_image(phys, x = 1.5, y = 100, scale = 11)
  
  # Create the final combined plot
  plot_summary <- ggarrange(gg_scale, gg_taxa, gg_col, gg_tox,
                             labels = c('(a)', '(b)', '(c)', '(d)'),
                             ncol = 2,
                             nrow = 2,
                             label.x = -0.03,
                             align = 'h')
  
  annotate_figure(plot_summary,
                  left = text_grob('no. effects', rot = 90, hjust = -0.2))
  
  ggsave('../figs/fig_summary.png', height = 10, width = 11)
  
  