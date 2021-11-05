#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(dplyr)
library(stringr)

###########
# GLOBALS #
###########

mo_gc_table <- snakemake@input[["mo_gc_table"]]
coverage_file <- snakemake@input[["coverage_file"]]

####################
## test normality ##
####################

gc_table <- fread(mo_gc_table)
coverage <- fread(coverage_file)

##gc depth table
depth_table <- coverage[,c(1,7)]
gc_depth_stats <- merge(gc_table, depth_table, by.x="#Name", by.y="#rname")

##remove other
gc_depth_stats <- subset(gc_depth_stats, !(plot_label=="Other contig"))
gc_depth_stats <- subset(gc_depth_stats, !(`#Name`=="scf1605"))

gc_depth_stats$plot_label <- str_replace_all(gc_depth_stats$plot_label, "BUSCO and viral contig", "BUSCO contig")


##test normality of meandepth
shapiro <- shapiro.test(gc_depth_stats$meandepth)
shapiro_chars <- capture.output(print(shapiro))
writeLines(shapiro_chars, con=file(snakemake@output[["shapiro_res"]]))

##cannot use parametric tests

##########################
## non-parametric tests ##
##########################

summary <- group_by(gc_depth_stats, plot_label) %>%
  summarise(
    count = n(),
    mean = mean(meandepth, na.rm = TRUE),
    sd = sd(meandepth, na.rm = TRUE),
    median = median(meandepth, na.rm = TRUE),
    IQR = IQR(meandepth, na.rm = TRUE)
  )
summary_chars <- capture.output(print(summary))
writeLines(summary_chars, con=file(snakemake@output[["summary_stats"]]))

##does meandepth content between any groups differ significantly
kruskal_res <- kruskal.test(meandepth ~ plot_label, data = gc_depth_stats)
kruskal_chars <- capture.output(print(kruskal_res))
writeLines(kruskal_chars, con=file(snakemake@output[["kruskal_res"]]))

##which groups differ significantly
pair_wilcox_res <- pairwise.wilcox.test(gc_depth_stats$meandepth, gc_depth_stats$plot_label,
                                        p.adjust.method = "BH")
wilcox_chars <- capture.output(print(pair_wilcox_res))
writeLines(wilcox_chars, con=file(snakemake@output[["wilcox_res"]]))

#write log
sessionInfo()