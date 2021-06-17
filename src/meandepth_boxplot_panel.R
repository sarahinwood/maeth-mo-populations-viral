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
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########
##depth files
st_depth_file_i44 <- snakemake@input[["st_depth_file_i44"]]
st_depth_file_i4 <- snakemake@input[["st_depth_file_i4"]]
st_depth_file_i60 <- snakemake@input[["st_depth_file_i60"]]
st_depth_file_i68 <- snakemake@input[["st_depth_file_i68"]]
st_depth_file_i76 <- snakemake@input[["st_depth_file_i76"]]

##scaffold labels
mh_gc_table <- snakemake@input[["mh_gc_table"]]

########
# MAIN #
########

## read in depth tables ##
st_depth_44 <- fread(st_depth_file_i44)
st_depth_4 <- fread(st_depth_file_i4)
st_depth_68 <- fread(st_depth_file_i68)
st_depth_76 <- fread(st_depth_file_i76)

##add sample label so tables can be joined
st_depth_44$Sample <- "44"
st_depth_4$Sample <- "4"
st_depth_68$Sample <- "68"
st_depth_76$Sample <- "76"
##join tables
full_depth_table <- rbind(st_depth_44, st_depth_4, st_depth_68, st_depth_76)

## scaffold ID table ##
full_scaffold_table <- fread(mh_gc_table, header=TRUE)
scaffold_table <- subset(full_scaffold_table, !(plot_label=="Other contig"))
scaffold_table$plot_label <- factor(scaffold_table$plot_label, levels=c("BUSCO contig", "BUSCO and viral contig", "Viral contig"))

## full table for plotting ##
st_depth_labels <- merge(full_depth_table, scaffold_table, by.x="#rname", by.y="#Name", all.y=TRUE)

##########################
##making panel box plot ##
##########################

##height/width in inches
pdf(snakemake@output[["boxplot_panel"]], height=7.5, width=10)
ggplot(st_depth_labels, aes(x=plot_label, y=meandepth, colour=plot_label))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  ylab("Mean sequencing depth")+
  scale_colour_viridis(discrete=TRUE, direction=-1)+
  coord_cartesian(ylim = c(0,40))+
  theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.title = element_blank())+
  facet_wrap(~Sample)
dev.off()

pdf(snakemake@output[["boxplot_panel_nofacet"]], height=7.5, width=10)
ggplot(st_depth_labels, aes(x=plot_label, y=meandepth, colour=plot_label))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  ylab("Mean sequencing depth")+
  scale_colour_viridis(discrete=TRUE, direction=-1)+
  coord_cartesian(ylim = c(0,40))+
  theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank())+
  facet_wrap(~Sample)
dev.off()

fwrite(st_depth_labels, snakemake@output[['depth_table']])

##add strip.background = element_blank(),
    strip.text.x = element_blank()
    ##to remove facet wrap titles

#write log
sessionInfo()
