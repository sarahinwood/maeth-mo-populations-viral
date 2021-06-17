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

samtools_depth <- snakemake@input[["samtools_depth"]]
##scaffold labels
scaffold_id_table <- snakemake@input[["scaffold_id_table"]]

########
# MAIN #
########

##scaffold ID table
full_scaffold_table <- fread(scaffold_id_table, header=TRUE)
scaffold_table <- subset(full_scaffold_table, !(plot_label=="Other contig"))

##depth table
st_depth_names <- c("#Name", "BP", "depth")
st_depth <- fread(samtools_depth, col.names=st_depth_names)

##merge with table of scaffold ids for plotting
st_depth_labels <- merge(st_depth, scaffold_table, by='#Name', all.y=TRUE)
##order legend labels
st_depth_labels$plot_label <- factor(st_depth_labels$plot_label, levels=c("BUSCO contig", "BUSCO and viral contig", "Viral contig"))

#####################
## making box plot ##
#####################

pdf(snakemake@output[["boxplot_y_zoom"]])
ggplot(st_depth_labels, aes(x=plot_label, y=depth, colour=plot_label))+
  geom_boxplot(outlier.shape=NA)+
  theme_bw(base_size=18)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE, direction=-1)+
  coord_cartesian(ylim = c(0,200))
dev.off()

jpeg(snakemake@output[["boxplot"]])
ggplot(st_depth_labels, aes(x=plot_label, y=depth, colour=plot_label))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE, direction=-1)
dev.off()

#write log
sessionInfo()
