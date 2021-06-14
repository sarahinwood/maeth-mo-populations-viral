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
st_depth_names <- c("#Name", "BP", "depth")
st_depth_44 <- fread(st_depth_file_i44, col.names=st_depth_names)
st_depth_4 <- fread(st_depth_file_i4, col.names=st_depth_names)
st_depth_60 <- fread(st_depth_file_i60, col.names=st_depth_names)
st_depth_68 <- fread(st_depth_file_i68, col.names=st_depth_names)
st_depth_76 <- fread(st_depth_file_i76, col.names=st_depth_names)

##add sample label so tables can be joined
st_depth_44$Sample <- "44"
st_depth_4$Sample <- "4"
st_depth_60$Sample <- "60"
st_depth_68$Sample <- "68"
st_depth_76$Sample <- "76"

full_depth_table <- rbind(st_depth_44, st_depth_4, st_depth_60, st_depth_68, st_depth_76)

## scaffold ID table ##
scaffold_table <- fread(mh_gc_table, header=TRUE)
scaffold_table$plot_group <- tstrsplit(scaffold_table$plot_label, " and", keep=c(1))
scaffold_table$scaffold_number <- tstrsplit(scaffold_table$'#Name', "_", keep=c(2))
scaffold_table$scaffold_number <- as.numeric(as.character(scaffold_table$scaffold_number))
setorder(scaffold_table, scaffold_number)

## full table for plotting ##
st_depth_labels <- merge(full_depth_table, scaffold_table, by="#Name", all.x=TRUE)
st_depth_boxpl <- subset(st_depth_labels, !(plot_group=="Other contig"))
##remove outlier contigs
st_depth_boxpl <- subset(st_depth_labels, !(`#Name` %in% mh_depth_outliers))
##order legend labels
st_depth_boxpl$plot_group <- factor(st_depth_boxpl$plot_group, levels=c("BUSCO contig", "BUSCO and viral contig", "Viral contig"))

##########################
##making panel box plot ##
##########################
##height/width in inches
pdf(snakemake@output[["boxplot_panel"]], height=7.5, width=10)
ggplot(st_depth_boxpl, aes(x=reorder(`#Name`, scaffold_number), y=depth, colour=plot_group))+
  geom_boxplot(outlier.shape=NA)+
  theme_bw(base_size=18)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun.y=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE)+
  coord_cartesian(ylim = c(0,35))+
  facet_wrap(~Sample)
dev.off()

#write log
sessionInfo()
