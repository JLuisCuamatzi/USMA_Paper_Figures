# script to plot the coverage
#### script for R version = r/3.6.1
# load libraries
library(data.table)
library(zoo)
library(ggplot2)
library(tidyverse)
library(dplyr)

## Arguments:
# COVERAGE_FILE=''
# WINDOW=''
# STEP=''
# NORM_TABLE=''
# 

# to indicate arguments in command line
args=(commandArgs(TRUE)) 

if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# indicate the number of threads for data.table library
if( getDTthreads() < 10 || getDTthreads() > 10){
  setDTthreads(threads = 10)
}else{print(paste("WARNING: data.table package is running with ", getDTthreads(), " threads.", sep=''))}

## INPUTS:
SIZE <- as.numeric(WINDOW)      # window size for coverage estimate SIZE <- as.numeric(1000) 
STEP <- as.numeric(STEP)        # step in coverage estimate   STEP <- as.numeric(1000)
# 1) row coverage file
print(COVERAGE_FILE)
df <- fread(COVERAGE_FILE)      # read the file


## ESTIMATE THE RAW COVERAGE IN NON-OVERLAPPING WINDOWS OF SIZE N

names(df) <- c("Chr","Position","depth") # rename the headers
sample.median <- median(df$depth) #estimate the global coverage
# calculate the median in windows of size n
depth2plot <- df[, .(window.start = rollapply(Position, width=SIZE, by=STEP, FUN=min, align="left"),
                     window.end = rollapply(Position, width=SIZE, by=STEP, FUN=max, align="left"),
                     coverage.median = rollapply(depth, width=SIZE, by=STEP, FUN=median, align="left")), .(Chr)] 

depth2plot$Chr <- as.numeric(gsub("USMA_521_v2_", "", depth2plot$Chr)) # remove the string "USMA_521_v2_" in the chromosome name


## ESTIMATE THE NORMALIZED COVERAGE BY WHOLE CHROMOSOME 
depth4chr <- df %>% group_by(Chr) %>% summarise(Median.cov.chr = median(depth))   # 1) median coverage by chr
depth4chr$Norm.global.coverage <- depth4chr$Median.cov.chr/sample.median          # 2) normalized coverage by chr
depth4chr$Chr <- as.numeric(gsub("USMA_521_v2_", "", depth4chr$Chr))              # 3) remove the string "USMA_521_v2_" in the chromosome name

## ADD THE MEDIAN COVERAGE BY CHR TO depth2plot object
depth2plot <- depth2plot %>% left_join(select(depth4chr, Chr, Median.cov.chr), by = c("Chr" = "Chr"))

## ESTIMATE THE NORMALIZED COVERAGE IN NON-OVERLAPPING WINDOWS
Norm.cov <- depth2plot %>% group_by(Chr) %>% mutate(Norm.cov = coverage.median/sample.median) %>% setDT()
Norm.cov$Sample <- rep(SAMPLE, length(Norm.cov$Chr)) # Norm.cov$Sample <- rep("2021EE01", length(Norm.cov$Chr))
Norm.cov.mitoch <- Norm.cov[Chr == 24] # individual df for mitoch (USMA_521_v2_24)


## WRITE TABLE WITH THE RAW AND NORMALIZED COVERAGE IN NON-OVERLAPPING WINDOWS
write.table(Norm.cov, NORM_TABLE, row.names = F, quote = F, sep = "\t")

if (file.exists(NORM_TABLE)){
  print("The output file was succesfully written!")
} else {
  print("Error in the write!")
}

## remove 24 to 28! Just keep nuclear chr
Norm.cov <- Norm.cov[Chr != 24]
Norm.cov <- Norm.cov[Chr != 25]
Norm.cov <- Norm.cov[Chr != 26]
Norm.cov <- Norm.cov[Chr != 27]
Norm.cov <- Norm.cov[Chr != 28]


## PLOTTING
# Raw coverage
plot <- ggplot(depth2plot, aes(x = window.end, y = coverage.median, colour = as.factor(Chr))) + 
  geom_point() + 
  facet_grid(.~ Chr, space = 'free_x', scales = "free_x" ) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,sample.median*3)) + 
  geom_hline(yintercept = sample.median, linetype = "dashed", color = "red", size = 2) + 
  labs(title = paste("Sample: ", SAMPLE, " (", SIZE, " non-overlapping bp)", sep = ""), 
       x = "Chromosome", 
       y = "Median coverage depth") + 
  theme(legend.position = "none", panel.spacing.x = grid::unit(0, "cm"), 
        panel.border = element_rect(colour = "grey", size = 0.1), panel.ontop = FALSE, 
        axis.title.x = element_text(face = "bold", color = "black", size = 22), 
        axis.title.y = element_text(face = "bold", color = "black", size = 22), 
        axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 18), 
        axis.ticks.x= element_blank(), 
        plot.title = element_text(face = "bold", color = "red", size = 24, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold")) # make the plot

ggsave(OUTPUT_PLOT_1, plot = plot, width = 20, height = 12, units = "in", dpi = 300) # save plot with raw coverage

# Normalized coverage
plot.2 <- ggplot(Norm.cov, aes(x = window.end, y = Norm.cov, color = as.factor(Chr))) +
  geom_area(aes(alpha = 0.2, fill = as.factor(Chr))) +
  scale_y_continuous(limits = c(0,5))+
  geom_hline(yintercept = 1) +
  facet_wrap(.~ Chr, scales = "free_x", strip.position = "bottom") +
  theme_bw() + 
  labs(title = paste("Sample: ", SAMPLE, " (", SIZE, " non-overlapping bp)", sep = ""), 
       x = "Chromosome", 
       y = "Median coverage depth") + 
  theme(legend.position = "none", panel.spacing.x = grid::unit(0, "cm"), 
        panel.border = element_rect(colour = "white", size = 0.1), panel.ontop = FALSE, 
        axis.title.x = element_text(face = "bold", color = "black", size = 22), 
        axis.title.y = element_text(face = "bold", color = "black", size = 22), 
        axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 18), 
        axis.ticks.x= element_blank(), 
        plot.title = element_text(face = "bold", color = "red", size = 24, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold")); plot.2

ggsave(OUTPUT_PLOT_2, plot = plot.2, width = 20, height = 12, units = "in", dpi = 300) # save the plot

# Normalized coverage in mitoch

plot.3 <- ggplot(Norm.cov.mitoch, aes(x = window.end, y = Norm.cov, color = as.factor(Chr))) +
  geom_area(aes(alpha = 0.2, fill = as.factor(Chr))) +
  scale_y_continuous()+
  geom_hline(yintercept = 1) +
  facet_wrap(.~ Chr, scales = "free_x", strip.position = "bottom") +
  theme_bw() + 
  labs(title = paste("Sample: ", SAMPLE, " (", SIZE, " non-overlapping bp)", sep = ""), 
       x = "Chromosome", 
       y = "Median coverage depth") + 
  theme(legend.position = "none", panel.spacing.x = grid::unit(0, "cm"), 
        panel.border = element_rect(colour = "white", size = 0.1), panel.ontop = FALSE, 
        axis.title.x = element_text(face = "bold", color = "black", size = 22), 
        axis.title.y = element_text(face = "bold", color = "black", size = 22), 
        axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 18), 
        axis.ticks.x= element_blank(), 
        plot.title = element_text(face = "bold", color = "red", size = 24, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold")); plot.2

ggsave(OUTPUT_PLOT_3, plot = plot.3, width = 20, height = 12, units = "in", dpi = 300)

rm(list = ls()) # clean the environment





