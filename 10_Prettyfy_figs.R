library(Seurat)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
library(scales)
library(stringr)
library(ggplot2)
library(ggrepel)
library(SpatialExperiment)
setwd(dir = "/Users/ASUS/Desktop/IJC/DUTRENEO/")
arm <- c("HotDUTRENEO","HotSTANDAR","ColdSTANDAR")
# Create the ggplot
# v_plot<- ggplot(data = big_df, aes(x = meandiff_CN, y = -log10(Adj.pval), col = diffexpressed, label = delabel)) +
#   geom_vline(xintercept = c(-0.2, 0.2), col = "black", linetype = "dashed") +
#   geom_point(size = 1.5, alpha = 0.7) +  # Set point size and transparency
#   scale_color_manual(values = c('Cr/Pr - Nr ∈ [-0.2,0.2]' = "grey", 'Cr/Pr - Nr ∈ (0.2, Inf)' = 'orange', 'Cr/Pr - Nr ∈ (-Inf,-0.2)' = 'darkolivegreen'),labels = c(expression(paste(bar('CR/PR'), ' - ',bar('NR') %in% '(0.2, Inf)')), expression(paste(bar('CR/PR'), ' - ',bar('NR') %in% '(-Inf,-0.2)')),expression(paste(bar('CR/PR'), ' - ',bar('NR') %in% '[-0.2,0.2]')))) + # Custom legend labels with mathematical symbols
#   theme_classic() +  # Use a minimal theme for a clean look
#   geom_label_repel(
#     aes(label = delabel), 
#     size = 3, 
#     show.legend  = F,
#     na.rm = TRUE, 
#     max.overlaps = Inf,  # Allow all labels to be plotted
#     force = 2,  # Increase force of repulsion for labels
#     force_pull = 0.5,  # Adjust force of attraction to data points
#     box.padding = 0.5,  # Adjust padding around label boxes
#     point.padding = 0.5,  # Adjust padding around data points
#     min.segment.length = 0,# Ensures all labels are connected with segments
#   ) +  
#   # Add dashed vertical lines
#   labs(
#     title = "HOT_DUTRENEO - Cancer Border",  # Add a dynamic title, replace with paste0(arm[b], " Arm - Cancer Border") if arm is defined
#     x = expression(paste(bar('CR/PR'), ' - ', bar('NR'))),  # Label x-axis
#     y = "-log10 Adjusted P-value",  # Label y-axis
#     color = ""  # Legend title
#   ) +
#   xlim(-0.4,0.4) + 
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Center and bold title
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size=12),# Increase axis title size
#     legend.position = "top",  # Move legend to top
#     legend.title = element_text(size = 14),  # Increase legend title size
#     legend.text = element_text(size = 9), # Increase legend text size
#     legend.margin = margin(t = -15)
#   ) +
#   guides(
#     col = guide_legend(
#       title.position = "top", 
#       nrow = 3,
#       override.aes = list(
#         shape = c(16, 16)  # Use colored dots as legend symbols
#       )
#     )  # Set legend title position, number of rows, and change legend symbols to colored dots
#   )
# v_plot
# ggsave(filename = paste0("Results_Correlation/DEFINITIVE_PLOTS/Volcano_",arm[b],".png"),plot = v_plot,width = 7,height = 7)

# HOT DUTRENEO
b=1
big_df <- readRDS(paste0("Results_Correlation/RDS/VolcanoDF_",arm[b],".RDS"))
big_df$diffexpressed <- "Cr/Pr - Nr ∈ [-0.2,0.2]"
#big_df <- unique(big_df[,c(11,13,14,15)])
big_df$diffexpressed[-log10(big_df$Adj.pval) > -log10(0.001) & big_df$meandiff_CN<(-0.2)] <- "Cr/Pr - Nr ∈ (-Inf,-0.2)"
big_df$diffexpressed[-log10(big_df$Adj.pval) > -log10(0.001) & big_df$meandiff_CN>(0.2)] <- "Cr/Pr - Nr ∈ (0.2, Inf)"
big_df$delabel <- as.factor(big_df$delabel)

sorted_C <- big_df[which(big_df$diffexpressed == "Cr/Pr - Nr ∈ (0.2, Inf)"),] %>% arrange(Adj.pval)
top_C <- sorted_C[1:10, ]
top_C <- top_C[which(!is.na(top_C$delabel)),]
sorted_N <- big_df[which(big_df$diffexpressed == "Cr/Pr - Nr ∈ (-Inf,-0.2)"),] %>% arrange(Adj.pval)
top_N <- sorted_N[1:10, ]
top_N <- top_N[which(!is.na(top_N$delabel)),]
top_10 <- rbind(top_N,top_C)

v_plot<- ggplot(data = big_df, aes(x = meandiff_CN, y = -log10(Adj.pval), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.2, 0.2), col = "black", linetype = "dashed") +
  geom_point(size = 1.5, alpha = 0.7) +  # Set point size and transparency
  scale_color_manual(values = c('Cr/Pr - Nr ∈ [-0.2,0.2]' = "grey",
                                'Cr/Pr - Nr ∈ (0.2, Inf)' = 'darkolivegreen',
                                'Cr/Pr - Nr ∈ (-Inf,-0.2)' = 'orange'
                                ),
                     labels = c("Higher in NR",
                                "Higher in CR/PR",
                                "")) + # Custom legend labels with mathematical symbols
  theme_classic() +  # Use a minimal theme for a clean look
  geom_label_repel(
    data = top_10,
    aes(label = delabel), 
    size = 3, 
    show.legend  = F,
    na.rm = TRUE, 
    max.overlaps = Inf,  # Allow all labels to be plotted
    force = 20,  # Increase force of repulsion for labels
    force_pull = 0.5,  # Adjust force of attraction to data points
    box.padding = 0.5,  # Adjust padding around label boxes
    point.padding = 0.2,  # Adjust padding around data points
    min.segment.length = 0,# Ensures all labels are connected with segments
  ) +  
  # Add dashed vertical lines
  labs(
    title = "HOT_DUTRENEO - Cancer Border",  # Add a dynamic title, replace with paste0(arm[b], " Arm - Cancer Border") if arm is defined
    x = expression(paste(bar('CR/PR'), ' - ', bar('NR'))),  # Label x-axis
    y = "-log10 Adjusted P-value",  # Label y-axis
    color = ""  # Legend title
  ) +
  xlim(-0.4,0.4) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Center and bold title
    axis.title = element_text(size = 14),
    axis.text = element_text(size=12),# Increase axis title size
    legend.position = "top",  # Move legend to top
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 14), # Increase legend text size
    legend.margin = margin(t = -15)
  ) +
  guides(
    col = guide_legend(
      title.position = "top", 
      nrow = 2,
      override.aes = list(
        shape = c(16, 16, NA),  # Use colored dots as legend symbols, set NA for the group to be removed
        linetype = c(1, 1, NA),  # Set NA for the group to be removed
        size = 5
      )
    )  # Set legend title position, number of rows, and change legend symbols to colored dots
  )
v_plot
ggsave(filename = paste0("Results_Correlation/DEFINITIVE_PLOTS/Volcano_Simplified_",arm[b],".png"),plot = v_plot,width = 7,height = 7)


# HOT STANDARD
b=2
big_df <- readRDS(paste0("Results_Correlation/RDS/VolcanoDF_",arm[b],".RDS"))
big_df$diffexpressed <- "Cr/Pr - Nr ∈ [-0.2,0.2]"
big_df$diffexpressed[-log10(big_df$Adj.pval) > -log10(0.001) & big_df$meandiff_CN<(-0.2)] <- "Cr/Pr - Nr ∈ (-Inf,-0.2)"
big_df$diffexpressed[-log10(big_df$Adj.pval) > -log10(0.001) & big_df$meandiff_CN>(0.2)] <- "Cr/Pr - Nr ∈ (0.2, Inf)"
big_df$delabel[which(big_df$diffexpressed!="Cr/Pr - Nr ∈ [-0.2,0.2]")] <- big_df$parejas[which(big_df$diffexpressed!="Cr/Pr - Nr ∈ [-0.2,0.2]")]
big_df$delabel <- as.factor(big_df$delabel)
big_df <- unique(big_df[,c(11,13,14,15)])



sorted_C <- big_df[which(big_df$diffexpressed == "Cr/Pr - Nr ∈ (0.2, Inf)"),] %>% arrange(Adj.pval)
top_C <- sorted_C[1:10, ]
top_C <- top_C[which(!is.na(top_C$delabel)),]
sorted_N <- big_df[which(big_df$diffexpressed == "Cr/Pr - Nr ∈ (-Inf,-0.2)"),] %>% arrange(Adj.pval)
top_N <- sorted_N[1:10, ]
top_N <- top_N[which(!is.na(top_N$delabel)),]
top_10 <- rbind(top_N,top_C)

v_plot<- ggplot(data = big_df, aes(x = meandiff_CN, y = -log10(Adj.pval), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.2, 0.2), col = "black", linetype = "dashed") +
  geom_point(size = 1.5, alpha = 0.7) +  # Set point size and transparency
  scale_color_manual(values = c('Cr/Pr - Nr ∈ [-0.2,0.2]' = "grey",
                                'Cr/Pr - Nr ∈ (0.2, Inf)' = 'darkolivegreen',
                                'Cr/Pr - Nr ∈ (-Inf,-0.2)' = 'orange'
  ),
  labels = c("Higher in NR",
             "Higher in CR/PR",
             "")) + # Custom legend labels with mathematical symbols
  theme_classic() +  # Use a minimal theme for a clean look
  geom_label_repel(
    data = top_10,
    aes(label = delabel), 
    size = 3, 
    show.legend  = F,
    na.rm = TRUE, 
    max.overlaps = Inf,  # Allow all labels to be plotted
    force = 20,  # Increase force of repulsion for labels
    force_pull = 0.5,  # Adjust force of attraction to data points
    box.padding = 0.5,  # Adjust padding around label boxes
    point.padding = 0.2,  # Adjust padding around data points
    min.segment.length = 0,# Ensures all labels are connected with segments
  ) +  
  # Add dashed vertical lines
  labs(
    title = "HOT_STANDARD - Cancer Border",  # Add a dynamic title, replace with paste0(arm[b], " Arm - Cancer Border") if arm is defined
    x = expression(paste(bar('CR/PR'), ' - ', bar('NR'))),  # Label x-axis
    y = "-log10 Adjusted P-value",  # Label y-axis
    color = ""  # Legend title
  ) +
  xlim(-0.4,0.4) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Center and bold title
    axis.title = element_text(size = 14),
    axis.text = element_text(size=12),# Increase axis title size
    legend.position = "top",  # Move legend to top
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 14), # Increase legend text size
    legend.margin = margin(t = -15)
  ) +
  guides(
    col = guide_legend(
      title.position = "top", 
      nrow = 2,
      override.aes = list(
        shape = c(16, 16, NA),  # Use colored dots as legend symbols, set NA for the group to be removed
        linetype = c(1, 1, NA),# Set NA for the group to be removed
        size = 5
      )
    )  # Set legend title position, number of rows, and change legend symbols to colored dots
  )
v_plot
ggsave(filename = paste0("Results_Correlation/DEFINITIVE_PLOTS/Volcano_Simplified_",arm[b],".png"),plot = v_plot,width = 7,height = 7)


b=3
big_df <- readRDS(paste0("Results_Correlation/RDS/VolcanoDF_",arm[b],".RDS"))
big_df$diffexpressed <- "Cr/Pr - Nr ∈ [-0.2,0.2]"
big_df$diffexpressed[-log10(big_df$Adj.pval) > -log10(0.001) & big_df$meandiff_CN<(-0.2)] <- "Cr/Pr - Nr ∈ (-Inf,-0.2)"
big_df$diffexpressed[-log10(big_df$Adj.pval) > -log10(0.001) & big_df$meandiff_CN>(0.2)] <- "Cr/Pr - Nr ∈ (0.2, Inf)"
big_df$delabel[which(big_df$diffexpressed!="Cr/Pr - Nr ∈ [-0.2,0.2]")] <- big_df$parejas[which(big_df$diffexpressed!="Cr/Pr - Nr ∈ [-0.2,0.2]")]
big_df$delabel <- as.factor(big_df$delabel)
big_df <- unique(big_df[,c(11,13,14,15)])


sorted_C <- big_df[which(big_df$diffexpressed == "Cr/Pr - Nr ∈ (0.2, Inf)"),] %>% arrange(Adj.pval)
top_C <- sorted_C[1:10, ]
top_C <- top_C[which(!is.na(top_C$delabel)),]
sorted_N <- big_df[which(big_df$diffexpressed == "Cr/Pr - Nr ∈ (-Inf,-0.2)"),] %>% arrange(Adj.pval)
top_N <- sorted_N[1:10, ]
top_N <- top_N[which(!is.na(top_N$delabel)),]
top_10 <- rbind(top_N,top_C)

# Cold Specific
v_plot<- ggplot(data = big_df, aes(x = meandiff_CN, y = -log10(Adj.pval), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.2, 0.2), col = "black", linetype = "dashed") +
  geom_point(size = 1.5, alpha = 0.7) +  # Set point size and transparency
  scale_color_manual(values = c('Cr/Pr - Nr ∈ [-0.2,0.2]' = "grey",
                                'Cr/Pr - Nr ∈ (-Inf,-0.2)' = 'orange',
                                'Cr/Pr - Nr ∈ (0.2, Inf)' = 'darkolivegreen'),
                     labels = c("Higher in NR",
                                "",
                                "Higher in CR/PR")) + # Custom legend labels with mathematical symbols
  theme_classic() +  # Use a minimal theme for a clean look
  geom_label_repel(
    data = top_10,
    aes(label = delabel), 
    size = 3, 
    show.legend  = F,
    na.rm = TRUE, 
    max.overlaps = Inf,  # Allow all labels to be plotted
    force = 50,  # Increase force of repulsion for labels
    force_pull = 0.5,  # Adjust force of attraction to data points
    box.padding = 0.5,  # Adjust padding around label boxes
    point.padding = 0.2,  # Adjust padding around data points
    min.segment.length = 0,# Ensures all labels are connected with segments
  ) +  
  # Add dashed vertical lines
  labs(
    title = "COLD_STANDARD - Cancer Border",  # Add a dynamic title, replace with paste0(arm[b], " Arm - Cancer Border") if arm is defined
    x = expression(paste(bar('CR/PR'), ' - ', bar('NR'))),  # Label x-axis
    y = "-log10 Adjusted P-value",  # Label y-axis
    color = ""  # Legend title
  ) +
  xlim(-0.5,0.5) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Center and bold title
    axis.title = element_text(size = 14),
    axis.text = element_text(size=12),# Increase axis title size
    legend.position = "top",  # Move legend to top
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 14), # Increase legend text size
    legend.margin = margin(t = -15)
  )+
  guides(
    col = guide_legend(
      title.position = "top", 
      nrow = 2,
      override.aes = list(
        shape = c(16, NA),  # Use colored dots as legend symbols, set NA for the group to be removed
        linetype = c(1,NA),  # Set NA for the group to be removed
        size = 5
      )
    )  # Set legend title position, number of rows, and change legend symbols to colored dots
  )
v_plot
ggsave(filename = paste0("Results_Correlation/DEFINITIVE_PLOTS/Volcano_Simplified_",arm[b],".png"),plot = v_plot,width = 7,height = 7)


a<- readRDS("RDS/Border_ID_Spots_HotStandar.RDS")
