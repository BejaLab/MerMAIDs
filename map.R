#install.packages(svglite)
suppressMessages(library(dplyr))
library(tidyr)
library(readr)
library(maps)
library(ggplot2)
library(ggthemes)


oga_abundance_matrix <- readr::read_tsv("./data/OGA_results/OGA_Data_type_1_plus_MerMAID_pro/data/abundance_matrix.csv", skip = 1)#, sep = "\t" )

df.tara_metadata <- readr::read_tsv("./data/OGA_results/OGA_Data_type_1_plus_MerMAID_pro/data/environmental_parameters.csv")#, sep = "\t")

df.oga <- oga_abundance_matrix[,-2] # remove the tax column

df.oga <- gather(df.oga,sample_ID, abundance, -Gene_ID)

MerMAID_ids <- read.table("./data/OGA_results/OGA_Data_type_1_plus_MerMAID_pro/homologs/filtered/MerMAID_ids.txt",col.names = c("Gene_ID"))
MerMAID_ids$rhodopsin_type <- "MerMAID"

df.oga <- left_join(df.oga, MerMAID_ids) # add Station info to the DF

df.oga$rhodopsin_type[is.na(df.oga$rhodopsin_type)] <- "non_MerMAID"

#sum abundance of all related MerMAIDs per sample
df.abundance <- aggregate(abundance ~ rhodopsin_type + sample_ID, df.oga, FUN=sum)

# calculate total rhodopsin abundance per sample and rename the column to total_rhod
total_rhod_per_sample <- aggregate(abundance ~ sample_ID, df.oga, FUN=sum)
colnames(total_rhod_per_sample)[2] <- "total_rhod"

# add total_rhod to df.abundance
df.abundance <- left_join(df.abundance, total_rhod_per_sample)

# calculate rhod_ratio % per sample
df.abundance$rhod_ratio <- (df.abundance$abundance/df.abundance$total_rhod)*100

# keep only MerMAID
df.abundance <- df.abundance %>% filter(rhodopsin_type == "MerMAID")

metadata <- df.tara_metadata  %>% select("sample_ID", "Station")
df.abundance <- inner_join(df.abundance, metadata) # add Station to df.abundance

df.abundance_per_station <- aggregate(rhod_ratio ~ Station, df.abundance, mean) #mean MerMAIDs ratio per station

metadata <- df.tara_metadata  %>%
  select("Station","latitude", "longitude")  %>%
  distinct(Station, .keep_all = TRUE)# subset metadata to station, lat and lon
df.abundance_per_station <- inner_join(df.abundance_per_station, metadata)

# Depth profile

df.detected_samples <-df.abundance %>% filter(rhod_ratio > 0.01) # keep only samples with rhod_ratio > 0.01%

nrow(df.detected_samples) # rhod_ratio > 0.01 %

#stats from those high abundance
summary(df.detected_samples)

### Extract depth
metadata <- df.tara_metadata  %>% select("sample_ID","depth","Depth (m)")

# add depth to dataframe
df.detected_samples <- inner_join(df.detected_samples, metadata)

# classify each sample into a depth_bin
df.detected_samples$depth_bin <- cut(df.detected_samples$`Depth (m)`, seq(0,1050,50), labels = seq(50,1050,50))
df.detected_samples$depth_bin <- as.numeric(paste(df.detected_samples$depth_bin))

bin_size <- aggregate(sample_ID ~ depth_bin, df.detected_samples, FUN = length)
colnames(bin_size)[2] <- "bin_size"
df.detected_samples <- left_join(df.detected_samples, bin_size)

write.csv(df.detected_samples, file = "./data/Depth_profile.csv")

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))


bp <- ggplot(df.detected_samples %>% filter(bin_size >= 5), aes(depth_bin,rhod_ratio)) + stat_boxplot(geom ='errorbar', aes(group = cut_width(depth_bin, 50))) +
  geom_boxplot(aes(group = cut_width(depth_bin, 50)), fill="#8ebdbe", color = "#4d6769",
               outlier.shape = NA,
               outlier.size = 2) + coord_flip() + main_theme + scale_x_reverse(breaks = seq(50,1000,50)) + scale_y_log10(breaks=c(.01,.1, 1, 3), labels=c("0.01","0.1", "1", "3"))#)

bp <- bp + labs(x = "Depth (m)", y="MerMAIDs/rhodopsins (%)")
bp <-bp + geom_jitter(data = df.detected_samples, width = 0.2, shape=21, fill="#E5EEEF", color="black")
bp
ggsave(plot = bp, filename = "./figures/Depth_profile.svg", width = 10, height = 10, units = "in")

gg <- ggplot()

wrld <- map_data("world")

xlims = c(-155, 70)
ylims = c(-70, 70)

p <- ggplot()
p <- p + theme(panel.background = element_rect(fill =NA),
               panel.border = element_rect(colour = "#000000",
                                           size = 1,
                                           linetype = "solid",
                                           fill = NA),
               axis.title = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_text(),
               axis.ticks.y = element_line(),
               legend.position="bottom",
               legend.background = element_rect(fill="white", colour = "black"),
               legend.key = element_rect(fill=NA))

#Draws the map and assigns background color for continents
p <-p + geom_polygon( data=wrld, aes(x=long, y=lat, group = group), colour="#b3b3b3", fill="#b3b3b3")

#Plots negative stations
neg_map <- p + geom_point( data=df.abundance_per_station %>% filter(rhod_ratio == 0),
                           shape = 4,
                           #                           colour = "#a7a7a7",
                           #                           fill = "#a7a7a7",
                           colour="blue",fill="blue",
                           size = 3,
                           aes(x=longitude, y=latitude)
)

#Add positive stations sized by abundance

p <- neg_map + geom_point( data=df.abundance_per_station %>%
                             filter(rhod_ratio > 0), 
                           shape=21, 
                           colour="#b21616", fill="#e84646",
                           aes(x=longitude, y=latitude, size=rhod_ratio)
)

p <- p + coord_fixed(xlim = xlims, ylim = ylims)
scale_values <- c(0.0000001, 0.01,0.1,0.2,0.4,0.6)
p <- p + scale_size(breaks = scale_values, labels=c("0", "0.01 %","0.1 %","0.2 %","0.4 %", "0.6%"), name = "Average MerMAIDs/rhodopsins per station (%)")# + sca(name = "New Legend Title")
p
ggsave(plot = p, filename = "./figures/map.svg", width = 10, height = 10, units = "in")
