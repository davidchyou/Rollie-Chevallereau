library(dplyr)
library(ggplot2)
library(gplots)

#Read in the raw data
pdata <- read.csv("PHAST_mismatch_count_nr_full_05_b.csv")

#Collecting frequencies for subfigures a, and compute the base data frame for subfigure b and c
pdata_simple <- data.frame(GENERA = as.character(pdata$GENERA),
                           TOTAL = pdata$TOTAL,
                           N_MM0 = pdata$N_MM0,
                           N_MM5 = pdata$N_MM7,
                           N_PHAST = pdata$N_PHAST,
                           N_PRIMED = pdata$N_PROPHAGE_MISMATCH,
                           N_NAIVE = pdata$ANY_MATCH - pdata$N_PROPHAGE_MISMATCH,
                           N_NO_TARGET = pdata$N_PHAST - pdata$ANY_MATCH,
                           N_SPACER = pdata$N_SPACER,
                           N_HIT = pdata$N_HIT,
                           N_GENOME = pdata$N_GENOME,
                           N_GENOME_PHAST = pdata$N_GENOME_PHAST,
                           TOTAL_SF = pdata$TOTAL_SHUFFLED)

pdata_simple_agg <- aggregate(TOTAL_SF~GENERA, data=pdata_simple, sum)

pdata_simple <- pdata_simple %>% group_by(GENERA) %>% slice(c(1))
pdata_simple <- as.data.frame(pdata_simple)
pdata_simple$TOTAL_SF <- pdata_simple_agg$TOTAL_SF / 10

pdata_simple_sum <- data.frame(GENERA = c("All"), 
                               TOTAL = sum(pdata_simple$TOTAL),
                               N_MM0 = sum(pdata_simple$N_MM0),
                               N_MM5 = sum(pdata_simple$N_MM5),
                               N_PHAST = sum(pdata_simple$N_PHAST),
                               N_PRIMED = sum(pdata_simple$N_PRIMED),
                               N_NAIVE = sum(pdata_simple$N_NAIVE),
                               N_NO_TARGET = sum(pdata_simple$N_NO_TARGET),
                               N_SPACER = sum(pdata_simple$N_SPACER),
                               N_HIT = sum(pdata_simple$N_HIT),
                               N_GENOME = sum(pdata_simple$N_GENOME),
                           	   N_GENOME_PHAST = sum(pdata_simple$N_GENOME_PHAST),
                               TOTAL_SF = sum(pdata_simple$TOTAL_SF))


#Collecting frequencies for subfigures b and c, and only include genera with >= 500 BLASTN hits

pdata <- pdata[pdata$TOTAL >= 500,]
pdata_simple <- pdata_simple[pdata_simple$TOTAL >= 500,]
genera_total <- paste(as.character(pdata$GENERA), " (n=", pdata$TOTAL, ")", sep="") 
pdata <- cbind(pdata, data.frame(GENERA_TOTAL = genera_total))

pdata_simple_sum_2 <- data.frame(GENERA = c("All"), 
                                 TOTAL = sum(pdata_simple$TOTAL),
                                 N_MM0 = sum(pdata_simple$N_MM0),
                                 N_MM5 = sum(pdata_simple$N_MM5),
                                 N_PHAST = sum(pdata_simple$N_PHAST),
                                 N_SPACER = sum(pdata_simple$N_SPACER),
                                 N_HIT = sum(pdata_simple$N_HIT),
                                 TOTAL_SF = sum(pdata_simple$TOTAL_SF))

pdata$GENERA <- factor(pdata$GENERA, levels = rev(unique(pdata$GENERA)))
pdata$GENERA_TOTAL <- factor(pdata$GENERA_TOTAL, levels = rev(unique(pdata$GENERA_TOTAL)))

#Rearranging all data frames for GGPLOTs

nn <- dim(pdata_simple)[1]
df_plot1 <- data.frame(COUNT = c(pdata_simple$N_SPACER, pdata_simple$N_HIT), TYPE = rep(c("Spacer", "BLASTN hit"), each = nn), GENERA = rep(pdata_simple$GENERA, 2))
df_plot1$GENERA <- factor(df_plot1$GENERA, levels = rev(unique(df_plot1$GENERA)))

df_plot2 <- data.frame(COUNT = c(pdata_simple$N_PHAST, pdata_simple$TOTAL_SF, pdata_simple$N_MM0, pdata_simple$N_MM5), 
                       TYPE = rep(c("Prophage", "Control", "0 mismatch", "1-5 misatches"), each = nn), GENERA = rep(pdata_simple$GENERA, 4))
df_plot2$GENERA <- factor(df_plot2$GENERA, levels = rev(unique(df_plot2$GENERA)))
df_plot2$TYPE <- factor(df_plot2$TYPE, levels = rev(c("Prophage", "Control", "0 mismatch", "1-5 misatches")))

freqs <- aggregate(COUNT~GENERA, data=df_plot2[! df_plot2$TYPE %in% c("Control", "Prophage"),], sum)
names(freqs) <- c("GENERA", "TOTAL")

df_plot2 <- merge(df_plot2, freqs, by = c("GENERA"))

df_plot3 <- data.frame(COUNT = c(pdata_simple_sum$N_PHAST, pdata_simple_sum$TOTAL_SF, pdata_simple_sum$N_MM0, pdata_simple_sum$N_MM5),
                       TYPE = c("Prophage", "Control", "0 mismatch", "1-5 misatches"))
df_plot3$TYPE <- factor(df_plot3$TYPE, levels = c("Prophage", "Control", "0 mismatch", "1-5 misatches"))

df_plot4 <- data.frame(COUNT = c(pdata_simple$N_NO_TARGET, pdata_simple$N_NAIVE, pdata_simple$N_PRIMED),
                       TYPE = rep(c("No targets at all", "Targets found but none has mismatches", "Targets with mismatches were found"), each = nn), GENERA = rep(pdata_simple$GENERA, 3))
df_plot4$GENERA <- factor(df_plot4$GENERA, levels = rev(unique(df_plot4$GENERA)))

pgp <- pdata_simple$N_GENOME_PHAST / pdata_simple$N_GENOME
sdgp <- sqrt(pgp * (1 - pgp) / pdata_simple$N_GENOME)
genera_total <- paste(as.character(pdata_simple$GENERA), " (n=", pdata_simple$N_GENOME, ")", sep="") 
df_plot5 <- data.frame(PROP_GENOME_PHAST = pgp, SD_GENOME_PHAST = sdgp, GENERA = pdata_simple$GENERA, GENERA_TOTAL = genera_total)
df_plot5$GENERA <- factor(df_plot5$GENERA, levels = rev(unique(df_plot5$GENERA)))
df_plot5$GENERA_TOTAL <- factor(df_plot5$GENERA_TOTAL, levels = rev(unique(df_plot5$GENERA_TOTAL)))

#Draw the figures in GGPLOTs

gg <- ggplot(pdata, aes(x=GENERA_TOTAL,y=factor(N_NT_MISMATCH_SPACER)))+coord_flip()+geom_tile(aes(fill=PROP))
gg <- gg + xlab("") + ylab("Number of mismatches in spacer") + scale_fill_gradient(high="red", low="white") + labs(fill="Proportion")

gg2 <- ggplot(df_plot1, aes(x=GENERA,y=COUNT,fill=TYPE)) + geom_bar(stat="identity",position=position_dodge()) + coord_flip() + labs(x = "", y ="Count", fill="")
gg3 <- ggplot(df_plot2, aes(x=GENERA,y=COUNT,fill=TYPE)) + geom_bar(stat="identity",position=position_dodge()) + coord_flip() + labs(x = "", y ="Count", fill="") 
gg3 <- gg3 + scale_fill_discrete(breaks = c("Prophage", "Control", "0 mismatch", "1-5 misatches")) 
gg3 <- gg3 + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))
gg4 <- ggplot(df_plot3, aes(x=TYPE,y=COUNT)) + geom_bar(stat="identity") + labs(x = "", y ="Count")
gg4 <- gg4 + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))

gg5 <- ggplot(pdata, aes(x=GENERA_TOTAL,y=factor(N_NT_MISMATCH_SPACER)))+coord_flip()+geom_tile(aes(fill=log(COUNT)))
gg5 <- gg5 + xlab("") + ylab("Number of mismatches in spacer") + scale_fill_gradient(high="red", low="white") + labs(fill="log(Count)")
gg5 <- gg5 + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))

gg6 <- ggplot(df_plot4, aes(x=GENERA,y=COUNT,fill=TYPE)) + geom_bar(stat="identity") + coord_flip() + labs(x = "", y ="Number of prophages", fill="")
gg6 <- gg6 + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))

gg7 <- ggplot(df_plot5, aes(x=GENERA_TOTAL,y=100*PROP_GENOME_PHAST)) + geom_bar(stat="identity",fill="grey50") + coord_flip() + labs(x = "", y ="% genomes")
gg7 <- gg7 + geom_errorbar(aes(ymin = 100*(PROP_GENOME_PHAST - 1.96 * SD_GENOME_PHAST), ymax = 100*(PROP_GENOME_PHAST + 1.96 * SD_GENOME_PHAST)), width = 0.2, position = position_dodge(0.9))
gg7 <- gg7 + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))

mtrx <- do.call(c, lapply(0:5, FUN = function(x){pdata[pdata$N_NT_MISMATCH_SPACER == x,]$PROP}))
mtrx <- matrix(mtrx, dim(pdata)[1]/6, 6)
rownames(mtrx) <- unique(as.character(pdata$GENERA_TOTAL))
colnames(mtrx) <- 0:5

hm <- heatmap.2(mtrx, 
                dendrogram = "row", 
                Colv=FALSE, 
                density.info = "none", 
                trace="none", 
                margins=c(6,15), 
                col=colorpanel(255, low="white", high="red"), 
                xlab = "Number of mismatches in spacer", 
                key=TRUE)

pdata_simple <- rbind(pdata_simple, pdata_simple_sum)

