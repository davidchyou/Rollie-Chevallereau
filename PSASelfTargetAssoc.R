library(ggplot2)
library(ape)

#Load up the raw data and create binary variables for presence of ACR and a complete CRISPR-Cas system.

df_has_acr <- read.csv("PSA_self_target_raw_data.csv")

df_has_acr <- df_has_acr[df_has_acr$HAS_ARRAY == 1,]

df_has_acr_cr <- df_has_acr[df_has_acr$HAS_CRISPR == 1,]
df_has_acr_no_acr <- df_has_acr[df_has_acr$HAS_ACR == 0,]

#Performing Fisher's exact test: +/- self-target vs +/- ACR, +/- self-target vs +/- comlete CRISPR, and tabulate the results

fm1 <- fisher.test(as.matrix(table(df_has_acr_cr$HAS_ST, df_has_acr_cr$HAS_ACR)))
lbl1 <- paste("Self-target, ACR with CRISPR\nOR=", round(fm1$estimate, 3), " P=", round(fm1$p.value, 3), sep = "")
if (fm1$p.value < 0.05) {
	lbl1 <- paste("Self-target, ACR with CRISPR\nOR=", round(fm1$estimate, 3), " P<0.05", sep = "")
}

df1 <- data.frame(LBL=c("ST-/ACR-", "ST+/ACR-", "ST-/ACR+", "ST+/ACR+"), 
                  FREQ=as.numeric(table(df_has_acr_cr$HAS_ST + 2*df_has_acr_cr$HAS_ACR)),
                  LBL2 = rep(lbl1, 4), PVAL = fm1$p.value)
                  
fm2 <- fisher.test(as.matrix(table(df_has_acr_no_acr$HAS_ST, df_has_acr_no_acr$HAS_CRISPR)))
lbl2 <- paste("Self-target, CRISPR without ACR\nOR=", round(fm2$estimate, 3), " P=", round(fm2$p.value, 3), sep = "")
if (fm2$p.value < 0.05) {
	lbl2 <- paste("Self-target, CRISPR, without ACR\nOR=", round(fm2$estimate, 3), " P<0.05", sep = "")
}

df2 <- data.frame(LBL=c("ST-/CRISPR-", "ST+/CRISPR-", "ST-/CRISPR+", "ST+/CRISPR+"), 
                  FREQ=as.numeric(table(df_has_acr_no_acr$HAS_ST + 2*df_has_acr_no_acr$HAS_CRISPR)),
                  LBL2 = rep(lbl2, 4), PVAL = fm2$p.value)
                  
fm3 <- fisher.test(as.matrix(table(df_has_acr_cr$HAS_ST_PRO, df_has_acr_cr$HAS_ACR)))
lbl3 <- paste("Self-target in prophages, ACR with CRISPR\nOR=", round(fm3$estimate, 3), " P=", round(fm3$p.value, 3), sep = "")
if (fm3$p.value < 0.05) {
	lbl3 <- paste("Self-target in prophages, ACR with CRISPR\nOR=", round(fm3$estimate, 3), " P<0.05", sep = "")
}

df3 <- data.frame(LBL=c("ST-/ACR-", "ST+/ACR-", "ST-/ACR+", "ST+/ACR+"), 
                  FREQ=as.numeric(table(df_has_acr_cr$HAS_ST_PRO + 2*df_has_acr_cr$HAS_ACR)),
                  LBL2 = rep(lbl3, 4), PVAL = fm3$p.value)
                  
fm4 <- fisher.test(as.matrix(table(df_has_acr_no_acr$HAS_ST_PRO, df_has_acr_no_acr$HAS_CRISPR)))
lbl4 <- paste("Self-target in prophages, CRISPR without ACR\nOR=", round(fm4$estimate, 3), " P=", round(fm4$p.value, 3), sep = "")
if (fm4$p.value < 0.05) {
	lbl4 <- paste("Self-target in prophages, CRISPR, without ACR\nOR=", round(fm4$estimate, 3), " P<0.05", sep = "")
}

df4 <- data.frame(LBL=c("ST-/CRISPR-", "ST+/CRISPR-", "ST-/CRISPR+", "ST+/CRISPR+"), 
                  FREQ=as.numeric(table(df_has_acr_no_acr$HAS_ST_PRO + 2*df_has_acr_no_acr$HAS_CRISPR)),
                  LBL2 = rep(lbl4, 4), PVAL = fm4$p.value)

#Plot the frequencies as bar-charts, with Fisher's exact test results displayed in the header.

dfout <- rbind(df1, df2)
dfout2 <- rbind(df3, df4)
dfout3 <- rbind(df1, df2, df3, df4)
pvs <- c(fm1$p.value, fm2$p.value, fm3$p.value, fm4$p.value)

gg <- ggplot(dfout, aes(x=LBL,y=FREQ)) + geom_bar(stat="identity")
gg <- gg + facet_grid(.~LBL2, scales="free") + xlab("") + ylab("Frequency") + theme(axis.text.x = element_text(angle=90,hjust=1.05,vjust=0.05))

gg2 <- ggplot(dfout2, aes(x=LBL,y=FREQ)) + geom_bar(stat="identity")
gg2 <- gg2 + facet_grid(.~LBL2, scales="free") + xlab("") + ylab("Frequency") + theme(axis.text.x = element_text(angle=90,hjust=1.05,vjust=0.05))

gg3 <- ggplot(dfout3, aes(x=LBL,y=FREQ)) + geom_bar(stat="identity")
gg3 <- gg3 + facet_grid(.~LBL2, scales="free") + xlab("") + ylab("Frequency") + theme(axis.text.x = element_text(angle=90,hjust=1.05,vjust=0.05))
