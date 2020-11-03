####################################################################################################
#This R-script belongs to:
#Tilte: Mycorrhizal fungi control value of phosphorus in trade symbiosis with host roots when exposed to abrupt 'crashes' and 'booms' of resource availability
#Authors: Anouk van 't Padje Gijsbert D.A. Werner, E. Toby Kiers
#Journal: New Phytologist
#DOI: 10.1111/nph.17055

#contact: anouk.vantpadje@wur.nl
#######################################################################################################
#load libraries
library(ART)
library(ordinal)
library(car)
library(RColorBrewer)
library(dplyr)
library(Rmisc)
library(MASS)
library(ggplot2)
library(multcompView)


########################################################################################################
#the error bar function
error.bar <- function(x, y, upper, lower=upper, length=0.03, ...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, col= "black", ...)
}

########################################################################################################
#open data
#setwd("")
means_per_petriplate_week2_3_4<- read.csv("datafile.csv")
means_per_petriplate_week2_3_4$treatment_week <- factor(means_per_petriplate_week2_3_4$treatment_week, levels=levels(means_per_petriplate_week2_3_4$treatment_week)[c(1,4,7, 2,5,8, 3,6,9)])
means_per_petriplate_week2_3_4$week <-as.factor(means_per_petriplate_week2_3_4$week)

intra_week2_3_4 <- read.csv("Intraradical_colonization.csv")

########################################################################################################

########################################################################################################
#Figure 1D: Quantum dot transfer per root
pdf("Figure_1D.pdf", width =3.43, height = 2)
op <- par(mar=c(1,1,1,1), xpd=NA)
means <- summarySE(means_per_petriplate_week2_3_4, measurevar = "core_QD.nmol..root.", groupvars = c("treatment_week"), na.rm= TRUE )

plot(means$core_QD.nmol..root.[c(1:3)], ylab = "Red quantum-dot-apatite \n per host root (nmol)", 
     xlab= "Days since event", cexaxis=0.5, pch=22, bg="#000000", cex=1, ylim=c(0.00,0.08), xlim=c(0.9,3), axes = FALSE)
points(c(1,2,3), means$core_QD.nmol..root.[c(4:6)], pch=23, cex=1, bg="#009E73")
points(c(1,2,3), means$core_QD.nmol..root.[c(7:9)], pch=24, cex=1, bg="#F0E442")

lines(means$core_QD.nmol..root.[c(1:3)], lty=1)
lines(means$core_QD.nmol..root.[c(4:6)], lty=2)
lines(means$core_QD.nmol..root.[c(7:9)], lty=4)

axis(side=1, at=c(1, 2, 3), label=c("7", "14", "21") , cex=0.75)
axis(side=2, at=c(0.0,0.01,0.02,0.03, 0.04, 0.05, 0.06, 0.07, 0.08), label=c("0.0", "0.01", "0.02", "0.03", "0.04",  "0.05", "0.06", "0.07", "0.08"), cex.axis=0.5, tck=c(-0.03),las=1)
axis(side=2, at=c(0.005,0.015,0.025,0.035, 0.045, 0.055, 0.065, 0.075), label=c("", "",  "", "",  "", "", "", ""), cex.axis=0.5, tck=c(-0.01))

arrows(x0=1,y0=means$core_QD.nmol..root.[1]-means$se[1], y1=means$core_QD.nmol..root.[1]+means$se[1], angle=90, code=3, length=0.03)
arrows(x0=2,y0=means$core_QD.nmol..root.[2]-means$se[2], y1=means$core_QD.nmol..root.[2]+means$se[2], angle=90, code=3, length=0.03)
arrows(x0=3,y0=means$core_QD.nmol..root.[3]-means$se[3], y1=means$core_QD.nmol..root.[3]+means$se[3], angle=90, code=3, length=0.03)

arrows(x0=1,y0=means$core_QD.nmol..root.[4]-means$se[4], y1=means$core_QD.nmol..root.[4]+means$se[4], angle=90, code=3, length=0.03)
arrows(x0=2,y0=means$core_QD.nmol..root.[5]-means$se[5], y1=means$core_QD.nmol..root.[5]+means$se[5], angle=90, code=3, length=0.03)
arrows(x0=3,y0=means$core_QD.nmol..root.[6]-means$se[6], y1=means$core_QD.nmol..root.[6]+means$se[6], angle=90, code=3, length=0.03)

arrows(x0=1,y0=means$core_QD.nmol..root.[7]-means$se[7], y1=means$core_QD.nmol..root.[7]+means$se[7], angle=90, code=3, length=0.03)
arrows(x0=2,y0=means$core_QD.nmol..root.[8]-means$se[8], y1=means$core_QD.nmol..root.[8]+means$se[8], angle=90, code=3, length=0.03)
arrows(x0=3,y0=means$core_QD.nmol..root.[9]-means$se[9], y1=means$core_QD.nmol..root.[9]+means$se[9], angle=90, code=3, length=0.03)


legend(0.9,0.09, legend=c("control", "crash", "boom"), lty=c(1,2,4), pch = c(22,23,24), pt.bg=c("#000000", "#009E73", "#F0E442"), cex=0.75, bty="n")

par(op)
dev.off()


mean_se(means_per_petriplate_week2_3_4$total_QD.nmol..root.[means_per_petriplate_week2_3_4$week=="4"]) #average uptake per root

mean_se(means_per_petriplate_week2_3_4$total_QD.nmol.mg.root.[means_per_petriplate_week2_3_4$week=="4"]) #average uptake per mg of root

#uptake per week
mean(means_per_petriplate_week2_3_4$total_QD.nmol..root.[means_per_petriplate_week2_3_4$week=="2"])
mean(means_per_petriplate_week2_3_4$total_QD.nmol..root.[means_per_petriplate_week2_3_4$week=="3"])
mean(means_per_petriplate_week2_3_4$total_QD.nmol..root.[means_per_petriplate_week2_3_4$week=="4"])

#red in root compartment
shapiro.test(means_per_petriplate_week2_3_4$core_QD.nmol..root.) # not normal
leveneTest(core_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
lm<-glm(core_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, type=2, test.statistic="F")
TukeyHSD(aov(lm))

#cyan IN root compartment (#STABLE compartment)
shapiro.test(means_per_petriplate_week2_3_4$stable_QD.nmol..root.)
leveneTest(stable_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4) # NO Homogeneity
lm<-glm(stable_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, type=2, test.statistic="F")
TukeyHSD(aov(lm))

#yellow IN root compartment (#treated compartment)
shapiro.test(means_per_petriplate_week2_3_4$treated_QD.nmol..root.)
leveneTest(treated_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
lm<-glm(treated_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, type=2, test.statistic="F")
TukeyHSD(aov(lm))





########################################################################################################
#Figure S1 TOTAL quantum-dot-apatite per root
pdf("FigureS1.pdf", width =5, height = 5)
op <- par(mar=c(5,5,2,2), xpd=NA)

means <- summarySE(means_per_petriplate_week2_3_4, measurevar = "total_QD.nmol..root.", groupvars = c("treatment_week"), na.rm= TRUE )

bar<-barplot(means$total_QD.nmol..root., 
             col = c("#000000", "#009E73", "#F0E442"), 
             ylab = "Total quantum-dot-apatite \n per host root (nmol)", xlab= "Days since event", 
             space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
             names.arg=c("", "7", "","", "14","","","21",""),
             ylim= c(0,max(0.2)), axes=FALSE)

axis(side=1, at=c(0, 3.8, 8.0, 11.9), label=c("", "", "", "") , cex=1)
axis(side=2, at=c(0,  0.05,  0.1, 0.150, 0.20), 
     label=c("0.0", "0.05", "0.10",  "0.15", "0.20"), cex.axis=1, tck=c(-0.03))
axis(side=2, at=c(0.025, 0.075, 0.125, 0.175), 
     label=c( "",  "","", ""), cex.axis=1, tck=c(-0.01))

error.bar(bar, means$total_QD.nmol..root., means$se)
legend(9,0.22, legend = c("control", "crash", "boom"), fill = c("#000000", "#009E73", "#F0E442"), box.col= NA, bg= NULL, bty="n" )


par(op)
dev.off()

#Stats
#TOTAL quantum-dot-apatite per root
shapiro.test(means_per_petriplate_week2_3_4$total_QD.nmol..root.)#normal
leveneTest(total_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
lm<-lm(total_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, type=2, test.statistic="F")
TukeyHSD(aov(lm))


########################################################################################################
#Figure S2: DRY ROOT MASS

pdf("Figure_S2.pdf", width =5, height = 5)
op <- par(mar=c(5,5,2,2), xpd=NA)

means <- summarySE(means_per_petriplate_week2_3_4, measurevar = "rootweight.mg.", groupvars = c("treatment_week"), na.rm= TRUE )

bar<-barplot(means$rootweight.mg., 
             col = c("#000000", "#009E73", "#F0E442"), 
             ylab = "Dry root weight (mg)", xlab= "Days since event", 
             space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
             names.arg=c("", "7", "","", "14","","","21",""),
             ylim= c(0,max(70)), axes=FALSE)

axis(side=1, at=c(0, 3.8, 8.0, 11.9), label=c("", "", "", "") , cex=1)
axis(side=2, at=c(0,  10,20,30,40,50,60,70), 
     label=c("0.0", "10", "20",  "30", "40", "50", "60", "70"), cex.axis=1, tck=c(-0.03))
axis(side=2, at=c(5,15,25,35,45,55,65), 
     label=c( "", "","", "","", "", ""), cex.axis=1, tck=c(-0.01))

error.bar(bar, means$rootweight.mg., means$se)
legend(9,80, legend = c("control", "crash", "boom"), fill = c("#000000", "#009E73", "#F0E442"), box.col= NA, bg= NULL, bty="n" )


par(op)
dev.off()

#statistics
shapiro.test(means_per_petriplate_week2_3_4$rootweight.mg.) # normal
leveneTest(rootweight.mg.~ treatment * week, data=means_per_petriplate_week2_3_4) # Homogeneity good
lm <- glm(rootweight.mg.~ treatment * week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, test.statistic="F")

TukeyHSD(aov(lm))



######################################################################################################
#Figure S3 Contribution per compartment at day 21
#create data frame
means_core <- summarySE(means_per_petriplate_week2_3_4, measurevar = "core_QD.nmol..root.", groupvars = c("treatment_week"), na.rm= TRUE )
means_stable <- summarySE(means_per_petriplate_week2_3_4, measurevar = "stable_QD.nmol..root.", groupvars = c("treatment_week"), na.rm= TRUE )
means_treated <- summarySE(means_per_petriplate_week2_3_4, measurevar = "treated_QD.nmol..root.", groupvars = c("treatment_week"), na.rm= TRUE )

counts <- means_core #root
counts$means_stable<-  means_stable$stable_QD.nmol..root.
counts$means_treated<-  means_treated$treated_QD.nmol..root.
row.names(counts) <- counts$treatment

counts$treatment_week <- NULL
counts$N <- NULL
counts$sd <- NULL
counts$se <- NULL
counts$ci <- NULL
counts <- t(counts)

#make figure
pdf("Figure_S3_all.pdf", width =5, height = 5)
op <- par(mar=c(5,5,2,2), xpd=NA)

bar<- barplot(counts[,c(1,4,7,2,5,8,3,6,9)],
              xlab="Treatment", ylim=c(0,0.18),
              ylab="Quantum-dot-apatite per host root (nmol)",
              col=c("red", "cyan", "yellow"), axes= FALSE,
              space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
              names.arg=c("", "control", "","", "crash","","","boom",""))

#error bars for root (bottem)
error.bar(bar, means_core$core_QD.nmol..root.[c(1,4,7,2,5,8,3,6,9)], means_core$se[c(1,4,7,2,5,8,3,6,9)])

# error bars for stable (middle)
error.bar(bar, means_stable$stable_QD.nmol..root.[c(1,4,7,2,5,8,3,6,9)]+
            means_core$core_QD.nmol..root.[c(1,4,7,2,5,8,3,6,9)], 
          means_stable$se[c(1,4,7,2,5,8,3,6,9)])

#error bars for treated (upper)
error.bar(bar, means_treated$treated_QD.nmol..root.[c(1,4,7,2,5,8,3,6,9)]+
            means_stable$stable_QD.nmol..root.[c(1,4,7,2,5,8,3,6,9)]+
            means_core$core_QD.nmol..root.[c(1,4,7,2,5,8,3,6,9)], means_treated$se[c(1,4,7,2,5,8,3,6,9)])


axis(side=1, at=c(0, 3.8, 8.0, 12), label=c( "", "", "",  ""), tick = NA)
axis(side=1, at=c(bar), label=c( "7", "14", "12",   "7", "14", "12",  "7", "14", "12"), tick = 0, padj=-1.3, cex.axis=0.7)
axis(side=2, at=c(0, 0.02,0.04,0.06,0.08, 0.10, 0.12, 0.14, 0.16, 0.18), 
     label=c("0", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12", "0.14", "0.16", "0.18"), cex.axis=0.7, tck=c(-0.03), las=1)
axis(side=2, at=c(0.01,0.03,0.05,0.07,0.09,0.11,0.13,0.15,0.17), 
     label=c("", "", "", "", "", "", "", "", ""), cex.axis=0.7, tck=c(-0.01), las=1)

legend(9, 0.18, legend=c("core", "stable", "manipulated"), fill=c("red", "cyan", "yellow"), cex=1, bty="n")

par(op)
dev.off()

# STATS
#red in root compartment
shapiro.test(means_per_petriplate_week2_3_4$core_QD.nmol..root.) # not normal
leveneTest(core_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
lm<-glm(core_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, type=2, test.statistic="F")

#cyan IN root compartment (#from STABLE compartment)
shapiro.test(means_per_petriplate_week2_3_4$stable_QD.nmol..root.)
leveneTest(stable_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4) # NO Homogeneity
lm<-glm(stable_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, type=2, test.statistic="F")

#yellow IN root compartment (# from treated compartment)
shapiro.test(means_per_petriplate_week2_3_4$treated_QD.nmol..root.)
leveneTest(treated_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
lm<-glm(treated_QD.nmol..root.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(lm))
Anova(lm, type=2, test.statistic="F")


########################################################################################################
#Figure 2A TOTAL  hypha per compartment(mg)
pdf("Figure3.pdf", width = 7.5, height = 5)
op <- par(mar=c(4,5,2,2), xpd=NA)

m <- rbind(c(0, 1, 0.55, 1), c(0, 0.4, 0, 0.55), c(0.3, 0.7, 0, 0.55), c(0.6, 1, 0, 0.55))
split.screen(m)

screen(1)

means <- summarySE(means_per_petriplate_week2_3_4, measurevar = "extraradical_hyphae_total.mg.", groupvars = c("treatment_week"), na.rm= TRUE )
bar<-barplot(means$extraradical_hyphae_total.mg., 
             col = c("#000000", "#009E73", "#F0E442"), ylab = "Total extraradical \n fungal biomass (mg)", cex.axis = 1, cex.lab= 0.8,
             xlab= "Days since event", space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
             names.arg=c("", "7", "","", "14","","","21",""),cex.names = 0.8,
             ylim= c(0,0.5), axes=FALSE)

error.bar(bar, means$extraradical_hyphae_total.mg., means$se)
legend(10,0.6, legend = c("control", "crash", "boom"), fill = c("#000000", "#009E73", "#F0E442"), box.col= NA, bg= NULL, cex=0.75 )

axis(side=1, at=c(0, 3.8, 8.0, 11.9), label=c("", "", "", "") , cex=0.5)
axis(side=2, at=c(0, 0.1, 0.20, 0.30, 0.40, 0.50), label=c("0.00", "0.10", "0.20",  "0.30",  "0.40", "0.50"), cex.axis=0.8,tck=c(-0.03), las=1)
axis(side=2, at=c(0.05, 0.15, 0.25, 0.35, 0.45), label=c("", "", "", "", ""), cex.axis=0.5, tck=c(-0.01))
text(0.1,0.6, "A", cex=1)

#FIGURE 2B-C-D
#TOTAL extraradical hypha per compartment
means <- summarySE(means_per_petriplate_week2_3_4, measurevar = "extraradical_hyphae_stable.mg.", groupvars = c("treatment_week"), na.rm= TRUE )
limits <- aes(ymax=means$biom + means$se, ymin=means$biom-means$se)

# # Fig 2B
#stable compartment
screen(2)
par(mar=c(4,5,2,2), xpd=NA)

bar<-barplot(means$extraradical_hyphae_stable.mg., 
             col = c("#000000", "#009E73", "#F0E442"), 
             col.axis= "black", col.main= "black", col.lab= "black",
             ylab = "Extraradical fungal biomass \n per compartment (mg)", xlab= "", 
             main="",
             space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
             names.arg=c("", "7", "","", "14","","","21",""),
             ylim= c(0,0.4), axes=FALSE, cex.lab=1, cex.axis=0.5, cex.lab=0.8)

error.bar(bar, means$extraradical_hyphae_stable.mg., means$se, length=0.03)

axis(side=1, at=c(0, 3.8, 8.0, 11.9), label=c("", "", "", "") , cex=0.5)
axis(side=1, at=c(1.9, 5.7, 9.5), label=c("7", "14", "21") , cex=0.5, tick=FALSE)
axis(side=2, at=c(0, 0.1, 0.20, 0.30, 0.40, 0.50), label=c("0.00", "0.10", "0.20",  "0.30",  "0.40", "0.50"), cex.axis=0.8,tck=c(-0.03), las=1)
axis(side=2, at=c(0.05, 0.15, 0.25, 0.35, 0.45), label=c("", "", "", "", ""), cex.axis=0.5, tck=c(-0.01))
text(0.1,0.6, "B", cex=1)
# # Fig 2C
#core compartment
screen(3)
means <- summarySE(means_per_petriplate_week2_3_4, measurevar = "extraradical_hyphae_core.mg.", groupvars = c("treatment_week"), na.rm= TRUE )

par(mar=c(4,5,2,2), xpd=NA)
bar<-barplot(means$extraradical_hyphae_core.mg., 
             col = c("#000000", "#009E73", "#F0E442"), 
             col.axis= "black", col.main= "black", col.lab= "black",
             ylab = "", xlab= "Days since event", 
             main="",
             space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
             ylim= c(0,0.4), axes=FALSE, cex.lab=1, cex.axis=0.5, cex.lab=0.8)

error.bar(bar, means$extraradical_hyphae_core.mg., means$se, length=0.03)

axis(side=1, at=c(0, 3.8, 8.0, 11.9), label=c("", "", "", "") , cex=0.5)
axis(side=1, at=c(1.9, 5.7, 9.5), label=c("7", "14", "21") , cex.lab=0.5, tick=FALSE, padj=0)
axis(side=2, at=c(0, 0.1, 0.20, 0.30, 0.40, 0.50), label=c("0.00", "0.10", "0.20",  "0.30",  "0.40", "0.50"), cex.axis=0.8,tck=c(-0.03), las=1)
axis(side=2, at=c(0.05, 0.15, 0.25, 0.35, 0.45), label=c("", "", "", "", ""), cex.axis=0.5, tck=c(-0.01))
text(0.1,0.6, "C", cex=1)

# # Fig 2D
#treated compartment
screen(4)
means <- summarySE(means_per_petriplate_week2_3_4, measurevar = "extraradical_hyphae_treated.mg.", groupvars = c("treatment_week"), na.rm= TRUE )

par(mar=c(4,5,2,2), xpd=NA)
bar<-barplot(means$extraradical_hyphae_treated.mg., 
             col = c("#000000", "#009E73", "#F0E442"), 
             col.axis= "black", col.main= "black", col.lab= "black",
             ylab = "", xlab= "", 
             main="",
             space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
             ylim= c(0,0.4), axes=FALSE, cex.lab=1.5, cex.axis=0.5, cex.lab=0.8)

error.bar(bar, means$extraradical_hyphae_treated.mg., means$se, length=0.03)

axis(side=1, at=c(0, 3.8, 8.0, 11.9), label=c("", "", "", "") , cex=0.5)
axis(side=1, at=c(1.9, 5.7, 9.5), label=c("7", "14", "21") , cex=0.5, tick=FALSE)
axis(side=2, at=c(0, 0.1, 0.20, 0.30, 0.40, 0.50), label=c("0.00", "0.10", "0.20",  "0.30",  "0.40", "0.50"), cex.axis=0.8,tck=c(-0.03), las=1)
axis(side=2, at=c(0.05, 0.15, 0.25, 0.35, 0.45), label=c("", "", "", "", ""), cex.axis=0.5, tck=c(-0.01))
text(0.1,0.6, "D", cex=1)

close.screen(all.screens = TRUE)
par(op)
dev.off()

# # STATS
# stats Total extra radiacal hyphal mass (fig 2A)
shapiro.test(means_per_petriplate_week2_3_4$extraradical_hyphae_total.mg.)
leveneTest(extraradical_hyphae_total.mg.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
total_biomass.aov <- glm(extraradical_hyphae_total.mg.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(total_biomass.aov))
Anova(total_biomass.aov, type=2, test.statistic = "F")
post <- TukeyHSD(aov(total_biomass.aov))
multcompLetters(post$`treatment:week`[,4])


#orange (572)= stable
shapiro.test(means_per_petriplate_week2_3_4$extraradical_hyphae_stable.mg.)
leveneTest(extraradical_hyphae_stable.mg.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
stable_biomass.aov <- glm(extraradical_hyphae_stable.mg.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(stable_biomass.aov))
Anova(stable_biomass.aov, type=2, test.statistic = "F")
post <- TukeyHSD(aov(stable_biomass.aov))
multcompLetters(post$`treatment:week`[,4])

#core
shapiro.test(means_per_petriplate_week2_3_4$extraradical_hyphae_core.mg.)
leveneTest(extraradical_hyphae_core.mg.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
core_biomass.aov <- glm(extraradical_hyphae_core.mg.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(core_biomass.aov))
Anova(core_biomass.aov, type=2, test.statistic = "F")
post <- TukeyHSD(aov(core_biomass.aov))
multcompLetters(post$`treatment:week`[,4])

#yellow (488)= treated
shapiro.test(means_per_petriplate_week2_3_4$extraradical_hyphae_treated.mg.)
leveneTest(extraradical_hyphae_treated.mg.~treatment*week, data=means_per_petriplate_week2_3_4) # Homogeneity good
treated_biomass.aov <- glm(extraradical_hyphae_treated.mg.~treatment*week, data=means_per_petriplate_week2_3_4)
qqnorm(resid(treated_biomass.aov))
Anova(treated_biomass.aov, test.statistic = "F")
pst <- TukeyHSD(aov(treated_biomass.aov))
multcompLetters(pst$`treatment:week`[,4])

#fold higher in core compartment
mean(means_per_petriplate_week2_3_4$extraradical_hyphae_core.mg.)/mean((mean(means_per_petriplate_week2_3_4$extraradical_hyphae_stable.mg.)+(means_per_petriplate_week2_3_4$extraradical_hyphae_treated.mg.)))
mean(means_per_petriplate_week2_3_4$extraradical_hyphae_stable.mg.)/mean(means_per_petriplate_week2_3_4$extraradical_hyphae_treated.mg.)



########################################################################################################
#FIGURE 3A intraradical hyphal abundance
means <- summarySE(intra_week2_3_4,measurevar = "intraradical_hyphae.copy..",  groupvars = c("treatment_week"),  na.rm= TRUE)

pdf("Figure4_swapped.pdf", width =5, height = 5)
op <- par(mar=c(5,5,1,1), xpd=NA)

bar<-barplot(means$intraradical_hyphae.copy.., 
             col = c("#000000", "#000000", "#000000", "#009E73", "#009E73","#009E73", "#F0E442", "#F0E442", "#F0E442"), 
             ylab = "Intraradical fungal abundace \n copy number (x 100000000)", 
             xlab= "Days since event", space= c(0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1),
             names.arg=c("7", "14", "21","7", "14","21","7","14","21"),
             ylim= c(0,800000000), axes=FALSE)


error.bar(bar, means$intraradical_hyphae.copy.., means$se)
legend(0, 800000000, legend = c("control", "crash", "boom"), fill = c("#000000", "#009E73", "#F0E442"), 
       box.col= NA, bg= NULL)

axis(side=1, at=c(0, 3.8, 8.0, 11.9), label=c("", "", "", "") , cex=1)
axis(side=2, at=c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), 
     label=c("0.00", "1.00", "2.00", "3.00", "4.00","5.00",  "6.00", "7.00", "8.00"), cex.axis=0.8, las=1, tck=c(-0.03))
axis(side=2, at=c(50000000, 150000000, 250000000, 350000000, 450000000, 550000000, 650000000, 750000000), 
     label=c("", "", "", "", "","",  "", ""), cex.axis=0.7, las=1, tck=c(-0.01))

par(op)
dev.off()

#stats
shapiro.test(intra_week2_3_4$intraradical_hyphae.copy..)
leveneTest(intraradical_hyphae.log.copy...~treatment*as.factor(week), data=intra_week2_3_4) # Homogeneity good
total_biomass.aov <- glm(intraradical_hyphae.log.copy...~treatment*as.factor(week), data=intra_week2_3_4)
qqnorm(resid(total_biomass.aov))
Anova(total_biomass.aov, type=2, test.statistic = "F")
pst <- TukeyHSD(aov(total_biomass.aov))
multcompLetters(pst$`treatment:as.factor(week)`[,4])


#comparison
1- (means$intraradical_hyphae.copy..[6]/means$intraradical_hyphae.copy..[9]) # CRASH vs EXTRA
1- (means$intraradical_hyphae.copy..[6]/means$intraradical_hyphae.copy..[3]) # CRASH vs CONTROL

0.33*means$intraradical_hyphae.copy..[9]


########################################################################################################
#FIGURE 3B: hyphal benefit: carbon to phosphorus
pdf("Figure_5.pdf", width =3.43, height = 2)
op <- par(mar=c(1,1,1,1), xpd=NA)

means <- summarySE(means_per_petriplate_week2_3_4,measurevar = "C.P_stable.treated.mg..mg.root.",  groupvars = c("treatment_week"),  na.rm= TRUE)

plot(means$C.P_stable.treated.mg..mg.root.[c(1:3)], ylim=c(6,8), 
     ylab = "Exchange rate \n (Extraradical fungal abundance/P transferred)", xlab="Days since event", cex=1, axes=FALSE, col= "white")

axis(side=1, at=c(1,2,3), label=c("7", "14", "21") , cex=0.5)
axis(side=2, at=c(6, 6.5, 7, 7.5, 8), label=c("6.0", "6.5", "7.0", "7.5", "8.0"), cex.axis=1, tck=c(-0.03), las = 0.5)
axis(side=2, at=c(6.25, 6.75, 7.25, 7.75), label=c("", "", "", ""), cex.axis=0.5, tck=c(-0.01), las=1)

points(means$C.P_stable.treated.mg..mg.root.[c(1:3)], pch=22, cex=1, bg= "#000000")#control
points(means$C.P_stable.treated.mg..mg.root.[c(4:6)], pch=23, cex=1, bg= "#009E73")#crash
points(means$C.P_stable.treated.mg..mg.root.[c(7:9)], pch=24, cex=1, bg= "#F0E442")#boom
lines(means$C.P_stable.treated.mg..mg.root.[c(4:6)], lty=2)
lines(means$C.P_stable.treated.mg..mg.root.[c(7:9)], lty=4)
lines(means$C2P[2:4], lty=1)

arrows(x0=1,y0=means$C.P_stable.treated.mg..mg.root.[1]-means$se[1], y1=means$C.P_stable.treated.mg..mg.root.[1]+means$se[1], angle=90, code=3, length=0.03)
arrows(x0=2,y0=means$C.P_stable.treated.mg..mg.root.[2]-means$se[2], y1=means$C.P_stable.treated.mg..mg.root.[2]+means$se[2], angle=90, code=3, length=0.03)
arrows(x0=3,y0=means$C.P_stable.treated.mg..mg.root.[3]-means$se[3], y1=means$C.P_stable.treated.mg..mg.root.[3]+means$se[3], angle=90, code=3, length=0.03)

arrows(x0=1,y0=means$C.P_stable.treated.mg..mg.root.[4]-means$se[4], y1=means$C.P_stable.treated.mg..mg.root.[4]+means$se[4], angle=90, code=3, length=0.03)
arrows(x0=2,y0=means$C.P_stable.treated.mg..mg.root.[5]-means$se[5], y1=means$C.P_stable.treated.mg..mg.root.[5]+means$se[5], angle=90, code=3, length=0.03)
arrows(x0=3,y0=means$C.P_stable.treated.mg..mg.root.[6]-means$se[6], y1=means$C.P_stable.treated.mg..mg.root.[6]+means$se[6], angle=90, code=3, length=0.03)

arrows(x0=1,y0=means$C.P_stable.treated.mg..mg.root.[7]-means$se[7], y1=means$C.P_stable.treated.mg..mg.root.[7]+means$se[7], angle=90, code=3, length=0.03)
arrows(x0=2,y0=means$C.P_stable.treated.mg..mg.root.[8]-means$se[8], y1=means$C.P_stable.treated.mg..mg.root.[8]+means$se[8], angle=90, code=3, length=0.03)
arrows(x0=3,y0=means$C.P_stable.treated.mg..mg.root.[9]-means$se[9], y1=means$C.P_stable.treated.mg..mg.root.[9]+means$se[9], angle=90, code=3, length=0.03)
text(3, 8, "*", cex=2)

legend(1,8, legend=c("control", "crash", "boom"), lty=c(1,2,4), pch = c(22,23,24), pt.bg=c("#000000", "#009E73", "#F0E442"), cex=0.5, bty="n")

close.screen(all.screens = TRUE)
par(op)
dev.off()

#time effect on both hyphal compartmentsshapiro.test(intra_week2_3_4$intraradical_hyphae.copy..)
#boom
shapiro.test(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment=="extra"])
leveneTest(C.P_stable.treated.mg..mg.root.~ week, data=means_per_petriplate_week2_3_4[means_per_petriplate_week2_3_4$treatment=="extra",])
lm <- glm(C.P_stable.treated.mg..mg.root.~ week, data=means_per_petriplate_week2_3_4[means_per_petriplate_week2_3_4$treatment=="extra",])
qqnorm(resid(lm))
Anova(lm,type= "2", test = "F")
TukeyHSD(aov(lm))


#control
shapiro.test(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment=="control"])
leveneTest(C.P_stable.treated.mg..mg.root.~ week, data=means_per_petriplate_week2_3_4[means_per_petriplate_week2_3_4$treatment=="control",])
lm <- lm(C.P_stable.treated.mg..mg.root.~ week, data=means_per_petriplate_week2_3_4[means_per_petriplate_week2_3_4$treatment=="control",])
qqnorm(resid(lm))
Anova(lm,type= "2", test = "F")
TukeyHSD(aov(lm))

#crash
shapiro.test(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment=="crash"])
leveneTest(C.P_stable.treated.mg..mg.root.~ week, data=means_per_petriplate_week2_3_4[means_per_petriplate_week2_3_4$treatment=="crash",])
lm <- lm(C.P_stable.treated.mg..mg.root. ~ week, data=means_per_petriplate_week2_3_4[means_per_petriplate_week2_3_4$treatment=="crash",])
qqnorm(resid(lm))
Anova(lm,type= "2", test = "F")
TukeyHSD(aov(lm))




#difference between each timepoint
wilcox.test(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="crash_2"],
            means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="extra_2"])

wilcox.test(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="crash_3"],
            means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="extra_3"])

wilcox.test(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="crash_4"],
            means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="extra_4"])

wilcox.test(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="control_4"],
            means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$treatment_week=="extra_4"])


#increase over time ?????
#from day 7 to 14
aveage_increase_7_14 <- mean(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$week=="3"]) 
- mean(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$week=="2"])

stdev_increase_7_17 <- sd(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$week=="3"]) 
+ sd(means_per_petriplate_week2_3_4$C.P_stable.treated.mg..mg.root.[means_per_petriplate_week2_3_4$week=="2"])


#from day 14 to 21
means$C.P_stable.treated.mg..mg.root.[c(7,8,9)] - means$C.P_stable.treated.mg..mg.root.[c(4,5,6)]
means_core$sd[c(7,8,9)]^2 + means_core$sd[c(4,5,6)]^2


