#!/usr/bin/env Rscript --vanilla

### This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
### This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
### You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

library(pbapply)
library(ggplot2)

## Define some functions

read_table_file <- function(fileName) {
    ## Function to read the table files. fileName can also be an URL
    table.df <- read.csv(fileName, sep="\t")
    colnames(table.df) = c("ProbeID", "Intensity", "DetectionPValue")
    table.df <- table.df[,!names(table.df) == "DetectionPValue"] # Actually, we don't need it
    table.df$Sample <- strsplit(basename(fileName), "_")[[1]][1]
    return(table.df)
}

getSymbols <- function(probes) {
    ## Get a vector of probes, and convert them into a vector of Symbols (genes). NAs can be returned
    
    symbols <- unlist(as.list(hgu133a.db::hgu133aSYMBOL[probes]))
    ## The above process fails to convert some important probes to genes, we need to manually do it.
    ##
    ## ARNTL
    ## BMAL1 a.k.a. ARNTL. That's the gene name used in the paper and the one I will favour.
    ## See here: http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/209824_s_at.html
    ## and here: http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/210971_s_at.html
    select <- probes == "209824_s_at" | probes == "210971_s_at"
    symbols[select] = "ARNTL"
    ## CSNK1E
    ## See here: http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/202332_at.html
    select <- probes == "202332_at"
    symbols[select] = "CSNK1E" 
    ## NR1D1
    ## It has two associated probes.
    ## See here: http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/31637_s_at.html
    ## and here: http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/204760_s_at.html
    select <- probes == "31637_s_at" | probes == "204760_s_at"
    symbols[select] = "NR1D1"

    return(symbols)
}

probeRLE <- function(probe.df) {
    ## Calculate the relative log expression (RLE) per probe (log2)
    log2medians <- by(probe.df,
                      probe.df$ProbeID,
                      function(a.df) median(log2(a.df$Intensity))) # calculate the median of the log2 intensity per probe
    rle <- sapply(1:nrow(probe.df),
                    function(i) log2(probe.df[i,"Intensity"]) - log2medians[probe.df[i,"ProbeID"]]) # calculate the deviation from it
    return(rle)
}

guessGender <- function(adonor, data) {
    ## Guess the donor's gender based on the PBMCs expression of 4 genes. See paper.
    a.df <- data[data$Donor == adonor & grepl("DDX3Y|KDM5D|USP9Y|RPS4Y1", data$Symbol), c("rle", "Sample", "Symbol")]
    a.m <- reshape2::acast(a.df, Sample ~ Symbol, value.var = "rle")
    a.pca <- prcomp(a.m)                                                                                   
    pc1 <- sum(a.pca[["rotation"]][,"PC1"])
    return(ifelse(pc1 < 0, "Female", "Male"))
}

## Fit a simple linear model
simplef <- function(gene, data) {
    ## simple model fitting function
    gene.df <- data[gene == data$Symbol,]
    gene.df$Donor <- droplevels(gene.df$Donor)
    control = lme4::lmerControl(check.conv.singular = "ignore") # otherwise we get lots of errors
    lme4::lmer(rle ~ 1 + (t | Gender) + (1 | Donor), data = gene.df, REML=FALSE, control = control)
}

## fit a cosine model
cosinef <- function(gene, data) {
    ## cosine model fitting function
    ## Quoting the paper:
    ##
    ## The general formula of the fitted model is given by:
    ##
    ## Yⱼᵢₖ = a + b cos(2 π tᵢₖ) + c sin(2 π tᵢₖ) + d (fixed covariates) + g(random intercepts) + εⱼᵢₖ
    ##
    ## where Yⱼᵢₖ represents the log2 expression of gene j for individual i recorded at time tᵢₖ,
    ## with tᵢₖ computed as the calendar day of the date of bleed divided by
    ## the total number of days within the equivalent year.
    ## The fixed covariates and random intercepts terms were data-set-specific.
    ## [..]
    ## The identity of each subject of the BABYDIET and of the asthma data sets were modelled as
    ## a random intercept in the corresponding models.
    gene.df <- data[gene == data$Symbol,]
    gene.df$Donor <- droplevels(gene.df$Donor)
    control = lme4::lmerControl(check.conv.singular = "ignore") # otherwise we get lots of warnings
    lme4::lmer(rle ~ 1 + cos(t) + sin(t) + (t | Gender) + (1 | Donor), data = gene.df, REML=FALSE, control = control)
}

## Predict a model
predict_log2I <- function(a_fit, a_time = seq(0,2*pi, by=0.1)) {
    ## Function that given a fit function and a time vector with values from 0 to 1, predicts a log2(I)
    x <- as.matrix(coef(a_fit)[["Donor"]][1, c("(Intercept)", "cos(t)", "sin(t)")])            # x = [a b c]
    y <- matrix(c(rep(1, times = length(a_time)),      
                  cos(a_time),                      #      ⸢      1       ...        1     ⸣
                  sin(a_time)),                     # y =  |cos(2 π t₁)   ...  cos(2 π tₘ) |
                nrow=3,                             #      ⸤sin(2 π t₁)   ...  sin(2 π tₘ) ⸥
                byrow=TRUE)
    I_hat <- (x %*% y)[1,] # Matrix multiplication x * y
    return(I_hat)
}

predict_seasonality <- function(data.df) {
    ## Detect genes with seasonality expression and return a data frame with the model prediction for a year

    geneNames <- levels(data.df$Symbol)
    ## Simple model
    print("Fitting a simple linear model for all genes")
    simple.fit.l <- pbsapply(geneNames, simplef, data = data.df )
    ## Cosine model
    print("Fitting a cosine linear model for all genes")
    cosine.fit.l <- pbsapply(geneNames, cosinef, data = data.df)
    ## Compare both models
    ## https://bookdown.org/ndphillips/YaRrr/comparing-regression-models-with-anova.html
    print("Comparing models with ANOVA")
    pvals = pbsapply(geneNames, function(g) {return(anova(simple.fit.l[[g]], cosine.fit.l[[g]])$"Pr(>Chisq)"[2])})
    ## Choosing significant genes for the cosine model
    alpha = 0.05
    new_alpha = alpha/length(pvals) # Bonferroni correction
    significant <- which(pvals < new_alpha)
    seasonality_genes <- geneNames[significant]
    ## Create a prediction with the significant genes
    time.v=seq(0,2*pi,by=.1) # Time vector
    predict.df <- do.call("rbind",
                          lapply(seasonality_genes,
                                 function(gene) data.frame(JDay  = (364*time.v/(2*pi))+1, # Julian day 1-365 https://en.wikipedia.org/wiki/Julian_day
                                                           I_hat = predict_log2I(cosine.fit.l[[gene]], a_time = time.v),
                                                           Symbol = gene)))
    predict.df$Symbol <- factor(predict.df$Symbol)
    return(predict.df)
}

guessSW <- function(data.df) {
    ## Given the data for a single gene, guess whether it's a winter or a summer gene. SW: Summer Winter
    ## Please no australian data for the guessing because they live upside down.
    ##
    ## The reason for doing this is:
    ## We have to guess because the reported genes in the supplement material are not labelled.
    ## We will do that by another fancy linear regression

    ## Our time variable is in radians, therefore ranges from 0 to 2*pi.
    ## Use it to define what values have been taken in winter/summer
    jday_to_rad <- function(jday) { return( (jday-1)*2*pi / (365-1) ) } # Julian day to radian form.

    ## Winter is Dec|Jan|Feb ❄️ 
    dec1st <- jday_to_rad(335)
    jan1st <- jday_to_rad(1)
    feb1st <- jday_to_rad(32)
    feb28th <- jday_to_rad(59)
    ## Summer is Jun|Jul|Aug ☀️
    jun1st <- jday_to_rad(152)
    jul1st <- jday_to_rad(184)
    aug1st <- jday_to_rad(213)
    aug31st <- jday_to_rad(243)

    ## cosine model for all the genes
    control = lme4::lmerControl(check.conv.singular = "ignore") # otherwise we get lots of warnings
    fit <- lme4::lmer(rle ~ 1 + (cos(t) + sin(t) | Symbol) + (t | GeoLocation) + (t | Gender) + (1 | Donor), data = data.df, control = control )
    ##   Here is the difference:__________________________

    ## Projecting the data into some selected dates
    sdates <- c(dec1st, jan1st, feb1st, feb28th, # This will do for winter
                jun1st, jul1st, aug1st, aug31st) # This for summer
    x <- as.matrix(coef(fit)[["Symbol"]][, c("(Intercept)", "cos(t)", "sin(t)")])            # x = [a b c]
    y <- matrix(c(rep(1, times = length(sdates)),      
                  cos(sdates),                      #      ⸢      1       ...        1     ⸣
                  sin(sdates)),                     # y =  |cos(2 π t₁)   ...  cos(2 π tₘ) |
                nrow=3,                             #      ⸤sin(2 π t₁)   ...  sin(2 π tₘ) ⸥
                byrow=TRUE)
    I_hat <- x %*% y # Matrix multiplication x * y
    rownames(I_hat) <- rownames(coef(fit)[["Symbol"]])
    ## Matrix I_hat should have this shape now:
    ##           date1 date2 date3 date4     => sdates
    ##    Î =  ⸢  .     .     .     .     ⸣ gene 1
    ##         |  .     .     .     .     | gene 2
    ##         ⸤  .     .     .     .     ⸥ gene 3

    ## Traversing the Î matrix now by row, should tell us what genes
    ## have a higher estimated intensity in winter/summer.
    n = length(sdates) # This has to be always even
    differences <- apply(I_hat, 1, function(x) { mean(x[1:n/2]) - mean(x[n/2:n]) }) # Difference in mean intensity winter-summer

    # data frame with our guess
    guess.df <- data.frame(row.names = names(differences),
                           SH = factor(sapply(differences, function(x) ifelse(x > 0, "winter", "summer")), # SH: Season High
                                       levels = c("winter", "summer")))
    return(guess.df)
}


### Asthma PBMc
### E-GEOD-19301

if (! file.exists("asthma.rds")) {
    ## Start from scratch

    ## First, load the Sample and Data Relationship Format (pheno data)
    sdrf.df <- read.csv("ftp://ftp.ebi.ac.uk/biostudies/nfs/E-GEOD-/301/E-GEOD-19301/Files/E-GEOD-19301.sdrf.txt", sep="\t", row.names=1)
    rownames(sdrf.df) <- gsub(" 1$", "", rownames(sdrf.df))
    sdrf.df$Donor <- as.factor(gsub("Donor: ", "", sapply(strsplit(sdrf.df$Comment..Sample_title.,";"), "[", 1)))
    sdrf.df$VisitType <- as.factor(gsub(" Visit Type: ", "", sapply(strsplit(sdrf.df$Comment..Sample_title.,";"), "[", 3)))
    
    ## Read the data directly from ArrayExpress
    ftp.df <- data.frame(URL = scan("E-GEOD-19301-ftplinks.txt", character()))
    row.names(ftp.df) <- gsub("_sample_table.txt$", "", basename(ftp.df$URL))
    ## Paper says: Only asthma patients defined as being in a quiet disease phase were included in our analyses (from 685 to 384 samples)
    select <- row.names(sdrf.df[grep("QUIET", sdrf.df$VisitType),])
    ## Load the selected samples tables into a data frame with probe intensity
    print("Loading asthma microarray data")
    table.l <- pblapply(ftp.df[select, "URL"], read_table_file)
    asthma.df <- do.call("rbind", table.l)
    asthma.df$Sample <- factor(asthma.df$Sample)
    asthma.df$ProbeID <- factor(asthma.df$ProbeID)
    rm(select) ; rm(table.l) ; rm(ftp.df)
    
    ## Map the probes to gene symbols
    ##
    ## probes can map to >= 0 genes. Let's used probes with the highest MAD across samples to select
    ## the representing probe.
    mad.df <- aggregate(Intensity ~ ProbeID, FUN = mad, data = asthma.df)
    mad.df$Symbol <- factor(getSymbols(as.character(mad.df$ProbeID)))
    mad.df <- mad.df[!is.na(mad.df$Symbol), ]
    mad.df <- mad.df[order(mad.df$Intensity, decreasing = TRUE), ]
    selectProbes <- mad.df[!duplicated(mad.df$Symbol), "ProbeID"]
    asthma.df <- asthma.df[asthma.df$ProbeID %in% selectProbes, ]
    asthma.df$ProbeID <- droplevels(asthma.df$ProbeID)
    ## Final and unique probe->symbol (gene) mapping
    asthma.df$Symbol <- factor(getSymbols(as.character(asthma.df$ProbeID)))
    rm(selectProbes) ; rm(mad.df)

    ## Adding some more information
    sampleNames <- as.character(asthma.df$Sample)
    asthma.df$Donor <- factor(sdrf.df[sampleNames, "Donor"])
    countries <- sdrf.df[sampleNames, "Characteristics.country."]
    asthma.df$GeoLocation <- factor(gsub("^GBR$", "GBR/IRL", gsub("^IRL$", "GBR/IRL", countries)))
    ## Date data
    sampleDates <- as.POSIXlt(sdrf.df[sampleNames, "FactorValue..SAMPLE.COLLECTION.DATE."])
    days_in_year <- rep(365, times = length(sampleDates))
    days_in_year[lubridate::leap_year(sampleDates)] <- 366
    asthma.df$t <- 2*pi*(lubridate::yday(sampleDates)-1)/(days_in_year-1) # time variable in radians
    rm(sampleNames) ; rm(countries) ; rm(sampleDates) ; rm(days_in_year) ; rm(sdrf.df)

    ## Calculate rle
    asthma.df$rle <- probeRLE(asthma.df)
    
    ## Guessing donors' gender. Not reported in the microarray metadata
    donors <- levels(asthma.df$Donor)
    gender.df <- data.frame(row.names = donors,
                            gender = sapply(donors, guessGender, data = asthma.df))
    asthma.df$Gender <- factor(gender.df[asthma.df$Donor, "gender"])
    rm(donors) ; rm(gender.df)
    
    ## Save the results
    saveRDS(asthma.df, "asthma.rds", compress = "xz")
    gc() # Free some memory
} else {
    # Read the precomputed data frame.
    asthma.df <- readRDS("asthma.rds")
}

print("Processing GBR|IRL samples")
if (! file.exists("asthma.GBR_IRL.rds")) {
    GBR_IRL.df <- predict_seasonality(asthma.df[asthma.df$GeoLocation == "GBR/IRL", ])
    saveRDS(GBR_IRL.df, "asthma.GBR_IRL.rds", compress = "xz")
} else {
    GBR_IRL.df <- readRDS("asthma.GBR_IRL.rds")
}
print("Processing AUS samples")
if (! file.exists("asthma.AUS.rds")) {
    AUS.df <- predict_seasonality(asthma.df[asthma.df$GeoLocation == "AUS", ])
    saveRDS(AUS.df, "asthma.AUS.rds", compress = "xz")
} else {
    AUS.df <- readRDS("asthma.AUS.rds")
}
print("Processing USA samples")
if (! file.exists("asthma.USA.rds")) {
    USA.df <- predict_seasonality(asthma.df[asthma.df$GeoLocation == "USA", ])
    saveRDS(USA.df, "asthma.USA.rds", compress = "xz")
} else {
    USA.df <- readRDS("asthma.USA.rds")
}

## Figure 3C plots
fig3c_plot <- function(data.df, title) {
    ## Create a ggplot object with some projection data
    
    ## This is a vector of when the months start in julian days
    months_jdays <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
    names(months_jdays) <- substr(month.name, 1, 3)

    ngenes_label = paste0(length(levels(data.df$Symbol)), " genes")
    plot <- ggplot(data.df, aes(x=JDay, y=I_hat)) +
        geom_line(aes(group=Symbol), linewidth = .05) +
        scale_x_continuous(limits = c(1, 365), expand = c(0, 0), breaks = months_jdays) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("Months") +
        ylab("Relative log2 (expression)") +
        ggpp::geom_text_npc(aes(npcx = "right", npcy = "top", label = ngenes_label)) +
        ggtitle(title)
    return(plot)
}

GBR_IRL.plot <- fig3c_plot(GBR_IRL.df, "United Kingdom and Ireland")
AUS.plot <- fig3c_plot(AUS.df, "Australia")
USA.plot <- fig3c_plot(USA.df, "United States")
fig3c.grob <- gridExtra::arrangeGrob(GBR_IRL.plot, AUS.plot, USA.plot,
                                     top=grid::textGrob("Seasonal genes in different populations of adult asthmatic patients"),
                                     nrow = 1)
ggsave("fig3c.png", plot = fig3c.grob, width=40, height=20, units="cm")


## Summer and winter genes
##
## Seasonal genes detected in the BABYDIET data, do they express the same pattern in the ashtma data? (figure 3D)
## We have the list of seasonal genes, but the summer and winter labels are not reported.
## We can try to guess them and see.

## Load the supplemental excel file from the paper. On the "Table 3" page we can find the seasonal genes in the BABYDIET dataset
bbdiet_seasonal <- readxl::read_excel("41467_2015_BFncomms8000_MOESM557_ESM.xls", sheet = "Table 3")
bbdiet_seasonal_genes <- dplyr::pull(bbdiet_seasonal,3)
bbdiet_seasonal_genes <- bbdiet_seasonal_genes[bbdiet_seasonal_genes != "NA" & bbdiet_seasonal_genes != "GENE SYMBOL" & !is.na(bbdiet_seasonal_genes)]
bbdiet_seasonal_genes <- toupper(bbdiet_seasonal_genes)
bbdiet_seasonal_genes <- unique(bbdiet_seasonal_genes) # 5087 Gene symbols

if (! file.exists("index.SW.rds")) {
    asthma.SW.df <- asthma.df[asthma.df$Symbol %in% bbdiet_seasonal_genes, ] # SW: Summer Winter
    asthma.SW.df$Symbol <- droplevels(asthma.SW.df$Symbol) # 3446 Gene Symbols
    ## TODO: We lost 1641 symbols, what happened to those? Is there a better way of finding them?
    index.SW.df <- guessSW(asthma.SW.df[asthma.SW.df$GeoLocation != "AUS", ])
    saveRDS(index.SW.df, "index.SW.rds", compress = "xz")
} else {
    index.SW.df <- readRDS("index.SW.rds")
}

## Fit a cosine model to the data
print("Fitting a cosine model in the asthma data for the summer/winter genes discovered in the babydiet data")
predict_sw <- function(data.df) {
    model.l <- pbsapply(levels(data.df$Symbol), cosinef, data = data.df)
    time.v=seq(0,2*pi,by=.1) # Time vector
    predict.df <- do.call("rbind",
                          lapply(names(model.l),
                                 function(gene) data.frame(JDay  = (364*time.v/(2*pi))+1,
                                                           I_hat = predict_log2I(model.l[[gene]],
                                                                                 a_time = time.v),
                                                           Symbol = gene)))
    return(predict.df)
}

print("GBR/IRL")
GBR_IRL.SW.df <- predict_sw(asthma.SW.df[asthma.SW.df$GeoLocation == "GBR/IRL", ])
GBR_IRL.SW.df$SH <- index.SW.df[GBR_IRL.SW.df$Symbol, "SH"]
print("AUS")
AUS.SW.df <- predict_sw(asthma.SW.df[asthma.SW.df$GeoLocation == "AUS", ])
AUS.SW.df$SH <- index.SW.df[AUS.SW.df$Symbol, "SH"]
print("USA")
USA.SW.df <- predict_sw(asthma.SW.df[asthma.SW.df$GeoLocation == "USA", ])
USA.SW.df$SH <- index.SW.df[USA.SW.df$Symbol, "SH"]

## Figure 3D plots
fig3d_plot <- function(data.df, title) {
    ## Create a ggplot object with some projection data
    
    ## This is a vector of when the months start in julian days
    months_jdays <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
    names(months_jdays) <- substr(month.name, 1, 3)

    plot <- ggplot(data.df, aes(x=JDay, y=I_hat)) +
        geom_line(aes(group=Symbol, colour = SH), linewidth = .05) +
        scale_x_continuous(limits = c(1, 365), expand = c(0, 0), breaks = months_jdays) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("Months") +
        ylab("Relative log2 (expression)") +
        scale_color_manual(values=c("green", "blue")) + # green for winter, blue for summer
        theme(legend.position = "none") +
        ggtitle(title)
    return(plot)
}

GBR_IRL.SW.plot <- fig3d_plot(GBR_IRL.SW.df, "United Kingdom and Ireland")
AUS.SW.plot <- fig3d_plot(AUS.SW.df, "Australia")
USA.SW.plot <- fig3d_plot(USA.SW.df, "United States")

fig3d.grob <- gridExtra::arrangeGrob(GBR_IRL.SW.plot, AUS.SW.plot, USA.SW.plot,
                                     top=gridtext::richtext_grob("<span style='color:blue'>Summer</span> and <span style='color:green'>winter</span> genes from BABYDIET in the adult asthmatic patients"),
                                    nrow = 1)
ggsave("fig3d.png", plot = fig3d.grob, width=40, height=20, units="cm")



## Correlation of inmune lineages characteristic genes with
## PBMCs seasonal genes (and no-seasonal)
##
## The paper doesn't clearly explain what PBMC data they used
## to calculate this correlation, but I bet it was the T1D data.
## Since, once again, we are missing the date data on T1D and BabyDiet data,
## we are going to do it with the asthma data.

## This list of 13 genes come from figure 4a (which we are trying to replicate)
##
## There is a problem, the SBK1 gene is not measured by any probe that I could find
## Alternatevely, we can use this tool to find a coexpressed gene in the inmune tissue:
## https://gccri.bishop-lab.uthscsa.edu/shiny/correlation-analyzer/
## And the gene with the highest correlation to SBK1 is CACNA1I.
## TBH, I am not completely confident that this gene swap is legit
## Any inmunologist in the room?
inmune_genes <- c("CD1E", "CRIP2", "CTLA4", "MS4A1", "RSAD2", "CACNA1I", "SRRM2", # SBK1 has been replace with CACNA1I
                  "CLIC3", "FCN1", "PDCD1LG2", "PEG10", "SPARC", "TNFRSF10C")
inmune_genes_labels <- c("DCs", "CD4+ T", "Activated CD4+ T", "B cells", "Activated NKs/DCs", "Activated NKs", "CD4+ T cells",
                         "NK cells",  "Monocytes Neutrophils", "Activated DC", "Activated B cells", "Macrophages", "Neutrophils")
names(inmune_genes_labels) <- inmune_genes

## Calculate the correlation matrix, this step takes quite some time. Better save the result.
print("Calculating the correlation of inmune lineages characteristics genes with the seasonal genes")
print("This step takes a lot of time unless it's saved already in disk"
if (! file.exists("asthma.cor.rds")) {
    asthma.cor.m <- HiClimR::fastCor(as.matrix(reshape2::dcast(asthma.df, Sample ~ Symbol, value.var = "rle")[,-1]))
    saveRDS(asthma.cor.m, "asthma.cor.rds", compress = "xz")
} else {
    asthma.cor.m <- readRDS("asthma.cor.rds")
}    
print("done!")
asthma.cor.df <- do.call("rbind",
                         lapply(inmune_genes,
                                function(ig) {
                                    cor.v <- asthma.cor.m[ig, -which(colnames(asthma.cor.m) == ig)]
                                    data.frame(InmuneGene = ig,
                                               Cor = cor.v,
                                               Symbol = names(cor.v))
                                }))
asthma.cor.df$InmuneGene <- factor(asthma.cor.df$InmuneGene, levels=inmune_genes)
seasonFlags <- asthma.cor.df$Symbol %in% bbdiet_seasonal_genes
asthma.cor.df$Seasonality <- factor(ifelse(seasonFlags, "Seasonal", "Non-Seasonal"), levels = c("Non-Seasonal", "Seasonal"))

## Let's divide the inmune genes in two groups, so we can separate the plot in two
group1 <- inmune_genes[1:7]
group2 <- inmune_genes[8:13]
graphLabels <- data.frame(InmuneGene = as.vector(group1),
                          Label = as.vector(inmune_genes_labels[group1]))
cor.1.plot <- ggplot(asthma.cor.df[asthma.cor.df$InmuneGene %in% group1, ], aes(x = Cor) ) +
    geom_histogram(aes(y = after_stat(density)), fill="black", bins=30) +
    facet_grid(InmuneGene ~ Seasonality) +
    geom_text(data = graphLabels, aes(label=Label, x = 0, y = 4)) +
    xlab("Spearman correlation") +
    ylab("Density")
graphLabels <- data.frame(InmuneGene = as.vector(group2),
                          Label = as.vector(inmune_genes_labels[group2]))
cor.2.plot <- ggplot(asthma.cor.df[asthma.cor.df$InmuneGene %in% group2, ], aes(x = Cor) ) +
    geom_histogram(aes(y = after_stat(density)), fill="black", bins=30) +
    facet_grid(InmuneGene ~ Seasonality) +
    geom_text(data = graphLabels, aes(label=Label, x = 0, y = 4)) +
    xlab("Spearman correlation") +
    ylab("Density")
fig4a.plot <- ggpubr::ggarrange(cor.1.plot, cor.2.plot, ncol=2)
ggsave("fig4a.png", plot = fig4a.plot, width=40, height=20, units="cm") # Unfortunately, opening the image you can see that there is no correlation in this data

