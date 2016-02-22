#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(argparse)

parser <- ArgumentParser(description="Plot and compare different population's MetaXcan against different population's predixcan")
parser$add_argument('--input_pattern',
                    help='input pattern for files',
                    default='resultsp/beta_%s_ref_%s_predixcan_zscores.csv')

parser$add_argument('--output',
                    help='path of predixcan results, without extension',
                    default='resultsp/simulated_grid.png')

arguments <- parser$parse_args(commandArgs(TRUE))

#STUDY = c("AFR", "EUR", "EAS", "SEL")
#REFERENCE = c("AFR", "EUR", "EAS", "SEL")
#
#STUDY_DISPLAY = c("AFR", "EUR", "EAS", "AFR+EUR+EAS")
#REFERENCE_DISPLAY = c("AFR", "EUR", "EAS", "AFR+EUR+EAS")

STUDY = c("AFR", "EUR", "EAS")
REFERENCE = c("AFR", "EUR", "EAS")
#
STUDY_DISPLAY = c("AFR", "EUR", "EAS")
REFERENCE_DISPLAY = c("AFR", "EUR", "EAS")


zscore_predixcan_grid_plot_object <- function(data) {
    support <- data %>% distinct(the_facet) %>% select(lin_eq, the_facet)


    p <- ggplot(data, aes(x=zscore, y=predixcan_result)) +
        theme(strip.text = element_text(size=20)) +
        theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) +
        geom_point(size = 2) +
        labs(y = "PrediXcan Association Result for Study Population", x = "MetaXcan Association Result for Study Population on Reference Population") +
        #scale_colour_gradient(limits=c(0,100), name="SNPs in\neach gene", low="blue", high="red", na.value="green") +
        facet_wrap(~the_facet,scales="fixed",ncol=length(REFERENCE)) +
        geom_abline(intercept=0, slope=1, colour="#335533") +
        geom_text(data=support, aes(x = -3, y = 2, label = lin_eq), parse = TRUE, colour = "#335533")
        #geom_text(aes(x = -3, y = 2, label = sprintf("%s",lin_eq)), parse = TRUE, colour = "#335533")

    return(p)
}

path = arguments$input_pattern

build_data <- function() {
    data <- data.frame()

    for (i in  1:length(STUDY)) {
        for (j in 1:length(REFERENCE)) {
            s <- STUDY[[i]]
            r <- REFERENCE[[j]]

            file_name <- sprintf(path, s, r)
            print(sprintf("loading %s", file_name))
            file_data <- read.csv(file_name)

            s_d <- STUDY_DISPLAY[[i]]
            r_d <- REFERENCE_DISPLAY[[j]]
            the_facet <- sprintf("%s GWAS study on %s reference", s_d, r_d)
            file_data$the_facet <- c(the_facet)

            maximum =max(file_data$n_snp, na.rm = TRUE)
            cut = 40
            file_data$snp_scale <- ifelse(file_data$n_snp < cut,
                as.integer(file_data$n_snp * 50.0 /cut),
                as.integer(50 + (file_data$n_snp - cut)*50.0/(maximum-cut)) )

            lm_eqn <- function(df){
                m <- lm(predixcan_result ~ zscore, df);
                eq <- substitute(italic(r)^2~"="~r2,
                    list(r2 = format(summary(m)$r.squared, digits = 3)))
                t <- as.character(as.expression(eq));
                return(t)
            }
            file_data$lin_eq <- lm_eqn(file_data)

            data <- rbind(data, file_data)
        }
    }
    return (data)
}

data <- build_data()
plot <- zscore_predixcan_grid_plot_object(data)

output = arguments$output
print("Saving image...")
png(filename=output, width=1200, height=1200)
print(plot)
dev.off()