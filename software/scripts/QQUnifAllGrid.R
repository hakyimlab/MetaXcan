library(argparse)
source("Plots.R")

parser <- ArgumentParser(description="Plot and compare different population's MetaXcan against different population's predixcan")
parser$add_argument('--input_pattern',
                    help='input pattern for files',
                    default='results_igrowth/beta_%s_ref_%s_predixcan_zscores.csv')

parser$add_argument('--title',
                    help='plot title',
                    default='Intrinsic Growth phenotype')

parser$add_argument('--output',
                    help='path of predixcan results, without extension',
                    default='results_igrowth/qq_igrowth_grid.png')

arguments <- parser$parse_args(commandArgs(TRUE))

path = arguments$input_pattern

STUDY = c("AFR", "EUR", "EAS", "SEL")
REFERENCE = c("AFR", "EUR", "EAS", "SEL")

STUDY_DISPLAY = c("AFR", "EUR", "EAS", "AFR+EUR+EAS")
REFERENCE_DISPLAY = c("AFR", "EUR", "EAS", "AFR+EUR+EAS")

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
            the_facet <- sprintf("%s study on %s summary statistics", s_d, r_d)

            p_val <- 2*pnorm(-abs(file_data$zscore))
            y <- -sort(log10(p_val))
            nn <- length(y)
            x <- -log10((1:nn)/(nn+1))

            c95 <- rep(0,nn)
            c05 <- rep(0,nn)
            for(k in 1:nn)
            {
                c95[k] <- qbeta(0.95,k,nn-k+1)
                c05[k] <- qbeta(0.05,k,nn-k+1)
            }
            y95 <- -log10(c95)
            y05 <- -log10(c05)
            b <- -log10(0.05/nn) #bonferroni

            data <- rbind(data, data.frame(x=x, y=y, y95=y95, y05=y05, b=b, the_facet=the_facet))
        }
    }
    return (data)
}

data <- build_data()

title <- arguments$title
p <- qq_unif_faceted(data, title, 4)

output = arguments$output
png(filename=output,width=1024,height=1280)
print(p)
dev.off()