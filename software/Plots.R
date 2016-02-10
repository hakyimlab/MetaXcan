library(ggplot2)
library(qqman)

process_zscore_file <- function(file_prefix) {
	file <- paste(file_prefix, ".csv", sep="")
	data <- read.csv(file)
    plot_values(file_prefix, data)
    plot_histogram(file_prefix, data)
}

plot_values <- function(file_prefix, data) {
	ordered_data <- data[order(data$zscore),]

    max_n = max(ordered_data$n, na.rm = TRUE)

	sel <- data.frame(Index = seq.int(nrow(ordered_data)), Values = ordered_data$zscore)
	sel$s_min <- ifelse(ordered_data$zscore > 0, 0, -(ordered_data$n*1.0/max_n))
	sel$s_max <- ifelse(ordered_data$zscore > 0, (ordered_data$n*1.0/max_n), 0)

	p1<-ggplot(sel,aes(x=Index,y=Values, ymin = s_min, ymax=s_max) ) +
				geom_pointrange(col='gray',alpha=0.7)+
				geom_point(aes(colour='black'))+
				coord_cartesian(ylim = c(-5, 5))

	image <- paste(file_prefix, ".png", sep="")
	png(filename=image,width=1024,height=768)
	print(p1)
	dev.off()

	return()
}

plot_histogram <- function(file_prefix, data) {
	image <- paste(file_prefix, "-histogram.png", sep="")
	png(filename=image,width=1024,height=768)
    hist(data$zscore, breaks=80)
	dev.off()

	return()
}

process_zscore_predixcan_file <- function(file_prefix) {
    file <- paste(file_prefix, ".csv", sep="")
    data <- read.csv(file)
    zscore_predixcan_plot(data, file_prefix)
    do_manhattan_plot_from_data(data, file_prefix)

    valid = c("VAR_g", "se_predixcan_beta") %in% names(data)
    if (all(valid)) {
        se_beta_error_proxy_plot(data, file_prefix)
    }
}

zscore_predixcan_plot_object <- function(data) {
    maximum =max(data$n_snp, na.rm = TRUE)
    cut = 40
    data$transform <- ifelse(data$n_snp < cut,
                as.integer(data$n_snp * 50.0 /cut),
                as.integer(50 + (data$n_snp - cut)*50.0/(maximum-cut)) )


    dummy = data[1,]

    lm_eqn <- function(df){
        m <- lm(predixcan_result ~ zscore, df);
        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
             list(a = format(coef(m)[1], digits = 2),
                  b = format(coef(m)[2], digits = 2),
                r2 = format(summary(m)$r.squared, digits = 3)))
        as.character(as.expression(eq));
    }

    text = lm_eqn(data)

    p <- ggplot(data, aes(x=zscore, y=predixcan_result, colour = transform)) +
        labs(y = "PrediXcan Association Result", x = "MetaXcan Association Result") +
        geom_abline(intercept=0, slope=1, colour="grey49") +
        geom_point(size = 1) +
        scale_colour_gradient(limits=c(0,100), low="blue", high="red", na.value="green") +
        geom_text(x = -2, y = 1, label = text, parse = TRUE, colour = "#335533", size = 8)

    return(p)
}

zscore_predixcan_plot <- function(data, file_prefix) {
    p <- zscore_predixcan_plot_object(data)

  	image <- paste(file_prefix, "-scatter.png", sep="")
	png(filename=image,width=1024,height=768)
	print(p)
	dev.off()
}

export_data <- function(input) {
    data <- read.csv(input)
    data$model_n[is.na(data$model_n)] <- 0
    data$zscore[is.infinite(data$zscore)] <- NA
    data <- data[complete.cases(data),]
    #data$zscore[is.na(data$zscore)] <- 0
    return(data)
}

do_manhattan_plot_from_data <- function(data, file_prefix) {
    c <- data[complete.cases(data),]
    if(!("pval" %in% colnames(c)))
    {
        c$pval <- 2*pnorm(-abs(c$zscore))
    }

    nn = length(data$zscore)
    p_b = 0.05/nn

    image <- paste(file_prefix, "-manhattan.png", sep="")
    png(filename=image,width=1024,height=768)
    manhattan(c, chr="chr", bp="base_position", p="pval", snp="gene_name",
            suggestiveline=-log10(p_b), genomewideline=FALSE, annotatePval=p_b, annotateTop=FALSE)
#    manhattan(c, chr="chr", bp="base_position", p="pval", snp="gene", annotatePval=1e-5, annotateTop=FALSE)
    dev.off()
}

se_beta_error_proxy_plot <- function(data, file_prefix) {
    data$true_value <- 1/data$se_predixcan_beta
    data$proxy <- sqrt(data$VAR_g)

    dummy = data[1,]

    lm_eqn <- function(df){
        m <- lm(true_value ~ proxy, df);
        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
             list(a = format(coef(m)[1], digits = 2),
                  b = format(coef(m)[2], digits = 2),
                r2 = format(summary(m)$r.squared, digits = 3)))
        as.character(as.expression(eq));
    }

    text = lm_eqn(data)

    p <- ggplot(data, aes(x=proxy, y=true_value), colour="#006699") +
        geom_point(size = 1) +
        stat_smooth(method = "lm", col = "red") +
        geom_text(x = 0.2, y = 60, label = text, parse = TRUE)


  	image <- paste(file_prefix, "-se_beta_var_g.png", sep="")
	png(filename=image,width=1024,height=768)
	print(p)
	dev.off()
}

zscore_qqunif <- function(input, output) {
    df <- read.csv(input)
    zscore_qqunif_from_data(df, output)
}

zscore_qqunif_from_data <- function(data, output) {
    df <- data
    if(!("pval" %in% colnames(df)))
    {
        df$pval <- 2*pnorm(-abs(df$zscore))
    }
    qqunif(output, p=df$pval)
}

qqunif <- function(file = NULL,p,BH=T,CI=T,...)
{
  #p <- p[!is.na(p)]
  #p <- p[!is.nan(p)]
  #nn = length(p)
  #xx =  -log10((1:nn)/(nn+1))
  yy <- -sort(log10(p))
  nn <- length(yy)
  xx <-  -log10((1:nn)/(nn+1))

  if (! is.null(file)) {
    png(filename=file,width=1024,height=768)
  }
  plot( xx,  yy,
    xlab=expression(Expected~~-log[10](italic(p))),
    ylab=expression(Observed~~-log[10](italic(p))),
    cex.lab=1.4,mgp=c(2,1,0),
    ... )
  abline(0,1,col='gray')
  if(BH)
  {
    abline(-log10(0.05),1, col='red',lty=1)
    abline(-log10(0.10),1, col='orange',lty=2)
    abline(-log10(0.25),1, col='yellow',lty=3)
    legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
             col=c('red','orange','yellow'),lty=1:3, cex=1)
    abline(h=-log10(0.05/nn)) ## bonferroni
  }
  if(CI)
  {
    ## create the confidence intervals
    c95 <- rep(0,nn)
    c05 <- rep(0,nn)
    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)
    ## this portion was posted by anonymous on
    ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html

    for(i in 1:nn)
    {
      c95[i] <- qbeta(0.95,i,nn-i+1)
      c05[i] <- qbeta(0.05,i,nn-i+1)
    }

    lines(xx,-log10(c95),col='gray')
    lines(xx,-log10(c05),col='gray')

    if (! is.null(file)) {
        dev.off()
    }
  }
}

qq_unif_faceted <- function(df, title=NULL, columns=2) {
    p <- ggplot( df, aes(x=x, y=y)) +
        theme_bw() +
        theme(plot.title = element_text(lineheight=3, face="bold")) +
        xlab(expression(Expected~~-log[10](italic(p)))) +
        ylab(expression(Observed~~-log[10](italic(p)))) +
        geom_point(colour = 'black', shape=1) +
        facet_wrap(~the_facet,scales="fixed",ncol=columns)

    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }

    p <- p + geom_abline(intercept=0, slope=1, colour='gray')

    p <- p +
        geom_abline(intercept=-log10(0.05), slope=1, aes(colour='FDR = 0.05'), linetype=1, show_guide = TRUE) +
        geom_abline(intercept=-log10(0.10), slope=1, aes(colour='FDR = 0.10'),linetype=2, show_guide = TRUE) +
        geom_abline(intercept=-log10(0.25), slope=1, aes(colour='FDR = 0.25'),linetype=3, show_guide = TRUE) +
        scale_colour_manual("",
                      labels = c('FDR = 0.05', 'FDR = 0.10', 'FDR = 0.25'),
                      breaks = c('FDR = 0.05', 'FDR = 0.10', 'FDR = 0.25'),
                      values = c("red", "orange", "yellow")) +
        theme(legend.position="bottom")
  #scale_colour_manual(name='', values=c('FDR = 0.05'='red', 'FDR = 0.10'='orange', 'FDR = 0.25'='yellow'))
#legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
#             col=c('red','orange','yellow'),lty=1:3, cex=1)
    p <- p + geom_abline(data=df,aes(intercept=b, slope=0), colour='black') ## bonferroni
#
    p <- p + geom_line(data=df, aes(x=x, y=y95), colour='gray')
    p <- p + geom_line(data=df, aes(x=x, y=y05), colour='gray')
    return(p)
}

build_qq_unif_grid <- function(files, facets, title, columns, output) {
    df <- data.frame(y=numeric(0), x=numeric(0), y95=numeric(), y05=numeric(), b=numeric(0), the_facet=character(0))
    for (i in 1:length(files)) {
        file = files[[i]]
        data <- read.csv(file)

        the_facet <- facets[[i]]
        p_val <- 2*pnorm(-abs(data$zscore))
        y <- -sort(log10(p_val))
        nn <- length(y)
        x <- -log10((1:nn)/(nn+1))

        c95 <- rep(0,nn)
        c05 <- rep(0,nn)
        for(j in 1:nn)
        {
            c95[j] <- qbeta(0.95,j,nn-j+1)
            c05[j] <- qbeta(0.05,j,nn-j+1)
        }
        y95 <- -log10(c95)
        y05 <- -log10(c05)
        b <- -log10(0.05/nn) #bonferroni

        df <- rbind(df, data.frame(x=x, y=y, y95=y95, y05=y05, b=b, the_facet=the_facet))
    }

    p <- qq_unif_faceted(df, title=title, columns=columns)
    #png(filename="results/igrowth_qqunif.png",width=768,height=800)
    png(filename=output,width=768,height=800)
    print(p)
    dev.off()
}
