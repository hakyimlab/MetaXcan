library(ggplot2)

beta_plot <- function(z, z_q, file_prefix) {
    lm_eqn <- function(x,y){
        m <- lm(y ~ x);
        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
             list(a = format(coef(m)[1], digits = 2),
                  b = format(coef(m)[2], digits = 2),
                r2 = format(summary(m)$r.squared, digits = 3)))
        as.character(as.expression(eq));
    }

    text = lm_eqn(z, z_q)

    df <- data.frame(z=z, z_q=z_q)
    p <- ggplot(df, aes(x=z, y=z_q)) +
        geom_point(size = 1) +
        stat_smooth(method = "lm", col = "#335533") +
        geom_text(x = -0.1, y = 5, label = text, parse = TRUE, colour = "#335533")

  	image <- paste(file_prefix, "-scatter.png", sep="")
	png(filename=image,width=1024,height=768)
	print(p)
	dev.off()
}

beta_z_file = gzfile("intermediate/beta_z/T1D_chr1.assoc.dosage.gz")
beta_z = read.delim(beta_z_file, sep=" ")

beta_z_q_file = gzfile("intermediate/beta_z_q/T1D_chr1.assoc.dosage.gz")
beta_z_q = read.delim(beta_z_q_file, sep=" ")

beta_plot(beta_z$beta_z, beta_z_q$beta_z, "results/beta")
