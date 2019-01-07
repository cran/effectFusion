
plotHPD <- function(mcmc, alpha = 0.05, title, labels) {
    
    mcmc <- as.matrix(mcmc)
    k <- ncol(mcmc)
    
    if (k > 1) {
        interval <- t(apply(mcmc, 2, function(x) hpdMCMC(x, alpha)))
        pm <- colMeans(mcmc)
        c <- 2:(k + 1)
        l <- interval[, 1]
        u <- interval[, 2]
        set <- data.frame(pm, c, l, u)
        bar_width <- 0.05
    } else {
        interval <- hpdMCMC(mcmc, alpha)
        pm <- mean(mcmc)
        c <- 2
        l <- interval[1]
        u <- interval[2]
        set <- data.frame(pm, c, l, u)
        bar_width <- 0
    }
    
    
    gg_style <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "grey10"), 
        panel.grid.major = ggplot2::element_line(colour = "grey90"), panel.grid.minor = ggplot2::element_line(colour = "grey90"), 
        axis.text = ggplot2::element_text(size = 10), axis.title.x = ggplot2::element_blank(), 
        axis.title.y = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 10, 
            margin = ggplot2::margin(10, 0, 10, 0)), axis.text.x=ggplot2::element_text(color = "black", size=8, angle=30, vjust=.8, hjust=0.8), 
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.9))
    
    
    p <- ggplot2::ggplot(data = set, ggplot2::aes(y = pm, x = c)) + ggplot2::geom_errorbar(data = set, 
        ggplot2::aes(ymin = l, ymax = u), width = bar_width, colour = "cyan3", size = 0.9) + 
        ggplot2::geom_point(colour = "darkcyan", size = 2.5) + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        colour = "darkgrey") + ggplot2::ggtitle(title) + ggplot2::scale_x_continuous(labels = labels, breaks = 2:(k+1)) + gg_style
    
    return(p)
    
}

