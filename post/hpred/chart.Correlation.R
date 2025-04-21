chart.Correlation.nostars <- function (R, histogram = TRUE, method = c("pearson", "kendall", 
                                                                       "spearman"), a=F, ...) 
{
  x = checkData(R, method = "matrix")
  if (missing(method)) method = method[1]

  
  panel.cor <- function(x, y, digits = 2, prefix = "r=", 
                        use = "pairwise.complete.obs", method = "pearson", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")

      
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    # Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
    #                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
    #                                                                           "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0) / 1)
    # text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  
  f <- function(t) dnorm(t, mean = mean(x), sd = sd.xts(x))
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, 
         axes = FALSE, main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 2)
    rug(x)
  }
  # p,ot dots, reference line, and smoothed line
  my_line <- function(x,y,...) {
    points(x,y, pch=1, cex= .25, ...)
    abline(a = 0,b = 1, lwd=2, ...) 
    lines(lowess(x,y), lwd=2, col='red') }
    
  if (histogram) 
    pairs(x, gap = 0, lower.panel = my_line, upper.panel = panel.cor, 
          diag.panel = hist.panel)
  else pairs(x, gap = 0, lower.panel = my_line, upper.panel = panel.cor)
}
