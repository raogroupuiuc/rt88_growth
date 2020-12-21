summarizeFit <- function(fit, calcFC = TRUE, method = "separate", adjust.method = "BH",
                         removeNames = c("Block", "Row", "Column", "ControlType", "Status"),
                         removeCoefs = NULL, addAnova = FALSE) {

#extract data in MArrayLM object to a data frame
fit <- as.data.frame(fit)

#determine which columns are coefficients, p-values and which ones 
#will be the adjusted p-values
coefs <- grep("coefficients", names(fit))
pvals <- grep("^p.value", names(fit))

if (!is.null(removeCoefs)){
    coefs <- coefs[-removeCoefs]
    pvals <- pvals[-removeCoefs]
}

if(addAnova) {
  coefs <- c(coefs, which(names(fit) == "F"))
  pvals <- c(pvals, which(names(fit) == "F.p.value"))
}

fdrvals <- (ncol(fit) + 1):(ncol(fit)+length(pvals))

#Back-calculate fold-changes, if desired
if (calcFC) {
  if (addAnova) {
    fit[,coefs[-length(coefs)]] <- 2^abs(fit[,coefs[-length(coefs)]]) * sign(fit[,coefs[-length(coefs)]])
  } else {
    fit[,coefs] <- 2^abs(fit[,coefs]) * sign(fit[,coefs])
  }
  
  names(fit) <- gsub("coefficients", "FC", names(fit))
}
else
    names(fit) <- gsub("coefficients", "logFC", names(fit))



#Do p-value adjustment - code adapted from decidedTests for separate and global

method <- match.arg(method, c("separate", "global"))
adjust.method <- match.arg(adjust.method, c("none", "bonferroni", 
                                            "holm", "BH", "fdr", "BY"))
if (adjust.method == "fdr") 
    adjust.method <- "BH"
switch(method, separate = {
    p <- fit[,pvals]
    for (j in 1:ncol(p)) {
        o <- !is.na(p[, j])
        p[o, j] <- p.adjust(p[o, j], method = adjust.method)
    }
    fit <- cbind(fit, p)
}, global = {
    p <- fit[,pvals]
    o <- !is.na(p)
    p[o] <- p.adjust(p[o], method = adjust.method)
    fit <- cbind(fit, p)
})

names(fit)[pvals] <- gsub("p.value", "rawP", names(fit)[pvals])
names(fit)[fdrvals] <- gsub("p.value", "FDR", names(fit)[fdrvals])

#Re-organize results and only include the columns starting with "genes",
# the Amean, the fold-changes, raw-pval and adjust p-values

out <- fit[,c(grep("^genes", names(fit)), grep("^Amean$", names(fit)), as.vector(t(matrix(c(coefs,pvals, fdrvals),length(coefs),3))))]

#Remove "genes" from beginning of names and remove the desired columns
names(out) <- gsub("genes.","", names(out), fixed = TRUE)
out <- out[, !names(out) %in% removeNames]

return(out)
}