#' Function override in of microbiomeseq to calculate pair-wise anovas allowing more options and Welch-ANOVA
#'
#' @param df dataframe
#' @param meta_table the metadata df
#' @param grouping_column the actual column to test variables of
#' @param pValueCutoff as it says (default 0.05)
#' @param anova.type Defaults to classic ANOVA ("anova"). Option to use Welch-ANOVA ("welch").
#'
#' @returns list
#'
#' @import rstatix
#' @export
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' @author Additional credits go the the original function perform_anova_diversity from [MicrobiomeSeq](https://github.com/umerijaz/microbiomeSeq) by Alfred Ssekagiri which I modified and extended with some convenience options (mainly) for myself.
#' @note For Welch-ANOVA we use the [rstatix package](https://www.rdocumentation.org/packages/rstatix).
#' 
#' Updates:
#'          2022-08-02 Alx added Welch-ANOVA option
#'          2019-11-01 Alx initial checkin
#' 

perform_anova <-  function (df, meta_table, grouping_column, pValueCutoff=0.05, anova.type="anova" )
  {
  require(rstatix) #package is required for the welch anova option.
    
  if( ! anova.type %in% c("anova","welch") ) {
    warning("Unsupported anova method specified. Should be one of 'anova' or 'welch'. Defaulting to classic ANOVA.\n")
  }  
  
  df$.group. <- meta_table[, grouping_column]  # Add grouping column to df

  # Perform actual pair-wise ANOVA or Welch-ANOVA
  if (anova.type == "welch") {
    cat("Performing Welch-ANOVA\n")
    pval <- aggregate(value ~ measure, data = df, function(x) {
      sprintf("%.2g", tryCatch(
        summary(welch_anova_test(x ~ df$.group.))$p[1],
        error = function(e) NULL
      ))
    })
  } else {
    cat("Performing classic ANOVA\n")
    pval <- aggregate(value ~ measure, data = df, function(x) {
      sprintf("%.2g", tryCatch(
        summary(aov(x ~ df$.group.))[[1]][["Pr(>F)"]][1],
        error = function(e) NULL
      ))
    })
  }  
  
  pval <- pval[!pval$pvalue == "",]
  pval <- pval[as.numeric(pval$pvalue) <= pValueCutoff,]
  pval$pvalue <- sapply(as.numeric(pval$pvalue), function(x) {
    as.character(cut(
      x,
      breaks = c(-Inf, 0.001, 0.01, pValueCutoff,
                 Inf),
      label = c("***", "**", "*", "")
    ))
  })
  
    if (length(unique(as.character(meta_table[, grouping_column]))) >
        2) {
      df$measure <- as.character(df$measure)
      if (dim(pval)[1] > 0) {
        for (i in seq(1:dim(pval)[1])) {
          df[df$measure == as.character(pval[i, measure]),
             "measure"] = paste(as.character(pval[i, measure]),
                                as.character(pval[i, pvalue]))
        }
      }
      df$measure <- as.factor(df$measure)
    }
    s <- combn(unique(as.character(df[, grouping_column])),
               2)
    df_pw <- NULL
    for (k in unique(as.character(df$measure))) {
      bas <- max(df[(df$measure == k), "value"])
      #use max 10% of barplot added for error bars if 10 bars are present
      inc <- 0.1 * 0.1 * bas
      bas <- bas + inc
      for (l in 1:dim(s)[2]) {
        tmp <- c(k, s[1, l], s[2, l], bas, paste(sprintf("%.2g",
                                                         tryCatch(
                                                           summary(aov(as.formula(
                                                             paste("value ~",
                                                                   grouping_column)
                                                           ), data = df[(df$measure ==
                                                                           k) &
                                                                          (df[, grouping_column] == s[1, l] | df[,
                                                                                                                 grouping_column] == s[2, l]),]))[[1]][["Pr(>F)"]][1],
                                                           error = function(e)
                                                             NULL
                                                         )), "", sep = ""))
        if (!is.na(as.numeric(tmp[length(tmp)]))) {
          if (as.numeric(tmp[length(tmp)]) < pValueCutoff) {
            if (is.null(df_pw)) {
              df_pw <- tmp
            }
            else {
              df_pw <- rbind(df_pw, tmp)
            }
            bas <- bas + inc
          }
        }
      }
    }
    if (!is.null(df_pw)) {
      df_pw <- data.frame(row.names = NULL, df_pw)
      names(df_pw) <- c("measure", "from", "to", "y", "p")
    }
    out <- list(df_pw = df_pw, df = df)
    return(out)
  }

