############################################################
##### General function for Table 1 Statistics          #####
############################################################

# Based on data.table syntax
# Requires input of a vector of continuous variables names (characters), 
#                   a vector of categorical variables names (characters),
#                   name of the data table (character),
#                   name of the categorical variable that you want to test the difference by (character)

table1 <- function(cnt_vars, cat_vars, dat, by_var, sig_digit = 3) {
        
        print(str_c("Summarize continuous variables in all subjects (dataset: ", dat, ")"))
        cnt_sum <- vector("list", length(cnt_vars))
        for (c in 1:length(cnt_vars)) {
                
                eval(parse(text = str_c("mean <- ", dat, "[, mean(", cnt_vars[c], ", na.rm = T)]")))
                eval(parse(text = str_c("sd <- ", dat, "[, sd(", cnt_vars[c], ", na.rm = T)]")))
                eval(parse(text = str_c("na.count <- ", dat, "[, sum(is.na(", cnt_vars[c], "))]")))
                tmp <- as.data.table(cbind(mean, sd, na.count))
                
                names(cnt_sum)[c] <- cnt_vars[c]
                cnt_sum[[c]] <- tmp
        }
        print(cnt_sum)
        
        cat("\n")
        
        print(str_c("Summarize categorical variables in all subjects (dataset: ", dat, ")"))
        cat_sum <- vector("list", length(cat_vars))
        for (c in 1:length(cat_vars)) {
                
                eval(parse(text = str_c("tmp <- ", dat, "[, .N, by = ", cat_vars[c], "][order(", 
                                        cat_vars[c], ")]")))
                tmp[, percent := round(N/sum(N), sig_digit)]
                
                names(cat_sum)[c] <- cat_vars[c]
                cat_sum[[c]] <- tmp
        }
        print(cat_sum)
        
        if (!is.na(by_var)) {
                cat("\n")
                print(str_c("Number of subejcts by ", by_var, " (dataset: ", dat, ")"))
                eval(parse(text = str_c(dat, "[, .N, by = ", by_var, "]"))) %>% print()
                cat("\n")
        
                print(str_c("Summarize continuous variables by ", by_var, " (dataset: ", dat, ")"))
                cnt_sum_byvar <- vector("list", length(cnt_vars))
                for (c in 1:length(cnt_vars)) {
                        
                        eval(parse(text = str_c("mean <- ", dat, "[, mean(", cnt_vars[c], ", na.rm = T), by = ", by_var, "]")))
                        setnames(mean, "V1", "mean")
                        eval(parse(text = str_c("sd <- ", dat, "[, sd(", cnt_vars[c], ", na.rm = T), by = ", by_var, "]")))
                        setnames(sd, "V1", "sd")
                        eval(parse(text = str_c("na.count <- ", dat, "[, sum(is.na(", cnt_vars[c], ")), by = ", by_var, "]")))
                        setnames(na.count, "V1", "na.count")
                        
                        tmp <- merge(mean, sd, by = by_var)
                        tmp <- merge(tmp, na.count, by = by_var)
                        
                        eval(parse(text = str_c("ttest.p <- ", dat, "[, t.test(", cnt_vars[c], " ~ ", by_var, ")]$p.value")))
                        tmp[1, ttest.p := round(ttest.p, sig_digit)]
                        
                        names(cnt_sum_byvar)[c] <- cnt_vars[c]
                        cnt_sum_byvar[[c]] <- tmp
                }
                print(cnt_sum_byvar)
                
                print(str_c("Summarize categorical variables by ", by_var, " (dataset: ", dat, ")"))
                cat_sum_byvar <- vector("list", length(cat_vars))
                for (c in 1:length(cat_vars)) {
                        
                        eval(parse(text = str_c("tmp <- ", dat, "[", cat_vars[c], " != '', .N, by = .(", by_var, ", ", 
                                                cat_vars[c], ")][, percent := round(N/sum(N), 3), by = ", by_var, "]")))
                        eval(parse(text = str_c("tmp <- tmp[order(", by_var, ", ", cat_vars[c], ")]")))
                        
                        chisq.p <- chisq.test(matrix(tmp[, N], nrow = (nrow(tmp)/2), ncol = 2))$p.value
                        suppressWarnings(tmp[1, chisq.p := round(chisq.p, sig_digit)])
                        
                        names(cat_sum_byvar)[c] <- cat_vars[c]
                        cat_sum_byvar[[c]] <- tmp
                }
                print(cat_sum_byvar)
        }
}
