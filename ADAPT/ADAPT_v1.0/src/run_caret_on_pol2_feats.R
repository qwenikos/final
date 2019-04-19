if(!require(gbm)){
    install.packages("gbm")
    library(gbm)
}

if(!require(caret)){
    install.packages("caret",version = "6.0-81",
                 repos = "http://cran.r-project.org", 
                 dependencies = c("Depends", "Imports", "Suggests"))
    library(caret)
}
/home/nikos/biothesis/github/TSS/11_PAPER_NEW/ADAPT_algorithm/ADAPT_v1.0/src/tools/R/caret_6.0-81.tar.gz

suppressMessages(library(gbm))
suppressMessages(library(caret))

args <- commandArgs(trailingOnly=TRUE)

options(max.print=999999)

super_model <- readRDS(args[1])
train_data <- read.csv(args[2], header=TRUE, sep=",")

final_preds <- predict(super_model,train_data[,2:14],type="prob")
print(final_preds)