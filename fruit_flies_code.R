#dryR: https://github.com/naef-lab/dryR
#discorhythm: https://www.bioconductor.org/packages/release/bioc/html/DiscoRhythm.html
#https://bioconductor.org/packages/release/bioc/vignettes/DiscoRhythm/inst/doc/disco_workflow_vignette.html#65_Oscillation_Detection
# library(devtools)
# install_github("matthewcarlucci/DiscoRhythm", build_vignettes=TRUE)
library(tidyverse)
library(DiscoRhythm)

test <- as.data.frame(t(read.csv("/Users/gregorymatthews/Dropbox/CDSC/fruit_flies/MX7590~2(data)_clean_greg.csv", header = FALSE)))
names(test) <- test[1,]
test <- test[-1,]

table(test$treatment)

test[,9:ncol(test)] <- apply(test[,9:ncol(test)],2,as.numeric)

test <- test %>% filter(treatment != "pool - pool") 
test <- test %>% mutate(replicate = substr(treatment,nchar(treatment),nchar(treatment)), 
                        trt = substr(treatment,1,nchar(treatment)-4),
                        index = rep(1:24,8),
                        time = rep(seq(0,46,2),8),
                        day = rep(rep(c(1,2),each = 12),8)) 
test$trt[test$trt == "clock856"] <- "clk856"
names(test) <- gsub(" ","",names(test))

ggplot(aes(x = time, y = zymostenol, color = factor(day)), data = test) + geom_line(alpha = 0.5) + geom_smooth(se = F) + facet_grid(replicate~trt)
ggplot(aes(x = time, y = xylulose, color = factor(day)), data = test) + geom_line(alpha = 0.5) + geom_smooth(se = F) + facet_grid(replicate~trt)


"lignocericacid"

greg <- test %>% filter(trt == "clk856") %>%  select(lignocericacid, time, day, replicate)
ggplot(aes(x = time, y = lignocericacid, color = factor(replicate)), data = greg) + geom_line() 


#Disco rhythm code
#We use mxsample as the sample ID.  
#the metabolomic is the "gene"
#Prefix Time*_Unique Id_Replicate Id

#First filter by treatment
disco_sub <- test %>% filter(trt == "clk856") %>% mutate(discoID = paste0("hr",time,"_",index,"_",replicate)) %>% select(c(9:584,590))
t_disco_sub <- as.data.frame(t(disco_sub[,-ncol(disco_sub)]))
names(t_disco_sub) <- disco_sub$discoID
t_disco_sub$IDs <- row.names(t_disco_sub)
t_disco_sub <- t_disco_sub %>% relocate(IDs)

"lignocericacid"


indata <- t_disco_sub
#Converting dfrom DF to SE object.
se <- discoDFtoSE(indata)

#Check input 
selectDataSE <- discoCheckInput(se)

#PCA
discoPCAres <- discoPCA(selectDataSE)

PeriodRes <- discoPeriodDetection(selectDataSE,
                                  timeType="linear",
                                  main_per=24)


OVpca <- discoPCA(FinalSE)
OVpcaSE <- discoDFtoSE(data.frame("PC"=1:ncol(OVpca$x),t(OVpca$x)),
                       colData(FinalSE))
knitr::kable(discoODAs(OVpcaSE,period = 24,method = "CS")$CS,
             format = "markdown")


mustache <- discoODAs(selectDataSE,
                         period=24,
                         method=c("JTK","CS"),
                         ncores=1,
                         circular_t=FALSE)
mustache$JTK %>% arrange(pvalue) %>% View()
mustache$CS %>% arrange(pvalue) %>% View()

plot(-log(sort(mustache$CS$pvalue)), ylim = c(0,10))
abline(h = -log(0.05/576), col = "red")
length(mustache$JTK$pvalue)

plot(-log(sort(runif(576))), ylim = c(0,10))
abline(h = -log(0.05/576), col = "red")



# library(SummarizedExperiment)
# Metadata <- colData(selectDataSE)
# knitr::kable(discoDesignSummary(Metadata),format = "markdown")










