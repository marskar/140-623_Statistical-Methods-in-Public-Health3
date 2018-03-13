#' ---
#' title: "140.623.01 - Statistical Methods in Public Health III"
#' subtitle: "Assignment 4: Survival in Framingham Heart Study"
#' author: "Martin Skarzynski"
#' institute: Johns Hopkins School of Public Health
#' date: March 13, 2018
#' ---

#' ## Learning Objectives:

#' Students who successfully complete this section will be able to:

#' - Analyze the relationship between grouped survival time data and baseline covariates of interest using log-linear Poisson regression models.
#' - Check the assumptions for Poisson regression and use other models (such as negative binomial) as appropriate.
#' - Summarize the findings in a brief fashion for public health readers.
#' - Document and archive the steps of the statistical analysis.
#' 

#' Data Set:
#' The Framingham Heart Study is a long term prospective study of the etiology of cardiovascular disease among a population of free living subjects in the community of Framingham, Massachusetts. Individuals were followed for 24 years. These data are binned into 5- year intervals (1876 days each) and stratified by gender and baseline current smoking, age category, BMI category, diabetes, and blood pressure medications (see Coding Description on the next page). The data are stored in the csv data set FraminghamPS4bin.csv which may be downloaded from the course website.
#' 

#' Methods:
#' Use the data set described above and the appropriate statistical analyses to address the specific learning objectives.  Hints: The hints shown below are based on a dataset with the name framData, read in with the following code. In the following list of commands, if you want to look at differences by other variables than drug, you should change the variable name! Create a new .R file to type/run your commands so that you will have a record of your analysis.

library(readr)
framData = read_csv("FraminghamPS4bin.csv")

#' a. Explore the data using descriptive statistics:

#'- table()
#'- prop.table()
#'- summary() etc

dim(framData)
str(framData)
summary(framData)
library(purrr, help)
map(framData, class)


#' b. Explore several Poisson regression models using these grouped survival data and select between models:

model1 = glm(D ~ gender, offset = log(Y), data =  framData, family=poisson(link="log"))
summary(model1)
AIC(model1)

#' c. Check the assumptions of your Poisson models; use other models as appropriate:

# Pearson chi-square goodness-of-fit test (like poisgof in Stata)
X2 = sum(residuals(model1, type = "pearson")^2); X2
df = model1$df.residual; df
pval = 1-pchisq(X2, df); pval

# Negative binomial regression
library(MASS)
model2 = glm.nb(D ~ gender + offset(log(Y)), data=framData)
summary(model2)
AIC(model2)

#' d. Save your R script file that documents and archives the steps of your statistical analysis.
#' This file will make your analysis “reproducible.”
#' e. Summarize your findings in a brief report (less than two pages with at most one table and
#' one figure) as if for a biomedical/public health journal.
#' 

#' A suggested format is:
#' 

#' - Introduction – a few sentences about the research question(s)
#' - Data description – simple tabulations describing individual characteristics
#' - Results from multiple models that address question(s) (e.g., bivariate and multivariable)
#' - Graphical display that presents evidence in the data relevant to your scientific question.
#' 

model3 = glm(D ~ gender + cursmoke + diabetes + bpmeds + bmicat + agecat, offset = log(Y), data =  framData, family=poisson(link="log"))
summary(model3)
AIC(model3)


library(broom, help)
tidy(model3)
coef(model3)
tidy(model3)$p.value

df <- data.frame(adj_RR = round(exp(coef(model3)), 4),
           lower_CI = round(exp(confint(model3)[,1]), 4),
           upper_CI = round(exp(confint(model3)[,2]), 4),
           p_value = round(tidy(model3)$p.value, 4))
rownames(df) <- rownames(confint(model3))
df

#install.packages("captioner")
library(captioner, help)
figs <- captioner(prefix="Figure")
tbls <- captioner(prefix="Table")
library(knitr)

#' ## Introduction 

#' The Framingham Heart Study is a prospective study that followed study participants for 24 years in an attempt to better understand the etiology of cardiovascular disease. The study population is the community of Framingham, Massachusetts. The data that were obtain from the study were binned into 1875-day intervals (roughly 5 years). Additionally, categorical variables were created from the ages and body mass indices (BMI) of study participants.

#' ## Data Description

#' The research question that I will try to answer in this report is whether there is a relationship between grouped survival time data and baseline covariates of interest in the Framingham Heart Study. To answer this question, I will use a binned version of the Framingham data set and log-linear Poisson regression models. I hypothesize that smoking status, BMI and agecat will have a strong effect on  the death rate. Specifically, I expect that male, obese, diabetic study participants that smoke and belong to the oldest age category will have a higher death rate. I will also assess whether anti-hypertensive medications provide any benefit.

#' 

#' ## Results

#' I calculated descriptive statistics and determined that the overall mean and median death rates are `r mean(framData$Rate)` and `r median(framData$Rate)`. Interestingly, the `Rate` variable is skewed highly to the right indicating that there are outliers with high death rates. The study population is `r sum(framData$gender)/length(framData$gender)*100`% female; out of the total `r length(framData$gender)` study participants, `r sum(framData$gender)` were women and `r length(framData$gender)-sum(framData$gender)` were men. Roughly `r sum(framData$bpmeds)/length(framData$bpmeds)*100`% of the study participants took blood pressure medication. All of the variables in my log-linear model (`r coef(model3)`) were statistically significant (based on an $\alpha$ value of 0.5). The results of the model of summarized in Table 1. The first column shows the death rate ratio, the second and third columns show the lower and upper confidence intervals (respectively) and the final column shows the first four subzero digits of the p-value. All of the death rate ratios were above 1 except for gender indication that being a current smoker, taking anti-hypertensive medication, being diabetic, being obese and being elderly were all associated with a higher death rate, while being female meant that participants were less likely to die in the Framingham study. The `diabetes` variable had the highest coefficient, but also the widest confidence interval. I then 

#' ## Conclusions

#' The drug tested in this study DCPA did not statistically significantly increase survival according to univariate or multivariable cox proportional hazards analyses. The conclusion I draw from this randomized trial is that DPCA is not an effective treatment for patients with primary biliary cirrhosis. Alarmingly, the drug appears to increase the risk of death for men and patients with least advanced disease stage as determined by histology. The analysis described herein also present the possibility that bilirubin could be a prognostic biomarker for primary biliary cirrhosis. This work is only the beginning and more precise answers to the research questions discussed in the introduction will require further inspection with models more precisely adapted to each research question. 

#' ## Graphical Display 
#' I decided to plot the log death rates per bin for every observation and color each of the variables of interest. The mean for each subgroup is shown as a horizontal bar. Interestingly, the diabetes (`diabetes=1`; bottom-left) and higher age categories (`agecat=3` and `4`; top-left) were consistently associated with a higher death rate. As for the bpmeds variable, I believe that the high coefficient associated with this variable would disappear if we controlled for blood pressure, as participants with the highest blood pressure would be most likely to be perscribed anti-hypertensive medicine and most likely to die of cardiovascular complications. In conclusion, this analysis presents a multivariate log-linear model and univariate plots that highlight a potentially important link between various variables and the death rate in the Framingham Heart Study. Among the variables studied (`r coef(model3)`) diabetes stood out as having the strong association (highest coefficient) with the death rate. Further research is needed to improve our understanding of the interactions and etiologies of diabetes and cardiovascular disease.

knitr::kable(df, format = "markdown")

#' `r tbls(name="adj_RR_tbl", caption = "Adjusted Hazard Ratio Estimates of Death obtained from Proportional Hazards Regression.")`

bins <- unique(framData$tbin)

par(mfrow=c(3,3), mar = c(0, 0, 0, 0), oma = c(4, 4, 0.1, 0.1))

plot(log(Rate) ~ jitter(tbin, 1),
     xaxt='n', ann=FALSE,
     data = framData, col = agecat)
ctgs <- unique(framData$agecat)
for(bin in bins){
    for(ctg in ctgs){
        avg <- log(mean(framData$Rate[framData$tbin == bin & framData$agecat == ctg]))
        lines(c(bin-250, bin+250), (c(avg, avg)), col = ctg, lwd = 3)
}}

legend("top",
       legend=paste0("agecat=", unique(framData$agecat)),
       pch = 1,
       col=1:length(unique(framData$agecat)),
       cex = 0.75)

plot(log(Rate) ~ jitter(tbin, 1),
     yaxt='n', xaxt='n', ann=FALSE,
     data = framData, col = bmicat)
ctgs <- unique(framData$bmicat)
for(bin in bins){
    for(ctg in ctgs){
        avg <- log(mean(framData$Rate[framData$tbin == bin & framData$bmicat == ctg]))
        lines(c(bin-250, bin+250), (c(avg, avg)), col = ctg, lwd = 3)
}}
legend("top",
       legend=paste0("bmicat=", unique(framData$bmicat)),
       pch = 1,
       col=1:length(unique(framData$bmicat)),
       cex = 0.75)


plot(log(Rate) ~ jitter(tbin, 1),
     yaxt='n', xaxt='n', ann=FALSE,
     data = framData, col = cursmoke + 1)
ctgs <- unique(framData$cursmoke)
for(bin in bins){
    for(ctg in ctgs){
        avg <- log(mean(framData$Rate[framData$tbin == bin & framData$cursmoke == ctg]))
        lines(c(bin-250, bin+250), (c(avg, avg)), col = ctg+1, lwd = 3)
}}
legend("top",
       legend=paste0("cursmoke=", unique(framData$cursmoke)),
       pch = 1,
       col=1:length(unique(framData$cursmoke)),
       cex = 0.75)

plot(log(Rate) ~ jitter(tbin, 1), data = framData, col = diabetes + 1)
ctgs <- unique(framData$diabetes)
for(bin in bins){
    for(ctg in ctgs){
        avg <- log(mean(framData$Rate[framData$tbin == bin & framData$diabetes == ctg]))
        lines(c(bin-250, bin+250), (c(avg, avg)), col = ctg+1, lwd = 3)
}}
legend("top",
       legend=paste0("diabetes=", unique(framData$diabetes)),
       pch = 1,
       col=1:length(unique(framData$diabetes)),
       cex = 0.75)


plot(log(Rate) ~ jitter(tbin, 1),
     yaxt='n', ann=FALSE,
     data = framData, col = bpmeds + 1)
ctgs <- unique(framData$bpmeds)
for(bin in bins){
    for(ctg in ctgs){
        avg <- log(mean(framData$Rate[framData$tbin == bin & framData$bpmeds == ctg]))
        lines(c(bin-250, bin+250), (c(avg, avg)), col = ctg+1, lwd = 3)
}}
legend("top",
       legend=paste0("bpmeds=", unique(framData$bpmeds)),
       pch = 1,
       col=1:length(unique(framData$bpmeds)),
       cex = 0.75)

plot(log(Rate) ~ jitter(tbin, 1),
     yaxt='n', ann=FALSE,
     data = framData, col = gender + 1)
ctgs <- unique(framData$gender)
for(bin in bins){
    for(ctg in ctgs){
        avg <- log(mean(framData$Rate[framData$tbin == bin & framData$gender == ctg]))
        lines(c(bin-250, bin+250), (c(avg, avg)), col = ctg+1, lwd = 3)
}}
legend("top",
       legend=paste0("gender=", unique(framData$gender)),
       pch = 1,
       col=1:length(unique(framData$gender)),
       cex = 0.75)

#' `r figs(name="km_plots_fig", caption = "Survival of Primary Biliary Cirrhosis patients treated with D-penicillimin (DPCA) or a placebo.")`