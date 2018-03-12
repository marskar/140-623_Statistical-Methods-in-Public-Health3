#' ---
#' title: "140.623.01 - Statistical Methods in Public Health III"
#' subtitle: "Assignment 2: Survival in Primary Biliary Cirrhosis"
#' author: "Martin Skarzynski"
#' institute: Johns Hopkins School of Public Health
#' date: March 8, 2018
#' ---

#' ## Learning Objectives:

#' Students who successfully complete this section will be able to:

#' - To evaluate whether the drug DPCA prolongs life in patients.
#' - To identify baseline characteristics of patients which predict longer survival.
#' - Analyze the survival time data (without grouping) by the Kaplan-Meier estimate of the
#' survival function, the log- rank statistic, and Cox proportional hazards model.
#' - Check the estimated model for its consistency with the observed data; in particular, check
#' the proportional hazards assumption using the complementary log-log plot of the
#' estimated survival function.
#' - Summarize the findings for public health readers and document and archive the steps of
#' the statistical analysis by creating a script file in R.

#' 
#' ## Data Set:

#' Between January 1974 and May 1984, a double-blinded randomized trial on patients with primary biliary cirrhosis (PBC) of the liver was conducted at the Mayo clinic. A total of 312 patients were randomized to either receive the drug D-penicillin (DPCA) or a placebo. Patients were followed until they died from PBC or until censoring, either because of administrative censoring (withdrawn alive at end of study), death not attributable to PBC, liver transplantation, or loss to follow-up. At baseline, a large number of clinical, biochemical, serological and histologic measurements were recorded on each patient. This data set is a subset of the original data, and includes information on each patient's time to death or censoring, treatment, age, gender, serum bilirubin, and histologic disease stage (1-4).

#' The variables included in this dataset include:

#'- case: unique patient ID number
#'- sex: 0 = male, 1 = female (coded as "Female" and "Male" in the csv file rather than 0/1)
#'- drug: 0 = placebo, 1 = DPCA
#'- bil : serum bilirubin in mg/dl
#'- survyr: time (in years) to death or censoring
#'- death: indicator = 1 if patient died, 0 if censored
#'- ageyr: age in years [continuous variable]
#'- histo: histologic disease stage (1 - 4) [categorical variable]
#'- agecat: age categories, coded as "< 45 yrs", "45 - 55 yrs", and ">= 55 yrs"

#' Also included in the data set for your possible use are the following indicator (dummy) variables:

#' 
#' Age Indicators (indicator versions of agecat):
#' - agegr_2: 1 if patient is 45-55 years old, 0 otherwise
#' - agegr_3: 1 if patient is >= 55 years old, 0 otherwise
#' Histologic Stage Indicators:
#' - hstage2: 1 if patient is in Stage 2, 0 otherwise
#' - hstage3: 1 if patient is in Stage 3, 0 otherwise
#' - hstage4: 1 if patient is in Stage 4, 0 otherwise
#' 

#' The data are stored in the csv data set pbctrial.csv, which may be
#' downloaded from the course website.

#' ## Methods:

#' Use the data set described above and the appropriate statistical analyses to address the specific
#' learning objectives listed on the first page.
#' Hints: The hints shown below are based on a dataset with the name pbcData, read in with the
#' following code. In the following list of commands, if you want to look at differences by other
#' variables than drug, you should change the variable name! Create a new .R file to type/run your
#' commands so that you will have a record of your analysis.

setwd("~/github/140-623_Statistical-Methods-in-Public-Health3")
library(readr)
pbcData = read_csv("pbctrial.csv")

#' a. Explore the data using descriptive statistics:

#'- table()
#'- prop.table()
#'- summary() etc

dim(pbcData)
str(pbcData)
summary(pbcData)
sum(pbcData$death)
sum(pbcData$death)/length(pbcData$death)*100
library(purrr, help)
map(pbcData, class)
pbcData$histo <- as.factor(pbcData$histo)
pbcData$agecat <- as.factor(pbcData$agecat)
map(pbcData, class)
round(prop.table(table(pbcData[c("death", "drug", "sex")])), 3)

#' b. Define a survival object, defining the time variable (survyr) and the event (death == 1).
#' To do this, you must first install and load the "survival" package:

# install.packages("survival")
library(survival)


## only run this the first time

pbcData$SurvObj = with(pbcData, Surv(survyr, death == 1))

#' c. Explore differences in time to death by different baseline variables using graphs and complementary log-log plots.

# estimate survival curves for entire sample
km.overall = survfit(SurvObj ~ 1, data = pbcData,
type="kaplan-meier", conf.type="log-log")
km.overall
summary(km.overall)
# estimate survival curves for drug group
km.drug = survfit(SurvObj ~ drug, data = pbcData,
type="kaplan-meier", conf.type="log-log")
km.drug
summary(km.drug)

# plot km curves
plot(km.overall)
plot(km.drug)
# log rank test for equality of survivor functions
survdiff(SurvObj ~ drug, data=pbcData)
# complimentary log-log plot
plot(km.drug, fun="cloglog", ylab="log(-log(Survival Probability)",
xlab="Analysis time (shown on log scale)")

#' d. Fit several Cox proportional hazards regression models to the ungrouped survival data:

model1 = coxph(SurvObj ~ drug, data = pbcData)
summary(model1)

model2 = coxph(SurvObj ~ sex + bil + histo, data = pbcData)
summary(model2)

#' e. Save your R script file that documents and archives the steps of your statistical analysis.
#' This file will make your analysis "reproducible."
#' f. Summarize your findings in a brief report (less than two pages with at most one table and
#' one figure) as if for a biomedical/public health journal.
#' A suggested format is:
#' - Introduction - a few sentences about the research question(s)
#' - Data description - simple tabulations describing patient characteristics
#' - Results from multiple models that address question(s) (e.g., bivariate and multivariable)
#' - Graphical display that presents evidence in the data relevant to your scientific question.
#' 

#' ## Introduction 

#' The research question that I will try to answer in this report is whether D-penicillin (DPCA) provided any benefit for the primary biliary cirrhosis (PBC) patient population as a whole (n=312) and for sub-groups based on sex, age and histologic disease stage in a double-blinded randomized trial conducted at the Mayo clinic between January 1974 and May 1984. I hypothesize that the drug effect (if any) will be diminished in the higher age categories and disease stages. In other words, I expect that there will be differences in drug reponse as measured by time to death between the 4 disease stages and 3 age categories, specifically that older patients and those with more advanced disease will be more difficult to treat, which will result in a shorter survival time. I will also assess whether serum bilirubin level is a prognostic marker and whether drug benefit will differ among men versus women.
#' 

#' ## Results

#' I calculated descriptive statistics and determined that the overall median survival time was around 5 years. As for patient characteristics, the representation across age categories and disease stages appears to spread relatively evenly. The `age` and `survyr` variables appear to be normally distributed with a slight rightward skew. Interestingly, bilirubin is skewed highly to the right (mean = 3.3 mg/dl, median = 1.4 mg/dl) indicating that there are outliers with high bilirubin values. The patient population is 88% female; out of the total 312 patients, 276 were women and only 36 were men. Ages of patients ranged from 26 to 78 years, with a median age of ~50 years. Roughly three-thirds of of the patients were in a histologic stage 3 or 4. Mortality was high during the study. In the data collected, approximately 40% (125 out of 312) of study participants died from primary biliary cirrhosis.

#' There was no statistically significant (using an $\alpha$ of .05) difference between patients in the placebo and drug groups. Overall, simple Cox proportional hazards regression analysis showed that the drug group had a 6% greater hazard of death than the placebo group. Multivariable Cox regression analysis that included sex, age categories, bilirubin levels, and histologic disease stage in the model showed a 12% greater hazard in the group assigned the drug, though this results was also not statistically significant. Table 1 summarizes the results of the multivariate Cox regression analysis. I used the Wald and likelihood ratio tests to assess the statistical significance of the variables in the multivariate Cox regression model. There was a statistically significant increase in hazard of death for males versus females (HR = 1.71, p = 0.027), and those in the highest age category versus the lowest age category (HR = 1.71, p=0.031) and most advanced (histo = 4) histologic disease stage versus the least advanced (histo = 1) disease stage (HR = 15.0, p = 0.008). There was also a 16% higher hazard of death with every unit (mg/dl) increase of serum bilirubin (p < 0.001).

#' I also calculated Kaplan-Meier estimates of sub-groups based on the variables in the multivariate Cox regression model and the `drug` variable. These analyses did not indicate that the drug might be beneficial to some types of patients. Kaplan-Meier estimates of the survivor functions for various sample sub-groupings were calculated. Simple Cox regression models were used to evaluate univariate associations between patient characteristics and survival. Serum bilirubin was the only continuous covariate in the regression models and I converted this into a categorical variable (binary) that was assigned a value of 1 if serum bilirunbin level was above the median and 0 if serum bilirunbin level was below the median. I plotted the Kaplan-Meier estimates against time. Shockingly, men taking the drug appear to have a much shorter survival time than men taking placebo or women in either treatment group(Figure 1 top-left). This may help to explain why males had a 71% greater hazard of death compared to otherwise similar females (p = 0.027) in multivariate Cox regression analysis. The drug also appeared to have a negative effect on survival in patients with earliest stage of disease (histo = 1) compared to later stages (Figure 1 bottom-left). The categorical variable I created using bilirubin levels appears to cleanly divide patients with the best and worst survival in both treatment groups (Figure 1 bottom-right).
#' 

#' ## Conclusions
#' The drug tested in this study DCPA did not statistically significantly increase survival according to univariate or multivariable cox proportional hazards analyses. The conclusion I draw from this randomized trial is that DPCA is not an effective treatment for patients with primary biliary cirrhosis. Alarmingly, the drug appears to increase the risk of death for men and patients with least advanced disease stage as determined by histology. The analysis described herein also present the possibility that bilirubin could be a prognostic biomarker for primary biliary cirrhosis. This work is only the beginning and more precise answers to the research questions discussed in the introduction will require further inspection with models more precisely adapted to each research question.
#' 

# First, I will produce a few simple summaries of drug response based on `sex`, `agecat` and `histo` variables. First some basic exploratory data analysis will let me know if I am on the right track with the variables I have chosen. If there is no difference between the median survival times of the groups I am interested in, it will be unlikely that I will see anything significant in my model. 
library(dplyr)
# install.packages("broom")
library(broom)
pbcData %>%
    group_by(sex, drug) %>%
    summarise(med_surv = median(survyr))
pbcData %>%
    group_by(agecat, drug) %>%
    summarise(med_surv = median(survyr))
pbcData %>%
    group_by(histo, drug) %>%
    summarise(med_surv = median(survyr))



# I decided to put all variables of interest into one model rather creating multiple models that address each of the above questions, because the instructions say to have at most one figure and one table. If any of the results are statistically significant, I can explore the question further with a more specific model in the future. 


cox_all_var = coxph(formula = SurvObj ~ drug + sex + bil + histo + agecat, data = pbcData)
all_var_summary <- summary(cox_all_var)
library(broom, help)
cox_all_var %>%
tidy()
cox_all_var %>%
summary()
coef(cox_all_var)
tidy(cox_all_var)["p.value"]
coef(cox_all_var) %>%
summary()

df <- data.frame(adj_HR = round(exp(coef(cox_all_var)), 3),
           lower_CI = round(exp(confint(cox_all_var)[,1]), 3),
           upper_CI = round(exp(confint(cox_all_var)[,2]), 3),
           p_value = round(tidy(cox_all_var)["p.value"], 3))
rownames(df) <- rownames(confint(cox_all_var))

#install.packages("captioner")
library(captioner, help)
figs <- captioner(prefix="Figure")
tbls <- captioner(prefix="Table")
library(knitr)
knitr::kable(df, format = "markdown")

#' `r tbls(name="adj_HR_tbl", caption = "Adjusted Hazard Ratio Estimates of Death obtained from Proportional Hazards Regression.")`

#'

# Plotting
par(mfrow=c(2,2), mar = c(0, 0, 0, 0), oma = c(4, 4, 0.1, 0.1))
palette()
# sexplot
km_sex = survfit(SurvObj ~ drug + sex, data = pbcData,
type="kaplan-meier", conf.type="log-log")
plot(km_sex, las = 1,
     xaxt='n', ann=FALSE,
     col = 1:8)
legend("bottomleft",
       legend=names(km_sex$strata),
       col=1:length(km_sex$strata),
       cex = 0.75,
       lty=c(1,1), # gives the legend appropriate symbols (lines) 
       lwd=c(2.5,2.5))

# ageplot
km_age = survfit(SurvObj ~ drug + agecat, data = pbcData,
type="kaplan-meier", conf.type="log-log")
plot(km_age,
     xaxt='n', yaxt='n', ann=FALSE,
     col = 1:8, xlab = "Time", ylab = "Survival")
legend("bottomleft",
       legend=names(km_age$strata),
       col=1:length(km_age$strata),
       cex = 0.65,
       lty=c(1,1), # gives the legend appropriate symbols (lines) 
       lwd=c(2.5,2.5))

# histoplot
km_histo = survfit(SurvObj ~ drug + histo, data = pbcData,
type="kaplan-meier", conf.type="log-log")
plot(km_histo, las = 1, col = 1:8)
legend("bottomleft",
       legend=names(km_histo$strata),
       col=1:length(km_histo$strata),
       cex = 0.6,
       lty=c(1,1), # gives the legend appropriate symbols (lines) 
       lwd=c(2.5,2.5))



# bilplot
# To make a similar plot with the `bil` variable, I will first create a new categorical (binary) variable called `bilcat`.
pbcData['bilcat'] <- ifelse(pbcData["bil"][[1]]>median(pbcData["bil"][[1]]), 1, 0)
head(pbcData)
km_bil = survfit(SurvObj ~ drug + bilcat, data = pbcData,
type="kaplan-meier", conf.type="log-log")
plot(km_bil,
     yaxt='n', ann=FALSE,
     cex.lab = 0.75,
     col = 1:8)
legend("bottomleft",
       legend=names(km_bil$strata),
       col=1:length(km_bil$strata),
       cex = 0.75,
       lty=c(1,1), # gives the legend appropriate symbols (lines) 
       lwd=c(2.5,2.5))
mtext("Time (years)", side = 1, outer = TRUE, cex = 1.15, line = 2.2, col = "black")
mtext("Survival", side = 2, outer = TRUE, cex = 1.15, line = 2.2, col = "black")
#' `r figs(name="km_plots_fig", caption = "Survival of Primary Biliary Cirrhosis patients treated with D-penicillimin (DPCA) or a placebo.")`

#' 

