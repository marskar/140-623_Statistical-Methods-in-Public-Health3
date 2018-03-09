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
library(readr)
pbcData = read_csv("pbctrial.csv")

#' a. Explore the data using descriptive statistics:
#'- table()
#'- prop.table()
#'- summary() etc

dim(pbcData)
str(pbcData)
summary(pbcData)
library(purrr, help)
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

model2 = coxph(SurvObj ~ sex + bil + as.factor(histo), data = pbcData)
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

#' Between January 1974 and May 1984 a double-blinded randomized trial on patients with primary biliary cirrhosis (PBC) of the liver was conducted at the Mayo clinic. A total of 312 patients were randomized to either receive the drug D-penicillimin (DPCA), or a placebo. Patients were followed until they died from PBC, or until censoring, either because of administrative censoring (withdrawn alive at end of study), death not attributable to PBC, liver transplantation, or loss to follow-up. At baseline clinical, biochemical, serological and histologic measurements were recorded on each patient. A sub-study was undertaken to test for increased survival amongst patients on the new treatment, and to investigate the association between survival and patients' age, gender, histologic stage of disease, and serum bilirubin level.

#' The research question that I will try to answer in this report is whether D-penicillin (DPCA), the drug tested in the PBC trial, provided any benefit for the patient population as a whole (n=312) and for sub-groups based on sex, age and disease stage. I hypothesize that the drug effect will not be different between the 3 age categories, but will depend on disease stage. In other wrods, I expect that there will be differences in time to death between the 4 disease stages, specifically that more advanced disease will be more difficult to treat, which will result in a shorter time to event. I will also assess whether bilirubin is a prognostic marker and whether drug benefit will differ among men versus women.
#' 

#' ## Data description 


#' There are a total of 312 patients and the median survival time was around 5 years. As for patient characteristics, the representation across age categories and disease stages appears to spread relatively evenly. The `age` and `survyr` variable appear to be normally distributed with a slight leftward skew. Interestingly, bilirubin is skewed highly to the left indicating that there are outliers with high bilirubin values.

#' 

#' ## Methods
#' Descriptive statistics were calculated to investigate sample characteristics. Kaplan-Meier estimates of the survivor functions for various sample sub-groupings were calculated. Simple Cox regression models were used to evaluate univariate associations between patient characteristics and survival. Multivariable Cox regression was used to examine the association between survival and multiple patient characteristics simultaneously. Serum bilirubin was the only continuous covariate in the regression models. Age was modeled as a categorical variable, based on tertiles in the sample, to allow for a non-linear relationship between age and the loghazard of death. Both Wald and likelihood ratio methods were used to test for the statistical significance of covariates in the final multiple proportional hazards model. Only predictors achieving statistical significance ($\alpha$ = .05) were included in the final multivariable model.
#' 

#' ## Study Enrollees
#' The sample consists of 312 patients with primary biliary cirrhosis enrolled from 1974 to 1984 at the Mayo Clinic in Rochester, MN. The sample is majority female (276 patients, 88%) with only 36 male patients (12%). The average patient age at enrollment was 50 years, and the sample age range was from 26 to 78 years. The majority (75%) of the patients were in a later stage of the disease (Histologic Stage 3 or 4) at the time of enrollment. Average serum bilirubin level among participants at time of enrollment was 3.3 mg/dl. At the time of this analysis, 125 patients (40%) had died from causes related to primary biliary cirrhosis. Results Patients in the drug group had 6% greater hazard ("risk") of death than those in the placebo group, but this result was not statistically significant (95% CI, -25% - 50%, p > .05).

#' 

#' Serum bilirubin level, patients age, and histologic stage of disease all had statistically significant ($\alpha$ = .05) positive univariate associations with the hazard of death. Males had 62% higher risk of death than females (95% CI 2% - 158%, p = .04). In a multivariable analysis, serum bilirubin level, gender, and histologic stage of disease were found to have statistically significant associations with patient survival. The hazard ratio associated with a 1 mg/dl increase in serum bilirubin level was 1.16 (95% CI 1.13 - 1.19), indicating that a patient's risk of death increases by 16% for each 1 mg/dl increase in serum bilirubin after adjustment for gender, disease stage and age. The hazard ratio of death for males relative to females was 1.70 (95% CI 1.05 - 2.74), indicating that males had a 70% increase in the hazard of death compared to otherwise similar females. Those patients in the highest stage (stage 4) of disease had greater than 14 times the adjusted risk (95% CI 2.01 - 107.40) of dying when compared to patients in the earliest stage (stage 1). Table 1 presents results from both the unadjusted and adjusted sets of analyses.

#' ## Graphical display 
#' I plotted cox model fitted values against ageyr and bil variables marking sex and disease stage (`histo`) with color and different symbols, respectively. 

#' ## Conclusions
#' DCPA was not found to be statistically significantly associated with increased survival in either univariate or multivariable analyses. As this was a randomized trial with 312 patients, we conclude that DPCA does not appear to be efficacious in the treatment of patients with primary biliary cirrhosis. While primary biliary cirrhosis is a disease that primarily affects females, the prognosis is significantly worse for males. Similarly, the risk of death is much worse for patients in later stages of the disease relative to those in the earlier stages. The results of this research suggest that improved screening techniques to identify the disease in affected patients early on, coupled with increased outreach to males at risk of developing PBC could result in a better overall prognosis for patients having this disease.
#' It is clear that all of the variables I picked are important in the final model although not all levels of the categorical variables were statistically significant. This work is only the beginning and more precise answers to the research questions discussed in the introduction will require further inspection with models more precisely adapted to each research question.
#' 

#' ## Results 
#' First, I will produce a few simple summaries of drug response based on `sex`, `agecat` and `histo` variables. First some basic exploratory data analysis will let me know if I am on the right track with the variables I have chosen. If there is no difference between the median survival times of the groups I am interested in, it will be unlikely that I will see anything significant in my model. From this initial analysis it looks like patients in the highest age category that were given placebo fare the worst. These results indicate that elderly patients my stand to benefit the most from taking the drug. Shockingly, men taking the drug appear to have a shorter survival time than with the drug, and do not survive as long as a women in general. Similarly, the drug appeared to have a negative effect on survival in patients with earliest stage of disease (histo = 1). Now I will take a similar approach to the data but using a cox proportional hazards model.
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

par(mfrow=c(2,2))


#' I decided to put all variables of interest into one model rather creating multiple models that address each of the above questions, because the instructions say to have at most one figure and one table. If any of the results are statistically significant, I can explore the question further with a more specific model in the future. 

cox_all_var = coxph(formula = SurvObj ~ drug + sex + bil + as.factor(histo) + as.factor(agecat), data = pbcData)
all_var_summary <- summary(cox_all_var)
library(broom, help)
cox_all_var %>%
tidy()
coef(cox_all_var) %>%
summary()
df <- data_frame(adj_HR = exp(coef(cox_all_var)),
           lower = confint(cox_all_var)[,1],
           upper = confint(cox_all_var)[,2])
library(knitr)
knitr::kable(df, format = "markdown")
rownames(confint(cox_all_var))
#' The results of the model 

#' Plotting
km_drug_sex = survfit(SurvObj ~ drug + sex, data = pbcData,
type="kaplan-meier", conf.type="log-log")
km_age = survfit(SurvObj ~ drug + as.factor(agecat), data = pbcData,
type="kaplan-meier", conf.type="log-log")
km_histo = survfit(SurvObj ~ drug + as.factor(histo), data = pbcData,
type="kaplan-meier", conf.type="log-log")

# histoplot
plot(km_histo, col = 1:8, xlab = "Time", ylab = "Survival")
legend("bottomleft",
       legend=names(km_histo$strata),
       col=1:8,
       cex=0.5,
       lty=c(1,1), # gives the legend appropriate symbols (lines) 
       lwd=c(2.5,2.5))


#' To make a similar plot with the `bil` variable, I will create a new categorical variable called `bilcat`.

# bilplot
pbcData['bilcat'] <- ifelse(pbcData["bil"][[1]]>median(pbcData["bil"][[1]]), 1, 0)
head(pbcData)
km_bil = survfit(SurvObj ~ drug + as.factor(bilcat), data = pbcData,
type="kaplan-meier", conf.type="log-log")
