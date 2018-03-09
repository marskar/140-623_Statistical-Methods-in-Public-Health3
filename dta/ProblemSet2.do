codebook
tab sex
tab drug
tab histo
tab ageyr
tab agecat
tab agegr_2
tab agegr_3
summarize

sts list
sts graph, by(sex) t1("Male vs Female") b2("Years (t)") l2("S(t) = Probability Survival > t")
sts graph, by(drug) t1("Placebo vs DPCA") b2("Years (t)") l2("S(t) = Probability Survival > t")
sts graph, by(histo) t1("Histologic Disease Stage") b2("Stages 1-4") l2("S(t) = Probability Survival > t")
sts graph, by(agecat) t1("Age Categories") b2("Categories >45, 45-55, >55") l2("S(t) = Probability Survival > t")
sts test sex
sts test drug
sts test histo
sts test agecat
stphplot, by(sex) t1("Male vs Female") b2("Years (t)") l2("S(t) - Probability Survival > t")
stphplot, by(drug) t1("Placebo vs DPCA") b2("Years (t)") l2("S(t) - Probability Survival > t")
stphplot, by(histo) t1("Histologic Disease Stage") b2("Years (t)") l2("S(t) - Probability Survival > t")
stphplot, by(agecat) t1("Age Categories") b2("Years (t)") l2("S(t) - Probability Survival > t")

stcox drug
est store A

stcox sex
est store B

stcox i.histo
est store C

stcox i.agecat
est store D

stcox bil
est store E

stcox drug sex 
est store F
stcox drug hstage2 hstage3 hstage4
est store G
stcox drug i.agecat
est store H
stcox drug bil
est store I
stcox drug sex hstage2 hstage3 hstage4
est store J
stcox drug sex i.agecat
est store K
stcox drug sex bil
est store L
stcox drug sex hstage2 hstage3 hstage4 i.agecat bil
est store M
stcox sex hstage2 hstage3 hstage4 i.agecat bil
est store N
est stats *
