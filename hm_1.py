""" Visiting the doctor
    ------------------
"""

# -----------------------------------------------------------------------------------
# 0. Import libraries
# -----------------------------------------------------------------------------------


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm


# -----------------------------------------------------------------------------------
# A. Import data and compute summary statistics for each variable including sample
#    count, mean, standard deviation, minimum, maximum, 25th 50th, and 75th quantiles.
# -----------------------------------------------------------------------------------


# Import dataset

data_original = pd.read_csv("PS1_insurance.csv")
data = data_original

# Fix variables

data["lottery"] = data["lottery"].replace(["Not selected", "Selected", ], [0, 1])
data["female"] = data["female"].replace(["0: Male", "1: Female"], [0, 1])

# Compute summary statistics

summary_stats = data.describe().round(2)
print(summary_stats)

# Export summary statistics to latex

print(summary_stats.to_latex(index=True))


# -----------------------------------------------------------------------------------
# B. Generate a new variable for each person's age (as of 2010), and a dummy variable
#    that equals one if the person had at least one doctor visit after the lottery.
#    Summarize these variables.
# -----------------------------------------------------------------------------------


# Generate variable for age

data["age"] = 2010 - data["birthyear"]

# Generate dummy variable for at least one doctor visit after lottery

data["doc_visit_dummy"] = (data["doc_visit_num"] > 0).astype(int)

# Compute summary statistics

summary_stats_2 = data[["age", "doc_visit_dummy"]].describe().round(2)
print(summary_stats_2)

# Export summary statistics to latex

print(summary_stats_2.to_latex(index=True))


# -----------------------------------------------------------------------------------
# C. To visualize the income distribution of the sample, plot a histogram with bin-
#     widths of $2500.
# -----------------------------------------------------------------------------------


# Generate the bins

bin_values = np.arange(np.min(data["hh_inc"]), np.max(data["hh_inc"]) + 2500, 2500)

# Plot the histogram

plt.figure()
plt.hist(data["hh_inc"], bins=bin_values, edgecolor='black', linewidth=1)
plt.title("Income distribution")
plt.xlabel("Income")
plt.ylabel("Absolute frequency")
plt.tight_layout()
plt.savefig('images/histogram_income_distribution.png')


# -----------------------------------------------------------------------------------
# D. Is the average person that was selected by the lottery different from those that
#    were not in terms of age? Or sex?
# -----------------------------------------------------------------------------------


# Separate selected and not selected groups

selected = data[data["lottery"] == 1]
not_selected = data[data["lottery"] == 0]

# Average age of selected and not selected groups

mean_age_selected = np.mean(selected["age"])
print("Mean age selected: ", np.round(mean_age_selected, 2))
mean_age_not_selected = np.mean(not_selected["age"])
print("Mean age not selected: ", np.round(mean_age_not_selected, 2))

# Average sex of selected and not selected groups

mean_sex_selected = np.mean(selected["female"])
print("Mean sex selected: ", np.round(mean_sex_selected, 2))
mean_sex_not_selected = np.mean(not_selected["female"])
print("Mean sex not selected: ", np.round(mean_sex_not_selected, 2))


# -----------------------------------------------------------------------------------
# E. If we are interested in the effect of this insurance availability on doctor
#    visits, what would be a potential concern about drawing conclusions from simply
#    comparing outcomes (i.e., number of doctor visits) of those above the income
#    threshold and those below?
# -----------------------------------------------------------------------------------


# See document


# -----------------------------------------------------------------------------------
# F. Regress the number of doctor visits on whether a person was selected by the
#    lottery. Interpret the results statistically and economically. Are these results
#    causal?
# -----------------------------------------------------------------------------------


# Regression of number of doctor visits on lottery selection

Y_1 = data["doc_visit_num"]

X_1 = data["lottery"]
X_1 = sm.add_constant(X_1)

ols_1 = sm.OLS(Y_1, X_1)
fit_1 = ols_1.fit()

summary_fit_1 = fit_1.summary()
print(summary_fit_1)

# Export regression to latex

print(summary_fit_1.as_latex())


# -----------------------------------------------------------------------------------
# G. Include age, sex, and household income to the right-hand side of the regression
#    above individually, and then all together. What happens to the coefficient on
#    lottery? Discuss.
# -----------------------------------------------------------------------------------


# Compute correlation matrix

correlation_matrix = data[["doc_visit_num", "lottery", "age", "female", "hh_inc"]].corr("pearson").round(2)
print(correlation_matrix)

# Export correlation matrix to latex

print(correlation_matrix.to_latex(index=True))

# Regression of number of doctor visits on lottery selection and age

Y_2 = data["doc_visit_num"]

X_2 = data[["lottery", "age"]]
X_2 = sm.add_constant(X_2)

ols_2 = sm.OLS(Y_2, X_2)
fit_2 = ols_2.fit()

summary_fit_2 = fit_2.summary()
print(summary_fit_2)

# Export regression to latex

print(summary_fit_2.as_latex())

# Regression of number of doctor visits on lottery selection and sex

Y_3 = data["doc_visit_num"]

X_3 = data[["lottery", "female"]]
X_3 = sm.add_constant(X_3)

ols_3 = sm.OLS(Y_3, X_3)
fit_3 = ols_3.fit()

summary_fit_3 = fit_3.summary()
print(summary_fit_3)

# Export regression to latex

print(summary_fit_3.as_latex())

# Regression of number of doctor visits on lottery selection and household income

Y_4 = data["doc_visit_num"]

X_4 = data[["lottery", "hh_inc"]]
X_4 = sm.add_constant(X_4)

ols_4 = sm.OLS(Y_4, X_4)
fit_4 = ols_4.fit()

summary_fit_4 = fit_4.summary()
print(summary_fit_4)

# Export regression to latex

print(summary_fit_4.as_latex())

# Regression of number of doctor visits on lottery selection, age, sex and household income

Y_5 = data["doc_visit_num"]

X_5 = data[["lottery", "age", "female", "hh_inc"]]
X_5 = sm.add_constant(X_5)

ols_5 = sm.OLS(Y_5, X_5)
fit_5 = ols_5.fit()

summary_fit_5 = fit_5.summary()
print(summary_fit_5)

# Export regression to latex

print(summary_fit_5.as_latex())


# -----------------------------------------------------------------------------------
# H. Statistically and economically interpret the other coefficient estimates,
#    including the constant.
# -----------------------------------------------------------------------------------


# See document


# -----------------------------------------------------------------------------------
# I. Conditional on the factors above, is the relationship between age and doctors
#    visits the same for men and women?
# -----------------------------------------------------------------------------------


# Regression of number of doctor visits on lottery selection, age household income for females

females = data[data["female"] == 1]

Y_6 = females["doc_visit_num"]

X_6 = females[["lottery", "age", "hh_inc"]]
X_6 = sm.add_constant(X_6)

ols_6 = sm.OLS(Y_6, X_6)
fit_6 = ols_6.fit()

summary_fit_6 = fit_6.summary()
print(summary_fit_6)

# Export regression to latex

print(summary_fit_6.as_latex())

# Regression of number of doctor visits on lottery selection, age household income for males

males = data[data["female"] == 0]

Y_7 = males["doc_visit_num"]

X_7 = males[["lottery", "age", "hh_inc"]]
X_7 = sm.add_constant(X_7)

ols_7 = sm.OLS(Y_7, X_7)
fit_7 = ols_7.fit()

summary_fit_7 = fit_7.summary()
print(summary_fit_7)

# Export regression to latex

print(summary_fit_7.as_latex())


# -----------------------------------------------------------------------------------
# J. Suppose you think that a linear relationship between age and doctors visits in
#    not reasonable. How does the number of doctors’ visits change with a 10% increase
#    in age.
# -----------------------------------------------------------------------------------


# Create the new variable age squared

data["age_squared"] = data["age"] ** 2

# Regression of number of doctor visits on age and age squared

Y_8 = data["doc_visit_num"]

X_8 = data[["lottery", "age", "age_squared"]]
X_8 = sm.add_constant(X_8)

ols_8 = sm.OLS(Y_8, X_8)
fit_8 = ols_8.fit()

summary_fit_8 = fit_8.summary()
print(summary_fit_8)

# Export regression to latex

print(summary_fit_8.as_latex())


# -----------------------------------------------------------------------------------
# K. Consider the specification: NumVisits = ß0 + ß1female + ß2 age + ß3 HHincome + u
#    Show that you recover the same estimate ß1 when you
#    i. estimate the full regression and
#    ii. regress the part of NumVisits unexplained by age and HHincome on the part of
#        female unexplained by age and HHincome.
# -----------------------------------------------------------------------------------


# Regression of number of doctor visits on sex, age and income

Y_9 = data["doc_visit_num"]

X_9 = data[["female", "age", "hh_inc"]]
X_9 = sm.add_constant(X_9)

ols_9 = sm.OLS(Y_9, X_9)
fit_9 = ols_9.fit()

summary_fit_9 = fit_9.summary()
print(summary_fit_9)

# Export regression to latex

print(summary_fit_9.as_latex())

# Regression of number of doctor visits on age and income

Y_10 = data["doc_visit_num"]

X_10 = data[["age", "hh_inc"]]
X_10 = sm.add_constant(X_10)

ols_10 = sm.OLS(Y_10, X_10)
fit_10 = ols_10.fit()

summary_fit_10 = fit_10.summary()
print(summary_fit_10)

# Export regression to latex

print(summary_fit_10.as_latex())

# Regression of sex on age and income

Y_11 = data["female"]

X_11 = data[["age", "hh_inc"]]
X_11 = sm.add_constant(X_11)

ols_11 = sm.OLS(Y_11, X_11)
fit_11 = ols_11.fit()

summary_fit_11 = fit_11.summary()
print(summary_fit_11)

# Export regression to latex

print(summary_fit_11.as_latex())

# Regression of number of doctor visits unexplained by age and income on sex unexplained by age and income

Y_12 = fit_10.resid

X_12 = fit_11.resid
X_12 = sm.add_constant(X_12)

ols_12 = sm.OLS(Y_12, X_12)
fit_12 = ols_12.fit()

summary_fit_12 = fit_12.summary()
print(summary_fit_12)

# Export regression to latex

print(summary_fit_12.as_latex())