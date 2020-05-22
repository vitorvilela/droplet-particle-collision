# Import section
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import read_csv
from pandas import set_option
from pandas import plotting
import math
from scipy import stats 

print('\nFactorial desing for dimensionless lamella area at t*=0.918')
filenames = ['factorialDesign.csv'] 

print('\nReynolds (a) and Weber (b) 2^k design')
names = ['replicate', 'reynolds', 'weber', 'time', 'area*']

print('\nk=2 factors')
k = 2

print('\nn=3 replicates')
n = 3

dataset_original = read_csv(filenames[0], names=names)
print(dataset_original)

# the symbols (1), a, b, and ab represent the total of the response observation at all n replicates taken at the treatment combination
list_o = []
list_o.append(dataset_original.at[0, 'area*'])
list_o.append(dataset_original.at[4, 'area*'])
list_o.append(dataset_original.at[8, 'area*'])
o = sum(list_o)
print('\nlist_o')
print(list_o)
print(o)

list_b = []
list_b.append(dataset_original.at[1, 'area*'])
list_b.append(dataset_original.at[5, 'area*'])
list_b.append(dataset_original.at[9, 'area*'])
b = sum(list_b)
print('\nlist_b')
print(list_b)
print(b)

list_a = []
list_a.append(dataset_original.at[2, 'area*'])
list_a.append(dataset_original.at[6, 'area*'])
list_a.append(dataset_original.at[10, 'area*'])
a = sum(list_a)
print('\nlist_a')
print(list_a)
print(a)

list_ab = []
list_ab.append(dataset_original.at[3, 'area*'])
list_ab.append(dataset_original.at[7, 'area*'])
list_ab.append(dataset_original.at[11, 'area*'])
ab = sum(list_ab)
print('\nlist_ab')
print(list_ab)
print(ab)

print('\nDataset list')
dataset_list = list_a + list_b + list_ab + list_o
print(dataset_list)


# The effects of interest in the 2^2 design are the main effects A and B and the two-factor interaction AB
print('\nEffect denominator (ED)')
ED = n*pow(2, k-1)
print(ED)

# To estimate the main effect of A, we would average the observations on the right side of the square where A is at the high level, and subtract from this the average of the observations on the left side of the square where A is at the low level
print('\nContrast of A')
CA = (a + ab - b - o)
print(CA)
print('Factor A main effect (Reynolds)')
A = CA / ED
print(A)

print('\nContrast of B')
CB = (b + ab - a - o)
print(CB)
print('Factor B main effect (Weber)')
B = CB / ED
print(B)

# The AB interaction is estimated by taking the difference in the diagonal averages
print('\nContrast of AB')
CAB = (ab + o - a - b)
print(CAB)
print('AB effect (Reynolds-Weber interaction)')
AB = CAB / ED
print(AB)




# Analyses of Variance
print('\nVariance denominator (VD)')
VD = n*pow(2, k)
print(VD)

print('\nSum of squares of A')
SSA = CA*CA / VD 
print(SSA)

print('\nSum of squares of B')
SSB = CB*CB / VD 
print(SSB)

print('\nSum of squares of AB')
SSAB = CAB*CAB / VD 
print(SSAB)



print('\nTotal sum of squares')
SST = (a + b + ab + o) / VD
print('Old SST:')
print(SST)

dataset_mean = np.mean(dataset_list)
SST = 0
for d in dataset_list:
  SST += (d - dataset_mean) * (d - dataset_mean)
print('New SST:')  
print(SST)

dataset_sum = sum(dataset_list)
SST = 0
for d in dataset_list:
  SST += d*d
SST = SST - (dataset_sum*dataset_sum)/VD  
print('And another way to compute SST:')  
print(SST)



print('\nError sum of squares')
SSE = SST - SSA - SSB - SSAB
print(SSE)



# One degree of freedom is associated with each effect (two levels minus one) 
# so that the mean squared error of each effect equals the sum of squares
print('\nEffects degrees of freedom')
DOF_EFFECT = 1
print(3*DOF_EFFECT)

print('\nMean square of effect A')
MSA = SSA / DOF_EFFECT
print(MSA)

print('\nMean square of effect B')
MSB = SSB / DOF_EFFECT
print(MSB)

print('\nMean square of effect AB')
MSAB = SSAB / DOF_EFFECT
print(MSAB)

# The analysis of variance is completed by computing the total sum of squares SST
# (with 4n-1 degrees of freedom) as usual, and obtaining the error sum of squares SSE
# (with 4(n-1) − degrees of freedom) by subtraction
print('\nTotal sum of squares degrees of freedom')
DOF_SST = 4*n-1
print(DOF_SST)

print('\nError sum of squares degrees of freedom')
DOF_SSE = 4*(n-1)
print(DOF_SSE)

print('\nMean square error')
MSE = SSE / DOF_SSE
print(MSE)

print('\nThe statement H0 :mu1 = m2 is called the null hypothesis and H1 :mu1 != m2 is called the alternative hypothesis. The alternative hypothesis specified here is called a two-sided alternative hypothesis because it would be true if m1 < m2 or if m1 > m2.')
















print('\nModels and Residual analysis')
# Y = B0 + B1.x1 + B2.x2 + B3.x1x2 + e

# The intercept B0 is the grand average of all 12 observations
print('\nConstant')
B0 = (a + b + ab + o)/(n*pow(2, k))
print(B0)

# The slope B1 is one-half the effect estimate of A
print('\nCoefficient of A')
B1 = A / 2
print(B1)

print('\nCoefficient of B')
B2 = B / 2
print(B2)

print('\nCoefficient of AB')
B3 = AB / 2
print(B3)

print('\nThe standard error of a coefficient')
SE = math.sqrt(MSE/VD)
print(SE)


# The t-test also tells you how significant the differences are.
# In other words it lets you know if those differences could have happened by chance

# A One sample t-test tests the mean of a single group against a known mean
# The t-statistic to test H0 : B = 0

# Using a significance level of 0.2 [80% confidence level] (0.1 for each direction in a two-tailed test).
# From table: t-student critical (8 dof) = 1.397

print('\nUsing a significance level of 0.1 [80% confidence level] (0.1 for the only direction in a single-tailed test).')
print('From table: t-student critical (8 dof) = 1.397')
print('which is lesser than 3.280 for factor B (Weber) : can state it is a significant coefficient')
print('but is greater than 1.329 for factor A (Reynolds); the same for AB 0.220 : cannot state it is a significant coefficient')

print('\nDegrees of freedom associated with mean square error')
print(DOF_SSE)

print('\nt-statistic for factor A (Reynolds)')
TB1 = abs(B1/SE)
print(TB1)

print('\nt-statistic for factor B (Weber)')
TB2 = abs(B2/SE)
print(TB2)

print('\nt-statistic for AB (Reynolds-Weber)')
TB3 = abs(B3/SE)
print(TB3)

# A p-value is the probability that the results from your sample data occurred by chance

# When using a two-tailed test, regardless of the direction of the relationship you hypothesize, you are testing for the possibility of the relationship in both directions.  For example, we may wish to compare the mean of a sample to a given value x using a t-test. Our null hypothesis is that the mean is equal to x. A two-tailed test will test both if the mean is significantly greater than x and if the mean significantly less than x. The mean is considered significantly different from x if the test statistic is in the top 2.5% or bottom 2.5% of its probability distribution, resulting in a p-value less than 0.05.

# Using a significance level of 0.2, we should expect a p-value (2*stats.t.cdf()) less than 0.2

# Our null hypothesis is that the mean is equal to x. A one-tailed test will test either if the mean is significantly greater than x or if the mean is significantly less than x, but not both. Then, depending on the chosen tail, the mean is significantly greater than or less than x if the test statistic is in the top 5% of its probability distribution or bottom 5% of its probability distribution, resulting in a p-value less than 0.05.  The one-tailed test provides more power to detect an effect in one direction by not testing the effect in the other direction.

# Our null hypothesis is that the coefficient B is greater than zero. 
# Using a significance level of 0.1 [80% confidence level] (0.1 for the only direction in a single-tailed test), we should expect a p-value (stats.t.cdf()) less than 0.1
print('\np-value for B1 (Reynolds)')
PB1 = 1 - stats.t.cdf(TB1, df=DOF_SSE)
print(PB1)

print('\np-value for B2 (Weber)')
PB2 = 1 - stats.t.cdf(TB2, df=DOF_SSE)
print(PB2)

print('\np-value for B3 (Reynolds-Weber)')
PB3 = 1 - stats.t.cdf(TB3, df=DOF_SSE)
print(PB3)


# Cross-validation
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_1samp.html#scipy.stats.ttest_1samp

print('\nttest_ind')
# Calculate the T-test for the means of two independent samples of scores.
# This is a two-sided test for the null hypothesis that 2 independent samples have identical average (expected) values. This test assumes that the populations have identical variances by default.
# We can use this test, if we observe two independent samples from the same or different population, e.g. exam scores of boys and girls or of two ethnic groups. The test measures whether the average (expected) value differs significantly across samples. If we observe a large p-value, for example larger than 0.05 or 0.1, then we cannot reject the null hypothesis of identical average scores. If the p-value is smaller than the threshold, e.g. 1%, 5% or 10%, then we reject the null hypothesis of equal averages.
t_re, p_re = stats.ttest_ind(np.asarray(list_a + list_ab), np.asarray(list_b + list_o))
print("t_re = " + str(t_re))
print("p_re = " + str(2*p_re))
t_we, p_we = stats.ttest_ind(np.asarray(list_b + list_ab), np.asarray(list_a + list_o))
print("t_we = " + str(t_we))
print("p_we = " + str(2*p_we))

print('\nttest_ind_from_stats')
# T-test for means of two independent samples from descriptive statistics.
# This is a two-sided test for the null hypothesis that two independent samples have identical average (expected) values.
result_we = stats.ttest_ind_from_stats(mean1=np.mean(np.asarray(list_b + list_ab)), std1=np.std(np.asarray(list_b + list_ab)), nobs1=6, mean2=np.mean(np.asarray(list_a + list_o)), std2=np.std(np.asarray(list_a + list_o)), nobs2=6)
print("t_we = " + str(result_we[0]))
print("p_we = " + str(2*result_we[1]))

print('\nttest_rel')
# Calculate the T-test on TWO RELATED samples of scores, a and b.
# This is a two-sided test for the null hypothesis that 2 related or repeated samples have identical average (expected) values.
# Examples for the use are scores of the same set of student in different exams, or repeated sampling from the same units. The test measures whether the average score differs significantly across samples (e.g. exams). 
t_re, p_re = stats.ttest_rel(np.asarray(list_a + list_ab), np.asarray(list_b + list_o))
print("t_re = " + str(t_re))
print("p_re = " + str(2*p_re))
t_we, p_we = stats.ttest_rel(np.asarray(list_b + list_ab), np.asarray(list_a + list_o))
print("t_we = " + str(t_we))
print("p_we = " + str(2*p_we))

print('\nttest_1samp')
# Calculate the T-test for the mean of ONE group of scores.
# This is a two-sided test for the null hypothesis that the expected value (mean) of a sample of independent observations a is equal to the given population mean, popmean.
popmean = np.mean(np.asarray(list_a + list_o))
t_we, p_we = stats.ttest_1samp(np.asarray(list_b + list_ab), popmean)
print("t_we = " + str(t_we))
print("p_we = " + str(2*p_we))



