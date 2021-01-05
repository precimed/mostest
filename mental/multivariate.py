import sys
import pandas as pd
import scipy as sp
import scipy.linalg as la
from scipy.stats import chi2
import numpy as np

def mahalanobis(x=None, data=None, cov=None):

# """Compute the Mahalanobis Distance between each row of x and the data  
#    x    : vector or matrix of data with, say, p columns.
#    data : ndarray of the distribution from which Mahalanobis distance of each observation of x is to be computed.
#    cov  : covariance matrix (p x p) of the distribution. If None, will be computed from data.
# """
    x_minus_mu = x - np.mean(data)
    if not cov:
        cov = np.cov(data.values.T)
        inv_covmat = sp.linalg.inv(cov)
        left_term = np.dot(x_minus_mu, inv_covmat)
        mahal = np.dot(left_term, x_minus_mu.T)
        return mahal.diagonal()

filepath = sys.argv[1]
df = pd.read_csv(filepath, sep='\t')
nf = len(df.columns)
fd = nf - 3
df['PVAL'] = -1
stepwise = 10000
loopsize = len(df)//stepwise+(len(df)%stepwise > 0)
for l in range(0, loopsize):
    print(l)
    start = l * stepwise
    end = start + stepwise
    if end > len(df): end = len(df)
    df_x = df.iloc[start:end, 3:nf]
    df_x['mahala'] = mahalanobis(x=df_x, data=df[df.columns.tolist()[3:nf]])

    # Compute the P-Values
    df_x['p_value'] = 1 - chi2.cdf(df_x['mahala'], fd)
    df.iloc[start:end,nf] = df_x['p_value']

df.to_csv(filepath.split('.')[0] + '_multi.csv', index=None, sep='\t', header=True)
