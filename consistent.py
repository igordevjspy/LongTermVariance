# This simple script aims at estimating long-term variance estimator using parametric and non-parametric methods. 
# The documentation is under proofreading and will be out soon
# More to come soon: 
    # - non parametric long run variance based on sample autocorrelation using different weighting functions
    # - Quadratic spectral Kerlal non-parametric long run variance
    # - Pre-whitened long-run variance estimator
    
import csv
import numpy 
from numpy import ndarray
from numpy import fft
import statsmodels
from statsmodels import tsa
import scikits.statsmodels.api as sm 
import scikits
import scipy
import math
import pandas as pd


########################################
## Parametric Long-Run variance Estimator
## based on an AR(p) specifications, with parameters estimate by OLS

data =pd.read_csv('inflation.csv', sep='\t')
inflation = data['inflation'].values
x = inflation


def AR_OLS(x,p):
    import pandas as pd
    import statsmodels.api as smf
    col = ['x_' + str(i) for i in range(p+1)]
    index = [i for i in range(len(x))]
    df = pd.DataFrame(columns=col, index=index)
    df[col[0]] = x
    for i in range(1,p+1):
        df[col[i]].values[i:len(x)] = x[0:len(x)-i] 
    df = df[df.index>p-1]
    print df
    df.index = [i for i in range(len(x)-6)]
    results = smf.OLS(df[col[0]],smf.add_constant(df[col[1:len(col)]])).fit()
    results.params
    fitted = results.predict()
    sigma_2 = sum( (df[col[0]].values - fitted)**2)/float(len(x))
    return sigma_2/float((1-sum(results.params[1:7]))**2)

## return the parametric long term variance based on a AR(p)
AR_OLS(x,6)

    
### 2. Use the parameters estimated by the best of a periodogranm to an implied AR(p) spectrum
global x, F_T, lambda_T, z

T = len(x)
l = 2/float(T)
x= inflation
F_T = [i + 1/float(2) for i in range(-T/2,T/2)] 
lambda_T = dot(2*float(math.pi)/T,F_T )

# Z^2 is our periodogram
def Z(x,l,T):    
    if l<0: 
        return  math.sqrt(2)*T**(-0.5)*sum([math.cos(2*math.pi*l*j/float(T))*x[j] for j in range(T)])
    if l>0:
        return math.sqrt(2)*T**(-0.5)*sum([math.sin(2*math.pi*l*j/float(T))*x[j] for j in range(T)])
    if l== 0:
        return T**(-0.5)*sum(x)

# our AR part   
def Phi(phi,x):
    return 1 - sum([phi[i]*x**(i+1) for i in range(len(phi))])
    
    

def fonction(phi):
    spectrum = [2*math.pi*Phi(phi,np.exp(-complex(0,1)*lamb))*Phi(phi,np.exp(complex(0,1)*lamb)) for lamb in lambda_T]
    f = [z[i]*spectrum[i] for i in range(len(F_T))]
    return real(sum(f)).tolist()
 

# z is our periodogram 
z = [Z(x,l,T)**2 for l in F_T]
init_guess = [0,0,0,0,0,0]

sol = scipy.optimize.minimize(fonction,init_guess)
params = sol.x








    
  
    
