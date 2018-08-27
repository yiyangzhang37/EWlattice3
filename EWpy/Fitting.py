import numpy as np
from scipy.stats import linregress

class PowerLawBase:
    '''
    Fitting for power-law-like relations.
    The relation can be transformed to a linear relation,
    With two fitting parameters: slope_ and intercept_.
    '''
    x_ = None
    y_ = None
    slope_ = None
    intercept_ = None
    stderr_ = None
    def __init__(self):
        return None
    
    def fit(self, x, y):
        self.x_ = np.array(x)
        self.y_ = np.array(y)
        lx, ly = self.data_transform()
        slope, intercept, rvalue, pvalue, stderr = linregress(lx, ly)
        self.slope_ = slope
        self.intercept_ = intercept
        self.stderr_ = stderr
    
    def data_transform(self):
        return None, None
    
    def predict(self, x):
        return None
    

class PowerLaw(PowerLawBase):
    '''y = a*x^k
    '''   
    def data_transform(self):
        return np.log(self.x_), np.log(self.y_)
        
    def predict(self, x):
        return np.exp(self.intercept_) * np.power(x, self.slope_)

class ExpPowerLaw(PowerLawBase):
    '''y = exp(-a*x^k)
    '''
    def data_transform(self):
        return np.log(self.x_), np.log(-np.log(self.y_))
    
    def predict(self, x):
        return np.exp(-np.exp(self.intercept_)*np.power(x, self.slope_))

class InverseExp(PowerLawBase):
    '''y = a*exp(k/x^n), where n is an input parameter.
    '''
    n_ = 1.0
    def __init__(self, n = 1.0):
        self.n_ = n
        
    def data_transform(self):
        return 1.0/np.power(self.x_, self.n_), np.log(self.y_)
    
    def predict(self, x):
        return np.exp(self.intercept_)*np.exp(self.slope_/np.power(x, self.n_))