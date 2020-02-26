"""
Mako Interview Question.

This programme takes an input file containing trading data, calculates the implied volatility for each trade and ouputs a new csv file.

# Requires: Python 3.6, numpy, pandas, scipy

Installation of NumPy, SciPy and Pandas using pip --> $ python3 -m pip install --user numpy scipy pandas

"""
import unittest
import numpy as np
import pandas as pd
from scipy.stats import norm



# !!!!!!!! (NOT NEEDED)
"""
UsePandas library to import data from file, check formatting.
The input must be a string --> "<filename>.csv"
""" 
def format_csv(input_csv):
    df = pd.read_csv(input_csv)
    if len(df[0]) != 9:
        raise Exception("Input csv file must contain 9 colums: 'ID', Underlying Type',\
                         'Underlying', 'Risk-free rate', 'Days To Expiry', 'Strike', \
                         'Option Type','Model Type, 'Market Price'. ")
    
    return df



"""
Alternative method for computing the cumulative distribution function N(x).
Currently we are using scipy.norm.cdf() for speed and accuracy. 
"""
def cum_dist(x):
    z = (39*x)/2*np.sqrt(2*np.pi) - 111./2*np.arctan(35*x/(111*np.sqrt(2*np.pi)))
    N = 0.5*np.tanh(z) + 0.5
    return N


"""
            MATHEMATICAL FUNCTIONS
S - 
K - 
r - 
q - 
sig - Market volatility
T - Maturity date
C - Computed call option market price 
P - Computed put option market price
"""
# !!!!!!!!
def BlkSc_call(S, K, r, q, sig, T):
    """
    This function computes the market price of a CALL option...ADD DESCRIPTION and REFERENCES
    """
    d_one = (np.log(S/K) + (r - q + 0.5*sigma*sigma) * T)/(sigma * np.sqrt(T))
    d_two = d_one - sigma*np.sqrt(T)
    # (optional) replace norm.cdf() with cum_dist()
    N_d1 = norm.cdf(d_one)
    N_d2 = norm.cdf(d_two)

    C = np.exp(-q*T)*S*N_d1 - np.exp(-r*T)*K*N_d2 # Need to check that the formula is correct
    return C

"""
Here we use the call option formula when defining the Black Scholes pricing function for a put option.  
"""
def BlkSc_put(S, K, r, q, sig, T):
    """
    This computes the market price of a PUT option...ADD DESCRIPTION and REFERENCES.
    Uses a formula relating call and put prices.
    """
    return BlkSc_call(S, K, r, q, sig, T) - S + np.exp(-r*T)*K # Need to check the variables are correct.

print(BlkSc_call(100, 100, 0.1, 0, 0.2, 1))

# !!!!!!!!


#!!!!!!!!
"""
Black Scholes Model volatility calculator.
"""
def BSMImVol(VARIABLES, option_type):
    """ 
    This function solves the inverse pricing problem for market volatility based on the Black Scholes model of option pricing.
    The method iterates using the bisection method until the volatitlity solution is within tolerance and converges linearly
    (we could use a faster iteration method, such as Newton's method, but let's keep it simple for now).
    """
    if option_type == "Call":
        if BlkSc() > SOMETHING:
            return np.float('nan')
        H = 2.
        L = 0.
        # Start bisection method
        while (H - L) > 10e-8: 
            if BlkSc_call(VARIABLES, (H+L)/2) > call:
                H = (H + L)/2.
            else:
                L = (H + L)/2.

        return BlkSc_call(VARIABLES, (H+L)/2)

    elif option_type =="Put":
        if BlkSc_put() > SOMETHING:
            return np.float('nan')
        H = 2.
        L = 0.
        # Start bisection method
        while (H - L) > 10e-8: 
            if BlkSc_put(VARIABLES, (H+L)/2) > call:
                H = (H + L)/2
            else:
                L = (H + L)/2

        return BlkSc_put(VARIABLES, (H+L)/2)

    else:
        return np.float('nan')


"""
Bacheleier Model volatility calculator.
"""
# !!!!!!!!!!
def BAImVol(VARIABLES, option_type):
    return 0




# !!!!!!!
def main(input_csv):
    """
    This function takes a csv file containing market trade data, computes implied 
    volatility and writes a new csv file.  The input must be a string --> "<filename>.csv".
    teh output will be output.csv.
    """
    data = pd.read_csv(input_csv) # Label rows and columns properly
    if len(data[0]) != 7: # You may not need this if your import_csv function works
        raise Exception("Database not formatted correctly.")

    print(data)
    quit()   
    output_csv = 0
    """Iterate over rows in input csv"""
    for i in range(len()):
        # Be mindful of the number of if branches you are using
        if UNDERLYING_TYPE[i] == "Stock":
            SPOT[i] = x
        elif UNDERLYING_TYPE[i] == "Future":
            SPOT[i] = y
        else: 
            spot = np.float('nan')

        _id = data[i].id
        strike = data[i].strike
        risk_free = data[i].risk_free
        expiry = data[i].expiry
        model_type = data[i].model_type
        price = data[i].price # May be better to use variables directly 
    
        IMPLIED_VOLATILITY[i] = 0 # Use iterative function above 

        #ouput_csv[i] = CREATE_CSV --> [_id, SPOT, strike, risk_free, expiry, option_type, model_type, IMPLIED_VOLATILITY, price]
    
    return output_csv

if __name__=='__main__':
    main("input")
