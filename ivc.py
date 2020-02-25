# This programme calculates the implied volatility for....

# Interview question for Mako

# Requires: Python 3.6, Pandas, numpy, scipy

import unittest
import numpy as np
import pandas as pd
from scipy.stats import norm




# !!!!!!!! (NOT NEEDED)
def import_csv():
    df = pd.read_csv()
    return df

# !!!!!!! (NOT NEEDED?)
def cum_dist(x):
    z = (39*x)/2*np.sqrt(2*np.pi) - 111./2*np.arctan(35*x/(111*np.sqrt(2*np.pi)))
    N = 0.5*np.tanh(z) + 0.5
    return N


# !!!!!!!!

"""
S - 
K - 
r - 
q - 
sig - Market volatility
T - Maturity date
C - Computed call option market price 
P - Computed put option market price
"""
def BlkSc_call(S, K, r, q, sig, T):
    """This function computes the market price of a CALL option...ADD DESCRIPTION and REFERENCES"""

    d_one = (np.log(S/K) + (r - q + 0.5*sigma*sigma) * T)/(sigma * np.sqrt(T))
    d_two = d_one - sigma*np.sqrt(T)
    # (optional) replace norm.cdf() with cum_dist()
    N_d1 = norm.cdf(d_one)
    N_d2 = norm.cdf(d_two)

    C = np.exp(-q*T)*S*N_d1 - np.exp(-r*T)*K*N_d2 # Need to check that the formula is correct
    return C

"""Here we use the call option formula when defining the Black Scholes pricing function for a put option.  
"""
def BlkSc_put(S, K, r, q, sig, T):
    """This computes the market price of a PUT option...ADD DESCRIPTION and REFERENCES.
        Uses a formula relating call and put prices."""
    return BlkSc_call(S, K, r, q, sig, T) - S + np.exp(-r*T)*K # Need to check the variables are correct.

print(BlkSc_call(100, 100, 0.1, 0, 0.2, 1))

# !!!!!!!!


#!!!!!!!!
def BSMImVol(VARIABLES, option_type):
    """ This function solves the inverse pricing problem for market volatility based on the Black Scholes model of option pricing.
        Iterate until solution is tolerance using bisection method
       (we could use a faster iteration method but let's keep it simple for now).
    """
    if option_type == "Call":
        if BlkSc() > SOMETHING:
            return np.float('nan')
        H = 2.
        L = 0.
        # Start bisection method
        while (H - L) > 10e-8: 
            if BlkSc_call(VARIABLES, (H+L)/2) > call:
                H = (H + L)/2
            else:
                L = (H + L)/2

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


# !!!!!!!!!!
def BAImVol(VARIABLES, option_type):
    return 0






# !!!!!!!
def main(input_csv):
    data = pd.read_csv(input_csv) # Label rows and columns properly

    
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
