# This programme calculates the implied volatility for....

# Interview question for Mako

# Requires: Python 3.6, Pandas, numpy, scipy

import unittest
import numpy as np
import pandas as pd
from scipy.stats import norm


pd.read_csv("input")

def import_csv():
    df = pd.read_csv()

    return df


def BlkSc_call(S, K, r, q, sig, T):
    """This computes the call..."""

    

    d_one = (np.log(S/K) + (r - q + 0.5*sigma*sigma) * T)/(sigma * np.sqrt(T))
    d_two = d_one - sigma*np.sqrt(T)
    N_d1 = norm.cdf(d_one)
    N_d2 = norm.cdf(d_two)

    return call

def BlkSc_put(S, K, r, q, sig, T):
    """This computes the put..."""

    return put

def BSMImVol():
    """Iterate until solution is .. tolerance using bisection method.\
       Evaluates the Black Scholes formula."""
    H = 2.
    L = 0.
    while (H - L) > 10e-8: 
    # Start bisection method
        if BlkSc(,,,(H+L)/2,) > call:
            H = (H + L)/2
        else:
            L = (H + L)/2

    return BlkSc(,,,(H+L)/2,)



def BAImVol():


def cum_dist(x):
    z = (39*x)/2*np.sqrt(2*np.pi) - 111./2*np.arctan(35*x/(111*np.sqrt(2*np.pi)))
    N = 0.5*np.tanh(z) + 0.5
    return N

def inverse_func():
    return inverse_func

def iterative_function():

    while tol > 10e-9:       
        v_2 = inverse_func(v_1)
        tol = v_2 - v_1
        v_1 = v_2


def main(input_csv):
    data = pd.read_csv(input_csv)

    # Create Empty row
    output_template = [EMPTY TEMPLATE]
    """Iterate over rows in input csv"""
    for i in range(len()):
        if UNDERLYING_TYPE[i] == "Call":
            spot[i] = x
        elif UNDERLYING_TYPE[i] == "Put":
            spot[i] = y
        else: 
            spot = np.float('nan')

        CREATE_CSV --> [_id, spot, strike, risk_free, years_to_expiry, option_type, model_type, implied_volatility, price]
    

    return output_csv

if __name__=='__main__':
    main("input")
