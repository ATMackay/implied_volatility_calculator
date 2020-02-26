"""
Mako Interview Question.

This programme takes an input file containing trading data, calculates the implied volatility for each 
trade and ouputs a new csv file.

# Requirements: Python 3.6 (or later), numpy, pandas, scipy

Installation of NumPy, SciPy and Pandas using pip --> $ python3 -m pip install --user numpy scipy pandas

"""

import numpy as np
import pandas as pd
from scipy.stats import norm



"""
----------------------------------------------------------------------------------------------------------
            Import CSV
----------------------------------------------------------------------------------------------------------
"""

def import_csv(input_csv):
    """
    Use Pandas library to import data from file, check formatting.
    The input must be a string --> "<filename>.csv"
    """ 
    df = pd.read_csv(input_csv, na_values = ['no info', '.'], index_column = 0)
    if len(df[0]) != 9:
        raise Exception("Input csv file must contain 9 colums: 'ID', Underlying Type',\
                         'Underlying', 'Risk-free rate', 'Days To Expiry', 'Strike', \
                         'Option Type','Model Type, 'Market Price'. ")
    
    return df


"""
------------------------------------------------------------------------------------------------------------
                         MATHEMATICAL FUNCTIONS
------------------------------------------------------------------------------------------------------------
    
    VARIABLES
    ---------
S - Market price of the underlying asset
K - Strike Price
r - Risk-free interest rate
q - Dividend rate
sig - Market volatility
T - Maturity date

BlkSc_call - Computed call option price (Black Scholes model)
BlkSc_put - Computed put option price

Bach_call - Computed call option price (Bachelor model)
Bach_put - Computed put option price
"""

def cum_dist(x):
    """
    Alternative method for computing the cumulative distribution function N(x).
    Currently we are using scipy.norm.cdf() for speed and accuracy. 
    """
    z = (39*x)/(2*np.sqrt(2*np.pi)) - 111./2*np.arctan(35*x/(111*np.sqrt(2*np.pi)))
    N = 0.5*np.tanh(z) + 0.5
    return N


def BlkSc_call(S, K, r, sig, T):
    """
    This function computes the market price of a CALL option (non-dividend paying) using the Black Scholes model
    Reference: 
    (optional) replace norm.cdf() with cum_dist() to avoid dependency on SciPy.
    """
    d_one = (np.log(S/K) + (r + 0.5*sig*sig) * T)/(sig * np.sqrt(T))
    d_two = d_one - sig*np.sqrt(T)

    return S*norm.cdf(d_one) - np.exp(-r*T)*K*norm.cdf(d_two)


def BlkSc_put(S, K, r, sig, T):
    """
    This computes the market price of a PUT option
    ...Reference: 
    Uses a formula relating call and put prices.
    """
    return BlkSc_call(S, K, r, sig, T) - S + np.exp(-r*T)*K # Need to check the variables are correct.


def Bach_call(S, K, r, sig, T):
    """
    This function computes the market price of a CALL option (non-dividend paying) using the Bachelier model.
    ... Formula Reference :  http://unriskinsight.blogspot.com/2013/10/black-vs-bachelier-revisited.html 
    (optional) replace norm.cdf() with cum_dist() to avoid dependency on SciPy.
    """

    d_one = (S - K)/(sig*np.sqrt(T))
    C = ( (S - K)*norm.cdf(d_one) + sig*np.sqrt(T)*norm.pdf(d_one) )*np.exp(-r*T)

    return C


# !!!!!!!!!!
def Bach_put(S, K, r, sig, T):
    """
    This computes the market price of a PUT option using the Bacehelier model.
    ...Reference for formula -  http://unriskinsight.blogspot.com/2013/10/black-vs-bachelier-revisited.html
    Uses a formula relating call and put prices.
    """

    d_one = (S - K)/(sig*np.sqrt(T))
    P = ( (K - S)*norm.cdf(-d_one) + sig*np.sqrt(T)*norm.pdf(d_one) )*np.exp(-r*T)

    return P

"""
# Check example call/put values (TESTING)
print("Black Scholes Call:", BlkSc_call(100, 90, 0.01, 0.1, 1))
print("Black Scholes Put:", BlkSc_put(100, 90, 0.01, 0.1, 1))
print("Bachelier Call", Bach_call(100, 90, 0.01, 0.1, 1))
print("Bachelier Put", Bach_put(100, 90, 0.01, 0.1, 1))
quit()
"""



"""
Black Scholes model volatility calculator.
"""
def BSMImVol(S, K, r, T, option_type):
    """ 
    This function solves the inverse pricing problem for market volatility based on the Black Scholes model of option pricing.
    The function iterates using the bisection method until the volatility solution is within tolerance and converges linearly
    (we could use a faster iteration method, such as Newton's method, but let's keep it simple for now).
    """
    if option_type == "Call":
        H = 5.
        L = 0.
        # Start bisection method
        while (H - L) > 10e-8: 
            if BlkSc_call(S, K, r, (H+L)/2, T) > call:
                H = (H + L)/2.
            else:
                L = (H + L)/2.

        return BlkSc_call(S, K, r, (H+L)/2, T)

    elif option_type == "Put":
        H = 5.
        L = 0.
        # Start bisection method
        while (H - L) > 10e-8: 
            if BlkSc_put(S, K, r, (H+L)/2, T) > call:
                H = (H + L)/2
            else:
                L = (H + L)/2

        return BlkSc_put(S, K, r, (H+L)/2, T)

    else:
        return np.float('nan')


"""
Bachelier model volatility calculator.
"""
def BAImVol(S, K, r, sig, T, option_type):
    """
    This function solves the inverse pricing problem for market volatility based on the Bachelier model of option pricing.
    The function iterates using the bisection method until the volatility solution is within tolerance and converges linearly.
    """
    if option_type == "Call":
        H = 5.
        L = 0.
        # Start bisection method
        while (H - L) > 10e-8: 
            if BlkSc_call(S, K, r, (H+L)/2, T) > call:
                H = (H + L)/2.
            else:
                L = (H + L)/2.

        return BlkSc_call(S, K, r, (H+L)/2, T)

    elif option_type == "Put":
        H = 5.
        L = 0.
        # Start bisection method
        while (H - L) > 10e-8: 
            if Bach_put(S, K, r, (H+L)/2, T) > call:
                H = (H + L)/2
            else:
                L = (H + L)/2

        return Bach_put(S, K, r, (H+L)/2, T)

    else:
        return np.float('nan')



"""
-----------------------------------------------------------------------------------------------------------------------------
                         MAIN FUNCTION
-----------------------------------------------------------------------------------------------------------------------------
"""

# !!!!!!!
def main(input_csv):
    """
    This function takes a csv file containing market trade data, computes implied 
    volatility and writes a new csv file.  The input must be a string --> "<filename>.csv".
    thE output will be output.csv.
    """
    # check first few rows of input csv (TESTING)
    input_data = import_csv(input_csv) 
    print(input_data.head())  

    """
    df_new = pd.DataFrame(
    {
        'ID': []
        'Strike': []
        'Risk-Free Rate': []
        'Years to Expiry': []
        'Option Type': []
        'Model type': []
        'Implied Volatility': []
        'Market Price': []   
    }   
    )        
    """
    new_data = []
    # Iterate over rows in input csv
    for i in range(data.shape[0]-1):
        
        if data['Underlying type'][i+1] == 'Stock':
            spot = np.float('nan')
            # Still need to find formula/definition SOLUTION
        elif data['Underlying type'][i+1] == 'Future':
            spot = data['Underlying'][i+1]
        else: 
            spot = np.float('nan')

        # May be better to use variables directly 
        _id = data['ID'][i+1]
        strike = data['Strike'][i+1]
        risk_free = data['Risk-Free Rate'][i+1]
        expiry = data['Days To Expiry'][i+1]
        option_type = data['Option Type'][i+1]
        model_type = data['Model Type'][i+1]
        price = data['Market Price'][i+1] 

        if model_type == 'BlackScholes':   
            if option_type == 'Call':
                implied_volatility = BSMImVol(price, strike, risk_free, expiry, 'Call') 
            elif option_type == 'Put':
                implied_volatility = BSMImVol(price, strike, risk_free, expiry, 'Put')
        elif model_type == 'Bachelier':
            if option_type == 'Call':
                implied_volatility = BaMImVol(price, strike, risk_free, expiry,'Call') 
            elif option_type == 'Put':
                implied_volatility = BaMImVol(price, strike, risk_free, expiry, 'Put')
        else:
            implied_volatility = np.float('nan')


        # Create dictionary to store new row (with computed implied volatility) and add to list
        new_data.append({'ID': _id, 'Spot': spot, 'Strike': strike, 'Risk-Free Rate': risk_free, 'Years to Expiry': expiry, \
                     'Option Type': option_type, 'Model type': model_type, 'Implied Volatility':implied_volatility,\
                     'Market Price': price})
 
    #Create dataframe use list of dicts
    df_new = DataFrame(new_data)
   
    return df_new.to_csv("output.csv")
.

if __name__=='__main__':
    main("input.csv")
