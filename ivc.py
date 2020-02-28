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
import sys

"""
------------------------------------------------------
            PROGRESS BAR
------------------------------------------------------
"""

def update_progress(job_title, progress):
    length = 50                               # Modify this to change the length of the bar
    block = int(round(length*progress))
    msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 4))
    if progress >= 1: msg += " DONE\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()


"""
------------------------------------------------------------------------------------------------------------
                                    MATHEMATICAL FUNCTIONS
------------------------------------------------------------------------------------------------------------
    
    VARIABLES
    ---------
S - Market price of the underlying asset
K - Strike Price
r - Risk-free interest rate
sig - Implied volatility
T - Maturity date

Y - Price of call/put option derived from Back Scholes/Bachelier formula
"""

"""
-------------------------------------------------------------------------------------------
             ITERATIVE ROOT FINDING 
-------------------------------------------------------------------------------------------
"""

def newton(f, df, *args):
    """
    This function approximates the solution of f(x;.) = Y using Newton's method,

    f : f is a multivariable function taking 5 arguments,
        we are seeking a solution f(S, K, r, x, T) = Y searching for a solution in x space,
        the remaining parameters are fixed,

    df : partial derivative of f w.r.t. x,
    args: --> (Y, S, K, r, T) where
                f(S, K, r, x, T) = Y,
    x0 : initial guess,
    tol : stopping criteria |f(x)| < tol,
    max_iter : Maximum number of iterations of Netwon's method.

    """
    # Hardcoded paramters --> initial guess, tolerance and maximum iterations
    xn, tol, max_iter = 1., 10e-8, 100
    Y, S, K, r, T = (s for s in args)
    for i in range(max_iter):
        fxn = f(S, K, r, xn, T) - Y
        if abs(fxn) < tol:
            return xn
        dfxn = df(S, K, r, xn, T)
        if dfxn == 0:
            return np.float('nan')
        xn = xn - fxn/dfxn   
    return np.float('nan')


"""
-----------------------------------------------------------------------------------------------------
             FINANCIAL MODELLING FUNCTIONS
-----------------------------------------------------------------------------------------------------

BSC_call - Black Scholes model call option price.
BSC_call_dsig - Black Scholes model call option derivative with respect to sigma (implied volatility).

BSC_put - Black Scholes model put option price.
BSC_put_dsig - Black Scholes model put option derivative with respect to sigma (implied volatility).

BAC_call - Bachelier model call option price.
BAC_call_dsig - Bachelier model call option derivative with respect to sigma (implied volatility).

BAC_put - Bachelier model put option price.
BAC_put_dsig - Bachelier model put option derivative with respect to sigma (implied volatility).

"""
def BSC_call(S, K, r, sig, T):
    """
    This function computes the market price of a CALL option (non-dividend paying) using the Black Scholes model
    Reference: https://www.journals.uchicago.edu/doi/10.1086/260062 (cited in Wiki)
    (optional) replace norm.cdf() with closed form approximation to avoid dependency on SciPy.
    """
    d_one = (np.log(S/K) + (r + 0.5*sig*sig) * T)/(sig * np.sqrt(T))
    d_two = d_one - sig * np.sqrt(T)

    return S*norm.cdf(d_one) - np.exp(-r*T)*K*norm.cdf(d_two)

def BSC_call_dsig(S, K, r, sig, T):
    """
    Partial derivative of BSC_call w.r.t. sigma (used in Newton's method).
    Derivative of norm.cdf(x) w.r.t. x is norm.pdf(x).
    """
    d_one = (np.log(S/K) + (r + 0.5*sig*sig) * T)/(sig * np.sqrt(T))
    d_one_dsig =  (0.5 * T * sig * sig - np.log(S/K) - r * T )/(sig * sig * np.sqrt(T))
    d_two = d_one - sig*np.sqrt(T)
    d_two_dsig = d_one_dsig - np.sqrt(T)
    return  S*norm.pdf(d_one)*d_one_dsig - np.exp(-r*T)*K*norm.pdf(d_two)*d_two_dsig


def BSC_put(S, K, r, sig, T):
    """
    This computes the market price of a PUT option using the Balck Scholes Model 
    ...Reference: https://www.journals.uchicago.edu/doi/10.1086/260062 (cited in Wiki)
    Uses a formula relating call and put prices.
    """
    return BSC_call(S, K, r, sig, T) - S + np.exp(-r*T)*K 

def BSC_put_dsig(S, K, r, sig, T):
    """
    Partial derivative of BSC_call w.r.t. sigma (used in Newton's method).
    Identical to the BSC_call derivative.
    Derivative of norm.cdf(x) w.r.t. x is norm.pdf(x).
    """
    return BSC_call_dsig(S, K, r, sig, T)

def BAC_call(S, K, r, sig, T):
    """
    This function computes the market price of a CALL option (non-dividend paying) using the Bachelier model.
    ... Formula Reference :  http://unriskinsight.blogspot.com/2013/10/black-vs-bachelier-revisited.html 
    (optional) replace norm.cdf() with closed form approximation to avoid dependency on SciPy.
    """
    d_one = (S - K)/(sig*np.sqrt(T))
    return ( (S - K)*norm.cdf(d_one) + sig*np.sqrt(T)*norm.pdf(d_one) )*np.exp(-r*T)

def BAC_call_dsig(S, K, r, sig, T):
    """
    Partial derivative of Bac_call w.r.t. sigma.
    Derivative of norm.cdf(x) w.r.t. x is norm.pdf(x).
    """
    d_one = (S - K)/(sig*np.sqrt(T))
    d_one_sig = (K - S)/(sig*sig*np.sqrt(T))
    pdf_dx = -d_one*norm.pdf(d_one)
    return ( (S - K)*norm.pdf(d_one)*d_one_sig + np.sqrt(T)*norm.pdf(d_one) + sig*np.sqrt(T)*pdf_dx*d_one_sig )*np.exp(-r*T)

def BAC_put(S, K, r, sig, T):
    """
    This function computes the market price of a PUT option using the Bacehelier model.
    ...Reference for formula -  http://unriskinsight.blogspot.com/2013/10/black-vs-bachelier-revisited.html
    Uses a formula relating call and put prices.
    """
    d_one = (S - K)/(sig*np.sqrt(T))
    return ( (K - S)*norm.cdf(-d_one) + sig*np.sqrt(T)*norm.pdf(d_one) )*np.exp(-r*T)

def BAC_put_dsig(S, K, r, sig, T):
    """
    Partial derivative of Bac_call w.r.t. sigma.
    Identical to the Bac_call derivative.
    """
    return BAC_call_dsig(S, K, r, sig, T)



"""
--------------------------------------------------------------------------------------------------------------
                 IMPLIED VOLATILITY SOLVERS
--------------------------------------------------------------------------------------------------------------
"""

"""
Black Scholes model volatility calculator.
"""
def BSC_imvol(Y, S, K, r, T, op_type):
    """ 
    This function solves the inverse pricing problem for market volatility based on the 
    Black Scholes model of option pricing. The function iterates over IV parameter
    space using Newton's method until the volatility solution is within tolerance (hard coded).
    """
    if op_type == 'Call':
        return newton(BSC_call, BSC_call_dsig, Y, S, K, r, T)

    elif op_type == 'Put':
        return newton(BSC_put, BSC_put_dsig, Y, S, K, r, T)

    else:
        return np.float('nan')


"""
Bachelier model volatility calculator.
"""
def BAC_imvol(Y, S, K, r, T, op_type):
    """
    This function solves the inverse pricing problem for market volatility based 
    on the Bachelier model of option pricing. The function iterates over IV paramter
    space using Newton's method until the volatility solution is within tolerance (hard coded).
    """
    if op_type == "Call":
        return newton(BAC_call, BAC_call_dsig, Y, S, K, r, T)
    elif op_type == "Put":
        return newton(BAC_put, BAC_put_dsig, Y, S, K, r, T)
    else:
        return np.float('nan')

"""
-----------------------------------------------------------------------------------------------------------------------------
                                        MAIN FUNCTION
-----------------------------------------------------------------------------------------------------------------------------
"""

def main(input_csv, chunk_size, progress_bar):
    """
    This function takes a csv file containing market trade data, computes implied 
    volatility and writes the data to a new csv file.  The input must be 
    (string, int, bool) --> ("<filename>.csv", 1000, True)
    Larger chunk size will use more memory.
    The output will be "output.csv".
    """
    # Import csv using Pandas, assign Nan values to empty cells
    input_data = pd.read_csv(input_csv, na_values = ['no info', '.'] , index_col = 0) 
    rows = input_data.shape[0]

    # Check that dataframe has the correct number of columns
    if input_data.shape[1] != 8:
        raise Exception("Input dataframe not formatted ocrrectly. Should be an array with 9 columns.")
    col_list = ['ID', 'Spot', 'Strike', 'Risk-Free Rate', 'Years to Expiry', 
                'Option Type', 'Model Type', 'Implied Volatility', 'Market Price']  
    # Initialize dataframe
    new_data = []
    df_new = pd.DataFrame(new_data, index = None, columns = col_list)
    df_new.to_csv("output.csv", mode = 'a', header = False)
    # Iterate over rows in input csv
    for i in range(rows): 
        if progress_bar == True:
            update_progress("Computing implied volatility", i/rows)
 
        underlyting_type = input_data['Underlying Type'][i]
        underlying = input_data['Underlying'][i]
        strike = input_data['Strike'][i]
        risk_free = input_data['Risk-Free Rate'][i]
        years_to_expiry = input_data['Days To Expiry'][i]/365.
        option_type = input_data['Option Type'][i]
        model_type = input_data['Model Type'][i]
        op_market_price = input_data['Market Price'][i] 
        
        if underlyting_type == 'Stock':
            # Still need to find formula/definition, not 100% sure about this
            spot = underlying
        elif underlyting_type == 'Future':
            spot = op_market_price
        else: 
            spot = np.float('nan')


        if model_type == 'BlackScholes':   
            implied_volatility = BSC_imvol(op_market_price, underlying, strike, risk_free, years_to_expiry, option_type) 
        elif model_type == 'Bachelier':
            implied_volatility = BAC_imvol(op_market_price, underlying, strike, risk_free, years_to_expiry, option_type) 
        else:
            implied_volatility = np.float('nan')

        # Create dictionary to store new row (with computed implied volatility) and add to list
        new_data.append({'ID': i, 'Spot': spot, 'Strike': strike, 'Risk-Free Rate': risk_free, 'Years to Expiry': years_to_expiry, \
                     'Option Type': option_type, 'Model Type': model_type, 'Implied Volatility': implied_volatility,\
                     'Market Price': op_market_price})
        # Append dataframe to csv file every chunk_size entries
        if (i + 1) % chunk_size == 0:
            #Create dataframe use list of dicts
            df_new = pd.DataFrame(new_data, index = None, columns = col_list)
            df_new.to_csv("output.csv", mode = 'a')
            # Clear list
            new_data = []

    df_new = pd.DataFrame(new_data, index = None, columns = col_list)
    df_new.to_csv("output.csv", mode = 'a')
   
    
if __name__ == '__main__':
    main("input.csv", chunk_size=5000, progress_bar = True)
