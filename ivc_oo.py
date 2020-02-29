"""
Mako Interview Question.

This programme takes an input file containing trading data, calculates the implied volatility for each 
trade and ouputs a new csv file.

# Requirements: Python 3.6 (or later), numpy, pandas, scipy

Installation of NumPy, SciPy and Pandas using pip --> $ python3 -m pip install --user numpy scipy panda.

This Python code is OO. All solver functions are contained in the class 'Volsolver' which is instantiated .

"""

import numpy as np
import pandas as pd
from scipy.stats import norm
import sys


class VolSolver:
    """
    Python class containing mathematical functions for modelling call and put option prices.


    ------------------------------------------------------------------------ 
        VARIABLES
    ------------------------------------------------------------------------
    Y - Price of call/put option derived from Back Scholes/Bachelier formula
    S - Market price of the underlying asset
    K - Strike Price
    r - Risk-free interest rate
    sig - Implied volatility
    T - Maturity date

    f - pricing function
    df - partial derivative of pricing function
    """
    def __init__(self, op_type, model_type, Y, S, K, r, sig, T):       
        self.op_type = op_type
        self.model_type = model_type
        self.Y = Y
        self.S = S
        self.K = K
        self.r = r
        self.sig = sig
        self.T = T
        if self.model_type == 'BlackScholes':
            if self.op_type == 'Call':
                self.f = VolSolver.BSC_call
                self.df = VolSolver.BSC_call_dsig
            elif self.op_type == 'Put':
                self.f = VolSolver.BSC_put
                self.df = VolSolver.BSC_put_dsig
        elif self.model_type == 'Bachelier':
            if self.op_type == 'Call':
                self.f = VolSolver.BAC_call
                self.df = VolSolver.BAC_call_dsig
            elif self.op_type == 'Put':
                self.f = VolSolver.BAC_put
                self.df = VolSolver.BAC_put_dsig

    """
    -------------------------------------------------------------------------------------------
                 IMPLIED VOLATILITY SOLVER (NEWTON'S METHOD)
    -------------------------------------------------------------------------------------------
    """

    def newton_solver(self):
        """
        This function approximates the solution of f(x;.) = Y using Newton's method,
        f : f is a multivariable function taking 5 arguments,
            we are seeking a solution f(S, K, r, x, T) = Y searching for a solution in x space,
            the remaining parameters are fixed,

        df : partial derivative of f w.r.t. x,
        args: (Y, S, K, r, T) where
                f(S, K, r, x, T) = Y,
        x0 : initial guess,
        tol : stopping criteria |f(x)| < tol,
        max_iter : Maximum number of iterations of Netwon's method.

        """
        # Hardcoded paramters --> tolerance and maximum iterations
        tol, max_iter = 10e-8, 100
        xn = self.sig
        for i in range(max_iter):
            fxn = self.f(self) - self.Y
            if abs(fxn) < tol:
                return xn
            dfxn = self.df(self)
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
    def BSC_call(self):
        """
        This function computes the market price of a CALL option (non-dividend paying) using the Black Scholes model
        Reference: https://www.journals.uchicago.edu/doi/10.1086/260062 (cited in Wiki)
        (optional) replace norm.cdf() with closed form approximation to avoid dependency on SciPy.
        """
        S, K, r, sig, T = self.S, self.K, self.r, self.sig, self.T
        d_one = (np.log(S/K) + (r + 0.5*sig*sig) * T)/(sig * np.sqrt(T))
        d_two = d_one - sig * np.sqrt(T)

        return S*norm.cdf(d_one) - np.exp(-r*T)*K*norm.cdf(d_two)

    def BSC_call_dsig(self):
        """
        Partial derivative of BSC_call w.r.t. sigma (used in Newton's method).
        Derivative of norm.cdf(x) w.r.t. x is norm.pdf(x).
        Derivatives are calculated using a chain and product rule.
        """
        S, K, r, sig, T = self.S, self.K, self.r, self.sig, self.T
        d_one = (np.log(S/K) + (r + 0.5*sig*sig) * T)/(sig * np.sqrt(T))
        d_one_dsig =  (0.5 * T * sig * sig - np.log(S/K) - r * T )/(sig * sig * np.sqrt(T))
        d_two = d_one - sig*np.sqrt(T)
        d_two_dsig = d_one_dsig - np.sqrt(T)
        return  S*norm.pdf(d_one)*d_one_dsig - np.exp(-r*T)*K*norm.pdf(d_two)*d_two_dsig


    def BSC_put(self):
        """
        This computes the market price of a PUT option using the Balck Scholes Model 
        ...Reference: https://www.journals.uchicago.edu/doi/10.1086/260062 (cited in Wiki)
        Uses a formula relating call and put prices.
        """
        S, K, r, sig, T = self.S, self.K, self.r, self.sig, self.T
        return VolSolver.BSC_call(self) - S + np.exp(-r*T)*K 

    def BSC_put_dsig(self):
        """
        Partial derivative of BSC_call w.r.t. sigma (used in Newton's method).
        Identical to the BSC_call derivative.
        Derivative of norm.cdf(x) w.r.t. x is norm.pdf(x).
        """
        return VolSolver.BSC_call_dsig(self)

    def BAC_call(self):
        """
        This function computes the market price of a CALL option (non-dividend paying) using the Bachelier model.
        ... Formula Reference :  http://unriskinsight.blogspot.com/2013/10/black-vs-bachelier-revisited.html 
        (optional) replace norm.cdf() with closed form approximation to avoid dependency on SciPy.
        """
        S, K, r, sig, T = self.S, self.K, self.r, self.sig, self.T
        d_one = (S - K)/(sig*np.sqrt(T))
        return ( (S - K)*norm.cdf(d_one) + sig*np.sqrt(T)*norm.pdf(d_one) )*np.exp(-r*T)

    def BAC_call_dsig(self):
        """
        Partial derivative of Bac_call w.r.t. sigma.
        Derivative of norm.cdf(x) w.r.t. x is norm.pdf(x).
        Derivatives are calculated using a chain and product rule.
        """
        S, K, r, sig, T = self.S, self.K, self.r, self.sig, self.T
        d_one = (S - K)/(sig*np.sqrt(T))
        d_one_sig = (K - S)/(sig*sig*np.sqrt(T))
        pdf_dx = -d_one*norm.pdf(d_one)
        return ( (S - K)*norm.pdf(d_one)*d_one_sig + np.sqrt(T)*norm.pdf(d_one) + sig*np.sqrt(T)*pdf_dx*d_one_sig )*np.exp(-r*T)

    def BAC_put(self):
        """
        This function computes the market price of a PUT option using the Bacehelier model.
        ...Reference for formula -  http://unriskinsight.blogspot.com/2013/10/black-vs-bachelier-revisited.html
        Uses a formula relating call and put prices.
        """
        S, K, r, sig, T = self.S, self.K, self.r, self.sig, self.T
        d_one = (S - K)/(sig*np.sqrt(T))
        return ( (K - S)*norm.cdf(-d_one) + sig*np.sqrt(T)*norm.pdf(d_one) )*np.exp(-r*T)

    def BAC_put_dsig(self):
        """
        Partial derivative of Bac_call w.r.t. sigma.
        Identical to the Bac_call derivative.
        """
        return VolSolver.BAC_call_dsig(self)



class OptionData(VolSolver):
    """
    Python class containing main routinefor computing option data.
    """

    def __init__(self, input_csv, chunk_size, progress_bar):
        #super().__init__()       
        self.input_csv = input_csv
        self.chunk_size = chunk_size
        self.progress_bar = progress_bar

    def start(self):
        """
        This function takes a csv file containing market trade data, computes implied 
        volatility and writes the data to a new csv file.  The input must be 
        (string, int, bool) --> e.g. ("<filename>.csv", 1000, True)
        Larger chunk size will use more memory.
        The output will be file --> "output.csv".
        """
        # Import csv using Pandas, assign Nan values to empty cells
        input_data = pd.read_csv(self.input_csv, na_values = ['no info', '.'] , index_col = 0) 
        rows = input_data.shape[0]
        # Initialize dataframe
        new_data = []
        df_new = pd.DataFrame(new_data, index = None, columns = col_list)
        df_new.to_csv("output.csv", header = True, index = False)

        # Volatility initial guess
        sig_init_guess = 0.5
        # Iterate over rows in input csv
        for i in range(rows): 
            if self.progress_bar == True:
                # Display progress in terminal
                progress = i/rows
                msg = "\r{0}: {1}%".format("Computing implied volatility", round(progress*100, 4))
                sys.stdout.write(msg)
                sys.stdout.flush()
     
            underlyting_type = input_data['Underlying Type'][i]
            underlying = input_data['Underlying'][i]
            strike = input_data['Strike'][i]
            risk_free = input_data['Risk-Free Rate'][i]
            years_to_expiry = input_data['Days To Expiry'][i]/365.
            option_type = input_data['Option Type'][i]
            model_type = input_data['Model Type'][i]
            op_market_price = input_data['Market Price'][i] 
            
            if underlyting_type == 'Stock':
                # Still need to find formula/definition, not 100% sure about this.
                spot = underlying
            elif underlyting_type == 'Future':
                spot = op_market_price
            else: 
                spot = np.float('nan')

            # Calculate impled volatility using Volsolver.newton_solver()
            implied_volatility = VolSolver(option_type, model_type, op_market_price, underlying, strike, risk_free, sig_init_guess, years_to_expiry).newton_solver()

            # Create dictionary to store new row (with computed implied volatility) and add to list
            new_data.append({'ID': i, 'Spot': spot, 'Strike': strike, 'Risk-Free Rate': risk_free, 'Years to Expiry': years_to_expiry, \
                         'Option Type': option_type, 'Model Type': model_type, 'Implied Volatility': implied_volatility,\
                         'Market Price': op_market_price})
            # Append dataframe to csv file every chunk_size entries
            if (i + 1) % self.chunk_size == 0:
                #Create dataframe use list of dicts
                df_new = pd.DataFrame(new_data, index = None, columns = col_list)
                df_new.to_csv("output.csv", mode = 'a', na_rep = 'NaN', header = False, index = False)
                # Clear list
                new_data = []

        df_new = pd.DataFrame(new_data, index = None, columns = col_list)
        df_new.to_csv("output.csv", mode = 'a', na_rep = 'NaN', header = False, index = False)

        

"""
-----------------------------------------------------------------------------------------------------------------------------
                                        MAIN FUNCTION
-----------------------------------------------------------------------------------------------------------------------------
"""

def main(input_csv, chunk_size, progress_bar):
    """
    This function takes a csv file containing market trade data, computes implied 
    volatility and writes the data to a new csv file.  The input must be 
    (string, int, bool) --> e.g. ("<filename>.csv", 1000, True)
    Larger chunk size will use more memory.
    The output will be file --> "output.csv".
    """
    # Import csv using Pandas, assign Nan values to empty cells
    input_data = pd.read_csv(input_csv, na_values = ['no info', '.'] , index_col = 0) 
    rows = input_data.shape[0]

    # Check that dataframe has the correct number of columns
    if input_data.shape[1] != 8:
        raise Exception("Input dataframe not formatted ocrrectly. Should be an array with 8 columns.")
    col_list = ['ID', 'Spot', 'Strike', 'Risk-Free Rate', 'Years to Expiry', 
                'Option Type', 'Model Type', 'Implied Volatility', 'Market Price']  
    # Initialize dataframe
    new_data = []
    df_new = pd.DataFrame(new_data, index = None, columns = col_list)
    df_new.to_csv("output.csv", header = True, index = False)

    # Volatility initial guess
    sig_init_guess = 0.5
    # Iterate over rows in input csv
    for i in range(rows): 
        if progress_bar == True:
            # Display progress in terminal
            progress = i/rows
            msg = "\r{0}: {1}%".format("Computing implied volatility", round(progress*100, 4))
            sys.stdout.write(msg)
            sys.stdout.flush()
 
        underlyting_type = input_data['Underlying Type'][i]
        underlying = input_data['Underlying'][i]
        strike = input_data['Strike'][i]
        risk_free = input_data['Risk-Free Rate'][i]
        years_to_expiry = input_data['Days To Expiry'][i]/365.
        option_type = input_data['Option Type'][i]
        model_type = input_data['Model Type'][i]
        op_market_price = input_data['Market Price'][i] 
        
        if underlyting_type == 'Stock':
            # Still need to find formula/definition, not 100% sure about this.
            spot = underlying
        elif underlyting_type == 'Future':
            spot = op_market_price
        else: 
            spot = np.float('nan')

        # Calculate impled volatility using Volsolver.newton_solver()
        implied_volatility = VolSolver(option_type, model_type, op_market_price, underlying, strike, risk_free, sig_init_guess, years_to_expiry).newton_solver()

        # Create dictionary to store new row (with computed implied volatility) and add to list
        new_data.append({'ID': i, 'Spot': spot, 'Strike': strike, 'Risk-Free Rate': risk_free, 'Years to Expiry': years_to_expiry, \
                     'Option Type': option_type, 'Model Type': model_type, 'Implied Volatility': implied_volatility,\
                     'Market Price': op_market_price})
        # Append dataframe to csv file every chunk_size entries
        if (i + 1) % chunk_size == 0:
            #Create dataframe use list of dicts
            df_new = pd.DataFrame(new_data, index = None, columns = col_list)
            df_new.to_csv("output.csv", mode = 'a', na_rep = 'NaN', header = False, index = False)
            # Clear list
            new_data = []

    df_new = pd.DataFrame(new_data, index = None, columns = col_list)
    df_new.to_csv("output.csv", mode = 'a', na_rep = 'NaN', header = False, index = False)
   
    
if __name__ == '__main__':
    OptionData("input.csv", chunk_size = 5000, progress_bar = True).start()
    #main("input.csv", chunk_size = 5000, progress_bar = True)
