"""
Mako Interview Question.

This performs unit tests for ivc.py

# Requirements: Python 3.6 (or later), unittest, ivc

# Inputs for functions have been randomly chosen...

"""

import unittest
import numpy as np
import ivc


# Example functions for testing Newton solver
def g(S, K, r, x, T):
    return S + K + r + x*x*x + np.sin(0.5*x) + T
def g_dash(S, K, r, x, T):
    # Derivative of g w.r.t. x
    return  3*x*x + 0.5*np.cos(0.5*x)


class Testivc(unittest.TestCase):

    def test_BSC(self):
        # Test Black Scholes model functions
        r, sig, T = 0.01, 0.1, 0.01 # Constant risk-free, volatility, expiry date
        S, K  = 1.1, 0.9 # Call inputs
        call_ans = 0.20008999550015005
        self.assertAlmostEqual(ivc.BSC_call(S, K, r, sig, T), call_ans)
        S, K = K, S # Put inputs
        put_ans = 0.19989000549981684
        self.assertAlmostEqual(ivc.BSC_put(S, K, r, sig, T), put_ans)
        

    def test_BAC(self):
        # Test Bachelier model functions
        r, sig, T = 0.01, 0.1, 0.01
        S, K = 1.1, 0.9              # Call inputs
        call_ans = 0.19998000099996674
        self.assertAlmostEqual(ivc.BAC_call(S, K, r, sig, T), call_ans)

        S, K = K, S              # Put inputs
        put_ans = 0.19998000099996674
        self.assertAlmostEqual(ivc.BAC_put(S, K, r, sig, T), put_ans)


    def test_BSCvol(self):
        # Test Black Scholes call option IV
        iv_actual = 0.5
        r, T = 0.01, 0.01
        call_ans, S, K = 0.10038900828135255, 1.0, 0.9 
        self.assertAlmostEqual(ivc.imvol('Call', 'BlackScholes', call_ans, S, K, r, T), iv_actual)
        
        # Test Black Scholes put option IV
        put_ans, S, K = 0.10020234665104932, 0.9, 1.0
        self.assertAlmostEqual(ivc.imvol('Put', 'BlackScholes', put_ans, S, K, r, T), iv_actual)

    def test_BACvol(self):
        # Test Bachelier  call option IV
        iv_actual = 0.5
        r, T = 0.01, 0.1
        call_ans, S, K = 0.12516446483551702, 1.0, 0.9
        self.assertAlmostEqual(ivc.imvol('Call', 'Bachelier', call_ans, S, K, r, T), iv_actual)

        # Test Bachelier put option IV
        put_ans, S, K = 0.12516446483551702, 0.90, 1.0  
        f, df = ivc.BAC_put, ivc.BAC_put_dsig     
        self.assertAlmostEqual(ivc.imvol('Put', 'Bachelier', put_ans, S, K, r, T), iv_actual)

    def test_newton(self):
        """
        Find solution to g(100, 100, 1, x, 1) = 250 using the Newton solver,
        g is the (transcendental) function defined above
        """   
        ans = 3.6095232769758385 # Wolfram Alpha 
        x0 = 10
        self.assertAlmostEqual(ivc.newton(g, g_dash, x0, 250, 100, 100, 1, 1), ans)


if __name__ == 'main':
    unittest.main()



