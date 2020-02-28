# This programme calculates the implied volatility for....

# Interview question for Mako

# Requires: Python 3.6, Pandas, numpy

import unittest
import numpy as np
import ivc

"""
print("Black Scholes Call:", ivc.BSC_call(100, 90, 0.01, 0.1, 0.1))
print("Black Scholes Put:", ivc.BSC_put(90, 100, 0.01, 0.1, 0.1))
print("Bachelier Call", ivc.BAC_call(100, 90, 0.01, 0.1, 0.1))
print("Bachelier Put", ivc.BAC_put(90, 100, 0.01, 0.1, 0.1))
quit()
"""

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
        S, K  = 100, 90 # Call inputs
        call_ans = 10.090254370495515
        self.assertAlmostEqual(ivc.BSC_call(S, K, r, sig, T), call_ans)
        S, K = 90, 100 # Put inputs
        put_ans = 9.900431312968252
        self.assertAlmostEqual(ivc.BSC_put(S, K, r, sig, T), put_ans)
        

    def test_BAC(self):
        # Test Bachelier model functions
        r, sig, T = 0.01, 0.1, 0.01
        S, K = 100, 90              # Call inputs
        call_ans = 9.99000499833375
        self.assertAlmostEqual(ivc.BAC_call(S, K, r, sig, T), call_ans)

        S, K = 90, 100              # Put inputs
        call_ans = 9.99000499833375
        self.assertAlmostEqual(ivc.BAC_put(S, K, r, sig, T), put_ans)


    def test_BSC_imvol(self):
        # Test Black Scholes call option IV
        iv_actual = 0.1
        r, T = 0.01, 0.01
        call_ans, S, K = 10.090254370495515, 100, 90 
        self.assertAlmostEqual(ivc.BAC_imvol(call_ans, S, K, r, T, 'Call'), iv_actual)
        
        # Test Black Scholes put option IV
        put_ans, S, K = 9.900431312968252, 90, 100
        self.assertAlmostEqual(ivc.BAC_imvol(put_ans, S, K, r, T, 'Put'), iv_actual)

    def test_BAC_imvol(self):
        # Test Bachelier  call option IV
        iv_actual = 0.1
        r, T = 0.01, 0.01
        call_ans, S, K = 9.99000499833375, 100, 90
        self.assertAlmostEqual(ivc.BAC_imvol(call_ans, S, K, r, T, 'Call'), iv_actual)
        # Test Bachelier put option IV
        put_ans, S, K = 9.99000499833375, 90, 100       
        self.assertAlmostEqual(ivc.BAC_imvol(put_ans, S, K, r, T, 'Put'), iv_actual)

    def test_newton(self):
        """
        Find solution to g(100, 100, 1, x, 1) = 250 using the Newton solver,
        g is the (transcendental) function defined above
        """   
        ans = 3.6095232769758385 # Wolfram Alpha 
        self.assertAlmostEqual(ivc.newton(g, g_dash, 250, 100, 100, 1, 1), ans)



if __name__ == 'main':
    unittest.main()
 



