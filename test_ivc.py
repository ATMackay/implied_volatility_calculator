# This programme calculates the implied volatility for....

# Interview question for Mako

# Requires: Python 3.6, Pandas, numpy

import unittest
import numpy as np
import ivc




# Check example call/put values (TESTING)
print("Black Scholes Call:", BSC_call(100, 90, 0.01, 0.1, 1))
print("Black Scholes Put:", BSC_put(100, 90, 0.01, 0.1, 1))
print("Bachelier Call", BAC_call(100, 90, 0.01, 0.1, 1))
print("Bachelier Put", BAC_put(100, 90, 0.01, 0.1, 1))


# Local test Newton's method function
def g(S, K, r, x, T):
    return S + K + r + x*x*x + np.sin(0.5*x) + T
def g_dash(S, K, r, x, T):
    # Derivative of g w.r.t. x
    return S + K + r + 3*x*x + np.cos(0.5*x) + T





class Testivc(unittest.testcase):

    def tes_BSC(self):
        S, K, r, sig, T = 100, 90, 0.01, 0.1, 1. # Call inputs
        call_ans = 0
        self.assertEqual(BSC_call(S, K, r, sig, T), ans)
        self.assertEqual(BSC_put(S, K, r, sig, T), ans)
        # Test Bachelier model functions

    def test_BAC(self):
        S, K, r, sig, T = 90, 100, 0.01, 0.1, 1.
        ans = 0
        self.assertEqual(BAC_call(S, K, r, sig, T), ans)
        self.assertEqual(BAC_put(S, K, r, sig, T), ans)
        # Test Black Scholes model functions

    def test_BSC_imvol(self):
        # Test Black Scholes call option IV
        self.assertEqual(BAC_imvol(inputs1, 'Call'), ans1)
        # Test Black Scholes put option IV
        self.assertEqual(BAC_imvol(inputs1, 'Put'), ans1)

    def test_BAC_imvol(self):
        # Test Bachelier  call option IV
        self.assertEqual(BAC_imvol(inputs1, 'Call'), ans1)
        # Test Bachelier put option IV
        self.assertEqual(BAC_imvol(inputs1, 'Put'), ans1)

    def test_newton(self):
        """
        Find solution to g(100, 100, 1, x, 1) = 250 using the Newton solver,
        g is the (transcendental) function defined above
        """   
        ans = 3.60952 # Wolfram Alpha 
        self.assertAlmostEqual(newton(g, g_dash, 250, 100, 100, 1, 1), ans)



if __name__ == 'main':
    unittest.main()
 



