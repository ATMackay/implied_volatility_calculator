"""
Mako Interview Question.

This performs unit tests for ivc_oo.py

Requirements: Python 3.6 (or later), unittest, numpy, ivc_oo

To run ---> $ python3 -m unittest test_ivc_oo.py

"""

import unittest
import numpy as np
import ivc_oo


class Testivc(unittest.TestCase):

    def test_BSC(self):
        # Test Black Scholes model functions for sample inputs
        r, sig, T = 0.01, 0.1, 0.01 # Constant risk-free, volatility, expiry date
        S, K  = 1.1, 0.9 # Call inputs
        call_ans = 0.20008999550015005
        self.assertAlmostEqual(ivc_oo.VolSolver('Call', 'BlackScholes', call_ans, S, K, r, sig, T).BSC_call(), call_ans)
        S, K = K, S # Put inputs
        put_ans = 0.19989000549981684
        self.assertAlmostEqual(ivc_oo.VolSolver('Put', 'BlackScholes', put_ans, S, K, r, sig, T).BSC_put(), put_ans)

        # Test Black Scholes model zero case
        r, sig, T = 0., 10e-10, 10e-10
        S, K = 1.0, 1.0 
        # Black Scholes model should predict a value of zero for these inputs (both call and put)
        self.assertAlmostEqual(ivc_oo.VolSolver('Call', 'BlackScholes', call_ans, S, K, r, sig, T).BSC_call(), 0)
        self.assertAlmostEqual(ivc_oo.VolSolver('Put', 'BlackScholes', put_ans, S, K, r, sig, T).BSC_put(), 0)
        

    def test_BAC(self):
        # Test Bachelier model functions for sample inputs
        r, sig, T = 0.01, 0.1, 0.01
        S, K = 1.1, 0.9             
        call_ans = 0.19998000099996674
        self.assertAlmostEqual(ivc_oo.VolSolver('Call', 'Bachelier', call_ans, S, K, r, sig, T).BAC_call(), call_ans)

        S, K = K, S # Put inputs
        put_ans = 0.19998000099996674
        self.assertAlmostEqual(ivc_oo.VolSolver('Put', 'Bachelier', put_ans, S, K, r, sig, T).BAC_put(), put_ans)

        # Test Bachelier model zero case
        r, sig, T = 0., 10e-10, 10e-10
        S, K = 1.0, 1.0 
        # Bachelier model should predict a value of zero for these inputs
        self.assertAlmostEqual(ivc_oo.VolSolver('Call', 'Bachelier', call_ans, S, K, r, sig, T).BAC_call(), 0)
        self.assertAlmostEqual(ivc_oo.VolSolver('Put', 'Bachelier', put_ans, S, K, r, sig, T).BAC_put(), 0)


    def test_BSCvol(self):
        # Test Black Scholes call option IV
        iv_actual = 0.5
        r, T = 0.01, 0.01
        call_ans, S, K = 0.10038900828135255, 1.0, 0.9
        # Use Brenner Subrahmanyam (1988) approxiamtion as initial guess
        sig_init_guess = np.sqrt(2*np.pi/T)*call_ans/S 
        self.assertAlmostEqual(ivc_oo.VolSolver('Call', 'BlackScholes', call_ans, S, K, r, sig_init_guess, T).newton_solver(), iv_actual)
        
        # Test Black Scholes put option IV
        put_ans, S, K = 0.10020234665104932, 0.9, 1.0
        sig_init_guess = np.sqrt(2*np.pi/T)*put_ans/S
        self.assertAlmostEqual(ivc_oo.VolSolver('Put', 'BlackScholes', put_ans, S, K, r, sig_init_guess, T).newton_solver(), iv_actual)

    def test_BACvol(self):
        # Test Bachelier  call option IV
        iv_actual = 0.5
        r, T = 0.01, 0.1
        call_ans, S, K = 0.12516446483551702, 1.0, 0.9
        sig_init_guess = np.sqrt(2*np.pi/T)*call_ans/S 
        self.assertAlmostEqual(ivc_oo.VolSolver('Call', 'Bachelier', call_ans, S, K, r, sig_init_guess, T).newton_solver(), iv_actual)

        # Test Bachelier put option IV
        put_ans, S, K = 0.12516446483551702, 0.90, 1.0     
        sig_init_guess = np.sqrt(2*np.pi/T)*put_ans/S  
        self.assertAlmostEqual(ivc_oo.VolSolver ('Put', 'Bachelier', put_ans, S, K, r, sig_init_guess, T).newton_solver(), iv_actual)


if __name__ == 'main':
    unittest.main()
 



