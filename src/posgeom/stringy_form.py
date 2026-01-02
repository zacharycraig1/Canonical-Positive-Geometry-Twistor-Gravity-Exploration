import numpy as np
from scipy.optimize import root
import itertools

class StringyCanonicalForm:
    def __init__(self, exponents, coefficients):
        self.exponents = np.array(exponents)
        self.coefficients = np.array(coefficients)
        self.n_terms, self.dim = self.exponents.shape
    def evaluate_F(self, x):
        log_x = np.log(np.complex128(x))
        terms = self.coefficients * np.exp(np.dot(self.exponents, log_x))
        return np.sum(terms)
    def get_saddle_equations(self, x, s):
        x = np.complex128(x)
        log_x = np.log(x)
        log_monomials = np.dot(self.exponents, log_x)
        u = self.coefficients * np.exp(log_monomials)
        eqs = np.dot(u, self.exponents) + s
        return eqs
    def compute_hessian_det(self, x, s):
        x = np.complex128(x)
        log_x = np.log(x)
        log_monomials = np.dot(self.exponents, log_x)
        u = self.coefficients * np.exp(log_monomials)
        weighted_V = u[:, None] * self.exponents
        H = np.dot(weighted_V.T, self.exponents)
        return np.linalg.det(H)
    def find_saddle_points(self, s, x0_list=None):
        if x0_list is None:
            x0_list = [np.ones(self.dim) * 0.1, np.ones(self.dim), np.ones(self.dim)*10]
        solutions = []
        seen = []
        for x0 in x0_list:
            def func(y):
                x = np.exp(y)
                return np.real(self.get_saddle_equations(x, s)) 
            res = root(lambda y: func(y), np.log(x0), method='hybr', tol=1e-8)
            if res.success:
                sol_y = res.x
                sol_x = np.exp(sol_y)
                is_new = True
                for seen_x in seen:
                    if np.allclose(sol_x, seen_x, rtol=1e-4):
                        is_new = False
                        break
                if is_new:
                    seen.append(sol_x)
                    solutions.append(sol_x)
        return solutions
    def compute_canonical_form(self, s):
        saddles = self.find_saddle_points(s)
        if not saddles:
            return 0
        real_saddles = [x for x in saddles if np.all(np.isreal(x)) and np.all(x > 0)]
        if real_saddles:
            x_star = real_saddles[0]
            detH = self.compute_hessian_det(x_star, s)
            return 1.0 / detH
        else:
            x_star = saddles[0]
            detH = self.compute_hessian_det(x_star, s)
            return 1.0 / detH

