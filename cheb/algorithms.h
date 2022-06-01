#pragma once 
#ifndef _INCLUDE_CHEB_HEADERS_INTERNAL
#error This file is an internal header and should not be included directly. Please include <cheb/cheb.h> instead.
#endif
#include <cheb/constants.h>

namespace cheb {
	// Given a set of values this function returns the fft of these values
	// Note that this expects a mirrored input, i.e. one that results in n/2 +1 real fft coefficients
	// i.e. for a given set of actual values fc this mirror input can be achieved as
	// concatenate(2. * slice(fc, 0, 1), concatenate(slice(fc, 1), slice(fc, fc.size() - 1, 0, -1)));
	// This fft outputs normalized coefficients such that ifft(fft(x)) = x
	svec fft(svec vals);
	// Inverse fft, see fft above for a description
	svec ifft(svec vals);
	// Function to find all roots based on the coefficients to a chebyshev polynomial ak with a tolerance
	// level of 1e2 eps (by default). This function internally recursively evaluates the domain on sub-domains
	// of length lower than 50 for computational efficiency. Note that the roots returned here might not be exact,
	// especially with further subdivisions, and it is recommend to polish the results later with newtonroots
	// See also:
	//  - Good [1961]: The colleague matrix, a Chebyshev analogue of the companion matrix
	//	- Boyd [2002]: Computing zeros on a real interval through Chebyshev expansion and polynomial rootfinding
	//	- Trefethen [2013]: Approximation Theory and Approximation Practice
	svec rootsunit(svec ak, scalar htol = 1e2 * eps, int32_t splitDegree = 50);
	// Polishes the initial guesses of roots, i.e. estimated using rootsunit, using a simple Newton-Iteration
	// of the given maximum number of iterations and to within a given tolerance. This process relies on the
	// derivative of the chebyshev form of a function and not the derivative of the original function.
	// This process is done in unit space and not in interval space
	svec newtonroots(cheb::IntervalFunction fun, svec rts, scalar tol = 2.0 * eps, uint32_t maxiter = 10);
	// Barycentric Lagrange evaluation of a polynomial ( https://doi.org/10.1137/S0036144502417715 )
	// - xx is an array of evaluation points
	// - fk are the function values at the nodes of the polynomial
	// - xk are the nodes of the polynomial
	// - vk are the barycentric weights at the nodes
	svec bary(svec xx, svec fk, svec xk, svec vk);
	inline scalar bary(scalar xx, svec fk, svec xk, svec vk) {
		return bary(svec{ xx }, fk, xk, vk)[0];
	}
	// Clenshaw's method to evaluate a polynomial ( https://doi.org/10.1090/S0025-5718-1955-0071856-0 )
	// - xx is an array of evaluation points
	// - ak is the set of chebyshev coefficients at the nodes
	svec clenshaw(svec xx, const cheb::svec&  ak);
	inline scalar clenshaw(scalar xx, const cheb::svec& ak) {
		return clenshaw(svec{ xx }, ak)[0];
	}
	// This function returns the number of coefficients required for a specific function to be accurate
	// within some tolerance threshold. This function is used, for example, to find the required degree of
	// a polynomial to represent a function as a polynomial is only accurate, within tol, if this method
	// returns a length less than the number of coefficients.
	// See also http : Chopping a Chebyshev series https://arxiv.org/pdf/1512.01803v1.pdf
	len standard_chop(svec coeffs, scalar tol = 1e2 * eps);
	// Determine the coefficients of an arbitrary function using an adaptive approach. This adaptive 
	// approch relies on standard_chop to return a length less than the current polynomial degree.
	// If the function is overly complex and would require a degree of more than 2^16 the process is stopped
	// and the polynomial is returned as is (but a warning is emitted)
	svec adaptive(vfunc fun, uint32_t maxpow2 = 16);
	// Multiplies two sets of chebyshev coefficients of first kind, of equal length, together
	svec coeffmult(svec fc, svec gc);
	// returns Barycentric weights for Chebyshev points of 2nd kind of length n
	svec barywts2(len n);
	// returns n Chebyshev points of the second-kind
	svec chebpts2(len n);
	// Map function values at Chebyshev points of 2nd kind to first kind Chebyshev polynomial coefficients
	svec vals2coeffs2(svec vals);
	// Inverse of vals2coeffs, i.e. x == coeffs2vals2(vals2coeffs2(x))
	svec coeffs2vals2(svec coeffs);
}