#pragma once 
#ifndef _INCLUDE_CHEB_HEADERS_INTERNAL
#error This file is an internal header and should not be included directly. Please include <cheb/cheb.h> instead.
#endif
#include <cheb/utilities.h>

namespace cheb {
	/* 
	This is the main internal class that handles all of the actual CPR functions. This class defaults to 
	an interval of [-1. to 1.] for function intervals and is an immutable class. Instances can be 
	created as:
		- empty using default constructor, default constructed IFs have no valid internal state
		- constant defined by a scalar c and an interval over which this constant applies
		- identity functions model f(x) = x over a provided interval 
		- arbitrary functions (either func or vfunc) with either adaptive polynomial degree or
		  defined using a provided degree as an integer
		- as copies of other functions directly using their polynomial coefficients

	*/
	class IntervalFunction {
	public:
		Interval minterval = { -1.,1. };
		svec _coeffs = {};
		// Default constructor creating an invalid Intervalfunction
		IntervalFunction() :_coeffs(svec{}), minterval(-1, 1) {};
		// Direct constructor expecting a set of polynomial coefficients 
		explicit IntervalFunction(svec coeffs, Interval _interval = Interval{ -1.,1. }) :_coeffs(coeffs), minterval(_interval) {}
		// Identity function constructor
		IntervalFunction(Interval interval) : IntervalFunction((vals2coeffs2({ interval.a, interval.b })), interval) {}
		// Constant function constructor
		IntervalFunction(scalar c, Interval interval = Interval{ -1., 1. }) : IntervalFunction(svec{ c }, interval) {}
		// Arbitrary function constructor with fixed polynomial degree
		IntervalFunction(vfunc f, len n, Interval interval = Interval{ -1.,1. });
		IntervalFunction(func f, len n, Interval interval = Interval{ -1.,1. }) : IntervalFunction(funcToVfunc(f), n, interval) {}
		// Arbitrary function constructor with adaptive polynomial degree
		IntervalFunction(vfunc f, Interval interval = Interval{ -1.,1. });
		IntervalFunction(func f, Interval interval = Interval{ -1.,1. }) : IntervalFunction(funcToVfunc(f), interval) {}
		
		// Helper utilities
		bool operator==(const IntervalFunction& rhs);

		// Call options
		// This class can be called using either operator() (which performs a mapping of the provided
		// x to the internal interval range) or using eval (which expects a parameter in the range
		// -1 to 1. The evaluation can be done using either:
		// - The Clenshaw algorithm https://doi.org/10.1090/S0025-5718-1955-0071856-0
		// - Barycentric Lagrange Interpolation https://doi.org/10.1137/S0036144502417715
		// The operator() variant can be called with either scalars or svecs, eval expects a svec.
		scalar operator()(scalar x, method how = defaultEvalScheme) const;
		svec operator()(svec x, method how = defaultEvalScheme) const;
		svec eval(svec x, method how = defaultEvalScheme) const;
		svec call_clenshaw(svec x) const;
		svec call_bary(svec x) const;

		// Property functions
		Interval interval() const; // Returns the interval over which this function is defined
		svec endvalues() const; // Returns the values of the function at the interval ends
		svec support() const; // Returns the interval as a svec
		scalar vscale() const; // Returns the vertical magnitude of the function

		// coefficient functions
		svec coeffs() const; // Returns the underlying polynomial coefficients
		svec values() const; // Returns the function values at the polynomial nodes
		len size() const; // Returns the number of coefficients
		bool isconst() const; // Checks if there is only 1 coefficient, i.e. constant function
		bool isempty() const; // Checks of there are no coefficients, i.e. invalid function

		// Restricts the current Interval function to a narrower sub interval
		IntervalFunction restricted(Interval subinterval) const;
		// Simplifies the coefficients by chopping off unnecessary terms
		IntervalFunction simplified() const;
		// Increases/Decreases the polynomial degree either by chopping or appending zeros
		IntervalFunction prolonged(int32_t n) const;

		// Find all roots of the current polynomial over the relevant interval
		// If polish is true the resulting roots (from the polynomial) will be 
		// "polished" using a simple newton iteration using the derivative of 
		// the polynomial, in order to increase the accuracy of the result
		svec roots(bool polish = true, int32_t splitDegree = 50)const;
		// Returns an IntervalFunction that describes the anti-derivative of this function
		// Defined over the same interval (polynomial has degree n + 1)
		IntervalFunction antiDerivative() const;
		// Returns an IntervalFunction that describes the derivative of this function
		// Defined over the same interval (polynomial has degree n - 1)
		IntervalFunction derivative()const;
		// Returns the definite integral of this function over the relevant interval
		scalar defIntegral()const;

		// Arithmetic operations to support function arithmetic
		// Supported operations are:
		// Binary operators (take either one or two Interval functions and zero or one scalar):
		//	[+, -, *, /]
		// Binary funtions:
		//  [pow]
		// Unary operators: 
		//  [-]
		// Unary functions:
		//	[abs, sqrt]
		//  [acos, acosh, asin, asinh, atan, atanh]
		//  [cos, cosh, sin, sinh, tan, tanh]
		//  [exp, exp2, expm1, log, log2, log10, log1p]
#define binaryFriend(op)\
friend IntervalFunction (::operator op)(const IntervalFunction& lhs, const IntervalFunction& rhs);\
friend IntervalFunction (::operator op)(const IntervalFunction& lhs, const cheb::scalar& rhs);\
friend IntervalFunction (::operator op)(const cheb::scalar& lhs, const IntervalFunction& rhs);

		binaryFriend(+);
		binaryFriend(-);
		binaryFriend(*);
		binaryFriend(/);

		friend IntervalFunction pow(const IntervalFunction& lhs, const IntervalFunction& rhs);
		friend IntervalFunction pow(const IntervalFunction& lhs, const cheb::scalar& rhs);
		friend IntervalFunction pow(const cheb::scalar& lhs, const IntervalFunction& rhs);

		inline IntervalFunction operator-() const {
			return IntervalFunction(evaluate(_coeffs, [](scalar val) {return -val; }), interval());
		}

#define addUfunc(op)\
	 inline IntervalFunction op() { \
	return IntervalFunction([*this](scalar x){return ::op(this->operator()(x));}, interval()); \
}
		addUfunc(abs);
		addUfunc(acos);
		addUfunc(acosh);
		addUfunc(asin);
		addUfunc(asinh);
		addUfunc(atan);
		addUfunc(atanh);
		addUfunc(cos);
		addUfunc(cosh);
		addUfunc(exp);
		addUfunc(exp2);
		addUfunc(expm1);
		addUfunc(log);
		addUfunc(log2);
		addUfunc(log10);
		addUfunc(log1p);
		addUfunc(sinh);
		addUfunc(sin);
		addUfunc(tan);
		addUfunc(tanh);
		addUfunc(sqrt);
	};
}