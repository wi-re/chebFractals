#pragma once
#include <cheb/constants.h>
#include <cheb/utilities.h>
#include <cheb/IntervalFunction.h>
#include <optional>
#include <variant>

namespace cheb {
	/*
	This is the main class of chebfun. This class represents a chebyshev approximation of
	an arbitrary function over a given domain (default -1 to 1). This class is immutable.
	Instances of this class can be created as:
		- empty using default constructor, default constructed isntances have no valid internal state
		- constant defined by a scalar c and a domain over which this constant applies
		- identity functions model f(x) = x over a provided domain 
		- arbitrary functions (either func or vfunc) with either adaptive polynomial degree or
		  defined using a provided degree as an integer or a list of integers that applies to each interval
		- as ensembles of vectors of IntervalFunctions

	*/
	class Function {
	public:
		std::vector<IntervalFunction> funs;
		std::pair<svec, svec> breakdata;

		// Default constructor creating an invalid Function
		Function();
		// Direct constructor expecting a set of IntervalFunctions
		Function(std::vector<IntervalFunction> _funs);
		// Identity function constructor the string based one is a convenience alternative
		Function(Domain d);
		Function(std::string s, Domain d = Domain{ {-1.,1.} });
		// Constant function constructor
		explicit Function(scalar c, Domain d = Domain{ {-1.,1.} });
		// Arbitrary function constructor with fixed polynomial degree (per function or globally)
		Function(vfunc f, std::variant<lvec, len> n, Domain d = Domain{ {-1.,1.} });
		Function(vfunc f, Domain d = Domain{ {-1.,1.} });
		// Arbitrary function constructor with adaptive polynomial degree
		Function(func fd, std::variant<lvec, len> n, Domain d = Domain{ {-1.,1.} });
		Function(func fd, Domain d = Domain{ {-1.,1.} });

		// Call options
		// Evaluate the argument on the underlying function, the method can be chosen as either
		// - The Clenshaw algorithm https://doi.org/10.1090/S0025-5718-1955-0071856-0
		// - Barycentric Lagrange Interpolation https://doi.org/10.1137/S0036144502417715
		// The operator() variant can be called with either scalars or svecs, eval expects a svec.
		svec operator()(svec x, method how = defaultEvalScheme) const;
		scalar operator()(scalar x, method how = defaultEvalScheme) const;

		// Property functions
		scalar hscale() const; // Returns the horizontal extent of the function
		scalar vscale() const; // Returns the vertical magnitude of the function
		Domain domain() const; // Constructs a domain object from the internal breakpoints
		Domain support()const; // Returns the overall domain extent as a single domain object

		// IF functions
		bool isconst() const; // Checks if there is only 1 coefficient, i.e. constant function
		bool isempty() const; // Checks of there are no coefficients, i.e. invalid function
		std::size_t size() const; // Returns the number of underlying IFs
		svec breakpoints()  const; // Return the breakpoints between individual IFs
		scalar dot(const Function& f) const; // Returns the definite integral of this * f

		// Combines the current domain with another domain (of same extent) in order to
		// determine a function that's defined on the  superset of both domains breakpoints
		Function breakWith(Domain targetdomain) const;
		// Simplifies the coefficients of underlying IFs by chopping off unnecessary terms
		Function simplified() const;
		// Restricts the domain of this funciton to a narrower domain. Simplify controls if the
		// resulting IFs should be simplified, individually, before returning them
		Function restricted(Domain subinterval, bool simplify = true) const;

		// Find all roots of the current polynomial over the relevant interval
		// If polish is true the resulting roots (from the polynomial) will be 
		// "polished" using a simple newton iteration using the derivative of 
		// the polynomial, in order to increase the accuracy of the result
		svec roots(bool polish = true, int32_t splitDegree = 50) const;
		// Returns an IntervalFunction that describes the anti-derivative of this function
		// Defined over the same Domain (polynomials have degree n + 1)
		Function antiDerivative() const;
		// Returns an IntervalFunction that describes the derivatives of this function
		// Defined over the same Domain (polynomials have degree n - 1)R
		Function derivative() const;
		// Returns the definite integral of this function over the relevant domain
		scalar definiteIntegral() const;

		// Arithmetic operations to support function arithmetic
		// Supported operations are:
		// Binary operators (take either one or two Interval functions and zero or one scalar):
		//	[+, -, *, /]
		// Binary funtions:
		//  [pow, minimum, maximum]
		// Unary operators: 
		//  [-]
		// Unary functions:
		//	[absolute, sqrt]
		//  [acos, acosh, asin, asinh, atan, atanh]
		//  [cos, cosh, sin, sinh, tan, tanh]
		//  [exp, exp2, expm1, log, log2, log10, log1p]
#define chebBin(op)\
	friend Function (::operator op)(const Function& lhs, const Function& rhs);\
	friend Function (::operator op)(const cheb::scalar& lhs, const Function& rhs);\
	friend Function (::operator op)(const Function& lhs, const cheb::scalar& rhs);

		chebBin(+);
		chebBin(-);
		chebBin(/ );
		chebBin(*);

		Function operator-();
		friend Function pow(const Function& lhs, const Function& rhs);
		friend Function pow(const cheb::scalar& lhs, const Function& rhs);
		friend Function pow(const Function& lhs, const cheb::scalar& rhs);

		Function absolute();

		Function maximum(const Function& other) const;
		Function maximum(const scalar& c) const;
		Function minimum(const Function& other) const;
		Function minimum(const scalar& c) const;
	};
#define unaryFn(f, ns)\
inline auto f(const cheb::Function& fn) {\
	return cheb::Function([fn](cheb::scalar x) {return ns::f(fn(x)); }, fn.domain());\
}
	unaryFn(acos, std);
	unaryFn(acosh, std);
	unaryFn(asin, std);
	unaryFn(asinh, std);
	unaryFn(atan, std);
	unaryFn(atanh, std);

	unaryFn(cos, std);
	unaryFn(cosh, std);
	unaryFn(sin, std);
	unaryFn(sinh, std);
	unaryFn(tan, std);
	unaryFn(tanh, std);

	unaryFn(exp, std);
	unaryFn(exp2, std);
	unaryFn(expm1, std);
	unaryFn(log, std);
	unaryFn(log2, std);
	unaryFn(log10, std);
	unaryFn(log1p, std);
	unaryFn(sqrt, std);

}

#define chebBinFnFwd(op)\
 inline cheb::Function operator op(const cheb::Function& lhs, const cheb::Function& rhs);\
 inline cheb::Function operator op(const cheb::Function& lhs, const cheb::scalar& rhs);\
 inline cheb::Function operator op(const cheb::scalar& lhs, const cheb::Function& rhs);

#define chebBinFn(op)\
 inline cheb::Function operator op(const cheb::Function& lhs, const cheb::Function& rhs) {\
	using namespace cheb;\
	if (lhs.isempty() || rhs.isempty())\
		return cheb::Function();\
	auto newdom = lhs.domain().united(rhs.domain());\
	auto chbfn1 = lhs.breakWith(newdom);\
	auto chbfn2 = rhs.breakWith(newdom);\
	std::vector<cheb::IntervalFunction> newfuns;\
	for (int32_t i = 0; i < chbfn1.funs.size(); ++i)\
		newfuns.push_back((chbfn1.funs[i] op chbfn2.funs[i]).simplified());\
	return cheb::Function(newfuns);\
}\
 inline cheb::Function operator op(const cheb::Function& lhs, const cheb::scalar& rhs) {\
	using namespace cheb;\
	if (lhs.isempty())\
		return cheb::Function();\
	auto chbfn1 = lhs;\
	std::vector<cheb::IntervalFunction> newfuns;\
	for (int32_t i = 0; i < chbfn1.funs.size(); ++i)\
		newfuns.push_back((chbfn1.funs[i] op rhs));\
	return cheb::Function(newfuns);\
}\
 inline cheb::Function operator op(const cheb::scalar& lhs, const cheb::Function& rhs) {\
	using namespace cheb;\
	if (rhs.isempty())\
		return cheb::Function();\
	auto chbfn1 = rhs;\
	std::vector<cheb::IntervalFunction> newfuns;\
	for (int32_t i = 0; i < chbfn1.funs.size(); ++i)\
		newfuns.push_back((lhs op chbfn1.funs[i]));\
	return cheb::Function(newfuns);\
}


	chebBinFn(+);
	chebBinFn(-);
	chebBinFn(*);
	chebBinFn(/ );

#undef chebBinFn
#undef unaryFn
#undef chebBin
