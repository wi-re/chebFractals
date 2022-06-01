#include <test.h>
#include <iostream>
#include <cheb/cheb.h>

#include <vector>
#include <iostream>
#include <test.h>

using namespace cheb;

void testInterval() {
	//setUp
	Interval i1{ -2,3 };
	Interval i2{ -2,3 };
	Interval i3{ -1,1 };
	Interval i4{ -1,2 };
	{
		// test init
		auto [a, b] = Interval();
		assert(a == -1 && b == 1);
	}
	{
		// init disallow
		assertThrows((Interval{ 2,0 }));
		assertThrows((Interval{ 0,0 }));
	}
	{
		// test eq
		assert(Interval() == Interval());
		assert(i1 == i2);
		assert(i2 == i1);
		assertFalse(i3 == i2);
		assertFalse(i2 == i3);
	}
	{
		// test neq
		assertFalse(Interval() != Interval());
		assertFalse(i1 != i2);
		assertFalse(i2 != i1);
		assert((i3 != i2));
		assert((i2 != i3));
	}
	{
		// Test contains
		assertTrue(i2.contains(i1));
		assertTrue(i1.contains(i3));
		assertTrue(i1.contains(i4));
		assertFalse(i3.contains(i1));
		assertFalse(i4.contains(i1));
		assertFalse(!i2.contains(i1));
		assertFalse(!i1.contains(i3));
		assertFalse(!i1.contains(i4));
		assertTrue(!i3.contains(i1));
		assertTrue(!i4.contains(i1));
	}
	{
		//test maps
		auto random = randn(1000);
		auto yy = evaluate(random, [](scalar v) {return -1. + 2. * v; });
		auto interval = Interval{ -2,3 };
		auto vals = interval.mapFromInterval(interval(yy)) - yy;
		assertLessEqual(infnorm(vals), eps);
	}
	{
		// test interior
		auto npts = 1000;
		auto x1 = linspace(-2, 3, npts);
		auto x2 = linspace(-3, -2, npts);
		auto x3 = linspace(3, 4, npts);
		auto x4 = linspace(5, 6, npts);
		auto interval = Interval{ -2,3 };
		assertEqual(evaluate(interval.isinterior(x1), [](auto v) {return v ? 1 : 0; }).sum(), npts - 2);
		assertEqual(evaluate(interval.isinterior(x2), [](auto v) {return v ? 1 : 0; }).sum(), 0);
		assertEqual(evaluate(interval.isinterior(x3), [](auto v) {return v ? 1 : 0; }).sum(), 0);
		assertEqual(evaluate(interval.isinterior(x4), [](auto v) {return v ? 1 : 0; }).sum(), 0);
	}
}
void testDomain() {
	{
		// test init
		Domain{ { -2,1} };
		Domain{ { -2,0,1} };
		Domain{ linspace(-10,10,51) };
	}
	{
		// test disallow
		assertThrows((Domain{ {2} }));
		assertThrows((Domain{ {1,-1} }));
		assertThrows((Domain{ {-1,0,0} }));
	}
	{
		auto dom_a = Domain{ {-2,1} };
		auto dom_b = Domain{ {-2,0,1} };
		auto dom_c = Domain{ {-1,0,1,2} };
		vscalar res_a{ -2,1 };
		vscalar res_b{ -2,0,1 };
		vscalar res_c{ -1,0,1,2 };
		assertTrue(all(evaluate(zip(dom_a, res_a), [](auto p) {auto [l, r] = p; return l == r; })));
		assertTrue(all(evaluate(zip(dom_b, res_b), [](auto p) {auto [l, r] = p; return l == r; })));
		assertTrue(all(evaluate(zip(dom_c, res_c), [](auto p) {auto [l, r] = p; return l == r; })));
	}
	{
		// test intervals
		auto dom_a = Domain{ {-2,1} };
		auto dom_b = Domain{ {-2,0,1} };
		auto dom_c = Domain{ {-1,0,1,2} };
		std::vector<std::pair<scalar, scalar>> res_a{ {-2,1} };
		std::vector<std::pair<scalar, scalar>> res_b{ {-2,0},{0,1} };
		std::vector<std::pair<scalar, scalar>> res_c{ {-1,0},{0,1}, {1,2} };
		assertTrue(all(evaluate(zip(dom_a.intervals(), res_a), [](auto p) {auto [it, dom] = p; return it.a == dom.first && it.b == dom.second; })));
		assertTrue(all(evaluate(zip(dom_b.intervals(), res_b), [](auto p) {auto [it, dom] = p; return it.a == dom.first && it.b == dom.second; })));
		assertTrue(all(evaluate(zip(dom_c.intervals(), res_c), [](auto p) {auto [it, dom] = p; return it.a == dom.first && it.b == dom.second; })));
	}
	{
		// test contains
		auto d1 = Domain{ {-2,0,1,3,4,5} };
		auto d2 = Domain{ {-1,2} };
		auto d3 = Domain{ linspace(-10,10,1000) };
		auto d4 = Domain{ {-1,0,1,2} };
		assertTrue(d1.contains(d2));
		assertTrue(d3.contains(d1));
		assertTrue(d3.contains(d2));
		assertTrue(d3.contains(d2));
		assertTrue(d4.contains(d2));
		assertTrue(d2.contains(d4));
		assertFalse(d2.contains(d1));
		assertFalse(d1.contains(d3));
		assertFalse(d2.contains(d3));
	}
	{
		//test contians close
		auto tol = .8 * HTOL;
		auto d1 = Domain{ {-1,2} };
		auto d2 = Domain{ {-1 - tol,2 + tol} };
		auto d3 = Domain{ {-1 - 2 * tol,2 + 4 * tol} };
		assertTrue(d2.contains(d1));
		assertTrue(d2.contains(d1));
		assertFalse(d1.contains(d3));
	}
	{
		// test eq
		auto d1 = Domain{ {-2,0,1,3,5} };
		auto d2 = Domain{ {-2,0,1,3,5} };
		auto d3 = Domain{ {-1,1} };
		assertEqual(d1, d2);
		assertNotEqual(d1, d3);
	}
	{
		// test eq close
		auto tol = .8 * HTOL;
		auto d4 = Domain{ {-2,0,1,3,5} };
		auto d5 = Domain{ {-2 * (1 + tol),0 - tol,1 + tol,3 * (1 + tol),5 * (1 + tol)} };
		auto d6 = Domain{ {-2 * (1 + 2 * tol),0 - 2 * tol,1 + 2 * tol, 3 * (1 + 2 * tol),5 * (1 - 2 * tol)} };
		assertEqual(d4, d5);
		assertNotEqual(d5, d6);
	}
	{
		// test neq
		auto d1 = Domain{ {-2,0,1,3,5} };
		auto d2 = Domain{ {-2,0,1,3,5} };
		auto d3 = Domain{ {-1,1} };
		assertFalse(d1 != d2);
		assertTrue(d1 != d3);
	}
	{
		//def test_from_Function(self) :
		//	ff = Function(lambda x : np.cos(x), np.linspace(-10, 10, 11))
		//	Domain.from_Function(ff)
	}
	{
		// break points in
		//def test_breakpoints_in(self) :
		//d1 = Domain([-1, 0, 1])
		//d2 = Domain([-2, 0.5, 1, 3])

		//result1 = d1.breakpoints_in(d2)
		//self.assertIsInstance(result1, np.ndarray)
		//self.assertTrue(result1.size, 3)
		//self.assertFalse(result1[0])
		//self.assertFalse(result1[1])
		//self.assertTrue(result1[2])

		//result2 = d2.breakpoints_in(d1)
		//self.assertIsInstance(result2, np.ndarray)
		//self.assertTrue(result2.size, 4)
		//self.assertFalse(result2[0])
		//self.assertFalse(result2[1])
		//self.assertTrue(result2[2])
		//self.assertFalse(result2[3])

		//self.assertTrue(d1.breakpoints_in(d1).all())
		//self.assertTrue(d2.breakpoints_in(d2).all())
		//self.assertFalse(d1.breakpoints_in(Domain([-5, 5])).any())
		//self.assertFalse(d2.breakpoints_in(Domain([-5, 5])).any())
	}
	{
		//def test_breakpoints_in_close(self) :
		//tol = .8 * HTOL
		//d1 = Domain([-1, 0, 1])
		//d2 = Domain([-2, 0 - tol, 1 + tol, 3])
		//result = d1.breakpoints_in(d2)
		//self.assertFalse(result[0])
		//self.assertTrue(result[1])
		//self.assertTrue(result[2])
	}
	{
		// support
		auto dom_a = Domain{ {-2,1} };
		auto dom_b = Domain{ {-2,0,1} };
		auto dom_c = Domain{ linspace(-10,10,51) };
		assertTrue(dom_a.support() == Interval(-2., 1.));
		assertTrue(dom_b.support() == Interval(-2., 1.));
		assertTrue(dom_c.support() == Interval(-10., 10.));
	}
	{
		// size
		auto dom_a = Domain{ {-2,1} };
		auto dom_b = Domain{ {-2,0,1} };
		auto dom_c = Domain{ linspace(-10,10,51) };
		assertEqual(dom_a.size(), 2);
		assertEqual(dom_b.size(), 3);
		assertEqual(dom_c.size(), 51);
	}
	{
		// merge
		auto dom_a = Domain{ {-2,-1,0,1} };
		auto dom_b = Domain{ {-1.5,-.5,.5} };
		assertEqual(dom_a.merged(dom_b), (Domain{ {-2,-1.5,-1,-.5,0,0.5,1} }));
	}
	{
		// restrict
		auto dom_a = Domain{ {-2,-1,0,1} };
		auto dom_b = Domain{ {-1.5,-.5,.5} };
		auto dom_c = Domain{ linspace(-2,1,16) };
		assertEqual(dom_a.restricted(dom_b), (Domain{ {-1.5,-1,-.5,0,0.5} }));
		assertEqual(dom_a.restricted(dom_c), dom_c);
		assertEqual(dom_a.restricted(dom_a), dom_a);
		assertEqual(dom_b.restricted(dom_b), dom_b);
		assertEqual(dom_c.restricted(dom_c), dom_c);
		//# tests to check if catch breakpoints that are different by eps
//# (linspace introduces these effects)
		auto dom_d = Domain{ linspace(-.4,.4,2) };
		assertEqual(dom_c.restricted(dom_d), (Domain{ {-.4,-.2,0,.2,.4} }));
	}
	{
		// test restrict raises
		auto dom_a = Domain{ {-2,-1,0,1} };
		auto dom_b = Domain{ {-1.5,-.5,.5} };
		auto dom_c = Domain{ linspace(-2,1,16) };
		assertThrows(dom_b.restricted(dom_a));
		assertThrows(dom_b.restricted(dom_c));
	}
	{
		// test union
		auto dom_a = Domain{ {-2,0,2} };
		auto dom_b = Domain{ {-2,-1,1,2} };
		assertNotEqual(dom_a.united(dom_b), dom_a);
		assertNotEqual(dom_a.united(dom_b), dom_b);
		assertEqual(dom_a.united(dom_b), (Domain{ {-2,-1,0,1,2} }));
		assertEqual(dom_b.united(dom_a), (Domain{ {-2,-1,0,1,2} }));
	}
	{
		auto tol = .8 * HTOL;
		auto dom_a = Domain{ {-2,0,2} };
		auto dom_b = Domain{ {-2 + tol,-1 + tol,1 + tol,2 + tol} };
		assertEqual(dom_a.united(dom_b), (Domain{ {-2,-1,0,1,2} }));
		assertEqual(dom_b.united(dom_a), (Domain{ {-2,-1,0,1,2} }));
	}
	{
		auto dom_a = Domain{ {-2,0} };
		auto dom_b = Domain{ {2,3} };
		assertThrows(dom_a.united(dom_b));
		assertThrows(dom_b.united(dom_a));
	}
}
void testChebTech() {

	// test Chebyshev points
	assertEqual(chebpts2(0).size(), 0);
	assertEqual(vals2coeffs2(vscalar{}).size(), 0);
	assertEqual(coeffs2vals2(vscalar{}).size(), 0);
	for (int32_t k = 0; k < 10; ++k)
		assertLessEqual(infnorm(vals2coeffs2({ (double)k }) - vscalar{ (double)k }), eps);
	for (int32_t k = 0; k < 10; ++k)
		assertLessEqual(infnorm(coeffs2vals2({ (double)k }) - vscalar{ (double)k }), eps);

	auto scaled_tol = [](auto n) {
		auto 		tol = n < 20 ? 5e1 * eps : pow(log((double)n), 2.5) * eps;
		return tol;
	};
	// ------------------------------------------------------------------------
		// Tests to verify the mutually inverse nature of vals2coeffsand coeffs2vals
	// ------------------------------------------------------------------------
	auto vals2coeff2valsTester = [&](auto n) {
		auto values = randn(n);
		auto coeffs = vals2coeffs2(values);
		auto _values = coeffs2vals2(coeffs);
		assertLessEqual(infnorm(values - _values), scaled_tol(n));
	};
	auto coeffs2vals2valsTester = [&](auto n) {
		auto values = randn(n);
		auto coeffs = coeffs2vals2(values);
		auto _values = vals2coeffs2(coeffs);
		assertLessEqual(infnorm(values - _values), scaled_tol(n));
	};
	//vals2coeff2valsTester(5);
	//vals2coeff2valsTester(pow(2, 18) + 1);

	//coeffs2vals2valsTester(5);
	//coeffs2vals2valsTester(pow(2, 18) + 1);

//# ------------------------------------------------------------------------
//	# Add second - kind Chebyshev points test cases to ChebyshevPoints
//# ------------------------------------------------------------------------
	assertLessEqual(infnorm(chebpts2(1) - vscalar{ 0. }), eps);
	assertLessEqual(infnorm(chebpts2(2) - vscalar{ -1.,1. }), eps);
	assertLessEqual(infnorm(chebpts2(3) - vscalar{ -1.,0.,1. }), eps);
	assertLessEqual(infnorm(chebpts2(4) - vscalar{ -1., -.5, .5, 1. }), eps);
	assertLessEqual(infnorm(chebpts2(5) - vscalar{ -1., -pow(2.,-.5), 0., pow(2.,-.5), 1. }), eps);

	//# check the output is of the correct length, the endpoint values are - 1
	//	# and 1, respectively, and that the sequence is monotonically increasing
	auto chebptsLenTester = [](auto k) {
		auto pts = chebpts2(k);
		assertEqual(pts.size(), k);
		assertEqual(pts[0], -1.);
		assertEqual(*(std::end(pts) - 1), 1.);
		assertTrue(all(diff(pts), [](auto s) {return s > 0; }));
	};
	chebptsLenTester(5);
	chebptsLenTester((len)pow(2, 18) + 1);

	//ClassUsage
	auto ff = IntervalFunction([](vscalar x) { return std::sin(30. * x); }, 100);
	auto xx = -1 + 2. * randn(100);
	//xx = linspace(-1, 1, 10);

	//# tests for emptiness of Chebtech2 objects
	{
		auto f1 = IntervalFunction();
		auto f2 = IntervalFunction({ 1 });
		assertTrue(f1.isempty());
		assertFalse(f2.isempty());
		assertFalse(!(f1.isempty()));
		assertTrue(!(f2.isempty()));
	}
	// 	# tests for constantness of Chebtech2 objects
	{
		auto f1 = IntervalFunction();
		auto f2 = IntervalFunction({ 1 });
		assertFalse(f1.isconst());
		assertTrue(f2.isconst());
		assertTrue(!(f1.isconst()));
		assertFalse(!(f2.isconst()));
	}
	//# check the size() method is working properly
	{
		auto cfs = randn(10);
		assertEqual(IntervalFunction().size(), 0);
		assertEqual(IntervalFunction({ 1. }).size(), 1);
		assertEqual(IntervalFunction(cfs).size(), cfs.size());
	}
	//# test the different permutations of self(xx, ..)
	{
		ff(xx);
		auto b = ff(xx, method::bary);
		auto c = ff(xx, method::clenshaw);
		assertLessEqual(infnorm(b - c), 5e1 * eps);
	}
	{
		// test prolong
		assertEqual(ff.prolonged(0).size(), 0);
		assertEqual(ff.prolonged(1).size(), 1);
		assertEqual(ff.prolonged(20).size(), 20);
		assertEqual(ff.prolonged(ff.size()).size(), ff.size());
		assertEqual(ff.prolonged(200).size(), 200);
	}
	{
		// test vscale empty
		auto gg = IntervalFunction();
		assertEqual(gg.vscale(), 0.);
	}
	{
		auto gg = ff.simplified();
		//# check that simplify is calling standard_chop underneath
		assertEqual(gg.size(), standard_chop(ff.coeffs()));
		assertEqual(infnorm(subVec(ff.coeffs(), gg.size()) - gg.coeffs()), 0);
	}
	{
		// """Unit-tests for construction of Chebtech2 objects"""
		// # test n = 0 case separately
		auto vals = randn(0);
		auto fun = IntervalFunction(vals2coeffs2(vals));
		auto cfs = vals2coeffs2(vals);
		assertTrue(fun.coeffs().size() == 0);
		assertTrue(cfs.size() == 0);
		//# now test the other cases
		for (int32_t n = 1; n < 10; ++n) {
			auto vals = randn(n);
			auto fun = IntervalFunction(vals2coeffs2(vals));
			auto cfs = vals2coeffs2(vals);
			assertEqual(infnorm(fun.coeffs() - cfs), 0.);
		}
	}
	{
		auto x = IntervalFunction(Interval{ -1.,1. });
		auto s = -1. + 2. * randn(10000);
		assertEqual(infnorm(s - x(s)), 0.);
	}
	{
		auto coeffs = randn(10);
		auto f = IntervalFunction(coeffs);
		assertLessEqual(infnorm(f.coeffs() - coeffs), eps);
	}
	{
		auto ff = IntervalFunction(1.);
		assertEqual(ff.size(), 1);
		assertTrue(ff.isconst());
		assertFalse(ff.isempty());
	}
	{
		auto ff = IntervalFunction();
		assertEqual(ff.size(), 0);
		assertFalse(ff.isconst());
		assertTrue(ff.isempty());
	}
}


void testChebTech_functions() {
	//# --------------------------------------
	//	#          vscale estimates
	//# --------------------------------------
	std::vector<std::tuple<func, int32_t, scalar>> functions{
		std::make_tuple([](scalar x) {return sin(4. * pi * x); }, 40, 1.),
		std::make_tuple([](scalar x) {return cos(x); }, 15, 1.),
		std::make_tuple([](scalar x) {return cos(4. * pi * x); }, 39, 1.),
		std::make_tuple([](scalar x) {return exp(cos(4. * pi * x)); }, 181, exp(1.)),
		std::make_tuple([](scalar x) {return cos(3244. * x); }, 3389, 1.),
		std::make_tuple([](scalar x) {return exp(x); }, 15, exp(1.)),
		std::make_tuple([](scalar x) {return 1e10 * exp(x); }, 15, 1e10 * exp(1.)),
		std::make_tuple([](scalar x) {return 0. * x + 1; },  1, 1.)
	};

	for (auto [fun, n, vscale] : functions) {
		auto ff = IntervalFunction(fun, n);
		assertLessEqual(std::abs(ff.vscale() - vscale), .1 * vscale);
	}
	// calculus operations
	{
		auto emptyfun = IntervalFunction();
		assertEqual(emptyfun.defIntegral(), 0);
		assertTrue(emptyfun.antiDerivative().isempty());
		assertTrue(emptyfun.derivative().isempty());
	}
	//# --------------------------------------
	//	#           definite integrals
	//# --------------------------------------
	std::vector<std::tuple<func, int32_t, scalar, scalar>> def_integrals{
		std::make_tuple([](scalar x) {return sin(x); }, 14, .0, eps),
		std::make_tuple([](scalar x) {return sin(4. * pi * x); }, 40, .0, 1e1 * eps),
		std::make_tuple([](scalar x) {return cos(x); }, 15, 1.682941969615793, 2 * eps),
		std::make_tuple([](scalar x) {return cos(4. * pi * x); }, 39, 0.,2 * eps),
		std::make_tuple([](scalar x) {return exp(cos(4. * pi * x)); }, 182, 2.5321317555040166711,10 * eps),
		std::make_tuple([](scalar x) {return cos(3244. * x); }, 3389, 5.879599674161602e-04,5e2 * eps),
		std::make_tuple([](scalar x) {return exp(x); }, 15, exp(1) - exp(-1),4. * eps),
		std::make_tuple([](scalar x) {return 1e10 * exp(x); }, 15, 1e10 * (exp(1) - exp(-1)),4e10 * eps),
		std::make_tuple([](scalar x) {return 0. * x + 1; },  1, 2.,eps)
	};
	int32_t i = 0;
	for (auto [fun, n, integral, tol] : def_integrals) {
		//auto ff2 = Chebtech2(fun, n);
		auto ff = IntervalFunction(fun, n);
		auto sum = ff.defIntegral();
		//auto sum2 = ff2.sum();
		std::cout << sum << " - " << integral << std::endl;
		assertLessEqual(std::abs(sum - integral), tol);
	}
	//# --------------------------------------
	//	#          indefinite integrals
	//# --------------------------------------
	std::vector<std::tuple<func, func, int32_t, scalar>> indef_integrals{
		//std::make_tuple([](scalar x) {return 0. * x + 1; },	[](scalar x) {return x; },									1, 1. * eps),
		std::make_tuple([](scalar x) {return x; },			[](scalar x) {return 1. / 2. * x * x; },					2, 2. * eps),
		std::make_tuple([](scalar x) {return x * x; },		[](scalar x) {return 1. / 3. * x * x * x; },				3, 2. * eps),
		std::make_tuple([](scalar x) {return x * x * x; },		[](scalar x) {return 1. / 4. * x * x * x * x; },			4, 2. * eps),
		std::make_tuple([](scalar x) {return x * x * x * x; },	[](scalar x) {return 1. / 5. * x * x * x * x * x; },		5, 2. * eps),
		std::make_tuple([](scalar x) {return x * x * x * x * x; },	[](scalar x) {return 1. / 6. * x * x * x * x * x * x; },	6, 4. * eps),
		std::make_tuple([](scalar x) {return sin(x); },		[](scalar x) {return -cos(x); },							16, 2. * eps),
		std::make_tuple([](scalar x) {return cos(3. * x); },	[](scalar x) {return 1. / 3. * sin(3. * x); },				23, 2. * eps),
		std::make_tuple([](scalar x) {return exp(x); },		[](scalar x) {return exp(x); },								16, 3. * eps),
		std::make_tuple([](scalar x) {return 1e10 * exp(x); },[](scalar x) {return 1e10 * exp(x); },						16, 1e10 * 3. * eps)
	};
	for (auto [fun, dfn, n, tol] : indef_integrals) {
		auto ff = IntervalFunction(fun, n);
		auto gg = IntervalFunction(dfn, n + 1);

		//auto points = chebpts2(n+1);
		//auto values = evaluate(points, dfn);
		//auto coeffs2 = vals2coeffs2(values);


		auto coeffs = gg.coeffs();
		coeffs[0] = coeffs[0] - dfn(-1.);
		//debugPrint(dfn(-1.));
		auto fcoeffs = ff.antiDerivative().coeffs();
		auto absdiff = infnorm(fcoeffs - coeffs);
		//debugPrintVec(fcoeffs);
		//debugPrintVec(coeffs);
		assertLessEqual(absdiff, tol);
	}

	//# --------------------------------------
	//	#            derivatives
	//# --------------------------------------
	std::vector<std::tuple<func, func, int32_t, scalar>> derivatives{
		std::make_tuple([](scalar x) {return 0. * x + 1; },	[](scalar x) {return 0.; },									1, 1. * eps),
		std::make_tuple([](scalar x) {return x; },			[](scalar x) {return 1; },					2, 2. * eps),
		std::make_tuple([](scalar x) {return x * x; },		[](scalar x) {return 2. * x; },				3, 2. * eps),
		std::make_tuple([](scalar x) {return x * x * x; },		[](scalar x) {return 3. * x * x; },			4, 2. * eps),
		std::make_tuple([](scalar x) {return x * x * x * x; },	[](scalar x) {return 4. * x * x * x; },		5, 3. * eps),
		std::make_tuple([](scalar x) {return x * x * x * x * x; },	[](scalar x) {return 5. * x * x * x * x; }, 6, 4. * eps),
		std::make_tuple([](scalar x) {return sin(x); },		[](scalar x) {return cos(x); },							16, 5e1 * eps),
		std::make_tuple([](scalar x) {return cos(3. * x); },	[](scalar x) {return -3. * sin(3. * x); },				23, 5e2 * eps),
		std::make_tuple([](scalar x) {return exp(x); },		[](scalar x) {return exp(x); },								16, 2e2 * eps),
		std::make_tuple([](scalar x) {return 1e10 * exp(x); },[](scalar x) {return 1e10 * exp(x); },						16, 1e10 * 2e2 * eps)
	};
	for (auto [fun, dfn, n, tol] : derivatives) {
		auto ff = IntervalFunction(fun, n);
		auto gg = IntervalFunction(dfn, std::max(n - 1, 1));
		auto coeffs = gg.coeffs();
		auto fcoeffs = ff.derivative().coeffs();
		auto absdiff = infnorm(fcoeffs - coeffs);
		//debugPrintVec(fcoeffs);
		//debugPrintVec(coeffs);
		//debugPrintVec(fcoeffs - coeffs);
		//debugPrint(absdiff)
		assertLessEqual(absdiff, tol);
	}
	// Constructor Tests
	std::vector<std::tuple<func, const char*, int32_t, bool>> testfunctions{
		std::make_tuple([](scalar x) {return x * x * x + x * x + x + 1.1; }, "poly3(x)",4,true),
		std::make_tuple([](scalar x) {return exp(x); }, "exp(x)",15,false),
		std::make_tuple([](scalar x) {return sin(x); }, "sin(x)",14,true),
		std::make_tuple([](scalar x) {return .2 + .1 * sin(x); }, "(.2+.1sin(x))",14,false),
		std::make_tuple([](scalar x) {return cos(20. * x); }, "cos(20x)",51,true),
		std::make_tuple([](scalar x) {return 1.; }, "constfun",1,false),
		std::make_tuple([](scalar x) {return 0.; }, "zerofun",1,true)
	};

	for (auto [fun, name, len, b] : testfunctions) {
		auto fa = IntervalFunction(fun);
		assertEqual(fa.size(), len);
		auto fl = IntervalFunction(fun, 50);
		assertEqual(fl.size(), 50);
		auto fl2 = IntervalFunction(fun, 500);
		assertEqual(fl2.size(), 500);
	}
}

void testArithmetic() {
	std::vector<std::tuple<func, const char*, int32_t, bool>> testfunctions{
		std::make_tuple([](scalar x) {return x * x * x + x * x + x + 1.1; }, "poly3(x)",4,true),
		std::make_tuple([](scalar x) {return exp(x); }, "exp(x)",15,false),
		std::make_tuple([](scalar x) {return sin(x); }, "sin(x)",14,true),
		std::make_tuple([](scalar x) {return .2 + .1 * sin(x); }, "(.2+.1sin(x))",14,false),
		std::make_tuple([](scalar x) {return cos(20. * x); }, "cos(20x)",51,true),
		std::make_tuple([](scalar x) {return 1.; }, "constfun",1,false),
		std::make_tuple([](scalar x) {return 0.; }, "zerofun",1,true)
	};

	auto xx = -1. + 2. * randn(10);
	auto emptyfun = IntervalFunction();

	for (auto [fun, name, n, b] : testfunctions) {
		//std::cout << "IntervalFunction : " << name << std::endl;
		auto chebtech = IntervalFunction(fun, n);
		assertTrue((emptyfun + chebtech).isempty());
		assertTrue((chebtech + emptyfun).isempty());
		std::vector consts{ -1.,1.,10., -1e5 };
		for (auto c : consts) {
			auto f = [c, fun](vscalar x) {return c + evaluate(x, fun); };
			auto f1 = c + chebtech;
			auto f2 = chebtech + c;
			auto tol = 5e1 * eps * abs(c);
			assertLessEqual(infnorm(f(xx) - f1(xx)), tol);
			assertLessEqual(infnorm(f(xx) - f2(xx)), tol);
		}

		assertTrue((emptyfun - chebtech).isempty());
		assertTrue((chebtech - emptyfun).isempty());
		for (auto c : consts) {
			auto f = [c, fun](vscalar x) {return c - evaluate(x, fun); };
			auto g = [c, fun](vscalar x) {return evaluate(x, fun) - c; };
			auto f1 = c - chebtech;
			auto f2 = chebtech - c;
			//debugPrintVec(g(xx));
			//debugPrintVec(f2(xx));
			auto tol = 5e1 * eps * abs(c);
			assertLessEqual(infnorm(f(xx) - f1(xx)), tol);
			assertLessEqual(infnorm(g(xx) - f2(xx)), tol);
		}
		assertTrue((emptyfun * chebtech).isempty());
		assertTrue((chebtech * emptyfun).isempty());
		for (auto c : consts) {
			auto f = [c, fun](vscalar x) {return c * evaluate(x, fun); };
			auto g = [c, fun](vscalar x) {return evaluate(x, fun) * c; };
			auto f1 = c * chebtech;
			auto f2 = chebtech * c;
			//debugPrintVec(g(xx));
			//debugPrintVec(f2(xx));
			auto tol = 5e1 * eps * abs(c);
			assertLessEqual(infnorm(f(xx) - f1(xx)), tol);
			assertLessEqual(infnorm(g(xx) - f2(xx)), tol);
		}
		assertTrue((emptyfun / chebtech).isempty());
		assertTrue((chebtech / emptyfun).isempty());
		for (auto c : consts) {
			auto g = [c, fun](vscalar x) {return evaluate(x, fun) / c; };
			auto f2 = chebtech / c;
			//debugPrintVec(g(xx));
			//debugPrintVec(f2(xx));
			auto tol = 5e1 * eps * abs(c);
			assertLessEqual(infnorm(g(xx) - f2(xx)), tol);
			if (!b) {
				auto f = [c, fun](vscalar x) {return c / evaluate(x, fun); };
				auto f1 = c / chebtech;
				assertLessEqual(infnorm(f(xx) - f1(xx)), tol);
			}
		}
		assertTrue(pow(emptyfun,(chebtech)).isempty());
		assertTrue(pow(chebtech,(emptyfun)).isempty());


		for (auto [fun2, name2, len2, hasRoots2] : testfunctions) {
			//std::cout << name2 << std::endl;
			auto ff = IntervalFunction(fun, n);
			auto gg = IntervalFunction(fun2, len2);
			auto vscl = std::max(ff.vscale(), gg.vscale());
			auto lscl = std::max(ff.size(), gg.size());
			auto tol = 3. * vscl * lscl * eps;

			{
				auto FG = ff + gg;
				auto fun = [ff, gg](vscalar x) {return ff(x) + gg(x); };
				assertLessEqual(infnorm(fun(xx) - FG(xx)), tol);
			}
			{
				auto FG = ff - gg;
				auto fun = [ff, gg](vscalar x) {return ff(x) - gg(x); };
				assertLessEqual(infnorm(fun(xx) - FG(xx)), tol);
			}
			{
				auto FG = ff * gg;
				auto fun = [ff, gg](vscalar x) {return ff(x) * gg(x); };
				auto res = fun(xx);
				auto resFG = FG(xx);
				auto diff = fun(xx) - FG(xx);
				auto err = infnorm(diff);
				assertLessEqual(infnorm(fun(xx) - FG(xx)), tol);
			}
			if (!hasRoots2) {
				auto FG = ff / gg;
				auto fun = [ff, gg](vscalar x) {return ff(x) / gg(x); };
				assertLessEqual(infnorm(fun(xx) - FG(xx)), tol);
			}
		}
	}
	std::vector consts2{ 1,2 };
	auto fun = IntervalFunction([](scalar x) {return std::sin(x); }, 15);
	for (auto c : consts2) {
		auto f = [c, fun](vscalar x) {return pow(evaluate(x, fun), c); };
		auto g = [c, fun](vscalar x) {return pow(c, evaluate(x, fun)); };
		auto f1 = pow(fun, c);
		auto f2 = pow(c, fun);
		//debugPrintVec(g(xx));
		//debugPrintVec(f2(xx));
		auto tol = 5e1 * eps * abs(c);
		assertLessEqual(infnorm(f(xx) - f1(xx)), tol);
		assertLessEqual(infnorm(g(xx) - f2(xx)), tol);
	}
	fun = IntervalFunction([](scalar x) {return std::exp(x); }, 15);
	for (auto c : consts2) {
		auto f = [c, fun](vscalar x) {return pow(evaluate(x, fun), c); };
		auto g = [c, fun](vscalar x) {return pow(c, evaluate(x, fun)); };
		auto f1 = pow(fun,c);
		auto f2 = pow(c,fun);
		//debugPrintVec(g(xx));
		//debugPrintVec(f2(xx));
		auto tol = 5e1 * eps * abs(c);
		assertLessEqual(infnorm(f(xx) - f1(xx)), tol);
		assertLessEqual(infnorm(g(xx) - f2(xx)), tol);
	}
}
void testRoots() {
	{
		auto ff = IntervalFunction();
		assertEqual(ff.roots().size(), 0);
	}
	{
		auto ff = IntervalFunction(0.);
		auto gg = IntervalFunction(1.);
		assertEqual(ff.roots().size(), 0);
		assertEqual(gg.roots().size(), 0);
	}
	{
		std::vector<std::tuple<func, vscalar, scalar>> rootstestuns{
			std::make_tuple([](scalar x) {return 3. * x + 2; },			vscalar{-2. / 3.}, 1. * eps),
			std::make_tuple([](scalar x) {return x * x; },				vscalar{0,0}, 1. * eps),
			std::make_tuple([](scalar x) {return x * x + .2 * x - .08; },vscalar{-.4,.2}, 1. * eps),
			std::make_tuple([](scalar x) {return sin(x); },				vscalar{0}, 1. * eps),
			std::make_tuple([](scalar x) {return cos(2 * pi * x); },	vscalar{-0.75, -0.25, 0.25, 0.75}, 1. * eps),
			std::make_tuple([](scalar x) {return std::sin(100. * pi * x); },	linspace(-1,1,201), 1. * eps),
			std::make_tuple([](scalar x) {return sin(5. * pi / 2. * x); },	vscalar{-.8, -.4, 0, .4, .8}, 1. * eps)
		};

		for (auto [fun, rootsv, tol] : rootstestuns) {
			auto ff = IntervalFunction(fun);
			auto roots = ff.roots();
			debugPrintVec(roots);
			debugPrintVec(roots - rootsv);
			assertLessEqual(infnorm(roots - rootsv), tol);
		}
	}
}

void testUtilities() {
	{
		//"""Tests for the chebpy.core.utilities check_funs method"""
		auto f = [](scalar x) {return std::exp(x); };
		auto fun0 = IntervalFunction(f, Interval(-1, 0));
		auto fun1 = IntervalFunction(f, Interval(0, 1));
		auto fun2 = IntervalFunction(f, Interval(-.5, .5));
		auto fun3 = IntervalFunction(f, Interval(2, 2.5));
		auto fun4 = IntervalFunction(f, Interval(-3, 02));
		std::vector funs_a{ fun1,fun0,fun2 };
		std::vector funs_b{ fun1, fun2 };
		std::vector funs_c{ fun0, fun3 };
		std::vector funs_d{ fun1, fun4 };
		{
			auto funs = check_funs({});
			assertTrue(funs.size() == 0);
		}
		{
			auto funs = check_funs({ fun0, fun1 });
			assertTrue(funs[0] == fun0);
			assertTrue(funs[1] == fun1);
		}
		{
			auto funs = check_funs({ fun1, fun0 });
			assertTrue(funs[0] == fun0);
			assertTrue(funs[1] == fun1);
		}
		{
			// test overlap
			assertThrows(check_funs(funs_a));
			assertThrows(check_funs(funs_b));
		}
		{
			// test gap
			assertThrows(check_funs(funs_c));
			assertThrows(check_funs(funs_d));
		}
	}
	{
		//# tests for the chebpy.core.utilities compute_breakdata function
		auto f = [](scalar x) {return std::exp(x); };
		auto fun0 = IntervalFunction(f, Interval(-1, 0));
		auto fun1 = IntervalFunction(f, Interval(0, 1));
		{
			auto breaks = compute_breakdata({});
			assertTrue(breaks.first.size() == 0);
		}
		{
			auto [x, y] = compute_breakdata({ fun0 });
			assertLessEqual(infnorm(x - vscalar{ -1,.0 }), eps);
			assertLessEqual(infnorm(y - vscalar{ std::exp(-1.),std::exp(0.) }), 2. * eps);
		}
		{
			auto [x, y] = compute_breakdata({ fun0 , fun1 });
			assertLessEqual(infnorm(x - vscalar{ -1,.0,1. }), eps);
			assertLessEqual(infnorm(y - vscalar{ std::exp(-1.),std::exp(0.),std::exp(1.) }), 2. * eps);
		}
	}
}

void testFunction_Construction() {
	auto f = [](scalar x) {return std::exp(x); };
	auto fun0 = IntervalFunction(f, Interval(-1, 0));
	auto fun1 = IntervalFunction(f, Interval(0, 1));
	auto fun2 = IntervalFunction(f, Interval(-.5, .5));
	auto fun3 = IntervalFunction(f, Interval(2, 2.5));
	auto fun4 = IntervalFunction(f, Interval(-3, -2));
	std::vector funs_a{ fun1,fun0,fun2 };
	std::vector funs_b{ fun1, fun2 };
	std::vector funs_c{ fun0, fun3 };
	std::vector funs_d{ fun1, fun4 };
	{
		// init pass
		Function a(std::vector{ fun0 });
		Function b(std::vector{ fun1 });
		Function c(std::vector{ fun2 });
		Function d(std::vector{ fun0, fun1 });
	}
	{
		// init fail
		assertThrows(Function a(funs_a));
		assertThrows(Function b(funs_b));
		assertThrows(Function c(funs_c));
		assertThrows(Function d(funs_d));
	}
	{
		// init empty
		auto emptyfun = Function();
		assertEqual(emptyfun.size(), 0);
	}
	{
		// init const
		assertTrue(Function(1, Domain{ {-1.,1.} }).isconst());
		assertTrue(Function(-10, Domain{ linspace(-1,1,11) }).isconst());
		assertTrue(Function(3, Domain{ {-2.,0.,1.} }).isconst());
		assertTrue(Function(3.14, Domain{ linspace(-100,-90,11) }).isconst());

		assertFalse(Function(std::vector{ fun0 }).isconst());
		assertFalse(Function(std::vector{ fun1 }).isconst());
		assertFalse(Function(std::vector{ fun2 }).isconst());
		assertFalse(Function(std::vector{ fun0,fun1 }).isconst());
	}
	{
		// test identity
		std::vector _doms{
			linspace(-1,1,2),
			linspace(-1,1,11),
			linspace(-10,17,351),
			linspace(-9.3, 3.2,22),
			linspace(2.5,144.3, 2112) };
		for (auto dom : _doms) {
			auto ff = Function(dom);
			auto s = ff.support();
			auto a = s[0], b = s[1];
			auto xx = linspace(a, b, 1001);
			auto tol = eps * ff.hscale();
			auto fx = ff(xx);
			auto diff = fx - xx;
			//debugPrintVec(fx);
			//debugPrintVec(xx);
			//debugPrintVec(diff);
			assertLessEqual(infnorm(fx - xx), tol);
		}
		auto ff = Function(Domain{ -1.,1. });
		auto s = ff.support();
		auto a = s[0], b = s[1];
		auto xx = linspace(a, b, 1001);
		auto tol = eps * ff.hscale();
		assertLessEqual(infnorm(ff(xx) - xx), tol);
	}
	{
		// test_initfun_adaptive_continuous_domain
		auto ff = Function(f, { -2, -1 });
		assertEqual(ff.funs.size(), 1);
		auto [keys, values] = ff.breakdata;
		assertEqual(keys[0], -2.);
		assertEqual(keys[1], -1.);
		assertLessEqual(std::abs(values[0] - f(-2.)), eps);
		assertLessEqual(std::abs(values[1] - f(-1.)), eps);
	}
	{
		// test_initfun_adaptive_piecewise_domain
		auto ff = Function(f, { -2, 0, 1 });
		assertEqual(ff.funs.size(), 2);
		auto [keys, values] = ff.breakdata;
		assertEqual(keys[0], -2.);
		assertEqual(keys[1], 0.);
		assertEqual(keys[2], 1.);
		assertLessEqual(std::abs(values[0] - f(-2.)), eps);
		assertLessEqual(std::abs(values[1] - f(0.)), eps);
		assertLessEqual(std::abs(values[2] - f(1.)), eps);
	}
	{
		// def test_initfun_adaptive_raises(self):
		assertThrows(Function(f, Domain{ -2. }));
	}
	{
		// def test_initfun_fixedlen_continuous_domain(self):
		auto ff = Function(f, 20, { -2, -1 });
		assertEqual(ff.funs.size(), 1);
		auto [keys, values] = ff.breakdata;
		assertEqual(keys[0], -2.);
		assertEqual(keys[1], -1.);
		assertLessEqual(std::abs(values[0] - f(-2.)), 3. * eps);
		assertLessEqual(std::abs(values[1] - f(-1.)), 3. * eps);
	}
	{
		// def test_initfun_fixedlen_piecewise_domain_0(self) :
		auto ff = Function(f, vlen{ 20,30 }, { -2, 0, 1 });
		assertEqual(ff.funs.size(), 2);
		auto [keys, values] = ff.breakdata;
		assertEqual(keys[0], -2.);
		assertEqual(keys[1], 0.);
		assertEqual(keys[2], 1.);
		assertLessEqual(std::abs(values[0] - f(-2.)), 3. * eps);
		assertLessEqual(std::abs(values[1] - f(0.)), 3. * eps);
		assertLessEqual(std::abs(values[2] - f(1.)), 6. * eps);
	}
	{
		// def test_initfun_adaptive_raises(self):
		assertThrows(Function(f, 10, { -2. }));
		assertThrows(Function(f, vlen{ 30,40 }, { -1.,1. }));
	}
}
void testFunction_Properties() {
	auto f0 = Function();
	auto f1 = Function([](scalar x) {return x * x; }, { -1.,1. });
	auto f2 = Function([](scalar x) {return x * x; }, { -1.,0.,1.,2. });
	{
		// breakpoints
		//debugPrintVec(f1.breakpoints());
		//debugPrintVec(f2.breakpoints());
		assertEqual(f0.breakpoints().size(), 0);
		assertTrue(all(f1.breakpoints() == vscalar{ -1., 1. }));
		assertTrue(all(f2.breakpoints() == vscalar{ -1.,0., 1.,2. }));
	}
	{
		// domain
		Domain d1{ -1., 1. };
		Domain d2{ -1.,0.,1.,2. };
		assertEqual(f0.domain().size(), 0);
		assertEqual(f1.domain(), d1);
		assertEqual(f2.domain(), d2);
	}
	{
		// hscale
		assertEqual(f0.hscale(), 0);
		assertEqual(f1.hscale(), 1);
		assertEqual(f2.hscale(), 2);
	}
	{
		// isempty
		assertTrue(f0.isempty());
		assertFalse(f1.isempty());
		assertFalse(f2.isempty());
	}
	{
		// is const
		assertFalse(f0.isconst());
		assertFalse(f1.isconst());
		assertFalse(f2.isconst());
		auto c1 = Function([](scalar x) {return 3; }, 1, { -2,-1,0,1,2,3 });
		auto c2 = Function([](scalar x) {return -1; }, 1, { -2,3 });
		assertTrue(c1.isconst());
		assertTrue(c2.isconst());
	}
	{
		// support
		assertEqual(f0.support().size(), 0);
		assertTrue((f1.support() == vscalar{ -1.,1. }));
		assertTrue((f2.support() == vscalar{ -1.,2. }));
	}
	{
		// vscale
		assertEqual(f0.vscale(), 0.);
		assertEqual(f1.vscale(), 1.);
		assertEqual(f2.vscale(), 4.);
	}
}
void testFunction_ClassUsage() {
	auto f0 = Function();
	auto f1 = Function([](scalar x) {return x * x; }, { -1.,1. });
	auto f2 = Function([](scalar x) {return x * x; }, { -1.,0.,1.,2. });
	{
		std::vector _doms{
			linspace(-1,1,2),
			linspace(-1,1,11),
			linspace(-9.3,3.2,22)
		};
		for (auto dom : _doms) {
			auto f = Function([](scalar s) {return std::sin(s); }, 1000, dom);
			auto x = Function(f.support());
			auto s = x.support();
			auto a = s[0], b = s[1];
			auto pts = linspace(a, b, 10001);
			auto tol = eps * f.hscale();
			assertLessEqual(infnorm(x(pts) - pts), tol);
		}
	}
	{
		// restrict
		std::vector<vscalar> doms{ {-4.,4,},{-4.,0.,4.},{-2.,-1.,.3,1.,2.5} };
		for (auto dom : doms) {
			auto ff = Function([](scalar s) {return std::cos(s); }, 25, dom);
			auto yy = linspace(dom[0], *(std::end(dom) - 1), 1001);
			std::vector subdoms{ Domain(yy), Domain(yy[std::slice(2,5,1)]), Domain(yy[std::slice(0,500,2)]) };
			for (auto subdom : subdoms) {
				auto xx = linspace(subdom.breakpoints[0], *(std::end(subdom.breakpoints) - 1), 1001);
				auto gg = ff.restricted(subdom, false);
				auto vscl = ff.vscale();
				auto hscl = ff.hscale();
				auto lscl = 0.;
				for (const auto& fun : ff.funs)
					lscl = std::max(lscl, (scalar)fun.size());
				auto tol = vscl * hscl * lscl * eps;

				auto ffx = ff(xx);
				auto ggx = gg(xx);
				assertLessEqual(infnorm(ffx - ggx), tol);
				assertGreaterEqual(gg.funs.size(), subdom.size() - 1);
				for (const auto& fun : gg.funs) {
					assertEqual(fun.size(), 25);
				}
			}
		}
	}
	{
		//restrict empty
		assertTrue(f0.restricted({ -1,1. }, false).isempty());
	}
	{
		// test simplify
		auto dom = linspace(-2, 1.5, 13);
		auto f = Function([](vscalar x) {return std::cos(x); }, 70, dom).simplified();
		auto g = Function([](vscalar x) {return std::cos(x); }, dom);
		assertEqual(f.domain(), g.domain());
		for (int32_t i = 0; i < f.funs.size(); ++i)
			assertLessEqual(f.funs[i].size() - g.funs[i].size(), 1);
	}
	{
		// simplify empty
		assertTrue(f0.simplified().isempty());
	}
	{
		// test restrict
		auto dom1 = Domain(linspace(-2., 1.5, 13));
		auto dom2 = Domain(linspace(-1.7, 0.93, 17));
		auto dom3 = dom1.merged(dom2).restricted(dom2);
		auto f = Function([](vscalar x) {return std::cos(x); }, dom1).restricted(dom2);
		auto g = Function([](vscalar x) {return std::cos(x); }, dom3);
		assertEqual(f.domain(), g.domain());
		for (int32_t i = 0; i < f.funs.size(); ++i)
			assertLessEqual(f.funs[i].size() - g.funs[i].size(), 4);
	}
	{
		//restrict empty
		assertTrue(f0.restricted({ -1,1. }).isempty());
	}

}
void testFunction_Evaluation() {
	auto f0 = Function();
	auto f1 = Function([](vscalar x) {return x * x; }, { -1.,1. });
	auto f2 = Function([](vscalar x) {return x * x; }, { -1.,0.,1.,2. });

	{
		// call empty 
		assertEqual(f0(linspace(-1, 1, 100)).size(), 0);
	}
	{
		// call empty array
		assertEqual(f0(vscalar{}).size(), 0);
		assertEqual(f1(vscalar{}).size(), 0);
		assertEqual(f2(vscalar{}).size(), 0);
	}
	{
		// breakpoints
		auto x1 = f1.breakpoints();
		auto x2 = f2.breakpoints();
		assertTrue(all(f1(x1) == vscalar{ 1,1 }));
		assertTrue(all(f2(x2) == vscalar{ 1,0,1,4 }));
	}
	{
		// outside
		auto x = linspace(-3, 3, 100);
		assertTrue(all(evaluate(f1(x), [](scalar x) {return std::isfinite(x); })));
		assertTrue(all(evaluate(f2(x), [](scalar x) {return std::isfinite(x); })));
	}
	{
		// general eval
		auto f = [](vscalar x) { return std::sin(4. * x) + std::exp(std::cos(14. * x)) - 1.4; };
		auto npts = 50000;
		vscalar dom1{ -1.,1. };
		vscalar dom2{ -1,0,1 };
		vscalar dom3{ -2,0.3,1.2 };
		auto ff1 = Function(f, dom1);
		auto ff2 = Function(f, dom2);
		auto ff3 = Function(f, dom3);
		auto x1 = linspace(dom1[0], dom1[1], npts);
		auto x2 = linspace(dom2[0], dom2[2], npts);
		auto x3 = linspace(dom3[0], dom3[2], npts);
		assertLessEqual(infnorm(f(x1) - ff1(x1)), 5e1 * eps);
		assertLessEqual(infnorm(f(x2) - ff2(x2)), 5e1 * eps);
		assertLessEqual(infnorm(f(x3) - ff3(x3)), 6e1 * eps);
	}
}
#include <iomanip>
void testFunction_Calculus() {
	auto f = [](vscalar x) { return std::sin(4. * x - 1.4);  };
	auto g = [](vscalar x) { return std::exp(x); };
	auto df = [](vscalar x) { return 4. * std::cos(4. * x - 1.4); };
	auto If = [](vscalar x) {return -.25 * std::cos(4. * x - 1.4); };

	auto f1 = Function(f, { -1,1 });
	auto f2 = Function(f, { -3,0,1 });
	auto f3 = Function(f, { -2,.3,1.2 });
	auto f4 = Function(f, linspace(-1, 1, 11));
	std::vector fs{ f1,f2,f3,f4 };

	auto g1 = Function(g, { -1,1 });
	auto g2 = Function(g, { -3,0,1 });
	auto g3 = Function(g, { -2,.3,1.2 });
	auto g4 = Function(g, linspace(-1, 1, 11));
	std::vector gs{ g1,g2,g3,g4 };

	{
		// test sum
		assertLessEqual(std::abs(f1.definiteIntegral() - 0.372895407327895), 4. * eps);
		assertLessEqual(std::abs(f2.definiteIntegral() - 0.382270459230604), 4. * eps);
		assertLessEqual(std::abs(f3.definiteIntegral() + 0.008223712363936), 4. * eps);
		assertLessEqual(std::abs(f4.definiteIntegral() - 0.372895407327895), 4. * eps);
	}
	{
		auto f2 = Function(df, { -3,0,1 });
		// test diff
		auto xx = linspace(-5, 5, 10000);
		for (const auto& fi : fs) {
			auto h = fi.support();
			vscalar x = xx[(xx > h[0]) && (xx < h[1])];
			auto dfy1 = f(x);
			auto fdy1 = fi(x);
			auto dfy = df(x);
			auto fg = fi.derivative();
			auto fdy = fg(x);
			auto diff = fdy - dfy;
			auto diffn = fdy1 - dfy1;
			assertLessEqual(infnorm(fdy1 - dfy1), 2e3 * eps);
			assertLessEqual(infnorm(fdy - dfy), 2e3 * eps);
		}
	}
	{
		// test cumsum
		auto xx = linspace(-5, 5, 10000);
		for (const auto& fi : fs) {
			auto ff = Function(If, fi.domain());
			auto h = fi.support();
			vscalar x = xx[(xx > h[0]) && (xx < h[1])];

			auto dfy1 = f(x);
			auto fdy1 = fi(x);
			auto dfy = If(x);
			auto fg = fi.antiDerivative();
			auto fdy = fg(x);
			auto diff = fdy - dfy;
			auto diffn = fdy1 - dfy1;

			auto fa = If({ h[0] })[0];
			auto f1 = fi.antiDerivative()(x);
			auto f2 = If(x);
			auto ddiff = f1 - f2 + fa;
			//debugPrintVec(ddiff);
			//debugPrint(infnorm(ddiff));

			std::cout << std::setprecision(17) << infnorm(ddiff) << std::endl;
			std::cout << std::setprecision(17) << 4 * eps << std::endl;


			assertLessEqual(infnorm( ddiff), 4 * eps);
		}
	}
	{
		// sum empty
		auto f = Function();
		assertEqual(f.definiteIntegral(), 0.);
	}
	{
		// cumsum empty
		auto f = Function();
		assertTrue(f.antiDerivative().isempty());
	}
	{
		// diff empty
		auto f = Function();
		assertTrue(f.derivative().isempty());
	}
	{
		// DOT MISSING
	}
}
void testFunction_Roots() {
	auto f1 = Function([](vscalar x) {return std::cos(4. * pi * x); }, linspace(-10., 10., 101));
	auto f2 = Function([](vscalar x) {return std::sin(2. * pi * x); }, linspace(-1., 1., 5));
	auto f3 = Function([](vscalar x) {return std::sin(4. * pi * x); }, linspace(-10., 10., 101));
	{
		// test empty
		auto rts = Function().roots();
		assertEqual(rts.size(), 0);
	}
	{
		// multiple pieces
		auto rts = f1.roots();
		assertEqual(rts.size(), 80);
		//assertLessEqual(infnorm(rts - linspace(-9.875, 10., 80)), 10. * eps);
	}
	{
		// repeated roots
		auto rts = f2.roots();
		assertEqual(rts.size(), 5);
		assertLessEqual(infnorm(rts - f2.breakpoints()), eps);
	}
	{
		// multiple pieces + repeated
		auto rts = f3.roots();
		debugPrintVec(rts);
		assertEqual(rts.size(), 81);
		//assertLessEqual(infnorm(rts - linspace(-10., 10., 80)), 10. * eps);
	}
}

void testFunction_Arithmetic() {
	std::vector<std::tuple<vfunc, const char*, int32_t, bool>> testfunctions{
		std::make_tuple([](vscalar x) {return x * x * x + x * x + x + 1.1; }, "poly3(x)",4,true),
		std::make_tuple([](vscalar x) {return exp(x); }, "exp(x)",15,false),
		std::make_tuple([](vscalar x) {return sin(x); }, "sin(x)",14,true),
		std::make_tuple([](vscalar x) {return .2 + .1 * sin(x); }, "(.2+.1sin(x))",14,false),
		std::make_tuple([](vscalar x) {return cos(20. * x); }, "cos(20x)",51,true),
		std::make_tuple([](vscalar x) {return vscalar(1.,x.size()); }, "constfun",1,false),
		std::make_tuple([](vscalar x) {return vscalar(0.,x.size()); }, "zerofun",1,true)
	};
	std::vector<std::pair<Domain, scalar>> testdomains{
		{ {-1., 1.}, 2. * eps},
		{ {-2., 1.}, 1. * eps},
		{ {-1., 2.}, 1. * eps},
		{ {-5., 9.}, 35. * eps}
	};
	std::vector < std::tuple<vfunc, len, const char*>> powtestfunctions{
		std::make_tuple([](vscalar x) {return std::exp(x); }, 15, "exp"),
		std::make_tuple([](vscalar x) {return std::sin(x); }, 15, "sin"),
		std::make_tuple([](vscalar x) {return 2. - x; }, 2, "linear")
	};
	std::vector<std::pair<Domain, scalar>> powtestdomains{
		{ {-.5, .9}, 1. * eps},
		{ {-1.2, 1.3}, 1. * eps},
		{ {-2.2, 1.9}, 1. * eps},
		{ {.4, 1.3}, 1. * eps}
	};

	std::vector cs{ -1.,1.,10.,-1e5 };
	auto emptyfun = Function();
	auto yy = linspace(-1., 1., 2000);
	bool testConstants = false;
	if (testConstants)
	{
		// test constants
		for (auto [f, name, n, hasRoots] : testfunctions) {
			std::cout << name << std::endl;
			for (auto [domain, tol] : testdomains) {
				auto a = domain[0];
				auto b = domain[1];
				auto ff = Function(f, linspace(a, b, 13));
				std::cout << "Testing operator + " << std::endl;
				assertTrue((emptyfun + ff).isempty());
				assertTrue((ff + emptyfun).isempty());
				for (auto c : cs) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 11));
					auto g = [c, f](vscalar x) {return c + f(x); };
					auto gg1 = c + ff;
					auto gg2 = ff + c;
					auto vscl = ff.vscale();
					auto hscl = ff.hscale();
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto tol = 2. * abs(c) * vscl * hscl * lscl * eps;

					auto gy = g(xx);
					auto gy1 = gg1(xx);
					auto gy2 = gg2(xx);
					debugPrint(infnorm(gy - gy1));
					debugPrint(infnorm(gy - gy2));
					debugPrint(tol);
					assertLessEqual(infnorm(gy - gy1), tol);
					assertLessEqual(infnorm(gy - gy2), tol);
				}
				std::cout << "Testing operator - " << std::endl;
				assertTrue((emptyfun - ff).isempty());
				assertTrue((ff - emptyfun).isempty());
				for (auto c : cs) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 11));
					auto g1 = [c, f](vscalar x) {return c - f(x); };
					auto g2 = [c, f](vscalar x) {return f(x) - c; };
					auto gg1 = c - ff;
					auto gg2 = ff - c;
					auto vscl = ff.vscale();
					auto hscl = ff.hscale();
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto tol = 2. * abs(c) * vscl * hscl * lscl * eps;

					auto g1y = g1(xx);
					auto g2y = g2(xx);
					auto gy1 = gg1(xx);
					auto gy2 = gg2(xx);
					debugPrint(infnorm(g1y - gy1));
					debugPrint(infnorm(g2y - gy2));
					debugPrint(tol);

					assertLessEqual(infnorm(g1y - gy1), tol);
					assertLessEqual(infnorm(g2y - gy2), tol);
				}
				std::cout << "Testing operator - (unary)" << std::endl;
				assertTrue((-emptyfun).isempty());
				for (auto c : cs) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 11));
					auto g1 = [c, f](vscalar x) {return -f(x); };
					auto gg1 = -ff;
					auto vscl = ff.vscale();
					auto hscl = ff.hscale();
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto tol = 2. * abs(c) * vscl * hscl * lscl * eps;

					auto g1y = g1(xx);
					auto gy1 = gg1(xx);
					debugPrint(infnorm(g1y - gy1));
					debugPrint(tol);

					assertLessEqual(infnorm(g1y - gy1), tol);
				}
				std::cout << "Testing operator * " << std::endl;
				assertTrue((emptyfun * ff).isempty());
				assertTrue((ff * emptyfun).isempty());
				for (auto c : cs) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 11));
					auto g1 = [c, f](vscalar x) {return c * f(x); };
					auto g2 = [c, f](vscalar x) {return f(x) * c; };
					auto gg1 = c * ff;
					auto gg2 = ff * c;
					auto vscl = ff.vscale();
					auto hscl = ff.hscale();
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto tol = 2. * abs(c) * vscl * hscl * lscl * eps;

					auto g1y = g1(xx);
					auto g2y = g2(xx);
					auto gy1 = gg1(xx);
					auto gy2 = gg2(xx);
					debugPrint(infnorm(g1y - gy1));
					debugPrint(infnorm(g2y - gy2));
					debugPrint(tol);

					assertLessEqual(infnorm(g1y - gy1), tol);
					assertLessEqual(infnorm(g2y - gy2), tol);
				}
				std::cout << "Testing operator / " << std::endl;
				assertTrue((emptyfun / ff).isempty());
				assertTrue((ff / emptyfun).isempty());
				for (auto c : cs) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 11));
					auto g1 = [c, f](vscalar x) {return c / f(x); };
					auto g2 = [c, f](vscalar x) {return f(x) / c; };
					auto gg2 = ff / c;
					auto vscl = ff.vscale();
					auto hscl = ff.hscale();
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto tol = 2. * abs(c) * vscl * hscl * lscl * eps;

					auto g2y = g2(xx);
					auto gy2 = gg2(xx);
					debugPrint(tol);
					debugPrint(infnorm(g2y - gy2));
					assertLessEqual(infnorm(g2y - gy2), tol);
					if (!hasRoots) {
						auto g1y = g1(xx);
						auto gg1 = c / ff;
						auto gy1 = gg1(xx);
						debugPrint(infnorm(g1y - gy1));
						assertLessEqual(infnorm(g1y - gy1), 4. * tol);
					}
				}
			}
		}
		for (auto [f, n, name] : powtestfunctions) {
			std::cout << name << std::endl;
			std::vector cs{ 1.,2.,3. };
			for (auto [domain, tol] : powtestdomains) {
				auto a = domain[0];
				auto b = domain[1];
				auto ff = Function(f, linspace(a, b, 13));
				std::cout << "Testing operator ** " << std::endl;
				assertTrue(pow(emptyfun, ff).isempty());
				assertTrue(pow(ff, emptyfun).isempty());
				for (auto c : cs) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 11));
					auto g1 = [c, f](vscalar x) {return std::pow(c, f(x)); };
					auto g2 = [c, f](vscalar x) {return std::pow(f(x), c); };
					auto gg1 = pow(c, ff);
					auto gg2 = pow(ff, c);
					auto vscl = gg2.vscale();
					auto hscl = gg2.hscale();
					auto lscl = 0.;
					for (auto& fun : gg2.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto tol = 2. * abs(c) * vscl * hscl * lscl * eps;

					auto g1y = g1(xx);
					auto g2y = g2(xx);
					auto gy1 = gg1(xx);
					auto gy2 = gg2(xx);

					assertLessEqual(infnorm(g2y - gy2), tol);
					assertLessEqual(infnorm(g1y - gy1), tol);
				}
			}
		}
	}
	{
		for (auto [f, name, n, hasRoots] : testfunctions) {
			std::cout << name << std::endl;
			for (auto [domain, tol] : testdomains) {
				auto a = domain[0];
				auto b = domain[1];
				auto ff = Function(f, linspace(a, b, 13));
				std::cout << "Testing operator + " << std::endl;
				assertTrue((emptyfun + ff).isempty());
				assertTrue((ff + emptyfun).isempty());
				for (auto [g, name2, n2, hasRoots2] : testfunctions) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 4));
					auto gg = Function(g, linspace(a, b, 8));
					auto fg = [f, g](vscalar x) {return f(x) + g(x); };
					auto FG = ff + gg;

					auto vscl = std::max(ff.vscale(), gg.vscale());
					auto hscl = std::max(ff.hscale(), gg.hscale());
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					for (auto& fun : gg.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto htol = 2. * vscl * hscl * lscl * tol;

					assertEqual(ff.funs.size(), 3);
					assertEqual(gg.funs.size(), 7);
					assertEqual(FG.funs.size(), 9);

					auto fgy = fg(xx);
					auto FGy = FG(xx);
					std::cout << std::setprecision(17) << infnorm(fgy - FGy) << std::endl;
					std::cout << std::setprecision(17) << htol << std::endl;
					assertLessEqual(infnorm(fgy - FGy), htol);
				}
				std::cout << "Testing operator - " << std::endl;
				assertTrue((emptyfun - ff).isempty());
				assertTrue((ff - emptyfun).isempty());
				for (auto [g, name2, n2, hasRoots2] : testfunctions) {
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 4));
					auto gg = Function(g, linspace(a, b, 8));
					auto fg = [f, g](vscalar x) {return f(x) - g(x); };
					auto FG = ff - gg;

					auto vscl = std::max(ff.vscale(), gg.vscale());
					auto hscl = std::max(ff.hscale(), gg.hscale());
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					for (auto& fun : gg.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto htol = 2. * vscl * hscl * lscl * tol;

					assertEqual(ff.funs.size(), 3);
					assertEqual(gg.funs.size(), 7);
					assertEqual(FG.funs.size(), 9);

					auto fgy = fg(xx);
					auto FGy = FG(xx);
					assertLessEqual(infnorm(fgy - FGy), htol);
				}
				std::cout << "Testing operator * " << std::endl;
				assertTrue((emptyfun * ff).isempty());
				assertTrue((ff * emptyfun).isempty());
				for (auto [g, name2, n2, hasRoots2] : testfunctions) {
					std::cout << name << " x " << name2 << " @ " << a << " : " << b << std::endl;
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 4));
					auto gg = Function(g, linspace(a, b, 8));
					auto fg = [f, g](vscalar x) {return f(x) * g(x); };
					auto FG = ff * gg;

					auto vscl = std::max(ff.vscale(), gg.vscale());
					auto hscl = std::max(ff.hscale(), gg.hscale());
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					for (auto& fun : gg.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto htol = 12. * vscl * hscl * lscl * tol;

					assertEqual(ff.funs.size(), 3);
					assertEqual(gg.funs.size(), 7);
					assertEqual(FG.funs.size(), 9);

					auto fgy = fg(xx);
					auto FGy = FG(xx);
					std::cout << std::setprecision(17) << infnorm(fgy - FGy) << std::endl;
					std::cout << std::setprecision(17) << htol << std::endl;
					assertLessEqual(infnorm(fgy - FGy), htol);
				}
				std::cout << "Testing operator / " << std::endl;
				assertTrue((emptyfun / ff).isempty());
				assertTrue((ff / emptyfun).isempty());
				for (auto [g, name2, n2, hasRoots2] : testfunctions) {
					if (hasRoots2) continue;
					std::cout << name << " / " << name2 << " @ " << a << " : " << b << std::endl;
					auto xx = linspace(a, b, 1001);
					auto ff = Function(f, linspace(a, b, 4));
					auto gg = Function(g, linspace(a, b, 8));
					auto fg = [f, g](vscalar x) {return f(x) / g(x); };
					auto FG = ff / gg;

					auto vscl = std::max(ff.vscale(), gg.vscale());
					auto hscl = std::max(ff.hscale(), gg.hscale());
					auto lscl = 0.;
					for (auto& fun : ff.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					for (auto& fun : gg.funs)
						lscl = std::max(lscl, (scalar)fun.size());
					auto htol = 9. * vscl * hscl * lscl * tol;

					assertEqual(ff.funs.size(), 3);
					assertEqual(gg.funs.size(), 7);
					assertEqual(FG.funs.size(), 9);

					auto fgy = fg(xx);
					auto FGy = FG(xx);
					std::cout << std::setprecision(17) << infnorm(fgy - FGy) << std::endl;
					std::cout << std::setprecision(17) << htol << std::endl;
					assertLessEqual(infnorm(fgy - FGy), htol);
				}
			}
		}

	}


}

void testFunction() {
	{
		auto x = Function(std::string("x"), linspace(-1.25, 3, 11));
		auto y = Function(2, linspace(-1.25, 3, 11));
		auto g = pow(x, y).minimum(1.5);
		auto t = linspace(-1.25, 3, 21);
		auto f = [](vscalar x) { return evaluate(std::pow(x, 2), [](scalar x) {return std::min(x, 1.5); }); };
		auto i = [](vscalar x) { return std::pow(x, 2); };
		debugPrintVec(f(t));
		debugPrintVec(g(t));
		debugPrint(infnorm(f((t)-g(t))));
		assertLessEqual(infnorm(f(t) - g(t)), 1e1 * eps);
	}
	{
		auto x = Function(std::string("x"), linspace(-1.25, 3, 11));
		auto y = Function(2, linspace(-1.25, 3, 11));
		auto g = pow(x, y).maximum(1.5);
		auto t = linspace(-1.25, 3, 21);
		auto f = [](vscalar x) { return evaluate(std::pow(x, 2), [](scalar x) {return std::max(x, 1.5); }); };
		auto i = [](vscalar x) { return std::pow(x, 2); };
		//debugPrintVec(f(t));
		//debugPrintVec(g(t));
		//debugPrint(infnorm(f(t)-g(t)));
		assertLessEqual(infnorm(f(t) - g(t)), 1e1 * eps);
	}
}
