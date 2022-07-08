#pragma once 
#ifndef _INCLUDE_CHEB_HEADERS_INTERNAL
#error This file is an internal header and should not be included directly. Please include <cheb/cheb.h> instead.
#endif
#include <cheb/cheb.h>

namespace sturm {
	struct SturmSequence {
		using scalar = cheb::scalar;
		std::vector<scalar> _data;
		std::vector<std::size_t> lengths;
		std::vector<std::pair<scalar, std::vector<scalar>>> quotientSeries;
		std::size_t degree;
		std::size_t terms;
		double dynamicRange;
		double thresholdA;
		double thresholdB;

		bool useQuotientSeries = false;

		bool potentialProblem = false;

		int32_t trimseq(scalar* seq, int32_t n, double l);
		void writeZSeries(scalar* in, int32_t n, scalar* out);
		int32_t zseries_div(scalar* z1, int32_t lc1, scalar* z2, int32_t lc2, scalar* z3);
		int32_t chebdiv(scalar* z1, int32_t lc1, scalar* z2, int32_t lc2, scalar* z3);

		std::pair<int32_t, std::pair<scalar, std::vector<scalar>>> zseries_div_withq(scalar* z1, int32_t lc1, scalar* z2, int32_t lc2, scalar* z3);
		std::pair<int32_t, std::pair<scalar, std::vector<scalar>>> chebdiv_withq(scalar* z1, int32_t lc1, scalar* z2, int32_t lc2, scalar* z3);
		void symmetrize(scalar* z, int32_t n);
		void setDegreeWithZseries(std::size_t rank, scalar* zs, int32_t len);
		void generateSequence();
		void generateDerivative();
		void diffTerm(const cheb::svec& polyNomial);
		static cheb::svec differentiate(const cheb::svec& polyNomial);
		void initData(const cheb::svec& polyNomial);
		scalar* getDegree(int32_t rank);
		void setDegreeWithZseries(std::size_t rank, const cheb::svec& zs);
		void setDegree(std::size_t rank, const cheb::svec& data);
		void printDegree(int32_t rank);


		int32_t numRootsRemainderSeries(std::pair<scalar, scalar> lr = std::make_pair(-1.0, 1.0));
		std::pair< int32_t, int32_t> evalSturmRemainderSeries(std::pair<scalar, scalar> lr = std::make_pair(-1.0, 1.0));
		std::array< int32_t, 32> evalSturmRemainderSeries(std::array<scalar, 32> vals);

		int32_t numRootsQuotientSeries(std::pair<scalar, scalar> lr = std::make_pair(-1.0, 1.0));
		std::pair< int32_t, int32_t> evalSturmQuotientSeries(std::pair<scalar, scalar> lr = std::make_pair(-1.0, 1.0));
		std::array< int32_t, 32> evalSturmQuotientSeries(std::array<scalar, 32> vals);


		int32_t evalSturmQuotientSeries(scalar vals);
		int32_t evalSturmRemainderSeries(scalar vals);


	public:
		SturmSequence(const cheb::svec& polyNomial, double tA = 6.0, double tB = 5e3, bool qs = false);
		std::pair< int32_t, int32_t> printSturm(std::pair<scalar, scalar> lr = std::make_pair(-1.0, 1.0));

		int32_t numRoots(std::pair<scalar, scalar> lr = std::make_pair(-1.0, 1.0));
		std::pair< int32_t, int32_t> evalSturm(std::pair<scalar, scalar> lr = std::make_pair(-1.0, 1.0));
		std::array< int32_t, 32> evalSturm(std::array<scalar, 32> vals);
		int32_t evalSturm(scalar vals);

		double firstRoot(double li = -1.0, double ri = 1.0);
		double firstIntervalRoot(double li = -1.0, double ri = 1.0);
		std::pair<double,double> firstAndSecondIntervalRoot(double li = -1.0, double ri = 1.0);
	};
}