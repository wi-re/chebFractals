#include <cheb/cheb.h>
#include <cheb/algorithms.h>
//#include <armadillo>
#include <numbers>
#define FFTW_DLL 
#include <fftw3.h>
#include <array>
#include <mutex>
#include <cstring>
#include <cfloat>

namespace sturm {

	auto clen(int32_t n) {
		return n - (n + 1) / 2 + 1;
	}
	auto clenshaw(std::pair<double, double> xx, const double* ak, int32_t n) {
		if (n == 0)
			return std::pair<double, double>(0.0, 0.0);
		if (n == 1)
			return std::pair<double, double>(ak[0], ak[0]);
		//if (any(ak, [](auto s) {return std::isnan(s); }))
		//    return std::nan("") * vscalar(1., xx.size());
		double bk1l = 0.0, bk1r = 0.0, bk2l = 0.0, bk2r = 0.0;
		auto l = xx.first * 2.0;
		auto r = xx.second * 2.0;
		for (int32_t k = (int32_t)n - 1; k > 1; k -= 2) {
			bk2l = std::fma(l, bk1l, ak[k]) - bk2l;
			bk1l = std::fma(l, bk2l, ak[k - 1]) - bk1l;
			bk2r = std::fma(r, bk1r, ak[k]) - bk2r;
			bk1r = std::fma(r, bk2r, ak[k - 1]) - bk1r;
			//bk2l = (l * bk1l + ak[k]) - bk2l;
			//bk1l = (l * bk2l + ak[k - 1]) - bk1l;
			//bk2r = (r * bk1r + ak[k]) - bk2r;
			//bk1r = (r * bk2r + ak[k - 1]) - bk1r;
		}
		if ((n - 1) % 2 == 1) {
			auto temp = bk1l;
			bk1l = std::fma(l, bk1l, ak[1]) - bk2l;
			//bk1l = (l * bk1l + ak[1]) - bk2l;
			bk2l = temp;
			temp = bk1r;
			bk1r = std::fma(r, bk1r, ak[1]) - bk2r;
			//bk1r = (r * bk1r + ak[1]) - bk2r;
			bk2r = temp;
		}
		l = std::fma(0.5 * l, bk1l, ak[0]) - bk2l;
		r = std::fma(0.5 * r, bk1r, ak[0]) - bk2r;
		//l = (0.5 * l * bk1l + ak[0]) - bk2l;
		//r = (0.5 * r * bk1r + ak[0]) - bk2r;

		return std::make_pair(l, r);
	}
	auto clenshaw(double x, const double* ak, int32_t n) {
		if (n == 0)
			return 0.0;
		if (n == 1)
			return ak[0];
		//if (any(ak, [](auto s) {return std::isnan(s); }))
		//    return std::nan("") * vscalar(1., xx.size());
		double bk1l = 0.0, bk1r = 0.0, bk2l = 0.0, bk2r = 0.0;
		auto l = x * 2.0;
		for (int32_t k = (int32_t)n - 1; k > 1; k -= 2) {
			bk2l = ak[k] + l * bk1l - bk2l;
			bk1l = ak[k - 1] + l * bk2l - bk1l;
		}
		if ((n - 1) % 2 == 1) {
			auto temp = bk1l;
			bk1l = ak[1] + l * bk1l - bk2l;
			bk2l = temp;
		}
		l = ak[0] + 0.5 * l * bk1l - bk2l;

		return l;
	}
	auto clenshaw(std::array<double, 32> xx, const double* ak, int32_t n) {
		if (n == 0) {
			std::array<double, 32> res;
			res.fill(0.0);
			return res;
		}
		if (n == 1) {
			std::array<double, 32> res;
			res.fill(ak[0]);
			return res;
		}
		//if (any(ak, [](auto s) {return std::isnan(s); }))
		//    return std::nan("") * vscalar(1., xx.size());
		std::array<double, 32> bk1, bk2, res;
		bk1.fill(0.0);
		bk2.fill(0.0);
		for (int32_t k = (int32_t)n - 1; k > 1; k -= 2) {
			for (int32_t i = 0; i < 32; ++i) {
				bk2[i] = ak[k] + 2.0 * xx[i] * bk1[i] - bk2[i];
				bk1[i] = ak[k - 1] + 2.0 * xx[i] * bk2[i] - bk1[i];
			}
		}
		if ((n - 1) % 2 == 1) {
			for (int32_t i = 0; i < 32; ++i) {
				auto temp = bk1[i];
				bk1[i] = ak[1] + 2.0 * xx[i] * bk1[i] - bk2[i];
				bk2[i] = temp;
			}
		}
		for (int32_t i = 0; i < 32; ++i) {
			res[i] = ak[0] + xx[i] * bk1[i] - bk2[i];
		}
		return res;
	}
	auto sign(double x) {
		if (x >= 0.0) return 1.0;
		if (x < 0.0) return -1.0;
		return 0.0;
	}
	auto sign(std::pair<double, double> x) {
		return std::make_pair(sign(x.first), sign(x.second));
	}
	auto sign(std::array<double, 32> x) {
		std::array<double, 32> v;
		for (int32_t i = 0; i < 32; ++i)
			v[i] = sign(x[i]);
		return v;
	}
	// 0x1.408c1ep+2 0x1.a18707p+0 0x1.a0a766p+2 0x1.0ae690p-2 0x1.ba5d89p+0 0x1.a3ec29p-10 0x1.5efc37p-9 -0x1.58f94dp-9 -0x1.a7742ap-8 0x1.0f33fep-9 0x1.0490dap-8 -0x1.42d539p-10 -0x1.9c6e11p-10 0x1.30463ep-11 0x1.7fd4a6p-12 -0x1.b07516p-13 0x1.d1c2c9p-16 0x1.62b5cap-15 -0x1.4f520ap-14 0x1.f7b810p-18 0x1.858b30p-15 -0x1.9afe04p-17 -0x1.0c9ea6p-16 0x1.b1735bp-18 0x1.406bf3p-19 -0x1.f24228p-20 0x1.1836a4p-20 -0x1.5c9060p-27 -0x1.f18265p-21 0x1.621471p-22 0x1.5f57e6p-22 -0x1.8642e7p-23 -0x1.27b2bep-26 0x1.29eb1ap-25 -0x1.b21afap-25 0x1.7b77bfp-26 0x1.12c80dp-25 -0x1.a1abfap-26 -0x1.1e349cp-27 0x1.9f3b66p-27 -0x1.00bfadp-29 -0x1.9a1c4dp-29 0x1.a86efbp-29 -0x1.39c72ap-31 -0x1.c27878p-30 0x1.14f638p-30 0x1.da7dacp-32 -0x1.2a7610p-31


	//namespace bmp = boost::multiprecision;

		int32_t SturmSequence::trimseq(cheb::scalar* seq, int32_t n, double l) {
			if (n == 0)
				return 0;

			int32_t i = n - 1;
			auto num = (n + 1) / 2;

			cheb::scalar mag = -DBL_MAX;
			for (int32_t j = 0; j < n / 2; ++j)
				mag = num == j ? std::max(mag, std::abs(seq[j])) : std::max(mag, std::abs(seq[j]));
			mag = std::max(std::abs(seq[n / 2]), mag);
			cheb::scalar e = std::nextafter(mag, DBL_MAX) - mag;

			//e *= 4.0;

			auto t = 2.0 * 1.0 * l * std::log(l + 1.0) * e * dynamicRange * dynamicRange;
			t = 0.0;
			//auto t = 2.0 * 10.0 * l * l * e;

			if (std::log10(std::abs(mag) / std::abs(seq[i])) > thresholdA / std::log10(l+1.0))
			for (; i > n / 2 - 1; i--) {
				if (std::abs(seq[i-1]) / std::abs(seq[i]) < thresholdB * (l+1.0)) {
				//if (2.0 * std::abs(seq[i]) > t) {
					//std::cout << "Trimming sequence at " << i << " with " << std::hexfloat << seq[i] << std::endl;
					break;
				}
			}
			auto d = n - i - 1;
			//if (d != 0)
			//	std::cout << "Skipped values at " << i << " [ " << std::hexfloat << seq[i + 1] * 2.0 << " / " << t << " ], in iteration " << l << std::endl;
			//else
			//	std::cout << "Skipped nothing with " << std::hexfloat << seq[n - 1] * 2.0 << "/ " << t << " in iteration " << l << std::endl;
			//std::cout << "Chopped " << d << " digits " << n << " -> " << i + 1 - d<< std::endl;
			for (int32_t j = 0; j < i + 1 - d; ++j) {
				seq[j] = seq[d + j];
			}

			return i + 1 - d;
		}
		//#define USE_HIGH_PRECISION
#ifndef USE_HIGH_PRECISION
		void SturmSequence::writeZSeries(cheb::scalar* in, int32_t n, cheb::scalar* out) {
			for (auto i = 0; i < n; i++) {
				out[n - 1 + i] = in[i] * 0.5;
				out[n - 1 - i] = in[i] * 0.5;
			}
			out[n - 1] = in[0];
			auto m = -DBL_MAX;
			for (int32_t i = 0; i < 2 * n - 1; ++i)
				m = std::max(std::abs(out[i]), m);
			auto f = std::abs(in[0]);
			//for (int32_t i = 0; i < 2*n - 1; ++i)
			//	out[i] /= m;
		};

		int32_t SturmSequence::zseries_div(cheb::scalar* z1, int32_t lc1, cheb::scalar* z2, int32_t lc2, cheb::scalar* z3) {
			if (lc2 == 1 || lc1 < lc2) {
				memset(z3, 0x00, sizeof(cheb::scalar) * lc1);
				for (int32_t i = 0; i < lc1; ++i)
					z3[i] = -z3[i];
				return lc1;
			}
			auto scl = z2[0];
			auto i = 0;
			auto j = lc1 - lc2;
			for (int32_t i = 0; i < lc1; ++i)
				z3[i] = 0.;
			while (i < j) {
				//auto r = (z3[i] + z1[i]) * scl;
				auto r = (z1[i]);
				for (int32_t k = 0; k < lc2; ++k) {
					//z3[i + k] = -r * z2[k] + z3[i + k];
					//z3[j + k] = -r * z2[k] + z3[j + k];
					//z3[i + k] = std::fma(-r, z2[k] * scl, z3[i + k]);
					//z3[j + k] = std::fma(-r, z2[k] * scl, z3[j + k]);
					z1[i + k] = std::fma(-r, z2[k] / scl, z1[i + k]);
					z1[j + k] = std::fma(-r, z2[k] / scl, z1[j + k]);
				}
				i++;
				j--;
			}
			//for (int32_t k = 0; k < lc1; ++k)
				//z1[k] += z3[k];
			auto r = z1[i];
			for (int32_t k = 0; k < lc2 - 2; ++k) {
				//z3[k] = r * z2[k + 1] - z1[i + 1 + k];
				z3[k] = std::fma(-r, z2[k + 1] / scl, z1[i + 1 + k]);
			}
			for (int32_t k = 0; k < lc2 - 2; ++k) {
				z3[k] = -z3[k];
			}
			//std::cout << lc1 << " : " << lc2 << " : " << lc2 - 2 << std::endl;
			return lc2 - 2;
		}

		int32_t SturmSequence::chebdiv(cheb::scalar* z1, int32_t lc1, cheb::scalar* z2, int32_t lc2, cheb::scalar* z3) {
			if (clen(lc1) < clen(lc2)) {
				memset(z3, 0x00, sizeof(cheb::scalar) * lc1);
				for (int32_t i = 0; i < lc1; ++i)
					z3[i] = -z3[i];
				return lc1;
			}
			if (clen(lc2) == 1) {
				z3[0] = 0.0;
				return 1;
			}
			return zseries_div(z1, lc1, z2, lc2, z3);
		}

		std::pair<int32_t, std::pair<cheb::scalar, std::vector<cheb::scalar>>> SturmSequence::zseries_div_withq(scalar* z1, int32_t lc1, scalar* z2, int32_t lc2, scalar* z3) {
			if (lc1 < lc2) {
				memset(z3, 0x00, sizeof(cheb::scalar) * lc1);
				auto mag = -DBL_MAX;
				for (int32_t i = 0; i < lc1; ++i) {
					z3[i] = -z3[i];
					mag = std::max(mag, std::abs(z3[i]));
				}
				return std::make_pair(lc1, std::make_pair(mag, std::vector<scalar>{0.0}));
			}
			if (lc2 == 1) {
				std::vector<scalar> q(z1, z1 + lc1);
				for (auto& v : q)
					v /= z2[0];
				z3[0] = 0.0;
				return std::make_pair(1, std::make_pair(1.0, q));
			}
			auto scl = z2[0];
			auto i = 0;
			auto j = lc1 - lc2;
			auto dlen = j;
			std::vector<scalar> quo(j + 1, 0.0);
			for (int32_t i = 0; i < lc1; ++i)
				z3[i] = 0.;
			while (i < j) {
				//auto r = (z3[i] + z1[i]) * scl;
				auto r = (z1[i]);
				quo[i] = z1[i];
				quo[dlen - i] = r;
				for (int32_t k = 0; k < lc2; ++k) {
					//z3[i + k] = -r * z2[k] + z3[i + k];
					//z3[j + k] = -r * z2[k] + z3[j + k];
					//z3[i + k] = std::fma(-r, z2[k] * scl, z3[i + k]);
					//z3[j + k] = std::fma(-r, z2[k] * scl, z3[j + k]);
					z1[i + k] = std::fma(-r, z2[k] / scl, z1[i + k]);
					z1[j + k] = std::fma(-r, z2[k] / scl, z1[j + k]);
				}
				i++;
				j--;
			}
			//for (int32_t k = 0; k < lc1; ++k)
				//z1[k] += z3[k];
			auto r = z1[i];
			quo[i] = r;
			for (auto& q : quo)
				q /= scl;
			for (int32_t k = 0; k < lc2 - 2; ++k) {
				//z3[k] = r * z2[k + 1] - z1[i + 1 + k];
				z3[k] = std::fma(-r, z2[k + 1] / scl, z1[i + 1 + k]);
			}
			for (int32_t k = 0; k < lc2 - 2; ++k) {
				z3[k] = -z3[k];
			}
			auto mag = -DBL_MAX;
			for (int32_t i = 0; i < lc2 - 2; ++i)
				mag = std::max(mag, std::abs(z3[i]));
			//std::cout << lc1 << " : " << lc2 << " : " << lc2 - 2 << std::endl;
			return std::make_pair(lc2 - 2,std::make_pair(mag,quo));
		}
		std::pair<int32_t, std::pair<cheb::scalar, std::vector<cheb::scalar>>> SturmSequence::chebdiv_withq(scalar* z1, int32_t lc1, scalar* z2, int32_t lc2, scalar* z3) {
			if (clen(lc1) < clen(lc2)) {
				memset(z3, 0x00, sizeof(cheb::scalar) * lc1);
				for (int32_t i = 0; i < lc1; ++i)
					z3[i] = -z3[i];
				return std::make_pair(lc1, std::make_pair(1.0, std::vector<scalar>{ 0.0 }));
			}
			if (clen(lc2) == 1) {
				z3[0] = 0.0;
				std::vector<scalar> q(z1, z1 + lc1);
				for (auto& v : q)
					v /= z2[lc2 - 1];

				return std::make_pair(1, std::make_pair(1.0, q));
			}
			auto [l,q] =  zseries_div_withq(z1, lc1, z2, lc2, z3);

			//std::cout << "------------------------------------------" << std::endl;
			//std::cout << "rem: [";
			//for (int32_t i = 0; i < l; ++i)
			//	std::cout << z3[i] << " ";
			//std::cout << "]\n";
			//std::cout << "quo: ";
			//cheb::printVec(q.second);
			//std::cout << "\n";
			auto n = (q.second.size() + 1) / 2;
			q.second = std::vector<scalar>(q.second.begin() + n - 1, q.second.end());
			for (int32_t i = 1; i < q.second.size(); ++i)
				q.second[i] *= 2.0;
			//std::cout << "quo: ";
			//cheb::printVec(q.second);
			//std::cout << "\n";
			//std::cout << "magnitude: " << q.first << std::endl;
			return std::make_pair(l, q);
		}


		void SturmSequence::symmetrize(cheb::scalar* z, int32_t n) {
			int32_t i = 0;
			int32_t j = n - 1;
			while (i < j)
				//z[i++] = z[j--];
				z[j--] = z[i++];
		}
		void SturmSequence::setDegreeWithZseries(std::size_t rank, cheb::scalar* zs, int32_t len) {
			auto n = (len + 1) / 2;
			auto l = len - n + 1;
			lengths[rank] = l;
			auto m = -DBL_MAX;
			for (std::size_t k = 0; k < len; ++k)
				m = std::max(m, std::abs(zs[k]));
			//auto r = std::abs(zs[0]);
			//for (std::size_t k = 0; k < len; ++k)
			//	zs[k] /= m;
			for (std::size_t k = 0; k < lengths[rank]; ++k)
				_data[rank * degree + k] = (k == 0 ? zs[n - 1 + k] : 2.0 * zs[n - 1 + k]);
			//printDegree(rank);
		}
		void SturmSequence::generateSequence() {
			auto z2 = new scalar[2 * degree - 1];
			auto z1 = new scalar[2 * degree - 1];
			auto z0 = new scalar[2 * degree - 1];
			writeZSeries(getDegree(0), (int32_t)degree, z2);
			writeZSeries(getDegree(1), (int32_t)degree - 1, z1);

			memset(z0, 0x00, degree * sizeof(scalar));

			auto n2 = 2 * (int32_t)degree - 1;
			auto n1 = 2 * (int32_t)(degree - 1) - 1;
			int32_t i = 1;
			double prior = getDegree(1)[lengths[1] - 1];
			while (true) {
				auto [n0,q] = chebdiv_withq(z2, n2, z1, n1, z0);
				//cheb::printVec(z0);
				i++;
				symmetrize(z0, n0);
				n0 = trimseq(z0, n0, (double)i - 2.0);


				if (clen(n0) == 0 || z0[0] == 0.0)
					break;
				//std::cout << i << " => " << n2 << " @ " << clen(n2) << " -> " << n1 << " @ " << clen(n1) <<  " -> " << n0 <<  " @ " <<clen(n0) << " -> " << z0[0] << std::endl;;
				auto m = 1.0;
				for (auto v : q.second)
					m = std::max(m, v);
				//for (std::size_t k = 0; k < n0; ++k)
					//m = std::max(m, std::abs(z0[k]));
				for (std::size_t k = 0; k < n0; ++k)
					z0[k] /= m;

				q.first = m;
				quotientSeries.push_back(q);

				setDegreeWithZseries(i, z0, n0);

				auto max = -DBL_MAX;
				for (int32_t j = 0; j < lengths[i]; ++j)
					max = std::max(max, std::abs(getDegree(i)[j]));
				auto ratio = std::abs(getDegree(i)[0]) / max;
				if (lengths[i] >= 2 && n1 >= 1) {
					//
					//auto ratio = prior / std::abs(getDegree(i)[lengths[i] - 1]);
					if(std::abs(getDegree(i)[lengths[i]-2]) / std::abs(getDegree(i)[lengths[i] - 1]) > 1e4)
					//if(std::log10(ratio)>3.0)
					//if (std::log10(ratio) < -6.0)
					{
						//printf("%d / %d -> %d: %g / %g -> %g => %g ==> %d\n", n2, n1, n0, std::abs(getDegree(i)[0]), max, ratio, std::log10(ratio), (std::log10(ratio) < -10.0 ? 1 : 0));
						potentialProblem = true;
					}
					prior = std::abs(getDegree(i)[lengths[i] - 1]);
				}
				std::swap(z2, z1);
				std::swap(z1, z0);
				std::swap(n2, n1);
				std::swap(n1, n0);
			}
			delete[] z2;
			delete[] z1;
			delete[] z0;
			terms = i;
		}
#else
		void SturmSequence::writeZSeries(scalar* in, int32_t n, scalarHighPrecision* out) {
			for (auto i = 0; i < n; i++) {
				out[n - 1 + i] = in[i] * 0.5;
				out[n - 1 - i] = in[i] * 0.5;
			}
			out[n - 1] = in[0];
			scalarHighPrecision m = -std::numeric_limits<scalarHighPrecision>::max();
			for (int32_t i = 0; i < 2 * n - 1; ++i) {
				scalarHighPrecision v = bmp::abs(out[i]);
				m = (scalarHighPrecision)bmp::max(v, m);
			}
			auto f = std::abs(in[0]);
			for (int32_t i = 0; i < 2 * n - 1; ++i)
				out[i] /= m;
		};

		int32_t SturmSequence::zseries_div(scalarHighPrecision* z1, int32_t lc1, scalarHighPrecision* z2, int32_t lc2, scalarHighPrecision* z3) {
			if (lc2 == 1 || lc1 < lc2) {
				memset(z3, 0x00, sizeof(scalarHighPrecision) * lc1);
				for (int32_t i = 0; i < lc1; ++i)
					z3[i] = -z3[i];
				return lc1;
			}
			scalarHighPrecision scl = ((scalarHighPrecision)1.) / z2[0];
			auto i = 0;
			auto j = lc1 - lc2;
			for (int32_t i = 0; i < lc1; ++i)
				z3[i] = 0.;
			while (i < j) {
				auto r = z1[i];
				for (int32_t k = 0; k < lc2; ++k) {
					z3[i + k] = -r * z2[k] * scl + z3[i + k];
					z3[j + k] = -r * z2[k] * scl + z3[j + k];
				}
				i++;
				j--;
			}
			for (int32_t i = 0; i < lc1; ++i)
				z1[i] += z3[i];
			auto r = z1[i];
			for (int32_t k = 0; k < lc2 - 2; ++k) {
				z3[k] = r * z2[k + 1] * scl - z1[i + 1 + k];
			}
			return lc2 - 2;
		}

		int32_t SturmSequence::chebdiv(scalarHighPrecision* z1, int32_t lc1, scalarHighPrecision* z2, int32_t lc2, scalarHighPrecision* z3) {
			if (clen(lc1) < clen(lc2)) {
				memset(z3, 0x00, sizeof(scalarHighPrecision) * lc1);
				for (int32_t i = 0; i < lc1; ++i)
					z3[i] = -z3[i];
				return lc1;
			}
			if (clen(lc2) == 1) {
				z3[0] = 0.0;
				return 1;
			}
			return zseries_div(z1, lc1, z2, lc2, z3);
		}


		void SturmSequence::symmetrize(scalarHighPrecision* z, int32_t n) {
			int32_t i = 0;
			int32_t j = n - 1;
			while (i < j)
				z[j--] = z[i++];
		}
		void SturmSequence::setDegreeWithZseries(std::size_t rank, scalarHighPrecision* zs, int32_t len) {
			auto n = (len + 1) / 2;
			auto l = len - n + 1;
			lengths[rank] = l;
			for (std::size_t k = 0; k < lengths[rank]; ++k) {
				scalarHighPrecision res = 0.;
				if (k == 0)
					res = zs[n - 1 + k];// / bmp::abs(zs[n - 1]);
				else
					res = ((scalarHighPrecision)2.0) * zs[n - 1 + k];// / bmp::abs(zs[n - 1]);
				_data[rank * degree + k] = (double)res;
			}
			scalarHighPrecision m = -std::numeric_limits< scalarHighPrecision>::max();
			for (std::size_t k = 0; k < len; ++k) {
				scalarHighPrecision v = (scalarHighPrecision)bmp::abs(zs[k]);
				m = (scalarHighPrecision)bmp::max(m, v);
			}
			auto r = bmp::abs(zs[0]);
			for (std::size_t k = 0; k < len; ++k)
				zs[k] = (zs[k] / m);
			//printDegree(rank);
		}
		void SturmSequence::generateSequence() {
			auto z2 = new scalarHighPrecision[2 * degree - 1];
			auto z1 = new scalarHighPrecision[2 * degree - 1];
			auto z0 = new scalarHighPrecision[2 * degree - 1];
			writeZSeries(getDegree(0), (int32_t)degree, z2);
			writeZSeries(getDegree(1), (int32_t)degree - 1, z1);

			memset(z0, 0x00, degree * sizeof(scalarHighPrecision));

			auto n2 = 2 * (int32_t)degree - 1;
			auto n1 = 2 * (int32_t)(degree - 1) - 1;
			int32_t i = 1;
			while (true) {
				auto n0 = chebdiv(z2, n2, z1, n1, z0);
				//cheb::printVec(z0);
				i++;
				symmetrize(z0, n0);
				//n0 = trimseq(z0, n0);

				if (clen(n0) == 0 || z0[0] == 0.0)
					break;
				//std::cout << i << " => " << n2 << " @ " << clen(n2) << " -> " << n1 << " @ " << clen(n1) <<  " -> " << n0 <<  " @ " <<clen(n0) << " -> " << z0[0] << std::endl;;
				setDegreeWithZseries(i, z0, n0);
				std::swap(z2, z1);
				std::swap(z1, z0);
				std::swap(n2, n1);
				std::swap(n1, n0);
			}
			delete[] z2;
			delete[] z1;
			delete[] z0;
			terms = i;
		}

#endif

		void SturmSequence::diffTerm(const cheb::svec& polyNomial) {
			if (degree == 1) {
				lengths[1] = 1;
				_data[1] = 0.0;
				return;
			}
			lengths[1] = degree - 1;
			auto n = degree;
			auto zk = getDegree(1);
			if (n - 1 > 0)
				zk[n - 2] = 2.0 * ((scalar)(n - 2) + 1.0) * polyNomial[n - 2 + 1];
			if (n - 1 > 1)
				zk[n - 3] = 2.0 * ((scalar)(n - 3) + 1.0) * polyNomial[n - 3 + 1];
			for (int32_t i = (int32_t)n - 4; i >= 0; i -= 2)
				zk[i] = 2.0 * ((scalar)(i)+1.0) * polyNomial[(std::size_t)i + 1] + zk[i + 2];
			for (int32_t i = (int32_t)n - 5; i >= 0; i -= 2)
				zk[i] = 2.0 * ((scalar)(i)+1.0) * polyNomial[(std::size_t)i + 1] + zk[i + 2];
			zk[0] = .5 * zk[0];
		}
		void SturmSequence::generateDerivative() {
			if (degree == 1) {
				lengths[1] = 1;
				_data[1] = 0.0;
				return;
			}
			lengths[1] = degree - 1;
			auto n = degree;
			auto zk = getDegree(1);
			if (n - 1 > 0)
				zk[n - 2] = 2.0 * ((scalar)(n - 2) + 1.0) * _data[n - 2 + 1];
			if (n - 1 > 1)
				zk[n - 3] = 2.0 * ((scalar)(n - 3) + 1.0) * _data[n - 3 + 1];
			for (int32_t i = (int32_t)n - 4; i >= 0; i -= 2)
				zk[i] = 2.0 * ((scalar)(i)+1.0) * _data[(std::size_t)i + 1] + zk[i + 2];
			for (int32_t i = (int32_t)n - 5; i >= 0; i -= 2)
				zk[i] = 2.0 * ((scalar)(i)+1.0) * _data[(std::size_t)i + 1] + zk[i + 2];
			zk[0] = .5 * zk[0];

			auto mmc = -DBL_MAX;
			for (int32_t i = 0; i < n - 1; i++)
				mmc = std::max(mmc, std::abs(zk[i]));
			for (int32_t i = 0; i < n - 1; i++)
				zk[i] /= mmc;

		}
		cheb::svec SturmSequence::differentiate(const cheb::svec& polyNomial) {
			int32_t degree = polyNomial.size();
			auto n = degree;
			cheb::svec zk(n - 1, 0.0);
			if (n - 1 > 0)
				zk[n - 2] = 2.0 * ((scalar)(n - 2) + 1.0) * polyNomial[n - 2 + 1];
			if (n - 1 > 1)
				zk[n - 3] = 2.0 * ((scalar)(n - 3) + 1.0) * polyNomial[n - 3 + 1];
			for (int32_t i = (int32_t)n - 4; i >= 0; i -= 2)
				zk[i] = std::fma(2.0 * ((scalar)(i)+1.0), polyNomial[(std::size_t)i + 1], zk[i + 2]);
			for (int32_t i = (int32_t)n - 5; i >= 0; i -= 2)
				zk[i] = std::fma(2.0 * ((scalar)(i)+1.0), polyNomial[(std::size_t)i + 1], zk[i + 2]);
			zk[0] = .5 * zk[0];
			return zk;
		}
		void SturmSequence::initData(const cheb::svec& polyNomial) {
			degree = polyNomial.size();
			if (degree == 0)
				throw std::invalid_argument("Invalid polynomial length of 0.");
			lengths = std::vector<std::size_t>(std::max(degree,(std::size_t) 2ull), 0);
			_data = cheb::svec(degree > 1 ? degree * degree : 2ull, 0.0);
			auto p = polyNomial;
			//for (auto& a : p)
			//	a /= polyNomial[0];
			setDegree(0, p);
			auto mag = -DBL_MAX;
			for (auto v : polyNomial)
				mag = std::max(std::abs(v), mag);
			for (int32_t i = 0; i < degree; ++i)
				_data[i] /= mag;

			auto ma = -DBL_MAX;
			auto mi = DBL_MAX;
			for (auto v : polyNomial) {
				ma = std::max(std::abs(v), ma);
				mi = std::min(std::abs(v), mi);
			}
			dynamicRange = std::log2(ma / mi);

			generateDerivative();
			//diffTerm(p);
			if (degree > 1)
				generateSequence();
		}
		cheb::scalar* SturmSequence::getDegree(int32_t rank) {
			return _data.data() + rank * degree;
		}
		void SturmSequence::setDegreeWithZseries(std::size_t rank, const cheb::svec& zs) {
			auto n = (zs.size() + 1) / 2;
			auto l = zs.size() - n + 1;
			lengths[rank] = l;
			for (std::size_t k = 0; k < lengths[rank]; ++k)
				_data[rank * degree + k] = (k == 0 ? zs[n - 1 + k] : 2.0 * zs[n - 1 + k]);
			//printDegree(rank);
		}
		void SturmSequence::setDegree(std::size_t rank, const cheb::svec& data) {
			lengths[rank] = data.size();
			for (std::size_t k = 0; k < lengths[rank]; ++k)
				_data[rank * degree + k] = data[k];
			//printDegree(rank);
		}
		void SturmSequence::printDegree(int32_t rank) {
			std::cout << "Rank [ " << rank << " ] has degree " << lengths[rank] << " -> ";
			auto vec = cheb::svec(getDegree(rank), getDegree(rank + 1));
			std::cout << "[ ";
			if (vec.size() > 0)
				for (int32_t i = 0; i < vec.size() - 1; ++i)
					std::cout << vec[i] << " ";
			if (vec.size() > 0)
				std::cout << (*(std::end(vec) - 1));
			std::cout << "]\n";
		}

		SturmSequence::SturmSequence(const cheb::svec& polyNomial, double tA, double tB, bool qs):thresholdA(tA), thresholdB(tB),useQuotientSeries(qs) {
			initData(polyNomial);
		}
		std::pair<int32_t, int32_t> SturmSequence::printSturm(std::pair<scalar, scalar> lr) {
			std::cout << "Printing Sturm Sequence of degree " << degree << std::endl;
			auto [xl, xr] = lr;
			std::cout << "Interval range: " << xl << " x " << xr << " -> " << xr - xl << std::endl;
			auto init = clenshaw(std::pair{ xl , xr }, getDegree(0), (int32_t)lengths[0]);

			auto [pl, pr] = sign(init);
			std::cout << "Initivial values: " << init.first << " x " << init.second << " -> " << pl << " x " << pr << std::endl;
			//std::cout << init << " -> " << prev << std::endl;
			auto ctrs = std::make_pair(0, 0);
			for (int32_t i = 1; i < terms; ++i) {
				std::cout << "\tDegree: " << i << " => Counters: " << ctrs.first << " x " << ctrs.second << " ";
				auto v = clenshaw(std::pair{ xl, xr }, getDegree(i), (int32_t)lengths[i]);
				auto [vl, vr] = sign(v);
				std::cout << "\tValues: " << std::hexfloat << std::setprecision(13) << v.first << " x " << std::hexfloat << std::setprecision(13) << v.second << " -> " << vl << " x " << vr << " ";
				ctrs.first += pl * vl < -DBL_EPSILON ? 1 : 0;
				ctrs.second += pr * vr < -DBL_EPSILON ? 1 : 0;
				std::cout << "\tSign Change Value: " << pl * vl << " x " << pr << " x " << vr << std::endl;
				pl = vl;
				pr = vr;
			}
			std::cout << "Sturm Sequence evaluation finished with " << ctrs.first << " x " << ctrs.second << " => roots: " << ctrs.second - ctrs.first << std::endl;
			return ctrs;

		}


		int32_t SturmSequence::numRootsRemainderSeries(std::pair<scalar, scalar> lr) {
			//return 0;
			auto [xl, xr] = lr;
			auto init = clenshaw(std::pair{ xl , xr }, getDegree(0), (int32_t)lengths[0]);
			auto [pl, pr] = sign(init);
			//std::cout << init << " -> " << prev << std::endl;
			auto ctrs = std::make_pair(0, 0);
			for (int32_t i = 1; i < terms; ++i) {
				auto v = clenshaw(std::pair{ xl, xr }, getDegree(i), (int32_t)lengths[i]);
				auto [vl, vr] = sign(v);
				ctrs.first += pl * vl < 0.0 ? 1 : 0;
				ctrs.second += pr * vr < 0.0 ? 1 : 0;
				pl = vl;
				pr = vr;
			}
			auto [l, r] = ctrs;
			return l - r;
		}
		int32_t SturmSequence::numRootsQuotientSeries(std::pair<scalar, scalar> lr) {
			auto [xl, xr] = lr;
			auto p2 = clenshaw(std::pair{ xl , xr }, getDegree(0), (int32_t)lengths[0]);
			auto p1 = clenshaw(std::pair{ xl , xr }, getDegree(1), (int32_t)lengths[1]);
			auto [p1l, p1r] = sign(p1);
			auto [p2l, p2r] = sign(p2);
			//std::cout << init << " -> " << prev << std::endl;
			auto ctrs = std::make_pair(p1l * p2l < 0 ? 1 : 0, p1r* p2r < 0 ? 1 : 0);
			for (int32_t i = 2; i < terms; ++i) {
				auto [m, q] = quotientSeries[i - 2];
				auto cls = clenshaw(std::pair{ xl, xr }, q.data(), q.size());
				auto v = std::make_pair(
					std::fma(cls.first, p1.first, -p2.first) / m,
					std::fma(cls.second, p1.second, -p2.second) / m
				);
				auto [vl, vr] = sign(v);
				ctrs.first += p1l * vl < -DBL_EPSILON ? 1 : 0;
				ctrs.second += p1r * vr < -DBL_EPSILON ? 1 : 0;
				p1l = vl;
				p1r = vr;
				p2 = p1;
				p1 = v;
			}
			return ctrs.first - ctrs.second;
		}
		int32_t SturmSequence::numRoots(std::pair<scalar, scalar> lr) {
			if (useQuotientSeries)
				return numRootsQuotientSeries(lr);
			return numRootsRemainderSeries(lr);
		}

		std::pair<int32_t, int32_t> SturmSequence::evalSturmRemainderSeries(std::pair<scalar, scalar> lr ) {
			//return 0;
			auto [xl, xr] = lr;
			auto init = clenshaw(std::pair{ xl , xr }, getDegree(0), (int32_t)lengths[0]);
			auto [pl, pr] = sign(init);
			//std::cout << init << " -> " << prev << std::endl;
			auto ctrs = std::make_pair(0, 0);
			for (int32_t i = 1; i < terms; ++i) {
				auto v = clenshaw(std::pair{ xl, xr }, getDegree(i), (int32_t)lengths[i]);
				auto [vl, vr] = sign(v);
				ctrs.first += pl * vl < -DBL_EPSILON ? 1 : 0;
				ctrs.second += pr * vr < -DBL_EPSILON ? 1 : 0;
				pl = vl;
				pr = vr;
			}
			return ctrs;
		}
		std::pair<int32_t, int32_t> SturmSequence::evalSturmQuotientSeries(std::pair<scalar, scalar> lr) {
			auto [xl, xr] = lr;
			auto p2 = clenshaw(std::pair{ xl , xr }, getDegree(0), (int32_t)lengths[0]);
			auto p1 = clenshaw(std::pair{ xl , xr }, getDegree(1), (int32_t)lengths[1]);
			auto [p1l, p1r] = sign(p1);
			auto [p2l, p2r] = sign(p2);
			//std::cout << init << " -> " << prev << std::endl;
			auto ctrs = std::make_pair(p1l * p2l < 0 ? 1 : 0, p1r* p2r < 0 ? 1 : 0);
			for (int32_t i = 2; i < terms; ++i) {
				auto [m, q] = quotientSeries[i - 2];
				auto cls = clenshaw(std::pair{ xl, xr }, q.data(), q.size());
				auto v = std::make_pair(
					std::fma(cls.first, p1.first, -p2.first) / m,
					std::fma(cls.second, p1.second, -p2.second) / m
				);
				auto [vl, vr] = sign(v);
				ctrs.first += p1l * vl < -DBL_EPSILON ? 1 : 0;
				ctrs.second += p1r * vr < -DBL_EPSILON ? 1 : 0;
				p1l = vl;
				p1r = vr;
				p2 = p1;
				p1 = v;
			}
			return ctrs;
		}
		std::pair<int32_t, int32_t> SturmSequence::evalSturm(std::pair<scalar, scalar> lr) {
			if (useQuotientSeries)
				return evalSturmQuotientSeries(lr);
			return evalSturmRemainderSeries(lr);
		}
		std::array<int32_t,32> SturmSequence::evalSturmRemainderSeries(std::array<scalar, 32> vals) {
			//return 0;
			auto init = clenshaw(vals, getDegree(0), (int32_t)lengths[0]);
			auto ps = sign(init);
			//std::cout << init << " -> " << prev << std::endl;
			std::array<int32_t, 32> ctrs;
			ctrs.fill(0);
			for (int32_t i = 1; i < terms; ++i) {
				auto v = clenshaw(vals, getDegree(i), (int32_t)lengths[i]);
				auto vs = sign(v);
				for (int32_t i = 0; i < 32; ++i)
					ctrs[i] += ps[i] * vs[i] < 0.0 ? 1 : 0;
				ps = vs;
			}
			return ctrs;
		}
		std::array<int32_t, 32> SturmSequence::evalSturmQuotientSeries(std::array<scalar, 32> vals) {
			//return 0;
			auto p2 = clenshaw(vals, getDegree(0), (int32_t)lengths[0]);
			auto p1 = clenshaw(vals, getDegree(1), (int32_t)lengths[1]);
			//std::cout << init << " -> " << prev << std::endl;
			std::array<int32_t, 32> ctrs;
			ctrs.fill(0);
			auto p2s = sign(p2);
			auto p1s = sign(p1);
			for (int32_t i = 0; i < 32; ++i)
				ctrs[i] = p1s[i] * p2s[i] < 0 ? 1 : 0;

			for (int32_t i = 2; i < terms; ++i) {
				auto [m, q] = quotientSeries[i - 2];
				auto cls = clenshaw(vals, q.data(), q.size());
				std::array<scalar,32> v;
				for (int32_t i = 0; i < 32; ++i)
					v[i] = std::fma(cls[i], p1[i], -p2[i])/m;
				auto vs = sign(v);
				for (int32_t i = 0; i < 32; ++i)
					ctrs[i] += vs[i] * p1s[i] < 0 ? 1 : 0;
				p1s = vs;
				p2 = p1;
				p1 = v;
			}
			return ctrs;
		}

		std::array<int32_t, 32> SturmSequence::evalSturm(std::array<scalar, 32> vals) {
			if (useQuotientSeries)
				return evalSturmQuotientSeries(vals);
			return evalSturmRemainderSeries(vals);
		}

		int32_t SturmSequence::evalSturmRemainderSeries(scalar xl) {
			//return 0;
			auto init = clenshaw(xl, getDegree(0), (int32_t)lengths[0]);
			auto pl = sign(init);
			//std::cout << init << " -> " << prev << std::endl;
			auto ctrs = 0;
			for (int32_t i = 1; i < terms; ++i) {
				auto v = clenshaw(xl, getDegree(i), (int32_t)lengths[i]);
				auto vl = sign(v);
				ctrs += pl * vl < -DBL_EPSILON ? 1 : 0;
				pl = vl;
			}
			return ctrs;
		}
		int32_t SturmSequence::evalSturmQuotientSeries(scalar xl) {
			auto p2 = clenshaw(xl, getDegree(0), (int32_t)lengths[0]);
			auto p1 = clenshaw(xl, getDegree(1), (int32_t)lengths[1]);
			auto p1l = sign(p1);
			auto p2l = sign(p2);
			//std::cout << init << " -> " << prev << std::endl;
			auto ctrs = p1l * p2l < 0 ? 1 : 0;
			for (int32_t i = 2; i < terms; ++i) {
				auto [m, q] = quotientSeries[i - 2];
				auto cls = clenshaw(xl, q.data(), q.size());
				auto v = std::fma(cls, p1, -p2) / m;
				auto vl = sign(v);
				ctrs += p1l * vl < -DBL_EPSILON ? 1 : 0;
				p1l = vl;
				p2 = p1;
				p1 = v;
			}
			return ctrs;
		}
		int32_t SturmSequence::evalSturm(scalar vals) {
			if (useQuotientSeries)
				return evalSturmQuotientSeries(vals);
			return evalSturmRemainderSeries(vals);
		}

		double SturmSequence::firstRoot(double lli, double rri) {
			double li = lli;
			double ri = rri;
			if (degree == 1) {
				if (getDegree(0)[0] == 0.)
					return li;
				return -2.;
			}
			if (degree == 2) {
				return -getDegree(0)[0] / getDegree(0)[1];
			}
			auto dfdx = getDegree(1);
			auto df2dx = differentiate(cheb::svec(dfdx, dfdx + lengths[1]));
			auto x = li;

			auto [l, r] = evalSturm(std::make_pair(li, ri));
			if (l == r) return DBL_MAX;
			//std::cout << "Interval range: " << l << " : " << r << std::endl;
			auto gamma = 8.0;

			// pre narrow interval
			std::array<double, 32> pcandidates;
			for (int32_t i = 0; i < 32; ++i)
				pcandidates[i] = li + (ri - li) * (double)i / 31.0;
			auto seqc = evalSturm(pcandidates);
			int32_t start = 0;
			for (int32_t i = 0; i < 31; ++i) {
				if (seqc[i] != seqc[(std::size_t) i + 1]) {
					start = i;
					x = li + (ri - li) * (double)i / 31.0;
					ri = li + (ri - li) * (double)(i + 1) / 31.0;
					li = x;
					break;
				}
			}
			int32_t j = 0;
			int32_t i = 0;
			double dx = 0.0;
			double y = 0.0;
			for (; j < 8; ) {
				auto fx = clenshaw(x, getDegree(0), (int32_t)lengths[0]);
				auto f1x = clenshaw(x, dfdx, (int32_t)lengths[1]);
				auto f2x = clenshaw(x, df2dx.data(), (int32_t)df2dx.size());
				if (std::abs(f1x) < DBL_EPSILON) {
					if (std::abs(fx) < DBL_EPSILON)
						break;
					else
						f1x = -1.0;
				}
#define HALLEY
#ifdef HALLEY
				auto hx = -fx / f1x / (1.0 - fx / f1x * f2x / (2.0 * f1x));
				auto nx = -f1x / f2x;
				auto gx = -f1x;
				dx = -hx;
				if (hx < 0. || hx != hx) {
					if (nx > 0. && nx == nx)
						dx = -nx;
					else {
						if (gx > 0.)
							dx = -gx;
						else
							dx = gx;
					}
				}
				dx *= gamma;
				dx = std::clamp(dx, li - ri, 0.);
#else
				auto dx = -gamma * fx / f1x;
				if (dx > 0.0) {
					dx = -dx;
				}
#endif
				//auto dx = gamma * fx / f1x;
				//std::cout << "###################################################" << std::endl;
				//std::cout << "x = " << x << std::endl;
				//std::cout << "f(" << x << ") = " << clenshaw(x, getDegree(0), lengths[0]) << std::endl;
				//std::cout << "f(" << x << ") = " << clenshaw(x, getDegree(1), lengths[0]) << std::endl;
				//std::cout << "f/dx = " << dx << std::endl;
				//std::cout << "num roots: " << numRoots(std::make_pair(x, 1.0)) << std::endl;
				//if (dx > 0.0) {
					//std::cout << "Derivative points in wrong direction modifying gradient" << std::endl;
				//	dx = -dx;
					//std::cout << "Derivative was " << -dx << " is " << dx << std::endl;
				//}

				auto [lc, rc] = evalSturm(std::make_pair(x - 0.5 * dx, x - dx));
				//std::cout << "Sturm at half and full width: " << lc << " : " << rc << std::endl;
				//if (rc < l)
				//	std::cout << "Skipping zeroes at full step" << std::endl;
				//if (lc < l)
				//	std::cout << "Skipping zeroes at half step" << std::endl;
				//if (rc == l - 1 && dfdx(x + dx) <= 0.0)
				//    std::cout << "Safe skip" << std::endl;

				if (rc < l || lc < l) {
					//std::cout << "Testing different step widths" << std::endl;
					auto h = lc < l ? 0.5 : 1.0;
					std::array<double, 32> candidates;
					for (int32_t i = 0; i < 32; ++i)
						candidates[i] = x - dx * ((double)((std::size_t)i + 1)) / 32.0 * h;
					auto seq = evalSturm(candidates);
					if (l - seq[0] > 0) {
						gamma /= 4.0;
						//j++;
						dx = dx / 4.0;
						if (x - dx == x)break;
						li = x;
						ri = x - dx * h / 32.0;
						continue;
						//std::cout << "Skipped 0 with initial guess! reducing step width to " << gamma << " from  " << gamma * 4.0 << std::endl;
						//continue;
					}
					//gamma = 8.0;
					//std::cout << "Candidate intervals: " << std::endl;
					//int32_t i = 0;
					//for (; i < 32; ++i) {
					//    //std::cout << "\t" << i << " -> " << " x = " << candidates[i] << " -> " << seq[i] << " => " << seq[i] - r << std::endl;
					//    if (i >= 1 && seq[i - 1] - seq[i] != 0) {
					//        //std::cout << "\t\tSkipped " << seq[i - 1] - seq[i] << " zeroes." << std::endl;
					//        x = candidates[i - 1];
					//		dx = dx * ((double)(i)) / 32.0 * h;
					//        break;
					//    }
					//}
					//std::cout << "Linear search found x = " << (candidates[i - 1] ) << " for i = " << i - 1 << std::endl;
					//std::cout << "Binary search: " << std::endl;
					int32_t i = 16;
					int32_t s = 0;
					int32_t e = 31;
					for (int32_t j = 0; j <= 4; j++) {
						auto m = s + (e - s) / 2;
						auto ds = l - seq[m];
						//std::cout << "\tRange: [ " << s << " | " << m << " | " << e << "]" << std::endl;
						//std::cout << "\tdifference: " << ds << " [ " << l << " : " << seq[m] << "] => branching " << (ds > 0 ? "left" : "right") << std::endl;
						if (ds > 0)
							e = m;
						else
							s = m + 1;
					}
					// std::cout << "Binary search found x = " << candidates[s - 1] << " for i = " << s - 1 << std::endl;
					dx = dx * ((double)(s)) / 32.0 * h;
					x = candidates[(std::size_t)s - 1];
					ri = candidates[(std::size_t)s];
				}
				else
					x -= dx;
				++i;
				y = clenshaw(x, getDegree(0), (int32_t)lengths[0]);
				if (std::abs(y) < 1e0 * DBL_EPSILON || i > 64 || x - dx == x) {
					//std::cout << "Stopping Sturm-Newton method after " << i << " iterations with f(" << x << ") = " << y << std::endl;
					break;
				}
				li = x;
			}
			
//i = 0;
//j = 12;
int32_t it = 0;
static std::mutex mutex;
std::vector<cheb::scalar> points;
if (true)
			for (;; ) {
				auto fx = clenshaw(x, getDegree(0), (int32_t)lengths[0]);
				auto f1x = clenshaw(x, dfdx, (int32_t)lengths[1]);
				auto f2x = clenshaw(x, df2dx.data(), (int32_t)df2dx.size());

#define HALLEY
#undef HALLEY
#ifdef HALLEY
				auto hx = -fx / f1x / (1.0 - fx / f1x * f2x / (2.0 * f1x));
				auto nx = -f1x / f2x;
				auto gx = -f1x;
				dx = -hx;
				if (hx != hx) {
					if (nx == nx)
						dx = -nx;
					else {
						dx = -gx;
					}
				}
#else
				auto dx =0.1 *  fx / f1x;
				
#endif
				
				++it;

				{
					//std::lock_guard lg(mutex);
					//std::cout << "Iteration: " << it << "x : " << x << " dx: " << dx << " @ " << clenshaw(x, getDegree(0), (int32_t)lengths[0]) << " -> " << clenshaw(x - dx, getDegree(0), (int32_t)lengths[0]) << std::endl;
				}
				if (dx != dx)
					x = std::nextafter(x, DBL_MAX);
				else
					x = x - dx;
				x = std::clamp(x, lli, rri);
				//y = clenshaw(x - dx, getDegree(0), (int32_t)lengths[0]);
				if (it > 16 || x - dx == x || std::find(std::begin(points), std::end(points),x) != std::end(points)) {
					//std::cout << "Stopping Sturm-Newton method after " << i << " iterations with f(" << x << ") = " << y << std::endl;
					break;
				}
				points.push_back(x);
				//x = x - dx;
				
			}

			return std::abs(y) >= 1e2 * DBL_EPSILON ? firstRoot(std::nextafter(li, DBL_MAX), 1.0) : (i >= 32 ? DBL_MAX : x);


		}

		double SturmSequence::firstIntervalRoot(double lli, double rri) {
			double li = lli;
			double ri = rri;
			if (degree == 1) {
				if (getDegree(0)[0] == 0.)
					return li;
				return -2.;
			}
			if (degree == 2) {
				return -getDegree(0)[0] / getDegree(0)[1];
			}
			auto dfdx = getDegree(1);
			auto df2dx = differentiate(cheb::svec(dfdx, dfdx + lengths[1]));
			auto x = li;

			auto [l, r] = evalSturm(std::make_pair(li, ri));
			if (l == r) return DBL_MAX;

			for (int32_t i = 0; i < 52; ++i) {
				auto mi = (li + ri) / 2.0;
				auto m = evalSturm(mi);
				if (m == l) {
					li = mi;
				}
				else {
					ri = mi;
				}
			}
			x = li;
			auto y = clenshaw(x, getDegree(0), (int32_t)lengths[0]);
			return x; // std::abs(y) >= 1e2 * DBL_EPSILON ? firstRoot(std::nextafter(li, DBL_MAX), 1.0) : x;
		}


			std::pair<double,double> SturmSequence::firstAndSecondIntervalRoot(double lli, double rri) {
				double li = lli;
				double ri = rri;
				if (degree == 1) {
					if (getDegree(0)[0] == 0.)
						return std::make_pair(li,ri);
					return  std::make_pair(-2.,2.);
				}
				if (degree == 2) {
					return  std::make_pair(-getDegree(0)[0] / getDegree(0)[1],DBL_MAX);
				}
				auto dfdx = getDegree(1);
				auto df2dx = differentiate(cheb::svec(dfdx, dfdx + lengths[1]));
				auto x = li;

				auto [l, r] = evalSturm(std::make_pair(li, ri));
				if (l == r) return std::make_pair(DBL_MAX, DBL_MAX);
				if (r == l - 1) {
					for (int32_t i = 0; i < 52; ++i) {
						auto mi = (li + ri) / 2.0;
						auto m = evalSturm(mi);
						if (m == l) {
							li = mi;
						}
						else {
							ri = mi;
						}
					}
					x = li;
					auto y = clenshaw(x, getDegree(0), (int32_t)lengths[0]);
					auto r = std::abs(y) >= 1e2 * DBL_EPSILON ? firstRoot(std::nextafter(li, DBL_MAX), 1.0) : x;
					return std::make_pair(r, DBL_MAX);
				}
				else {
					auto l1 = li, r1 = ri;
					for (int32_t i = 0; i < 52; ++i) {
						auto mi = (l1 + r1) / 2.0;
						auto m = evalSturm(mi);
						if (m == l) {
							l1 = mi;
						}
						else {
							r1 = mi;
						}
					}
					double x1 = l1;
					double y1 = clenshaw(x1, getDegree(0), (int32_t)lengths[0]);
					double rt1 = std::abs(y1) >= 1e2 * DBL_EPSILON ? firstRoot(std::nextafter(l1, DBL_MAX), 1.0) : x1;
					x1 = rt1;
					double x2 = 2.0;
					int32_t ii = 0;
					auto l2 = l1;
					do  {
						++ii;
						//auto l2 = li, 
							auto r2 = ri;

						for (int32_t i = 0; i < 52; ++i) {
							auto mi = (l2 + r2) / 2.0;
							auto m = evalSturm(mi);
							if (m < r) return std::make_pair(x1, DBL_MAX);
							if (m >= l - ii) {
								l2 = mi;
							}
							else {
								r2 = mi;
							}
						}
						x2 = l2;
					} while (x1 >= x2);
					double y2 = clenshaw(x2, getDegree(0), (int32_t)lengths[0]);
					double rt2 = std::abs(y2) >= 1e2 * DBL_EPSILON ? firstRoot(std::nextafter(x2, DBL_MAX), 1.0) : x2;
					x2 = rt2;
					return std::make_pair(x1, x2);
				}
			}
}