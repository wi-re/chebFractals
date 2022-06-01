#include <cheb/cheb.h>
#include <cheb/utilities.h>
#include <cheb/Function.h>



namespace cheb {
	svec internal::merge_duplicates(const svec& arr, const svec& tols) {
		svec out;
		for (int32_t i = 0; i < arr.size() - 1; ++i)
			if (abs(arr[i] - arr[i + 1]) > tols[i])
				out.push_back(arr[i]);
		out.push_back(*(std::end(arr) - 1));
		return (out);
	}
	svec internal::merge_duplicates(const svec& arr, scalar tols) {
		return merge_duplicates(arr, svec(tols, arr.size()));
	}
	std::vector<int32_t> internal::sortindex(const std::vector<Interval>& intervals) {
		auto subintervals = intervals;
		svec leftbreakpts;
		std::vector<int32_t> indices;
		std::vector<std::pair<int32_t, scalar>> idx;
		int32_t i = 0;
		for (auto s : subintervals) {
			leftbreakpts.push_back(s.a);
			idx.push_back(std::make_pair(i++, s.a));
		}
		std::sort(idx.begin(), idx.end(), [](const auto& lhs, const auto& rhs) {return lhs.second < rhs.second; });
		int32_t j = 0;
		for (auto i : idx) {
			subintervals[j++] = intervals[i.first];
			indices.push_back(i.first);
		}
		for (int32_t k = 0; k < subintervals.size() - 1; ++k) {
			if ((subintervals[k + 1].a - subintervals[k].b) < 0) {
				throw std::invalid_argument("Interval Overlap");
			}
			if ((subintervals[k + 1].a - subintervals[k].b) > 0) {
				throw std::invalid_argument("IntervalGap");
			}
		}
		return indices;
	}
	std::vector<IntervalFunction> check_funs(const std::vector<IntervalFunction>& funs) {
		if (funs.size() == 0)
			return funs;
		std::vector<Interval> intervals;
		for (const auto& f : funs)
			intervals.push_back(f.interval());
		auto idx = internal::sortindex(intervals);
		std::vector<IntervalFunction> sortedfuns(funs.size());
		for (int32_t i = 0; i < funs.size(); ++i)
			sortedfuns[i] = funs[idx[i]];
		return sortedfuns;
	}
	std::pair<svec, svec> compute_breakdata(const std::vector<IntervalFunction>& funs) {
		if (funs.size() == 0)
			return std::pair<svec, svec>{};
		svec points, values;
		for (const auto& fn : funs) {
			auto h = fn.support();
			auto f = fn.endvalues();
			points.push_back(h[0]);
			points.push_back(h[1]);
			values.push_back(f[0]);
			values.push_back(f[1]);
		}
		auto xl = *std::begin(points);
		auto xr = *std::rbegin(points);
		auto yl = *std::begin(values);
		auto yr = *std::rbegin(values);
		auto n = (points.size() - 2);
		auto no = n / 2 + 2;
		svec x(no);
		svec y(no);
		*std::begin(x) = xl;
		*std::rbegin(x) = xr;
		*std::begin(y) = yl;
		*std::rbegin(y) = yr;
		if (n > 0) {
			auto xx = svec(std::begin(points) + 1, std::end(points) - 1);
			auto yy = svec(std::begin(values) + 1, std::end(values) - 1);
			for (int32_t i = 0; i < xx.size(); i += 2) {
				x[i / 2 + 1] = .5 * (xx[i] + xx[i + 1]);
				y[i / 2 + 1] = .5 * (yy[i] + yy[i + 1]);
			}
		}
		return std::make_pair((x), (y));
	}
}