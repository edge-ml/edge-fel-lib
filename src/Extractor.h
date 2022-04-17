//#pragma once

#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include <string>
#include <vector>
#include <algorithm>
#include "libmfcc.h"
#include "Cmplx.h"

namespace ex {

	class Extractor
	{
	public:

		float mean(std::vector<float>&);
		float mean_abs_dev(std::vector<float>&, float);
		float mean_geometric_abs(std::vector<float>&);
		float median(std::vector<float>);
		float median_abs_changes(std::vector<float>&);
		float median_changes(std::vector<float>&);
		float median_abs_dev(std::vector<float>&, float);
		float std_dev(std::vector<float>&, float);
		float avg_dev(std::vector<float>&, float);
		float var(std::vector<float>&, float);
		float abs_energy(std::vector<float>&);
		float kurtosis(std::vector<float>&, float, float);
		float skewness(std::vector<float>&, float, float);
		float zero_cross(std::vector<float>&);

		float max(std::vector<float>&);
		float abs_max(std::vector<float>&);
		float min(std::vector<float>&);
		float last_location_of_max(std::vector<float>&, float);
		float last_location_of_min(std::vector<float>&, float);
		float first_location_of_max(std::vector<float>&, float);
		float first_location_of_min(std::vector<float>&, float);
		float mean_n_abs_max(std::vector<float>, int);

		float mean_abs_changes(std::vector<float>&);
		float mean_changes(std::vector<float>&);
		float abs_sum_of_changes(std::vector<float>&);
		float change_quantile(std::vector<float>, float, float, int);

		float sum(std::vector<float>&);
		float range_count(std::vector<float>&, float, float);
		float non_zero_count(std::vector<float>&);
		float count_above(std::vector<float>&, float);
		float count_above_mean(std::vector<float>&, float);
		float count_below(std::vector<float>&, float);
		float count_below_mean(std::vector<float>&, float);
		float root_mean_square(std::vector<float>&, float);
		float quantile(std::vector<float>, float);
		float interquartile_range(std::vector<float>&);
		float negative_turnings(std::vector<float>&);
		float positive_turnings(std::vector<float>&);

		float autocorrelation(std::vector<float>&, int, float, float);
		float mfcc(std::vector<co::cd>&, int, int, int);
		std::vector<co::cd> fft(std::vector<float>&);
		std::vector<float> lpc(std::vector<float>&, int);
		std::vector<float> lpcc(std::vector<float>&, int);
		
	
	private:
		unsigned int bitReverse(unsigned int, int);
		float call_by_reference_median(std::vector<float>&);
	};

}
#endif
