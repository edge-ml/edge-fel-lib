//#pragma once

#ifndef EXTRACTIONHANDLER_H
#define EXTRACTIONHANDLER_H

#include "Extractor.h"
#include "ExtractionDelegate.h"
#include <cmath>

namespace eh {

	class ExtractionHandler
	{

	public:
		
		float handle_mean(std::string, std::vector<float>&);
		float handle_mean_abs_dev(std::string, std::vector<float>&);
		float handle_mean_geometric_abs(std::string, std::vector<float>&);
		float handle_median(std::string, std::vector<float>&);
		float handle_median_abs_changes(std::string, std::vector<float>&);
		float handle_median_changes(std::string, std::vector<float>&);
		float handle_median_abs_dev(std::string, std::vector<float>&);
		float handle_std_dev(std::string, std::vector<float>&);
		float handle_avg_dev(std::string, std::vector<float>&);
		float handle_var(std::string, std::vector<float>&);
		float handle_abs_energy(std::string, std::vector<float>&);
		float handle_kurtosis(std::string, std::vector<float>&);
		float handle_skewness(std::string, std::vector<float>&);
		float handle_zero_cross(std::string, std::vector<float>&);
		
		float handle_max(std::string, std::vector<float>&);
		float handle_abs_max(std::string, std::vector<float>&);
		float handle_min(std::string, std::vector<float>&);
		float handle_last_location_of_max(std::string, std::vector<float>&);
		float handle_last_location_of_min(std::string, std::vector<float>&);
		float handle_first_location_of_max(std::string, std::vector<float>&);
		float handle_first_location_of_min(std::string, std::vector<float>&);
		float handle_mean_n_abs_max(std::string, std::vector<float>&, float);

		float handle_mean_abs_changes(std::string, std::vector<float>&);
		float handle_mean_changes(std::string, std::vector<float>&);
		float handle_abs_sum_of_changes(std::string, std::vector<float>&);
		float handle_change_quantile(std::string, std::vector<float>&, float, float, float);

		float handle_sum(std::string, std::vector<float>&);
		float handle_range_count(std::string, std::vector<float>&, float, float);
		float handle_non_zero_count(std::string, std::vector<float>&);
		float handle_count_above(std::string, std::vector<float>&, float);
		float handle_count_above_mean(std::string, std::vector<float>&);
		float handle_count_below(std::string, std::vector<float>&, float);
		float handle_count_below_mean(std::string, std::vector<float>&);
		float handle_root_mean_square(std::string, std::vector<float>&);
		float handle_quantile(std::string, std::vector<float>&, float);
		float handle_interquartile_range(std::string, std::vector<float>&);
		float handle_negative_turnings(std::string, std::vector<float>&);
		float handle_positive_turnings(std::string, std::vector<float>&);

		float handle_autocorrelation(std::string, std::vector<float>&, float);
		std::vector<float> handle_mfcc(std::string, std::vector<float>&, float, float, float);
		std::vector<co::cd> handle_fft(std::string, std::vector<float>&);
		std::vector<float> handle_lpc(std::string, std::vector<float>&, float, float);
		std::vector<float> handle_lpcc(std::string, std::vector<float>&, float, float, float);

	};

}

#endif
