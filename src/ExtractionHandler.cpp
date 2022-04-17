#include "ExtractionHandler.h"

using namespace ed;
using namespace ex;
using namespace eh;
using namespace std;
using namespace co;

Extractor extractor;

//Handler function for mean, caches value
float ExtractionHandler::handle_mean(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}

	float value = extractor.mean(values);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for mean_abs_dev
float ExtractionHandler::handle_mean_abs_dev(string feature, vector<float>& values) {
	float mean = handle_mean("mean", values);
	return extractor.mean_abs_dev(values, mean);
}

//Handler function for mean_geometric_abs
float ExtractionHandler::handle_mean_geometric_abs(string feature, vector<float>& values) {
	return extractor.mean_geometric_abs(values);
}

//Handler function for median, caches value
float ExtractionHandler::handle_median(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}

	float value = extractor.median(values);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for median_abs_diff
float ExtractionHandler::handle_median_abs_changes(string feature, vector<float>& values) {
	return extractor.median_abs_changes(values);
}

//Handler function for median_diff
float ExtractionHandler::handle_median_changes(string feature, vector<float>& values) {
	return extractor.median_changes(values);
}

//Handler function for median_abs_dev, calls handle_median
float ExtractionHandler::handle_median_abs_dev(string feature, vector<float>& values) {
	float median = handle_median("median", values);
	return extractor.median_abs_dev(values, median);
}

//Handler function for std_dev, calls handle_var, caches value
float ExtractionHandler::handle_std_dev(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}

	float var = handle_var("var", values);
	float value = extractor.std_dev(values, var);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for avg_dev, calls handle_mean
float ExtractionHandler::handle_avg_dev(string feature, vector<float>& values) {
	float mean = handle_mean("mean", values);
	return extractor.avg_dev(values, mean);
}

//Handler function for var, calls handle_mean, caches value
float ExtractionHandler::handle_var(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}

	float mean = handle_mean("mean", values);
	float value = extractor.var(values, mean);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for abs_energy
float ExtractionHandler::handle_abs_energy(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}
	float value = extractor.abs_energy(values);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for kurtosis, calls handle_mean and handle_std_dev
float ExtractionHandler::handle_kurtosis(string feature, vector<float>& values) {
	float mean = handle_mean("mean", values);
	float std_dev = handle_std_dev("std_dev", values);
	return extractor.kurtosis(values, mean, std_dev);
}

//Handler function for skewness, calls handle_mean and handle_std_dev
float ExtractionHandler::handle_skewness(string feature, vector<float>& values) {
	float mean = handle_mean("mean", values);
	float std_dev = handle_std_dev("std_dev", values);
	return extractor.skewness(values, mean, std_dev);
}

//Handler function for zero_cross
float ExtractionHandler::handle_zero_cross(string feature, vector<float>& values) {
	return extractor.zero_cross(values);
}

//Handler function for max, caches value
float ExtractionHandler::handle_max(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}

	float value = extractor.max(values);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for abs_max
float ExtractionHandler::handle_abs_max(string feature, vector<float>& values) {
	return extractor.abs_max(values);
}

//Handler function for min, caches value
float ExtractionHandler::handle_min(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}

	float value = extractor.min(values);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for last location of maximum, calls handle_max
float ExtractionHandler::handle_last_location_of_max(string feature, vector<float>& values) {
	float max = handle_max("max", values);
	return extractor.last_location_of_max(values, max);
}

//Handler function for last location of minimum, calls handle_min
float ExtractionHandler::handle_last_location_of_min(string feature, vector<float>& values) {
	float min = handle_min("min", values);
	return extractor.last_location_of_min(values, min);
}

//Handler function for first location of maximum, calls handle_max
float ExtractionHandler::handle_first_location_of_max(string feature, vector<float>& values) {
	float max = handle_max("max", values);
	return extractor.first_location_of_max(values, max);
}

//Handler function for first location of minimum, calls handle_min
float ExtractionHandler::handle_first_location_of_min(string feature, vector<float>& values) {
	float min = handle_min("min", values);
	return extractor.first_location_of_min(values, min);
}

//Handler function for mean_n_abs_max, takes n as input
float ExtractionHandler::handle_mean_n_abs_max(string feature, vector<float>& values, float n) {
	return extractor.mean_n_abs_max(values, n);
}

//Handler function for mean_abs_change
float ExtractionHandler::handle_mean_abs_changes(string feature, vector<float>& values) {
	return extractor.mean_abs_changes(values);
}

//Handler function for mean_change
float ExtractionHandler::handle_mean_changes(string feature, vector<float>& values) {
	return extractor.mean_changes(values);
}

//Handler function for abs_sum_of_changes
float ExtractionHandler::handle_abs_sum_of_changes(string feature, vector<float>& values) {
	return extractor.abs_sum_of_changes(values);
}


//Handler function for change_quantiles, takes lower and upper quantile and a aggregation function as input
float ExtractionHandler::handle_change_quantile(string feature, vector<float>& values, float lower, float upper, float aggr) {
	return extractor.change_quantile(values, lower, upper, aggr);
}

//Handler function for sum, caches value
float ExtractionHandler::handle_sum(string feature, vector<float>& values) {
	if (ExtractionDelegate::doCache && ExtractionDelegate::calculated.count(feature)) {
		return ExtractionDelegate::calculated.at(feature);
	}

	float value = extractor.sum(values);
	ExtractionDelegate::checkAndInsert(feature, value);
	return value;
}

//Handler function for range_count, takes lower and upper bound for values as input
float ExtractionHandler::handle_range_count(string feature, vector<float>& values, float lower, float upper) {
	return extractor.range_count(values, lower, upper);
}

//Handler function for non_zero_count
float ExtractionHandler::handle_non_zero_count(string feature, vector<float>& values) {
	return extractor.non_zero_count(values);
}

//Handler function for count_above, takes lower bound x as input
float ExtractionHandler::handle_count_above(string feature, vector<float>& values, float x) {
	return extractor.count_above(values, x);
}

//Handler function for count_above_mean, calls handle_mean
float ExtractionHandler::handle_count_above_mean(string feature, vector<float>& values) {
	float mean = handle_mean("mean", values);
	return extractor.count_above(values, mean);
}

//Handler function for count_below, takes upper bound x as input
float ExtractionHandler::handle_count_below(string feature, vector<float>& values, float x) {
	return extractor.count_below(values, x);
}

//Handler function for count_below_mean, calls handle_mean
float ExtractionHandler::handle_count_below_mean(string feature, vector<float>& values) {
	float mean = handle_mean("mean", values);
	return extractor.count_below_mean(values, mean);
}

//Handler function for root_mean_square
float ExtractionHandler::handle_root_mean_square(string feature, vector<float>& values) {
	float energy = handle_abs_energy("abs_energy", values);
	return extractor.root_mean_square(values, energy);
}

//Handler function for quantiles, takes the quantile 0 <= q <= 1 as parameter
float ExtractionHandler::handle_quantile(string feature, vector<float>& values, float q) {
	return extractor.quantile(values, q);
}

//Handler function for interquartile_range
float ExtractionHandler::handle_interquartile_range(string feature, vector<float>& values) {
	return extractor.interquartile_range(values);
}

//Handler function for negative_turnings
float ExtractionHandler::handle_negative_turnings(string feature, vector<float>& values) {
	return extractor.negative_turnings(values);
}

//Handler function for positive_turnings
float ExtractionHandler::handle_positive_turnings(string feature, vector<float>& values) {
	return extractor.positive_turnings(values);
}

//Handler function for autocorrelation, takes the lag as input and calls handle_mean and handle_var
float ExtractionHandler::handle_autocorrelation(string feature, vector<float>& values, float lag) {
	float mean = handle_mean("mean", values);
	float var = handle_var("var", values);
	return extractor.autocorrelation(values, lag, mean, var);
}

//Handler function for fft, returns a vector of imaginary floats
vector<cd> ExtractionHandler::handle_fft(string feature, vector<float>& values) {
	return extractor.fft(values);
}

//Handler function for mfcc, takes samplingRate, numFilters and the coefficient number as parameter, calls handle_fft and transform result to real numbers
vector<float> ExtractionHandler::handle_mfcc(string feature, vector<float>& values, float sampling_rate, float num_filter, float m) {
	vector<cd> spectralData = handle_fft("fft", values);
	vector<float> coeffs;
	coeffs.reserve(m);
	for (int i = 1; i <= m; i++) {
		coeffs.push_back(extractor.mfcc(spectralData, sampling_rate, num_filter, i));
	}

	return coeffs;
}

//Handler function for lpc, creates a vector of n autocorrelations, takes amount of autocorrelations n as param
vector<float> ExtractionHandler::handle_lpc(string feature, vector<float>& values, float auto_n, float lpc_n) {
	int m = auto_n;
	vector<float> autoc;
	autoc.reserve(auto_n);

	while (m--)
	{
		float corr = 0;
		for (int i = 0; i < auto_n - m; i++)
		{
			corr += values[i] * values[i + m];
		}
		autoc.insert(autoc.begin(), corr/ auto_n);
	}
	return extractor.lpc(autoc, lpc_n);
}


//Handler function for lpcc, calls lcp, takes length of cepstrum as additional param
vector<float> ExtractionHandler::handle_lpcc(string feature, vector<float>& values, float auto_n, float lpc_n, float cep_length) {
	vector<float> lpc_coeffs = handle_lpc("lpc", values, auto_n, lpc_n);
	return extractor.lpcc(lpc_coeffs, cep_length);
}
