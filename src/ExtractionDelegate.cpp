#include "ExtractionDelegate.h"
#include <iostream>

using namespace ed;
using namespace eh;
using namespace std;
using namespace co;

bool ExtractionDelegate::doCache = false;
std::map<string, float> ExtractionDelegate::calculated;
ExtractionHandler handler;
vector<string> ExtractionDelegate::parameterFeatures = { "mean_n_abs_max",  "change_quantile", "range_count", "count_above", "count_below", "quantile", "autocorrelation" };

//Extracts the requested feature from the data, doesnt include features which return a list
float ExtractionDelegate::extractOne(string feature, vector<float>& values, std::map<string, float>& params) {
	auto iter = ExtractionDelegate::handlers.find(feature);
	if (!(iter == ExtractionDelegate::handlers.end())) {
		return (handler.*(iter->second))(feature, values);
	}
	else if (feature == "mean_n_abs_max" && params.find("mean_n_abs_max_n") != params.end()) {
		return handler.handle_mean_n_abs_max(feature, values, params.at("mean_n_abs_max_n"));
	}
	else if (feature == "change_quantile" && params.find("change_quantile_lower") != params.end() && params.find("change_quantile_upper") != params.end() && params.find("change_quantile_aggr") != params.end()) {
		return handler.handle_change_quantile(feature, values, params.at("change_quantile_lower"), params.at("change_quantile_upper"), params.at("change_quantile_aggr"));
	}
	else if (feature == "range_count" && params.find("range_count_lower") != params.end() && params.find("range_count_upper") != params.end()) {
		return handler.handle_range_count(feature, values, params.at("range_count_lower"), params.at("range_count_upper"));
	}
	else if (feature == "count_above" && params.find("count_above_x") != params.end()) {
		return handler.handle_count_above(feature, values, params.at("count_above_x"));
	}
	else if (feature == "count_below" && params.find("count_below_x") != params.end()) {
		return handler.handle_count_below(feature, values, params.at("count_below_x"));
	}
	else if (feature == "quantile" && params.find("quantile_q") != params.end()) {
		return handler.handle_quantile(feature, values, params.at("quantile_q"));
	}
	else if (feature == "autocorrelation" && params.find("autocorrelation_lag") != params.end()) {
		return handler.handle_autocorrelation(feature, values, params.at("autocorrelation_lag"));
	}
	else {
		return 0.0;
	}
}


//Extracts the requested feature from the data which return a list of floats. Fft is not included here!
vector<float> ExtractionDelegate::extractOneVectorial(string feature, vector<float>& values, std::map<string, float>& params) {
	if (feature == "mfcc" && params.find("mfcc_sampling_rate") != params.end() && params.find("mfcc_num_filter") != params.end() && params.find("mfcc_m") != params.end()) {
		return handler.handle_mfcc(feature, values, params.at("mfcc_sampling_rate"), params.at("mfcc_num_filter"), params.at("mfcc_m"));
	}
	else if (feature == "lpc" && params.find("lpc_auto_n") != params.end() && params.find("lpc_n") != params.end()) {
		return handler.handle_lpc(feature, values, params.at("lpc_auto_n"), params.at("lpc_n"));
	}
	else if (feature == "lpcc" && params.find("lpcc_auto_n") != params.end() && params.find("lpcc_n") != params.end() && params.find("lpcc_cep_length") != params.end()) {
		return handler.handle_lpcc(feature, values, params.at("lpcc_auto_n"), params.at("lpcc_n"), params.at("lpcc_cep_length"));
	}
	else {
		return {};
	}
}

//Extracts the requested features from the data, doesnt include features which return a list
std::map<string, float> ExtractionDelegate::extractSome(vector<string>& features, vector<float>& values, std::map<string, float>& params) {
	std::map<string, float> results;
	for (string feature : features) {
		float value = extractOne(feature, values, params);
		results.emplace(feature, value);
	}

	return results;
}

//Extracts all available features from the data, doesnt include features which return a list
std::map<string, float> ExtractionDelegate::extractAll(vector<float>& values, std::map<string, float>& params) {
	std::map<string, float> results;
	for (auto iter : ExtractionDelegate::handlers) {
		float value = (handler.*(iter.second))(iter.first, values);
		results.emplace(iter.first, value);
	}
	std::map<string, float> parameterFeatures = extractSome(ExtractionDelegate::parameterFeatures, values, params);
	results.insert(parameterFeatures.begin(), parameterFeatures.end());
	return results;
}

//Extracts the spectrum of the times series with fft
vector<cd> ExtractionDelegate::extractSpectrum(vector<float>& values) {
	return handler.handle_fft("fft", values);
}

//Helper method for fail-fast cache check and insertion
void ExtractionDelegate::checkAndInsert(string feature, float value) {
	if (doCache && !calculated.count(feature)) {
		calculated.emplace(feature, value);
	}
}

//Initializes the handler map 
void ExtractionDelegate::fillHandlerMap() {
	ExtractionDelegate::handlers.emplace("mean", &ExtractionHandler::handle_mean);
	ExtractionDelegate::handlers.emplace("mean_abs_dev", &ExtractionHandler::handle_mean_abs_dev);
	ExtractionDelegate::handlers.emplace("mean_geometric_abs", &ExtractionHandler::handle_mean_geometric_abs);
	ExtractionDelegate::handlers.emplace("median", &ExtractionHandler::handle_median);
	ExtractionDelegate::handlers.emplace("median_abs_changes", &ExtractionHandler::handle_median_abs_changes);
	ExtractionDelegate::handlers.emplace("median_changes", &ExtractionHandler::handle_median_changes);
	ExtractionDelegate::handlers.emplace("median_abs_dev", &ExtractionHandler::handle_median_abs_dev);
	ExtractionDelegate::handlers.emplace("std_dev", &ExtractionHandler::handle_std_dev);
	ExtractionDelegate::handlers.emplace("avg_dev", &ExtractionHandler::handle_avg_dev);
	ExtractionDelegate::handlers.emplace("var", &ExtractionHandler::handle_var);
	ExtractionDelegate::handlers.emplace("abs_energy", &ExtractionHandler::handle_abs_energy);
	ExtractionDelegate::handlers.emplace("kurtosis", &ExtractionHandler::handle_kurtosis);
	ExtractionDelegate::handlers.emplace("skewness", &ExtractionHandler::handle_skewness);
	ExtractionDelegate::handlers.emplace("zero_cross", &ExtractionHandler::handle_zero_cross);

	ExtractionDelegate::handlers.emplace("max", &ExtractionHandler::handle_max);
	ExtractionDelegate::handlers.emplace("abs_max", &ExtractionHandler::handle_abs_max);
	ExtractionDelegate::handlers.emplace("min", &ExtractionHandler::handle_min);
	ExtractionDelegate::handlers.emplace("last_location_of_max", &ExtractionHandler::handle_last_location_of_max);
	ExtractionDelegate::handlers.emplace("last_location_of_min", &ExtractionHandler::handle_last_location_of_min);
	ExtractionDelegate::handlers.emplace("first_location_of_max", &ExtractionHandler::handle_first_location_of_max);
	ExtractionDelegate::handlers.emplace("first_location_of_min", &ExtractionHandler::handle_first_location_of_min);

	ExtractionDelegate::handlers.emplace("mean_abs_changes", &ExtractionHandler::handle_mean_abs_changes);
	ExtractionDelegate::handlers.emplace("mean_changes", &ExtractionHandler::handle_mean_changes);
	ExtractionDelegate::handlers.emplace("abs_sum_of_changes", &ExtractionHandler::handle_abs_sum_of_changes);

	ExtractionDelegate::handlers.emplace("sum", &ExtractionHandler::handle_sum);
	ExtractionDelegate::handlers.emplace("non_zero_count", &ExtractionHandler::handle_non_zero_count);
	ExtractionDelegate::handlers.emplace("count_above_mean", &ExtractionHandler::handle_count_above_mean);
	ExtractionDelegate::handlers.emplace("count_below_mean", &ExtractionHandler::handle_count_below_mean);
	ExtractionDelegate::handlers.emplace("root_mean_square", &ExtractionHandler::handle_root_mean_square);
	ExtractionDelegate::handlers.emplace("interquartile_range", &ExtractionHandler::handle_interquartile_range);
	ExtractionDelegate::handlers.emplace("negative_turnings", &ExtractionHandler::handle_negative_turnings);
	ExtractionDelegate::handlers.emplace("positive_turnings", &ExtractionHandler::handle_positive_turnings);

}
