#include "Extractor.h"
#include <cmath>
#include <iostream>

using namespace ex;
using namespace std;
using namespace co;

//Returns the mean of the values
float Extractor::mean(vector<float>& values) {
  float sum = 0;
  for (auto& value : values) {
    sum += value;
  }

  return sum / values.size();
}

//Returns the mean of absolute deviation of a value from the mean
float Extractor::mean_abs_dev(vector<float>& values, float my_mean) {
  float sum = 0;
  for (auto& value : values) {
    sum += abs(value - my_mean);
  }

  return sum / values.size();
}

//Returns the geometric mean of the (absolute) values
float Extractor::mean_geometric_abs(vector<float>& values) {
  /*for (auto& value : values) {
  	if (value != 0) {
  		value = log(abs(value));
  	}
    }
    float log_mean = mean(values);
    return exp(log_mean);
  */

  /* for (auto& value : values) {
     if (value != 0) {
        value = pow(abs(value), (1.0f / values.size()));
      }
    }
    float log_mean = mean(values);
    return exp(log_mean);*/


  float product = 1;
  float temp = 1;
  for (int i = 0; i < values.size(); i++) {
    if (values.at(i) != 0) {
      if (i % 30 == 0 && i != 0) {
        product *= pow(temp, (1.0f / values.size()));
        temp = 1;
      }
      temp *= abs(values.at(i));
    } else {
      return 0.0;
    }
  }
  product *= pow(temp, (1.0f / (values.size())));

  return product;
  /*
    float product = 1;
    for (auto& value : values) {
  	if (value != 0) {
  		product = product * abs(value);
  	}
    }
    return pow(product, (1.0f / values.size()));
  */
}

//Helper method for median inside extractors
float Extractor::call_by_reference_median(vector<float>& values) {
  sort(values.begin(), values.end());
  int size = values.size();
  if (size % 2 == 0) {
    int mid = size / 2 - 1;
  
    return (values.at(mid) + values.at(mid + 1)) / 2;
  }

  return values.at(size / 2);
}

//Returns the median of the values
float Extractor::median(vector<float> values) {
  sort(values.begin(), values.end());
  int size = values.size();
  if (size % 2 == 0) {
    int mid = size / 2 - 1;
  
    return (values.at(mid) + values.at(mid + 1)) / 2;
  }

  return values.at(size / 2);
}

//Returns the median of all absolute changes
float Extractor::median_abs_changes(vector<float>& values) {
  vector<float> abs_changes;
  abs_changes.reserve(values.size());
  for (size_t i = 0; i < values.size() - 1; i++) {
    float diff = abs(values.at(i) - values.at(i + 1));
    abs_changes.push_back(diff);
  }

  return call_by_reference_median(abs_changes);
}

//Returns the median of all changes
float Extractor::median_changes(vector<float>& values) {
  vector<float> changes;
  changes.reserve(values.size());
  for (size_t i = 0; i < values.size() - 1; i++) {
    float diff = values.at(i + 1) - values.at(i);
    changes.push_back(diff);
  }

  return call_by_reference_median(changes);
}

//Returns the median of all absolute deviations to the median
float Extractor::median_abs_dev(vector<float>& values, float my_median) {
  vector<float> abs_dev;
  abs_dev.reserve(values.size());
  for (auto& value : values) {
    float dev = abs(value - my_median);
    abs_dev.push_back(dev);
  }

  return call_by_reference_median(abs_dev);
}

//Returns the standard deviation, which is the square root of the variance
float Extractor::std_dev(vector<float>& values, float my_var) {

  return sqrt(my_var);
}

//Returns the average absolute deviation from mean
float Extractor::avg_dev(vector<float>& values, float my_mean) {
  float abs_dev = 0;
  for (auto& value : values) {
    abs_dev += abs(value - my_mean);
  }

  return abs_dev / values.size();
}

//Returns the average quadratic deviation from mean
float Extractor::var(vector<float>& values, float my_mean) {
  float q_dev = 0;
  for (auto& value : values) {
    q_dev += (value - my_mean) * (value - my_mean);
    //q_dev += pow(value - my_mean, 2);
  }

  return q_dev / values.size();
}

//Returns the sum of the squared values
float Extractor::abs_energy(vector<float>& values) {
  float sum = 0;
  for (auto& value : values) {
    sum += (value * value);
    //sum += pow(value, 2);
  }

  return sum;
}

//Returns the mean of standardised values to the power of 4
float Extractor::kurtosis(vector<float>& values, float mean, float std_dev) {
  float sum = 0;
  for (auto& value : values) {
    float norm = (value - mean) / std_dev;
    sum += (norm * norm * norm * norm);
    //sum += pow(norm,4);
  }

  return sum / values.size();
}

//Returns the mean of standardised values to the power of 3
float Extractor::skewness(vector<float>& values, float mean, float std_dev) {
  float sum = 0;
  for (auto& value : values) {
    float norm = (value - mean) / std_dev;
    sum += (norm * norm * norm);
    //sum += pow(norm,3);
  }

  return sum / values.size();
}

//Returns the percentage of zero crossings of two consecutive values
float Extractor::zero_cross(vector<float>& values) {
  float crosses = 0;
  for (size_t i = 0; i < values.size() - 1; i++) {
    if (values.at(i) * values.at(i + 1) < 0) {
      crosses++;
    }
  }

  return crosses / values.size();
}

//Returns the maximum of all values
float Extractor::max(vector<float>& values) {
  float max = values.at(0);
  for (auto& value : values) {
    if (value > max) {
      max = value;
    }
  }

  return max;
}

//Returns the absolute maximum of all values
float Extractor::abs_max(vector<float>& values) {
  float max = abs(values.at(0));
  for (auto& value : values) {
    if (abs(value) > max) {
      max = abs(value);
    }
  }

  return max;
}

//Returns the minimum of all values
float Extractor::min(vector<float>& values) {
  float min = values.at(0);
  for (auto& value : values) {
    if (value < min) {
      min = value;
    }
  }

  return min;
}

//Returns the last location of the maximum value
float Extractor::last_location_of_max(vector<float>& values, float max) {
  size_t index = values.size() - 1;
  for (int i = values.size() - 2; i >= 0; i--) {
    if (values.at(i) == max) {
      index = i;
    }
  }

  return index;
}

//Returns the last location of the minimum value
float Extractor::last_location_of_min(vector<float>& values, float min) {
  size_t index = values.size() - 1;
  for (int i = values.size() - 2; i >= 0; i--) {
    if (values.at(i) == min) {
      index = i;
    }
  }

  return index;
}

//Returns the first location of the maximum value
float Extractor::first_location_of_max(vector<float>& values, float max) {
  float index = 0;
  for (int i = 1; i < values.size(); i++) {
    if (abs(values.at(i) - max) < numeric_limits<float>::epsilon()) {
      return i;
    }
  }

  return index;
}

//Returns the first location of the minimum value
float Extractor::first_location_of_min(vector<float>& values, float min) {
  float index = 0;
  for (int i = 1; i < values.size(); i++) {
    if (abs(values.at(i) - min) < numeric_limits<float>::epsilon()) {
      return i;
    }
  }

  return index;
}

//Returns the average of the n largest absolute values
float Extractor::mean_n_abs_max(vector<float> values, int n) {
  if (n > values.size()) {
    return 0.0;
  }
  for (auto& value : values) {
    value = abs(value);
  }
  sort(values.begin(), values.end());
  float sum = 0;
  for (int i = values.size() - 1; i > values.size() - 1 - n; i--) {
    sum += abs(values.at(i));
  }

  return sum / n;
}


//Returns the mean of the absolute differences of consecutive values
float Extractor::mean_abs_changes(vector<float>& values) {
  float sum_changes = 0.0;
  for (size_t i = 0; i < values.size() - 1; i++) {
    float diff = abs(values.at(i) - values.at(i + 1));
    sum_changes += diff;
  }

  return sum_changes / values.size();
}

//Returns the mean of differences of consecutive values
float Extractor::mean_changes(vector<float>& values) {
  float sum_changes = 0.0;
  for (size_t i = 0; i < values.size() - 1; i++) {
    float diff = values.at(i + 1) - values.at(i);
    sum_changes += diff;
  }

  return sum_changes / values.size();
}

//Returns the sum of the absolute differences of consecutive values
float Extractor::abs_sum_of_changes(vector<float>& values) {
  float sum = 0.0;
  for (size_t i = 0; i < values.size() - 1; i++) {
    float diff = abs(values.at(i) - values.at(i + 1));
    sum += diff;
  }

  return sum;
}

//Calculates the absolute differences of all values between the two quantiles and applies the aggregation function
float Extractor::change_quantile(vector<float> values, float lower, float upper, int aggr) {
  sort(values.begin(), values.end());

  for (size_t i = lower; i <= upper; i++) {
    float diff = abs(values.at(i) - values.at(i + 1));
    values.at(i) = diff;
  }
  if (aggr == 0) {
    return sum(values);
  }
  else if (aggr == 1) {
    return mean(values);
  }
  else if (aggr == 2) {
    return call_by_reference_median(values);
  }
  else if (aggr == 3) {
    return var(values, mean(values));
  }
  else if (aggr == 4) {
    return std_dev(values, var(values, mean(values)));
  }
  else {
    return 0;
  }
}

//Calculates the sum of all values
float Extractor::sum(vector<float>& values) {
  float sum = 0;
  for (auto value : values) {
    sum += value;
  }

  return sum;
}

//Calculates the sum of all values with a value between lower and upper
float Extractor::range_count(vector<float>& values, float lower, float upper) {
  float sum = 0;
  for (auto& value : values) {
    if (value > lower && value < upper) {
      sum += value;
    }
  }

  return sum;
}

//Returns the number of non-zero values
float Extractor::non_zero_count(vector<float>& values) {
  float amount = 0;
  for (auto& value : values) {
    if (value != 0) {
      amount++;
    }
  }

  return amount;
}

//Returns the percentage of values greater than x
float Extractor::count_above(vector<float>& values, float x) {
  float amount = 0;
  for (auto& value : values) {
    if (value > x) {
      amount++;
    }
  }

  return amount / values.size();
}

//Returns the percentage of values greater than mean
float Extractor::count_above_mean(vector<float>& values, float mean) {
  return count_above(values, mean);
}

//Returns the percentage of values lower than x
float Extractor::count_below(vector<float>& values, float x) {
  float amount = 0;
  for (auto& value : values) {
    if (value < x) {
      amount++;
    }
  }

  return amount / values.size();
}

//Returns the percentage of values lower than mean
float Extractor::count_below_mean(vector<float>& values, float mean) {
  return count_below(values, mean);
}

//Returns the root of the mean of all squares values
float Extractor::root_mean_square(vector<float>& values, float energy) {

  return sqrt(energy / values.size());
}

//Returns the value which is greater than q*n percent of all values
float Extractor::quantile(vector<float> values, float q) {
  sort(values.begin(), values.end());

  float nq = q * values.size();
  int i = floor(nq);
  if (abs(nq - i) > 0.00000001) {
    i = ceil(nq);
  
    return values.at(i - 1);
  }
  else {
  
    return (values.at(i - 1) + values.at(i)) / 2;
  }
}

//Returns the difference of the 3/4 and 1/4 quantiles
float Extractor::interquartile_range(vector<float>& values) {
  float upper = quantile(values, 0.75);
  float lower = quantile(values, 0.25);

  return upper - lower;
}

//Calculates the amount of local minimums
float Extractor::negative_turnings(vector<float>& values) {
  float amount = 0;
  for (int i = 0; i < values.size() - 2; i++) {
    if (values.at(i) > values.at(i + 1) && values.at(i + 1) < values.at(i + 2)) {
      amount++;
    }
  }

  return amount;
}

//Calculates the amount of local maximas
float Extractor::positive_turnings(vector<float>& values) {
  float amount = 0;
  for (int i = 0; i < values.size() - 2; i++) {
    if (values.at(i) < values.at(i + 1) && values.at(i + 1) > values.at(i + 2)) {
      amount++;
    }
  }

  return amount;
}

//Calculates a mfcc coefficient, implementation from
//https://github.com/jsawruk/libmfcc
//Copyright (c) 2010 Jeremy Sawruk
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
//to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
float Extractor::mfcc(vector<cd>& values, int samplingRate, int numFilters, int m) {
  return getCoefficient(values, samplingRate, numFilters, values.size(), m);
}

//Calculates an estimation of the autocorrelation for a specified lag
float Extractor::autocorrelation(vector<float>& values, int lag, float mean, float var) {
  float sum = 0;
  for (size_t i = 0; i < values.size() - lag; i++) {
    sum = (values.at(i) - mean) * (values.at(i + lag) - mean);
  }
  float corr = sum / ((values.size() - (float) lag) * var * var);

  return corr;
}


//Interative cooley-turkey fft algorithm, adapted from
//https://www.geeksforgeeks.org/iterative-fast-fourier-transformation-polynomial-multiplication/
vector<cd> Extractor::fft(std::vector<float>& values) {

  int n = values.size();
  float power = log2(values.size());
  int power_round = ceil(power);

  //Pad vector if necessary
  if (abs(power - power_round) > 0.0001) {
    n = pow(2, power_round);
  }
  vector<cd> res;
  res.reserve(n);

  for (int i = 0; i < n; i++) {
    cd c(0, 0);
    res.push_back(c);
  }

  // bit reversal of the given array
  for (unsigned int i = 0; i < n; i++) {
    int rev = bitReverse(i, power_round);
    cd c;
    if (rev < values.size()) {
      c.real = values[rev];
    }
    res.at(i) = c;
  }


  // j is iota
  const cd J(0, 1);
  for (int s = 1; s <= power_round; ++s) {
    int m = 1 << s; // 2 power s
    int m2 = m >> 1; // m2 = m/2 -1
    cd w(1, 0);

    // principle root of nth complex
    // root of unity.
    cd wm = e(mul(J, pi / m2));
    for (int j = 0; j < m2; ++j) {
      for (int k = j; k < n; k += m) {
        // t = twiddle factor
        cd t = mul(w, res[k + m2]);
        cd u = res[k];

        // similar calculating y[k]
        res[k] = add(u, t);

        // similar calculating y[k+n/2]
        res[k + m2] = sub(u, t);
      }
      w = mul(w, wm);
    }
  }

  return res;
}

// Utility function for reversing the bits of given index x
unsigned int Extractor::bitReverse(unsigned int x, int log2n) {
  int n = 0;
  for (int i = 0; i < log2n; i++)
  {
    n <<= 1;
    n |= (x & 1);
    x >>= 1;
  }
  return n;
}

//Calculates n lpc coefficients, implementation from
//https://github.com/jamiebullock/LibXtract
//Copyright (C) 2012 Jamie Bullock
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
//to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
vector<float> Extractor::lpc(vector<float>& autoc, int n) {
  vector<float> lpc(n - 1, 0);
  float error = autoc[0];
  if (error == 0.0) {
    return lpc;
  }

  for (int i = 0; i < n - 1; i++) {
    /* Sum up this iteration's reflection coefficient. */
    float r = -autoc[i + 1];
    for (int j = 0; j < i; j++) {
      r -= lpc[j] * autoc[i - j];
    }
    r /= error;

    /* Update LPC coefficients and total error. */
    lpc[i] = r;
    int j;
    for (j = 0; j < i / 2; j++)
    {
      float tmp = lpc[j];
      lpc[j] = r * lpc[i - 1 - j];
      lpc[i - 1 - j] += r * tmp;
    }
    if (i % 2) {
      lpc[j] += lpc[j] * r;
    }
    error *= 1 - r * r;
  }

  return lpc;
}

//Calculates a cepstrum length of lpc coefficients, implementation from
//https://github.com/jamiebullock/LibXtract
//Copyright (C) 2012 Jamie Bullock
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
//to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

vector<float> Extractor::lpcc(vector<float>& lpc_coeffs, int cep_length) {
  int order = lpc_coeffs.size() - 1; /* Eventually change this to Q = 3/2 p as suggested in Rabiner */
  vector<float> lpcc(cep_length, 0);

  for (int n = 1; n <= order && n <= cep_length; n++) {
    float sum = 0.0;
    for (int k = 1; k < n; k++) {
      sum += k * lpcc[k - 1] * lpc_coeffs[n - k];
    }
    lpcc[n - 1] = lpc_coeffs[n] + sum / n;
  }

  /* be wary of these interpolated values */
  for (int n = order + 1; n <= cep_length; n++) {
    float sum = 0.0;
    for (int k = n - (order - 1); k < n; k++)
      sum += k * lpcc[k - 1] * lpc_coeffs[n - k];
    lpcc[n - 1] = sum / n;
  }

  return lpcc;
}
