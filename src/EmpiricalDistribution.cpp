/*
Date:23 April 2018
Note:This implementation is from Sailfish. We modified few parts for our purposes.
*/
/**
 *  This implementation for keeping an empirical distribution and
 *  querying the PDF and CDF is taken (and ever-so-slightly-modified) from
 *  https://raw.githubusercontent.com/dcjones/isolator/master/src/emp_dist.cpp.
 *  The original author is Daniel C. Jones; NOT me (Rob Patro).
 **/

#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>
#include <limits>
#include <random>
//for Circall
 #include <cmath>
 //for Circall

#include "EmpiricalDistribution.hpp"

EmpiricalDistribution::EmpiricalDistribution(EmpiricalDistribution& other)
    : pdfvals(other.pdfvals) , cdfvals(other.cdfvals) , med(other.med),
      minVal(other.minVal), maxVal(other.maxVal) { isValid_.store(other.isValid_.load()); }


EmpiricalDistribution::EmpiricalDistribution() {}

EmpiricalDistribution::EmpiricalDistribution(
        const std::vector<uint32_t>& vals,
        const std::vector<uint32_t>& lens) {
    buildDistribution(vals, lens);
}


void EmpiricalDistribution::buildDistribution(
    	const std::vector<uint32_t>& vals,
        const std::vector<uint32_t>& lens) {
    assert(vals.size() == lens.size());
    auto n = vals.size();

    minVal = std::numeric_limits<uint32_t>::max();
    maxVal = 0;
    double valsum = 0;

    for (size_t i = 0; i < n; ++i) {
        minVal = (vals[i] < minVal) ? vals[i] : minVal;
        maxVal = (vals[i] > maxVal) ? vals[i] : maxVal;
        valsum += lens[i];
    }

    double cumpr = 0.0;
    unsigned int lastval = 0, maxval = 1;
    for (; lastval < n; ++lastval) {
        cumpr += lens[lastval] / valsum;
        maxval = vals[lastval];
        if (cumpr > 1.0 - 1e-6) {
            break;
        }
    }

    auto upperBound = lastval + 1;
    pdfvals.resize(upperBound);
    valsum = 0.0;
    for (unsigned int i = 0; i < upperBound; ++i) {
        valsum += lens[i];
    }

    for (unsigned int val = 0, i = 0; val < upperBound; ) {
        if (val == vals[i]) {
            pdfvals[val] = lens[i] / valsum;
            ++val;
            ++i;
        }
        else if (val < vals[i]) {
            pdfvals[val] = 0.0;
            ++val;
        }
    }
    
    cdfvals.resize(upperBound);
    cdfvals[0] = pdfvals[0];
    mu = 0.0;
    for (unsigned int val = 1; val < upperBound; ++val) {
        cdfvals[val] = cdfvals[val - 1] + pdfvals[val];
        mu += val * pdfvals[val];
    }

    //for Circall
    // compute sd
    sdval = 0.0;
    for (unsigned int i = 0; i < upperBound; ++i) {     
          sdval += (i-mu)*(i-mu)*pdfvals[i];
    }
    sdval = std::sqrt(sdval);
    //for Circall

    // compute median
    size_t i = 0, j = n - 1;
    unsigned int u = lens[0], v = lens[n - 1];
    while (i < j) {
        if (u <= v) {
            v -= u;
            u = lens[++i];
        }
        else {
            u -= v;
            v = lens[--j];
        }
    }
    med = vals[i];
    isValid_ = true;
}

uint32_t EmpiricalDistribution::minValue() const {
    return minVal;
}


uint32_t EmpiricalDistribution::maxValue() const {
    return maxVal;
}

bool EmpiricalDistribution::valid() const {
    return (isValid_ and (pdfvals.size() > 0));
}

float EmpiricalDistribution::median() const
{
    if (pdfvals.size() == 0) return NAN;
    else return med;
}

float EmpiricalDistribution::mean () const
{
    if (pdfvals.size() == 0) return NAN;
    else return mu;
}

//for Circall
float EmpiricalDistribution::sd () const
{
    if (pdfvals.size() == 0) return NAN;
    else return sdval;
}

//for Circall

float EmpiricalDistribution::pdf(unsigned int x) const
{
    return x < pdfvals.size() ? pdfvals[x] : 0.0;
}


float EmpiricalDistribution::cdf(unsigned int x) const
{
    return x < cdfvals.size() ? cdfvals[x] : 1.0;
}

std::vector<int32_t> EmpiricalDistribution::realize(uint32_t numSamp) const {
  // start at 0 instead of minVal
  size_t distSize = maxVal + 1;
  std::vector<double> paddedPDF(distSize, 0.0);
  for (size_t i = 0; i <= maxVal; ++i) {
    paddedPDF[i] = pdf(i);
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int32_t> d(paddedPDF.begin(), paddedPDF.end());

  std::vector<int32_t> samples(distSize, 0);
  for (size_t i = 0; i < numSamp; ++i) {
    ++samples[d(gen)];
  }

  return samples;
}


