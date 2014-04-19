#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>

#include "Spectrum.hpp"
#include "utils.hpp"

constexpr double EPS = 1e-6;

Spectrum::Spectrum(
    double first_frequency,
    double step_frequency,
    const std::vector<double> &values):
  values_(values),
  first_freq_(first_frequency),
  step_freq_(step_frequency) {}

Spectrum Spectrum::fromCsvFile(const std::string &path) {
  std::ifstream file(path);
  std::vector<double> frequency, value;

  std::string line;
  while (file >> line) {
    if (!line.empty()) {
      std::vector<std::string> fields = utils::split(line, ",");
      if (fields.size() != 2) {
        throw "Invalid spectrum file format near '" + line + "'";  
      }
      frequency.push_back(std::stod(fields[0]));
      value.push_back(std::stod(fields[1]));
    }
  }
  file.close();

  double first = frequency[0];
  double step = frequency[1] - frequency[0];
  for (int i = 1; i < static_cast<int>(frequency.size()); ++i) {
    if (abs(frequency[i] - frequency[i - 1] - step) > EPS) {
      throw "Invalid spectrum file: inconsistent step";
    }
  }

  return Spectrum(first, step, value);
}

void Spectrum::adjustTo(const Spectrum &a) {
  if (abs(step_freq_ - a.step_freq_) > EPS) {
    throw "Cannot adjust spectrums with different steps";
  }
  double step = step_freq_;

  double first_distance = (first_freq_ - a.first_freq_) / step;
  double tmp;
  if (modf(first_distance, &tmp) > EPS) {
    throw "Cannot adjust spectrums with different shifts";
  }

  double last_freq = first_freq_ + step * size();
  double a_last_freq = a.first_freq_ + step * a.size();
  double last_distance = (a_last_freq - last_freq) / step;

  int add_to_begin = round(std::max(0.0, first_distance));
  int add_to_end = round(std::max(0.0, last_distance));

  if (add_to_begin > 0) {
    first_freq_ = a.first_freq_;
    values_.insert(values_.begin(), add_to_begin, 0.0);
  }
  if (add_to_end > 0) {
    values_.insert(values_.end(), add_to_end, 0.0);
  }
}

std::vector<Spectrum> Spectrum::manyFromCsvFile(const std::string &path) {
  std::ifstream file(path);
  std::vector<double> frequency;
  std::vector<std::vector<double>> values;

  std::string line;
  while (file >> line) {
    if (!line.empty()) {
      std::vector<std::string> fields = utils::split(line, ",");
      if (!values.empty() && values.size() + 1 != fields.size()) {
        throw "Invalid spectrum file format near '" + line + "'";
      }
      frequency.push_back(std::stod(fields[0]));
      if (values.empty()) {
        values.resize(static_cast<int>(fields.size()) - 1);
      }
      for (int i = 1; i < static_cast<int>(fields.size()); ++i) {
        values[i - 1].push_back(std::stod(fields[i]));
      }
    }
  }
  file.close();

  double first = frequency[0];
  double step = frequency[1] - frequency[0];
  for (int i = 1; i < static_cast<int>(frequency.size()); ++i) {
    if (abs(frequency[i] - frequency[i - 1] - step) > EPS) {
      throw "Invalid spectrum file: inconsistent step";
    }
  }

  std::vector<Spectrum> result;
  for (auto value : values) {
    result.push_back(Spectrum(first, step, value));
  }
  return result;
}

int Spectrum::size() const {
  return values_.size();
}

double Spectrum::sum() const {
  return std::accumulate(values_.begin(), values_.end(), 0.0);
}

Spectrum Spectrum::normalized() const {
  return *this / sum();
}

Spectrum Spectrum::operator -() const {
  std::vector<double> v(values_.size());
  std::transform(values_.begin(), values_.end(), v.begin(), [](double x) { return -x; });
  return Spectrum(first_freq_, step_freq_, v);
}

Spectrum Spectrum::operator +() const {
  return *this;
}

double Spectrum::operator [](int i) const {
  return values_[i];
}

bool operator ==(const Spectrum &a, const Spectrum &b) {
  if (a.values_.size() != b.values_.size() ||
      a.first_freq_ != b.first_freq_ ||
      a.step_freq_ != b.step_freq_) {
    return false;
  }

  for (int i = 0; i < static_cast<int>(a.values_.size()); ++i) {
    if (abs(a.values_[i] - b.values_[i]) > EPS) {
      return false;
    }
  }
  
  return true;
}

bool operator !=(const Spectrum &a, const Spectrum &b) {
  return !(a == b);
}

Spectrum operator +(const Spectrum &a, const Spectrum &b) {
  std::vector<double> v(a.values_);
  for (int i = 0; i < static_cast<int>(v.size()); ++i) {
    v[i] += b.values_[i];
  }
  return Spectrum(a.first_freq_, a.step_freq_, v);
}

Spectrum operator +(double a, const Spectrum &b) {
  std::vector<double> v(b.values_.size());
  std::transform(b.values_.begin(), b.values_.end(), v.begin(), [a](double x) { return a + x; });
  return Spectrum(b.first_freq_, b.step_freq_, v);
}

Spectrum operator +(const Spectrum &a, double b) {
  std::vector<double> v(a.values_.size());
  std::transform(a.values_.begin(), a.values_.end(), v.begin(), [b](double x) { return x + b; });
  return Spectrum(a.first_freq_, a.step_freq_, v);
}

Spectrum operator -(const Spectrum &a, const Spectrum &b) {
  std::vector<double> v(a.values_);
  for (int i = 0; i < static_cast<int>(v.size()); ++i) {
    v[i] -= b.values_[i];
  }
  return Spectrum(a.first_freq_, a.step_freq_, v);
}

Spectrum operator -(double a, const Spectrum &b) {
  std::vector<double> v(b.values_.size());
  std::transform(b.values_.begin(), b.values_.end(), v.begin(), [a](double x) { return a - x; });
  return Spectrum(b.first_freq_, b.step_freq_, v);
}

Spectrum operator -(const Spectrum &a, double b) {
  std::vector<double> v(a.values_.size());
  std::transform(a.values_.begin(), a.values_.end(), v.begin(), [b](double x) { return x - b; });
  return Spectrum(a.first_freq_, a.step_freq_, v);
}

Spectrum operator *(double a, const Spectrum &b) {
  std::vector<double> v(b.values_.size());
  std::transform(b.values_.begin(), b.values_.end(), v.begin(), [a](double x) { return a * x; });
  return Spectrum(b.first_freq_, b.step_freq_, v);
}

Spectrum operator *(const Spectrum &a, double b) {
  std::vector<double> v(a.values_.size());
  std::transform(a.values_.begin(), a.values_.end(), v.begin(), [b](double x) { return x * b; });
  return Spectrum(a.first_freq_, a.step_freq_, v);
}

Spectrum operator *(const Spectrum &a, const Spectrum &b) {
  std::vector<double> v(a.values_);
  for (int i = 0; i < static_cast<int>(v.size()); ++i) {
    v[i] *= b.values_[i];
  }
  return Spectrum(a.first_freq_, a.step_freq_, v);
}

Spectrum operator /(const Spectrum &a, double b) {
  std::vector<double> v(a.values_.size());
  std::transform(a.values_.begin(), a.values_.end(), v.begin(), [b](double x) { return x / b; });
  return Spectrum(a.first_freq_, a.step_freq_, v);
}

Spectrum operator /(double a, const Spectrum &b) {
  std::vector<double> v(b.values_.size());
  std::transform(b.values_.begin(), b.values_.end(), v.begin(), [a](double x) { return a / x; });
  return Spectrum(b.first_freq_, b.step_freq_, v);
}

Spectrum operator /(const Spectrum &a, const Spectrum &b) {
  std::vector<double> v(a.values_);
  for (int i = 0; i < static_cast<int>(v.size()); ++i) {
    v[i] /= b.values_[i];
  }
  return Spectrum(a.first_freq_, a.step_freq_, v);
}


std::ostream & operator <<(std::ostream &stream, const Spectrum &t) {
  const int MAX_HEIGHT = 10;

  double mx = 0;
  for (int i = 0; i < t.size(); ++i) {
    mx = std::max(mx, t[i]);
    mx = std::max(mx, -t[i]);
  }

  int lower = 0, upper = 0;
  std::vector<int> height(t.size());
  for (int i = 0; i < t.size(); ++i) {
    height[i] = round(t[i] / mx * MAX_HEIGHT);
    lower = std::min(lower, height[i]);
    upper = std::max(upper, height[i]);
  }

  stream << "Specter: from " << t.first_freq_
    << " to " << t.first_freq_ + t.size() * t.step_freq_
    << " with step " << t.step_freq_ << std::endl;

  for (int i = upper; i > 0; --i) {
    for (auto h : height) {
      stream << ((h >= i) ? "#" : " ");
    }
    stream << std::endl;
  }

  for (auto h : height) {
    stream << "-";
  }
  stream << std::endl;

  for (int i = -1; i >= lower; --i) {
    for (auto h : height) {
      stream << ((h >= i) ? "#" : " ");
    }
    stream << std::endl;
  }

  return stream;
}
