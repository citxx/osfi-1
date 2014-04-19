#pragma once

#include <string>
#include <vector>

class Spectrum {
 public:
  explicit Spectrum(
      double first_frequency,
      double step_frequency,
      const std::vector<double> &values);

  static Spectrum fromCsvFile(const std::string &path);
  static std::vector<Spectrum> manyFromCsvFile(const std::string &path);

  void adjustTo(const Spectrum &a);

  int size() const;
  double sum() const;
  Spectrum normalized() const;
  Spectrum operator -() const;
  Spectrum operator +() const;
  double operator [](int i) const;

  friend bool operator ==(const Spectrum &a, const Spectrum &b);
  friend bool operator !=(const Spectrum &a, const Spectrum &b);
 
  friend Spectrum operator +(double a, const Spectrum &b);
  friend Spectrum operator +(const Spectrum &a, double b);
  friend Spectrum operator +(const Spectrum &a, const Spectrum &b);
  friend Spectrum operator -(double a, const Spectrum &b);
  friend Spectrum operator -(const Spectrum &a, double b);
  friend Spectrum operator -(const Spectrum &a, const Spectrum &b);
 
  friend Spectrum operator *(double a, const Spectrum &b);
  friend Spectrum operator *(const Spectrum &a, double b);
  friend Spectrum operator *(const Spectrum &a, const Spectrum &b);
  friend Spectrum operator /(const Spectrum &a, double b);
  friend Spectrum operator /(double a, const Spectrum &b);
  friend Spectrum operator /(const Spectrum &a, const Spectrum &b);

  friend std::ostream & operator <<(std::ostream &stream, const Spectrum &t);

 private:
  std::vector<double> values_;
  double first_freq_, step_freq_;
};

bool operator ==(const Spectrum &a, const Spectrum &b);
bool operator !=(const Spectrum &a, const Spectrum &b);

Spectrum operator +(const Spectrum &a, const Spectrum &b);
Spectrum operator -(const Spectrum &a, const Spectrum &b);

Spectrum operator *(double a, const Spectrum &b);
Spectrum operator *(const Spectrum &a, double b);
Spectrum operator *(const Spectrum &a, const Spectrum &b);
Spectrum operator /(const Spectrum &a, double b);
Spectrum operator /(double a, const Spectrum &b);
Spectrum operator /(const Spectrum &a, const Spectrum &b);

std::ostream & operator <<(std::ostream &stream, const Spectrum &t);
