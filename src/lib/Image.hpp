#pragma once

#include <string>
#include <vector>

class Image {
 public:
  explicit Image(int width, int height);

  static Image fromHDR(const std::string &path);

  void setPixel(int x, int y, double r, double g, double b);
  void savePPM(const std::string &path) const;
  void saveHDR(const std::string &path) const;

  void toneMapping();

 private:
  int width_, height_;
  std::vector<float> data_;
};
