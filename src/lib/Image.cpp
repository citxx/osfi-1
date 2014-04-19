#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

#include "rgbe.h"
#include "Image.hpp"

Image::Image(int width, int height):
  width_(width),
  height_(height),
  data_(3 * width * height, 0.0) {}

Image Image::fromHDR(const std::string &path) {
  FILE *file = fopen(path.c_str(), "rb");

  int width, height;
  RGBE_ReadHeader(file, &width, &height, nullptr);

  Image img(width, height);

  RGBE_ReadPixels_RLE(file, const_cast<float*>(img.data_.data()), 3, width, height);

  fclose(file);

  return img;
}

void Image::setPixel(int x, int y, double r, double g, double b) {
  int pixel_pos = 3 * (width_ * y + x);
  data_[pixel_pos] = r;
  data_[pixel_pos + 1] = g;
  data_[pixel_pos + 2] = b;
}

void Image::savePPM(const std::string &path) const {
  std::ofstream file(path);

  file << "P3" << std::endl;
  file << width_ << " " << height_ << std::endl;
  file << "255" << std::endl;
  float mn = *std::min_element(data_.begin(), data_.end());
  float mx = *std::max_element(data_.begin(), data_.end());
  for (auto value : data_) {
    file << static_cast<int>(round(255 * pow((value - mn) / (mx - mn), 1 / 2.2))) << std::endl;
  }

  file.close();
}

void Image::saveHDR(const std::string &path) const {
  FILE *file = fopen(path.c_str(), "wb");

  RGBE_WriteHeader(file, width_, height_, nullptr);
  RGBE_WritePixels_RLE(file, const_cast<float*>(data_.data()), 3, width_, height_);

  fclose(file);
}

void Image::toneMapping() {
  double b = 0.85;
  double mx = *std::max_element(data_.begin(), data_.end());
  for (float &v : data_) {
    v = log(v + 1) / log(mx + 1) / log(2 + 8 * pow(v / mx, log(b) / log(0.5)));
  }
}
