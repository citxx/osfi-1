#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "lib/Image.hpp"
#include "lib/Spectrum.hpp"
#include "lib/Vector3D.hpp"

const double PI = 3.14159265358979;
const int R_N = 50;
const int PHI_N = 60;

void printUsage(const std::string &name) {
  std::cout << "Usage: " << name << " <input file> <output file>" << std::endl;
}

bool parseArguments(
    int argc, char *argv[],
    std::string *input_path,
    std::string *output_path) {

  if (argc != 3) {
    return false;
  }

  *input_path = argv[1];
  *output_path = argv[2];

  return true;
}

int main(int argc, char *argv[]) {
  std::string input_path, output_path;
  bool correct_args = parseArguments(
      argc, argv,
      &input_path,
      &output_path);
  if (!correct_args) {
    std::cerr << "Invalid arguments" << std::endl;
    printUsage(argv[0]);
    return -1;
  }

  Image img = Image::fromHDR(input_path);
  img.toneMapping();
  img.savePPM(output_path);

  return 0;
}
