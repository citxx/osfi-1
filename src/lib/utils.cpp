#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "utils.hpp"

namespace utils {

std::vector<std::string> split(const std::string &s, const std::string &seps) {
  std::vector<std::string> parts;

  int pos = 0;
  for (int st = s.find_first_not_of(seps); st != std::string::npos; st = s.find_first_not_of(seps, pos)) {
    pos = s.find_first_of(seps, st);
    if (pos == std::string::npos) {
      pos = s.size();
    }
    parts.push_back(s.substr(st, pos - st));
  }
  return parts;
}

}
