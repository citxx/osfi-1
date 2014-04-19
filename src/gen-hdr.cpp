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
  std::cout << "Usage: " << name << " <options>" << std::endl
    << std::endl
    << "\t-N <N>                   set the width and height of generated image to <N>, default to 100" << std::endl
    << "\t-h, --height <height>    set z position of light source to <height>, default to 1.0 m" << std::endl
    << "\t-S, --area <area>        set area of the light source to <area>, default to 1.0 m^2" << std::endl
    << "\t-t, --theta <angle>      set the angle the light source rotated by to <angle>, default to 0.0" << std::endl
    << "\t-LSPD <path>             set the file with light source intencity spectrum, required" << std::endl
    << "\t-P, --power <power>      set the total power of the light source, default to 1.0 Watt" << std::endl
    << "\t-ASPD <path>             set the file with material reflection spectrum, required" << std::endl
    << "\t-LumSPD <path>           set the file with luminocity function, required" << std::endl
    << "\t-XYZSPD <path>           set the file with CIE XYZ transofrmation functions, required" << std::endl
    << "\t-o, --output <path>      set the name of generated image, default to 'default.ppm'" << std::endl;
}

bool parseArguments(
    int argc,
    char *argv[],
    int *N,
    double *height,
    double *area,  // Площадь
    double *theta,
    std::string *lspd,
    double *power,  // Мощность
    std::string *aspd,
    std::string *lum_spd,
    std::string *xyz_path,
    std::string *image_path) {
  bool is_lspd_set = false, is_aspd_set = false, is_lum_spd_set = false, is_xyz_spd_set = false;

  // Default values
  *N = 100;
  *height = 1.0;
  *area = 1.0;
  *theta = 0.0;
  *power = 1.0;
  *image_path = "default.ppm";

  if (argc % 2 == 0) {
    std::cerr << "Invalid number of arguments" << std::endl;
    return false;
  }

  for (int i = 1; i < argc; i += 2) {
    std::string arg = argv[i];
    if (arg == "-h" || arg == "--height") {
      *height = std::stod(argv[i + 1]);  // TODO: catch
    } else if (arg == "-S" || arg == "--area") {
      *area = std::stod(argv[i + 1]);  // TODO: catch
    } else if (arg == "-N") {
      *area = std::stoi(argv[i + 1]);  // TODO: catch
    } else if (arg == "-t" || arg == "--theta") {
      *theta = std::stod(argv[i + 1]) * PI / 180; // TODO: catch
    } else if (arg == "-P" || arg == "--power") {
      *power = std::stod(argv[i + 1]); // TODO: catch
    } else if (arg == "-LSPD") {
      is_lspd_set = true;
      *lspd = argv[i + 1];
    } else if (arg == "-ASPD") {
      is_aspd_set = true;
      *aspd = argv[i + 1];
    } else if (arg == "-LumSPD") {
      is_lum_spd_set = true;
      *lum_spd = argv[i + 1];
    } else if (arg == "-XYZSPD") {
      is_xyz_spd_set = true;
      *xyz_path = argv[i + 1];
    } else if (arg == "-o" || arg == "--output") {
      *image_path = argv[i + 1];
    } else {
      std::cerr << "Unexpected argument " << arg << std::endl;
      return false;
    }
  }
  return is_lspd_set && is_aspd_set && is_lum_spd_set && is_xyz_spd_set;
}

void getRgb(double x, double y, double z, double *r, double *g, double *b) {
  *r =  3.2404542 * x - 1.5371385 * y - 0.4985314 * z;
  *g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z;
  *b =  0.0556434 * x - 0.2040259 * y + 1.0572252 * z;
  if (*r < 0.0) *r = 0.0;
  if (*g < 0.0) *g = 0.0;
  if (*b < 0.0) *b = 0.0;
}

int main(int argc, char *argv[]) {
  int N;
  double height, area, theta, power;
  std::string lspd_path, aspd_path, lum_spd_path, image_path, xyz_path;
  bool correct_args = parseArguments(
      argc, argv,
      &N,
      &height,
      &area,
      &theta,
      &lspd_path,
      &power,
      &aspd_path,
      &lum_spd_path,
      &xyz_path,
      &image_path);
  if (!correct_args) {
    std::cerr << "Invalid arguments" << std::endl;
    printUsage(argv[0]);
    return -1;
  }

  Spectrum light_spectrum = Spectrum::fromCsvFile(lspd_path).normalized();
  Spectrum absorption_spectrum = Spectrum::fromCsvFile(aspd_path);
  Spectrum luminance_function = Spectrum::fromCsvFile(lum_spd_path);

  light_spectrum.adjustTo(absorption_spectrum);
  absorption_spectrum.adjustTo(light_spectrum);
  luminance_function.adjustTo(light_spectrum);
  std::cout << light_spectrum << absorption_spectrum << luminance_function;

  Spectrum reflected = luminance_function * light_spectrum * (1 - absorption_spectrum);
  std::cout << reflected;

  std::vector<Spectrum> xyz = Spectrum::manyFromCsvFile(xyz_path);
  Spectrum reflected_x = xyz[0] * reflected;
  Spectrum reflected_y = xyz[1] * reflected;
  Spectrum reflected_z = xyz[2] * reflected;
  std::cout << reflected_x << reflected_y << reflected_z;

  double radiance = power / (2 * PI * area);
  double radius = sqrt(area / PI);
  Vector3D normal(sin(theta), 0.0, -cos(theta));
  Vector3D r_u(0.0, 1.0, 0.0);
  Vector3D r_v = Vector3D::cross(normal, r_u).normalized();
  Vector3D z(0.0, 0.0, 1.0);

  std::cout << "normal = " << normal << std::endl;
  std::cout << "r_u = " << r_u << std::endl;
  std::cout << "r_v = " << r_v << std::endl;
  std::cout << "z = " << z << std::endl;

  double D_R = radius / R_N;
  double D_PHI = 2 * PI / PHI_N;

  Image img(N, N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      double x = static_cast<double>(i) / (N - 1) * 2 * height - height;
      double y = static_cast<double>(j) / (N - 1) * 2 * height - height;
      Vector3D v(x, y, -height);
      
      double irradiance = 0.0;
      for (int r_i = 0; r_i < R_N; ++r_i) {
        for (int phi_i = 0; phi_i < PHI_N; ++phi_i) {
          Vector3D r = r_i * D_R * (cos(phi_i * D_PHI) * r_u + sin(phi_i * D_PHI) * r_v);
          Vector3D q = v - r;

          double q2 = Vector3D::dot(q, q);
          double d_irradiance =
            radiance * fabs(Vector3D::dot(normal, q)) / q2 * Vector3D::dot(z, -q) / q2 * r_i * D_R * D_R * D_PHI;
          //if (d_irradiance < 0.000001 && x > 0) {
            //std::cout << Vector3D::dot(normal, q) << " " << fabs(Vector3D::dot(normal, q)) << " " << normal << " " << q << std::endl;
          //}
          irradiance += d_irradiance;
        }
      }
      //std::cout << v << " " << irradiance << std::endl;

      double cie_x = (irradiance * reflected_x).sum();
      double cie_y = (irradiance * reflected_y).sum();
      double cie_z = (irradiance * reflected_z).sum();
      double r, g, b;
      getRgb(cie_x, cie_y, cie_z, &r, &g, &b);
      img.setPixel(i, j, r, g, b);
    }
  }
  img.saveHDR(image_path);
}
