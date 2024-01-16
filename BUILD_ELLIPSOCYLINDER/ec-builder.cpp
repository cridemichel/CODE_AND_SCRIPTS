#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cassert>

namespace numbers
{
  inline constexpr double pi = 3.14159265358979323846;
}

template <typename T>
struct Vec3
{
  T x;
  T y;
  T z;
};

template <typename T>
Vec3<T> operator+(const Vec3<T>& left, const Vec3<T>& right)
{
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

template <typename T>
Vec3<T> operator-(const Vec3<T>& left, const Vec3<T>& right)
{
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vec3<T>& v)
{
  os << v.x << " " << v.y << " " << v.z;
  return os;
}

using Vec3d = Vec3<double>;
using Vec3u = Vec3<unsigned int>;

struct Model
{
  std::vector<Vec3d> positions;
  std::vector<Vec3d> normals;
  std::vector<Vec3u> faces;
  double volume;
};

namespace
{
  double dot(const Vec3d& left, const Vec3d& right)
  {
    return left.x * right.x + left.y * right.y + left.z * right.z;
  }

  [[nodiscard]] Vec3d normalize(const Vec3d& v)
  {
    const auto a = 1.0 / std::sqrt(dot(v, v));
    return {v.x * a, v.y * a, v.z * a};
  }

  Vec3d cross(const Vec3d& left, const Vec3d& right)
  {
    Vec3d result;
    result.x = left.y * right.z - left.z * right.y;
    result.y = left.z * right.x - left.x * right.z;
    result.z = left.x * right.y - left.y * right.x;
    return result;
  }
}

Model generate_ellipsocylinder(
    double total_length,
    double cylinder_lenght,
    double diameter,
    size_t stack_count = 4,  // number of polar divisions (poles excluded)
    size_t sector_count = 24 // number of azimuthal divisions
    );
Model generate_ellipsoid(
    const Vec3d& aspect,
    size_t stack_count,
    size_t sector_count);

void print_positions(const Model& model, const std::string& filename);
void print_positions_and_normals(const Model& model, const std::string& filename);
void print_obj(const Model& model);

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    std::cerr << "[ERROR]: usage [bin] total_length cylinder_lenght diameter" << std::endl;
    exit(1);
  }

  const auto e = generate_ellipsocylinder(
      std::stod(argv[1]),
      std::stod(argv[2]),
      std::stod(argv[3]), 
      4,
      24);
  
  print_positions(e, "vertices.dat");
  print_positions_and_normals(e, "vertnorm_opengl.dat");
  print_obj(e);

  return 0;
}

Model generate_ellipsocylinder(
    double total_length,
    double cylinder_lenght,
    double diameter,
    size_t stack_count,
    size_t sector_count)
{
  Model e;

  // ellipsoid semi-axes
  const auto aspect = Vec3d{(total_length - cylinder_lenght) / 2.0, diameter / 2.0, diameter / 2.0};

  // vertices generation
  e.positions = generate_ellipsoid(aspect, stack_count, sector_count).positions;
  const auto half_count = static_cast<size_t>(std::floor(static_cast<double>(e.positions.size()) / 2.0));
  for (size_t i = 0; i < half_count; ++i)
  {
    auto& p = e.positions[i];
    p.x -= cylinder_lenght / 2.0;
  }
  for (size_t i = half_count; i < e.positions.size(); ++i)
  {
    auto& p = e.positions[i];
    p.x += cylinder_lenght / 2.0;
  }
  std::cout << "Positions count: " << e.positions.size() << "\n";

  // indices generation
  for (size_t i = 0; i < sector_count; ++i)
  {
    const unsigned int i1 = 1 + (i + 0);
    const unsigned int i2 = 0;
    const unsigned int i3 = 1 + (i + 1) % sector_count;
    e.faces.push_back({i1, i2, i3});
  }
  for (size_t i = 0; i < stack_count - 1; ++i)
  {
    for (size_t j = 0; j < sector_count; ++j)
    {
      const unsigned int i1 = 1 + (i + 0) * sector_count + (j + 0);
      const unsigned int i2 = 1 + (i + 0) * sector_count + (j + 1) % sector_count;
      const unsigned int i3 = 1 + (i + 1) * sector_count + (j + 0);
      const unsigned int i4 = 1 + (i + 1) * sector_count + (j + 1) % sector_count;
      e.faces.push_back({i1, i4, i3});
      e.faces.push_back({i1, i2, i4});
    }
  }
  for (size_t i = 0; i < sector_count; ++i)
  {
    const unsigned int i1 = 1 + (stack_count - 1) * sector_count + (i + 0);
    const unsigned int i2 = 1 + (stack_count - 1) * sector_count + (i + 1) % sector_count;
    const unsigned int i3 = 1 + (stack_count - 0) * sector_count;
    e.faces.push_back({i1, i2, i3});
  }
  std::cout << "Faces count: " << e.faces.size() << "\n";

  // normals generation
  for (const auto& face : e.faces)
  {
    const auto& p1 = e.positions[face.x];
    const auto& p2 = e.positions[face.y];
    const auto& p3 = e.positions[face.z];
    const auto normal = normalize(cross(p2 - p1, p3 - p1));
    e.normals.push_back(normal);
  }
  std::cout << "Normals count: " << e.normals.size() << "\n";

  // volume computation
  e.volume = 0;
  for (const auto& face : e.faces)
  {
    const auto& p1 = e.positions[face.x];
    const auto& p2 = e.positions[face.y];
    const auto& p3 = e.positions[face.z];
    e.volume += std::abs(dot(p1, cross(p2, p3))) / 6.0;
  }
  std::cout << "Volume: " << e.volume << "\n";

  return e;
}

void print_positions(const Model& e, const std::string& filename)
{
  std::ofstream file(filename, std::ios::out);
  if (file.fail())
  {
    std::cerr << "[ERROR]: could not open file " << filename << std::endl;
    exit(1);
  }

  file << e.volume << "\n"; 
  for (const auto& p : e.positions)
  {
    file << p << "\n";
  }

  std::cout << "Written to: " << filename << "\n";
}

void print_positions_and_normals(const Model& e, const std::string& filename)
{
  std::ofstream file(filename, std::ios::out);
  if (file.fail())
  {
    std::cerr << "[ERROR]: could not open file " << filename << std::endl;
    exit(1);
  }

  for (size_t i = 0; i < e.faces.size(); ++i)
  {
    const auto& face = e.faces[i];
    file << e.positions[face.x] << " " << e.normals[i] << "\n";
    file << e.positions[face.y] << " " << e.normals[i] << "\n";
    file << e.positions[face.z] << " " << e.normals[i] << "\n";
  }
  
  std::cout << "Written to: " << filename << "\n";
}

void print_obj(const Model& e)
{
  std::ofstream file("ellipsocylinder.obj", std::ios::out);
  if (file.fail())
  {
    std::cerr << "[ERROR]: could not open file ellipsocylinder.obj" << std::endl;
    exit(1);
  }

  for (const auto& p : e.positions)
  {
    file << "v " << p.x << " " << p.y << " " << p.z << "\n";
  }

  for (const auto& n : e.normals)
  {
    file << "vn " << n.x << " " << n.y << " " << n.z << "\n";
  }

  for (const auto& face : e.faces)
  {
    static int normal_index = 1;
    file << "f ";
    file << face.x + 1 << "//" << normal_index << " ";
    file << face.y + 1 << "//" << normal_index << " ";
    file << face.z + 1 << "//" << normal_index << " ";
    file << "\n";
    normal_index++;
  }

  std::cout << "Written to: ellipsocylinder.obj" << "\n";
}

Model generate_ellipsoid(
    const Vec3d& aspect,
    size_t stack_count,
    size_t sector_count)
{
  Model se{};

  auto uv = [](const Vec3d& aspect, double theta, double phi)
  {
    Vec3d e{};
    e.x = aspect.x * std::cos(theta);
    e.y = aspect.y * std::sin(theta) * std::cos(phi);
    e.z = aspect.z * std::sin(theta) * std::sin(phi);
    return e;
  };

  if (stack_count % 2 != 0)
  {
    std::cerr << "[ERROR]: stack count must be even (pole excluded)" << std::endl;
    exit(1);
  }

  const auto half_stack_count = static_cast<size_t>(std::floor(static_cast<double>(stack_count) / 2.0));

  se.positions.push_back({-aspect.x, 0.0, 0.0});
  for (size_t i = 1; i <= half_stack_count; ++i)
  {
    const auto theta = numbers::pi - numbers::pi * static_cast<double>(i) / static_cast<double>(stack_count);
    for (size_t j = 0; j < sector_count; ++j)
    {
      const auto phi = 2.0 * numbers::pi * static_cast<double>(j) / static_cast<double>(sector_count);
      const auto [x, y, z] = uv(aspect, theta, phi);
      se.positions.push_back({x, y, z});
    }
  }
  for (size_t i = half_stack_count; i < stack_count; ++i)
  {
    const auto theta = numbers::pi - numbers::pi * static_cast<double>(i) / static_cast<double>(stack_count);
    for (size_t j = 0; j < sector_count; ++j)
    {
      const auto phi = 2.0 * numbers::pi * static_cast<double>(j) / static_cast<double>(sector_count);
      const auto [x, y, z] = uv(aspect, theta, phi);
      se.positions.push_back({x, y, z});
    }
  }
  se.positions.push_back({aspect.x, 0.0, 0.0});

  return se;
}
