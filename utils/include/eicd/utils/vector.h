#ifndef EICD_UTILS_VECTOR_HH
#define EICD_UTILS_VECTOR_HH

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <string>
#include <vector>

#include <Math/Vector4D.h>

#include <edm4hep/Vector3f.h>

#include <eicd/ReconstructedParticleCollection.h>
#include <eicd/ReconstructedParticleData.h>
#include <eicd/TrackParametersCollection.h>

namespace eicd::utils {

// inline edm4hep::Vector2f VectorFromPolar(const double r, const double theta)
// {
//  return {r * sin(theta), r * cos(theta)};
//}
inline edm4hep::Vector3f VectorFromSpherical(const double r, const double theta,
                                             const double phi) {
  using FloatType = decltype(edm4hep::Vector3f().x);
  const FloatType x = r * sin(theta) * cos(phi);
  const FloatType y = r * sin(theta) * sin(phi);
  const FloatType z = r * cos(theta);
  return {x, y, z};
}
template <class V> concept HasX = requires(V v) { v.x; };
template <class V> concept HasY = requires(V v) { v.y; };
template <class V> concept HasZ = requires(V v) { v.z; };
template <class V> concept Vector2D = HasX<V>&& HasY<V> && !HasZ<V>;
template <class V> concept Vector3D = HasX<V>&& HasY<V>&& HasZ<V>;
template <class V> concept VectorND = Vector2D<V> || Vector3D<V>;

template <Vector2D V> inline double mag(const V& v) {
  return std::hypot(v.x, v.y);
}
template <Vector3D V> inline double mag(const V& v) {
  return std::hypot(v.x, v.y, v.z);
}
template <Vector3D V> inline double theta(const V& v) {
  return std::atan2(std::hypot(v.x, v.y), v.z);
}
template <Vector3D V> inline double phi(const V& v) {
  return std::atan2(v.y, v.x);
}
template <Vector3D V> inline double eta(const V& v) {
  return -std::log(std::tan(0.5 * theta(v)));
}
} // namespace eicd::utils
template <eicd::utils::Vector2D V>
inline V operator+(const V& v1, const V& v2) {
  return {v1.x + v2.x, v1.y + v2.y};
}
template <eicd::utils::Vector3D V>
inline V operator+(const V& v1, const V& v2) {
  return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}
template <eicd::utils::Vector2D V>
inline double operator*(const V& v1, const V& v2) {
  return v1.x * v2.x + v1.y * v2.y;
}
template <eicd::utils::Vector3D V>
inline double operator*(const V& v1, const V& v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
template <eicd::utils::Vector2D V>
inline V operator*(const double d, const V& v) {
  return {d * v.x, d * v.y};
}
template <eicd::utils::Vector3D V>
inline V operator*(const double d, const V& v) {
  return {d * v.x, d * v.y, d * v.z};
}
template <eicd::utils::Vector2D V>
inline V operator*(const V& v, const double d) {
  return {d * v.x, d * v.y};
}
template <eicd::utils::Vector3D V>
inline V operator*(const V& v, const double d) {
  return {d * v.x, d * v.y, d * v.z};
}
template <eicd::utils::VectorND V>
inline V operator-(const V& v1, const V& v2) {
  return v1 + (-1. * v2);
}
template <eicd::utils::VectorND V>
inline double operator/(const V& v, const double d) {
  return (1. / d) * v;
}
#endif
