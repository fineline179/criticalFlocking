#pragma once
#include <cmath>

namespace sim {

// Minimal double-precision 3-vector for the simulation core (which has no Cinder dependency,
// so it can also run headless).
struct Vec3 {
    double x = 0, y = 0, z = 0;

    Vec3() = default;
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3  operator+(const Vec3 &o) const { return { x + o.x, y + o.y, z + o.z }; }
    Vec3  operator-(const Vec3 &o) const { return { x - o.x, y - o.y, z - o.z }; }
    Vec3  operator-() const { return { -x, -y, -z }; }
    Vec3  operator*(double s) const { return { x * s, y * s, z * s }; }
    Vec3  operator/(double s) const { return { x / s, y / s, z / s }; }
    Vec3 &operator+=(const Vec3 &o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3 &operator-=(const Vec3 &o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3 &operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

inline Vec3   operator*(double s, const Vec3 &v) { return v * s; }
inline double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3   cross(const Vec3 &a, const Vec3 &b)
{
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
inline double length2(const Vec3 &a) { return dot(a, a); }
inline double length(const Vec3 &a) { return std::sqrt(dot(a, a)); }
inline Vec3   normalize(const Vec3 &a)
{
    double l = length(a);
    return l > 0 ? a / l : Vec3(1, 0, 0);
}
// component of v perpendicular to the unit vector n
inline Vec3 perp(const Vec3 &v, const Vec3 &n) { return v - n * dot(v, n); }

// Rotate v about the unit axis k by angle a (Rodrigues' formula)
inline Vec3 rotate(const Vec3 &v, const Vec3 &k, double a)
{
    double c = std::cos(a), s = std::sin(a);
    return v * c + cross(k, v) * s + k * (dot(k, v) * (1 - c));
}

} // namespace sim
