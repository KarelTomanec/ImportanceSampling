#pragma once
const double epsilon = 1e-9; // Small value
// Vector 3D
struct vec3
{
	double x, y, z;

	vec3() { x = y = z = 0; }

	vec3(double x0, double y0, double z0 = 0) {
		x = x0; y = y0; z = z0;
	}

	vec3(double x0) {
		x = x0; y = x0; z = x0;
	}

	vec3 operator -() const {
		return vec3(-x, -y, -z);
	}

	vec3 operator*(double a)  const { return vec3(x * a, y * a, z * a); }
	vec3 operator*(const vec3 r)  const { return vec3(x * r.x, y * r.y, z * r.z); }
	vec3 operator/(const double r)  const {
		if (fabs(r) > epsilon)
			return vec3(x / r, y / r, z / r);
		else
			return vec3(0, 0, 0);
	}
	vec3 operator+(const vec3& v)  const { return vec3(x + v.x, y + v.y, z + v.z); }
	vec3 operator-(const vec3& v)  const { return vec3(x - v.x, y - v.y, z - v.z); }
	void operator+=(const vec3& v) { x += v.x, y += v.y, z += v.z; }
	void operator*=(double a) { x *= a, y *= a, z *= a; }
	double length()  const { return sqrt(x * x + y * y + z * z); }
	vec3 normalize() const {
		double l = length();
		if (l > epsilon)
			return (*this) / l;
		else
			return vec3(0, 0, 0);
	}
	double average() { return (x + y + z) / 3; }

	double luminance() {
		return x * 0.2126 + y * 0.7152 + z * 0.0722;
	}
};

inline vec3 operator *(const vec3& v, float a) {

	return vec3(a * v.x, a * v.y, a * v.z);
}

inline vec3 operator *(float a, const vec3& v) {
	return vec3(a * v.x, a * v.y, a * v.z);
}

// dot product of two vectors
double dot(const vec3& v1, const vec3& v2)
{
	return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// cross product of two vectors
vec3 cross(const vec3& v1, const vec3& v2)
{
	return vec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

vec3 reflect(const vec3& n, const vec3& v) {
	return 2.0 * dot(v, n) * n - v;
}


using color = vec3;
using point3 = vec3;