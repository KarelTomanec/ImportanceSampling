#pragma once
#include <stdio.h>
#include "material.h"

class Intersectable;

// Structure to store the result of ray tracing
struct Hit {
    double t;
    vec3 position;
    vec3 normal;
    Material* material;
    Intersectable* object;
    Hit() { t = -1; }
};

// Abstract 3D object 
class Intersectable {
public:
    Material* material;
    double power;
    virtual Hit intersect(const Ray& ray) = 0;
    virtual double pointSampleProb(double totalPower) {
        printf("Point sample on table\n");
        return 0;
    }
};

// Rectangle 2D in 3D space
class Rect :
    public Intersectable
{
    // anchor point, normal, 
    vec3 r0, normal, right, forward;
    double width, height; // size
public:
    Rect(vec3 _r0, vec3 _r1, vec3 _r2,
        double _width, double _height, Material* mat) {
        r0 = _r0;
        vec3 L = _r1 - r0;
        vec3 V = _r2 - r0;
        // compute normal
        normal = (L.normalize() + V.normalize()).normalize();
        material = mat;
        power = 0; // default - does not emit light
        width = _width; height = _height;
        // recompute directions to get rectangle
        right = cross(vec3(0, 0, 1), normal).normalize();
        forward = cross(normal, right).normalize();
    }

    // Compute intersection between a ray and the rectangle
    Hit intersect(const Ray& ray) {
        Hit hit;
        double denom = dot(normal, ray.dir);
        if (fabs(denom) > epsilon) {
            hit.t = dot(normal, r0 - ray.start) / denom;
            if (hit.t < 0) return hit;
            hit.position = ray.start + ray.dir * hit.t;
            double x = dot(hit.position - r0, right);
            double y = dot(hit.position - r0, forward);
            if (fabs(x) > width || fabs(y) > height) {
                hit.t = -1;
                return hit;
            }
            hit.normal = normal;
            hit.position = ray.start + ray.dir * hit.t;
            hit.material = material;
            hit.object = this;
        }
        return hit;
    }
};

// Sphere used as light source
struct Sphere :
    public Intersectable
{
    vec3 center;
    double  radius;

    Sphere(const vec3& cent, double rad, Material* mat, bool light) {
        const double targetPower = 60;
        center = cent; radius = rad;
        material = mat;
        if (light) {
            power = material->Le.average() * (4 * radius * radius * M_PI) * M_PI;
        }
        else {
            power = 1.0f;
        }
        material->Le = material->Le * (targetPower / power);
        power = targetPower;
    }

    Hit intersect(const Ray& r) {
        Hit hit;
        vec3 dist = r.start - center;
        double b = dot(dist, r.dir) * 2.0;
        double a = dot(r.dir, r.dir);
        double c = dot(dist, dist) - radius * radius;
        double discr = b * b - 4.0 * a * c;
        if (discr < 0) return hit;
        double sqrt_discr = sqrt(discr);
        double t1 = (-b + sqrt_discr) / 2.0 / a;
        double t2 = (-b - sqrt_discr) / 2.0 / a;
        if (t1 <= 0 && t2 <= 0) return hit;
        if (t1 <= 0 && t2 > 0)
            hit.t = t2;
        else
            if (t2 <= 0 && t1 > 0)
                hit.t = t1;
            else
                if (t1 < t2)
                    hit.t = t1;
                else
                    hit.t = t2;
        hit.position = r.start + r.dir * hit.t;
        hit.normal = (hit.position - center) / radius;
        hit.material = material;
        hit.object = this;
        return hit;
    }

    // find a random point with uniform distribution on that half sphere, which can be visible
    void sampleUniformPoint(const vec3& illuminatedPoint, vec3& point, vec3& normal) {
        do {
            // uniform in a cube of edge size 2
            normal = vec3(drandom() * 2 - 1, drandom() * 2 - 1, drandom() * 2 - 1);
            if (dot(illuminatedPoint - center, normal) < 0) continue;	// ignore surely non visible points
        } while (dot(normal, normal) > 1);	// finish if the point is in the unit sphere
        normal = normal.normalize();	// project points onto the surface of the unit sphere
        point = center + normal * radius;	// project onto the real sphere
    }

    double pointSampleProb(double totalPower) {
        return power / totalPower / (4 * radius * radius * M_PI);
    }
};