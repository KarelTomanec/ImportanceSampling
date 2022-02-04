#pragma once
#include <math.h>
#include "ray.h"

// Definition of the camera
class Camera
{
    // center of projection and orthogonal basis of the camera
    vec3 eye, lookat, right, up;
public:
    void set(const vec3& _eye, const vec3& _lookat, const vec3& _vup, double fov) {
        eye = _eye;
        lookat = _lookat;
        vec3 w = eye - lookat;
        double f = w.length();
        //right = cross(_vup, w).normalize() * f * tan(fov / 2);
        right = cross(_vup, w).normalize() * f * tan(fov / 2);
        up = cross(w, right).normalize() * f * tan(fov / 2);
    }
    Ray getRay(int X, int Y) { // X,Y - pixel coordinates, compute a primary ray
        vec3 dir = lookat +
            right * (2.0 * (X + 0.5) / screenWidth - 1) +
            up * (2.0 * (Y + 0.5) / screenHeight - 1) - eye;
        return Ray(eye, dir.normalize());
    }
};