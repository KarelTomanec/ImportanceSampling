#pragma once
#include "vec3.h"

// Structure for a ray
struct Ray {
    vec3 start, dir;
    Ray(const vec3& _start, const vec3& _dir) { start = _start; dir = _dir.normalize(); }
};