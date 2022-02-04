
#pragma once
#include "image.h"
#include <random>

class Light {
public:

    Light() = default;

    vec3 normal;

    virtual color sample_l(const point3& p, vec3* wi, float* pdf) = 0;

    virtual float pdf(const point3& p, const vec3& w, double totalPower) const = 0;

    virtual color get_illumination(const vec3& v) const = 0;

};

// The light source represented by a sphere
class SphereLight : public Light{
public:
    Sphere* sphere;
    vec3 p;
    SphereLight(Sphere* _sphere, vec3 _point, vec3 _normal) {
        sphere = _sphere, p = _point; normal = _normal;
    }

    color sample_l(const point3 & p, vec3* wi, float* pdf) override;

    float pdf(const point3& p, const vec3& w, double totalPower) const override;

    color get_illumination(const vec3& v) const override;

};


color SphereLight::sample_l(const point3& p, vec3* wi, float* pdf) {
    return sphere->material->Le;
}

float SphereLight::pdf(const point3& p, const vec3& w, double totalPower) const {
    return sphere->pointSampleProb(totalPower);
}

color SphereLight::get_illumination(const vec3& v) const {
    return sphere->material->Le;
}


class InfiniteAreaLight : public Light {
private:
    Distribution2D* distribution;

public:

    std::unique_ptr<img> radianceMap;

    InfiniteAreaLight() = default;

    InfiniteAreaLight(const char* filename) {


        normal = vec3(0, 0, 0);
        radianceMap = std::make_unique<img>(filename);
        int nu = radianceMap->width;
        int nv = radianceMap->height;
        float *im = new float[nu * nv];
       
        for (int v = 0; v < nv; ++v) {
            for (int u = 0; u < nu; ++u)
            {
                float luminance = (radianceMap->lookup(u, v)).luminance();
                im[v* nu + u] = luminance;
            }
        }

        // init sampling pdfs for infinite area light

        float* func = new float[nu * nv];


        for (int v = 0; v < nv; ++v)
        {
            float sinTheta = sinf(M_PI * float(v + 0.5f) / float(nv));
            // compute sampling distribution fo column u
            for (int u = 0; u < nu; u++)
            {
                func[v * nu + u] = im[v * nu + u] * sinTheta;
            }
        }

        // compute sampling distribution for columns of image

        distribution = new Distribution2D(func, nu, nv);

        delete[] im;
        delete[] func;
    }

    color sample_l(const point3& p, vec3* wi, float* pdf);

    color sample_l_debug() const;

    color get_illumination(const vec3& v) const;

    color get_illumination_environment(const vec3& v) const;

    float pdf(const point3& p, const vec3& w, double totalPower) const;
    

};

color InfiniteAreaLight::sample_l(const point3& p, vec3* wi, float* pdf) {

    const float u1 = drandom();
    const float u2 = drandom();

    float uv[2];
    float mapPdf;
    distribution->SampleContinuous(u1, u2, uv, &mapPdf);
    if (mapPdf == 0.f) return color(0, 0, 0);

    float theta = uv[1] * M_PI;
    float phi = uv[0] * 2.f * M_PI + M_PI / 2.f;
    //float phi = 2.0f * M_PI * u1;
    //float theta = acosf(1 - 2.0f * u2);

    float costheta = cosf(theta);
    float sintheta = sinf(theta);
    float sinphi = sinf(phi);
    float cosphi = cosf(phi);

    //*wi = vec3(sintheta * cosphi, sintheta * sinphi, costheta); // pbrt
    *wi = vec3(sintheta * cosphi, costheta, sintheta * sinphi); // funkcni
    //*wi = vec3(sintheta * sinphi, costheta, sintheta * cosphi); // z tabule

    normal = ((*wi) * (-1.0f)).normalize();
    
    *pdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    //*pdf = 1.0f / (4.0f * M_PI);
    if (sintheta == 0.f) *pdf = 0.f;


    //radianceMap->write(uv[0], uv[1], color(1, 0, 0));
    color L = radianceMap->lookupUV(uv[0], uv[1]);

    //std::cout << L.x << " " << L.y << " " << L.z << std::endl;
    return L;
}


color InfiniteAreaLight::get_illumination(const vec3& dir) const {

    float u = SphericalPhi(dir) * INV_2_PI;
    float v = SphericalTheta(dir) * INV_PI;

    //radianceMap->write(u, v, color(0, 1, 0)); // debug
    color c = radianceMap->lookupUV(u, v);
    return c;
}

color InfiniteAreaLight::get_illumination_environment(const vec3& dir) const {

    float u = SphericalPhi(dir) * INV_2_PI;
    float v = SphericalTheta(dir) * INV_PI;

    //radianceMap->write(u, v, color(0, 0, 1)); // debug
    return radianceMap->lookupUV(u, v);
}

float InfiniteAreaLight::pdf(const point3& p, const vec3& w, double totalPower) const {
    vec3 wi = w;
    float theta = SphericalTheta(wi);
    float phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float pd = distribution->Pdf(phi * INV_2_PI, theta * INV_PI) /
            (2.f * M_PI * M_PI * sintheta);

    //return 1.0f / (4.0f * M_PI);
    return pd;

}



