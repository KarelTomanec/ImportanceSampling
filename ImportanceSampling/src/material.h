#pragma once
#include <math.h>
#include "ray.h"
#include "utility.h"

// The definition of material surface (BRDF + emission)
class Material
{
public:
    vec3 Le; // the emmited power
    vec3 diffuseAlbedo; // albedo for diffuse component
    vec3 specularAlbedo; // albedo for specular component
    double  shininess;

    Material() { shininess = 0; Le = vec3(0); }


    // Evaluate the BRDF given normal, view direction (outgoing) and light direction (incoming)
    vec3 BRDF(const vec3& N, const vec3& V, const vec3& L) {
        vec3 brdf(0, 0, 0);
        double cosThetaL = dot(N, L);
        double cosThetaV = dot(N, V);

        if (cosThetaL <= epsilon || cosThetaV <= epsilon) return brdf;

        brdf = diffuseAlbedo / M_PI; // diffuse part

        vec3 R = reflect(N, L);
        double cosPhi = dot(V, R);

        if (cosPhi <= 0)
            return brdf; // farther by PI/2 from reflected direction
          // max-Phong specular BRDF: symmetric and energy conserving 
        return brdf + specularAlbedo * ((shininess + 1.0) / 2.0 / M_PI * pow(cosPhi, shininess) / fmax(cosThetaL, cosThetaV));
    }
    // BRDF.cos(theta) importance sampling for input normal, outgoing direction
    bool sampleDirection(const vec3& N, const vec3& V,
        vec3& L) { // output - the incoming light direction
    // To be implemented during exercise 1


        L = vec3(0, 0, 0);
        const float u_m = drandom();
        const float v_m = drandom();

        const float roulette = drandom();

        if (roulette < diffuseAlbedo.average()) {
            // diffuse sample
            const double x = sqrt(1.0 - u_m) * cos(2.0 * M_PI * v_m);
            const double y = sqrt(1.0 - u_m) * sin(2.0 * M_PI * v_m);
            const double z = sqrt(u_m);

            const vec3 k = N;
            const vec3 w = vec3(2.0 * drandom() - 1.0, 2.0 * drandom() - 1.0, 2.0 * drandom() - 1.0).normalize(); // random unit vector w != N
            const vec3 i = cross(N, w).normalize();
            const vec3 j = cross(i, k);


            L = i * x + j * y + k * z;

            if (dot(N, L) < 0) {
                return false;
            }
            return true;
        }
        else if (roulette < diffuseAlbedo.average() + specularAlbedo.average()) {
            // specular sample
            const double x = sqrt(1.0 - pow(u_m, 2.0 / (shininess + 1.0))) * cos(2.0 * M_PI * v_m);
            const double y = sqrt(1.0 - pow(u_m, 2.0 / (shininess + 1.0))) * sin(2.0 * M_PI * v_m);
            const double z = pow(u_m, 1.0 / (shininess + 1.0));

            const vec3 k = V;
            const vec3 i = cross(V, N).normalize();
            const vec3 j = cross(i, k);

            const vec3 R = i * x + j * y + k * z;

            L = (N * dot(N, R) * 2.0 - R).normalize();

            if (dot(N, L) < 0) {
                return false;
            }

            return true;
        }
        else {
            // the contribution is zero
            return false;
        }

    }
    // Evaluate the probability given input normal, view (outgoing) direction and incoming light direction
    double sampleProb(const vec3& N, const vec3& V, const vec3& L) {
        // To be implemented during exercise 1

        return diffuseAlbedo.average() * dot(N, L) / M_PI + specularAlbedo.average() * ((shininess + 1.0) / (2.0 * M_PI)) * pow((fmax(dot(V, reflect(N, L)), 0.0)), shininess);

    }
};

// Material used for light source
class LightMaterial : public Material {
public:
    LightMaterial(vec3 _Le) { Le = _Le; }
};

// Material used for objects, given how much is reflective/shiny
class TableMaterial : public Material {
public:
    TableMaterial(double shine) {
        shininess = shine;
        diffuseAlbedo = vec3(0.5);
        specularAlbedo = vec3(0.5);
        //specularAlbedo = vec3(0.05, 0.05, 0.05);
        //specularAlbedo = vec3(0.5, 0.5, 0.5);
        //diffuseAlbedo = vec3(0, 0, 0);
        //specularAlbedo = vec3(1.0, 1.0, 1.0);
    }
};