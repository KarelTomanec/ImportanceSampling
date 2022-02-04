#pragma once
#include <vector>
#include <math.h>
#include "camera.h"
#include "light.h"

int nIterations = 1; // how many iterations to render
int nTotalSamples = 32; // samples in one render iteration - should be even number

// Which sampling method should be used
enum Method { BRDF, LIGHT_SOURCE, INF_AREA, MULTIPLE } method;

// The scene definition with main rendering method
class Scene {
private:
    std::vector<Intersectable*> objects;
    double totalPower;
    int nLightSamples, nBRDFSamples;
    Camera camera;
    InfiniteAreaLight infAreaLight;

public:


    void build() {
        // Create a simple scene
        vec3 eyePos(0, 6, 18); // camera center
        vec3 lightCenterPos(0, 4, -6); // first light source

        // Create geometry - 4 rectangles
        objects.push_back(new Rect(vec3(0, -4, +2), eyePos, lightCenterPos, 8, 1, new TableMaterial(500)));
        objects.push_back(new Rect(vec3(0, -3.5, -2), eyePos, lightCenterPos, 8, 1, new TableMaterial(1000)));
        objects.push_back(new Rect(vec3(0, -2.5, -6), eyePos, lightCenterPos, 8, 1, new TableMaterial(5000)));
        objects.push_back(new Rect(vec3(0, -1, -10), eyePos, lightCenterPos, 8, 1, new TableMaterial(10000)));

        //TableMaterial* mat = new TableMaterial(5000);
        //objects.push_back(new Sphere(vec3(0, -2.5, -6), 4, mat, false));
        // Create 4 light sources
        objects.push_back(new Sphere(lightCenterPos + vec3(-4.5, 0, 0), 0.07, new LightMaterial(vec3(4, 2, 1)), true));
        objects.push_back(new Sphere(lightCenterPos + vec3(-1.5, 0, 0), 0.16, new LightMaterial(vec3(2, 4, 1)), true));
        objects.push_back(new Sphere(lightCenterPos + vec3(1.5, 0, 0), 0.4, new LightMaterial(vec3(2, 1, 4)), true));
        objects.push_back(new Sphere(lightCenterPos + vec3(4.5, 0, 0), 1, new LightMaterial(vec3(4, 1, 2)), true));

        // Set the camera
        camera.set(eyePos, vec3(0, 0, 0), vec3(0, 1, 0), 35.0 * M_PI / 180.0);

        totalPower = 0;
        for (int i = 0; i < objects.size(); i++) {
            totalPower += objects[i]->power; //  hit.t < 0 if no intersection
        }


        // Load infinite area light source image
        infAreaLight = InfiniteAreaLight("../Resources/EnvironmentMaps/raw023.hdr");
    }

    // Set the weight for the sampling method
    void setWeight(double wval) {
        for (int Y = 0; Y < screenHeight; Y++)
            for (int X = 0; X < screenWidth; X++)
                weight[Y * screenWidth + X] = wval;
    }

    // Render the scene
    void render() {
        // Total number of samples per pixel is: nIterators*nTotalSamples
        srand(1);
        char buffer[100];
        FILE* errorFile = 0;

        switch (method) {
        case BRDF:
            nBRDFSamples = nTotalSamples;
            nLightSamples = 0;
            errorFile = fopen("BRDF.txt", "w");
            setWeight(0.0);
            break;
        case LIGHT_SOURCE:
            nBRDFSamples = 0;
            nLightSamples = nTotalSamples;
            errorFile = fopen("light.txt", "w");
            setWeight(1.0);
            break;
        case INF_AREA:
            nBRDFSamples = 0;
            nLightSamples = nTotalSamples;
            errorFile = fopen("light.txt", "w");
            setWeight(1.0);
            break;
        case MULTIPLE:
            nBRDFSamples = nTotalSamples / 2;
            nLightSamples = nTotalSamples / 2;
            errorFile = fopen("light.txt", "w");
            setWeight(0.5);
            break;
        } // switch

        double cost = 0;
        bool debug = false;
        // How many iterations
        for (int iIter = 1; iIter <= nIterations; iIter++) {
            double error = 0;
            for (int Y = 0; Y < screenHeight; Y++) { // for all rows
                printf("%d\r", Y);
                //#pragma omp parallel for schedule(dynamic)
                for (int X = 0; X < screenWidth; X++) { // for all pixels in a row
                    //if (debug) { // debug particular pixel x,y, coordinates from pfsv (pfstools)
                    //    X = 307;
                    //    Y = screenHeight - 517;
                    //}

                    nLightSamples = (int)(weight[Y * screenWidth + X] * nTotalSamples + 0.5);
                    nBRDFSamples = nTotalSamples - nLightSamples;
                    cost += nBRDFSamples * costBRDF + nLightSamples * costLight;

                    // For a primary ray at pixel (X,Y) compute the color
                    vec3 color = trace(camera.getRay(X, Y));
                    double w = 1.0 / iIter; // the same weight for all samples for computing mean incrementally
                    image[Y * screenWidth + X] = color * w + image[Y * screenWidth + X] * (1.0 - w);

                    w = 1.0 / sqrt(iIter); // emphasize later samples
                    vec3 diff = reference[Y * screenWidth + X] - image[Y * screenWidth + X];
                    error += dot(diff, diff);
                    if (debug)
                        break;
                } // for X
                if (debug)
                    break;
            } // for Y
            double eff = 100000.0 * nIterations * nTotalSamples * screenWidth * screenHeight / error / cost;
            printf("Iter: %d, Error: %4.2f, Efficiency: %f, Relative Efficiency: %f\n",
                iIter, sqrt(error), eff, eff / referenceEfficiency);
            fprintf(errorFile, "%d, %f\n", iIter * nTotalSamples, sqrt(error));
        } // for iTer
        fclose(errorFile);

        infAreaLight.radianceMap->save("debug.hdr"); // debug
    } // render

    // Compute intersection between a rady and primitive
    Hit firstIntersect(const Ray& ray, Intersectable* skip) {
        Hit bestHit;
        for (int i = 0; i < objects.size(); i++) {
            if (objects[i] == skip) continue;
            //if (method == INF_AREA) {
                //if (i >= 4 && i <= 7) {
                //    continue;
                //}
            //}
            Hit hit = objects[i]->intersect(ray); //  hit.t < 0 if no intersection
            if (hit.t > epsilon) {
                if (bestHit.t < 0 || hit.t < bestHit.t) bestHit = hit;
            }
        }
        return bestHit;
    }

    // Sample the light source from all the light sources in the scene
    Light* sampleLightSource(const vec3& illuminatedPoint, vec3& direction, float& pdf) // the 3D point on an object
    {
        //if (method == LIGHT_SOURCE || method == BRDF || method == MULTIPLE) {
            while (true) { // if no light source is selected due to floating point inaccuracies, repeat
                double threshold = totalPower * drandom();
                double running = 0;
                for (int i = 0; i < objects.size(); i++) {
                    running += objects[i]->power; // select light source with the probability proportional to its power
                    if (running > threshold) {
                        Sphere* sphere = (Sphere*)objects[i];
                        vec3 point, normal;
                        // select a point on the visible half of the light source
                        ((Sphere*)objects[i])->sampleUniformPoint(illuminatedPoint, point, normal);
                        direction = point - illuminatedPoint;
                        return new SphereLight(sphere, point, normal);
                    } // if
                } // for i
            } // for ever
        //}
        //else {
        //    infAreaLight.sample_l(illuminatedPoint, &direction, &pdf);
        //    return  &infAreaLight;
        //}
    }

    // Trace a primary ray towards the scene
    vec3 trace(const Ray& r) {
        // error measures for the two combined techniques: used for adaptation
        Hit hit = firstIntersect(r, NULL);	// find visible point
        if (hit.t < 0) {
            //return infAreaLight.get_illumination_environment(r.dir.normalize());
            return vec3(0, 0, 0);
        }
        // The energy emanated from the material
        vec3 radianceEmitted = hit.material->Le;
        if (hit.material->diffuseAlbedo.average() < epsilon &&
            hit.material->specularAlbedo.average() < epsilon)
            return radianceEmitted; // if albedo is low, no energy can be reefleted
          // Compute the contribution of reflected lgiht
        vec3 radianceBRDFSampling(0, 0, 0), radianceLightSourceSampling(0, 0, 0);
        vec3 inDir = -r.dir;	// incident direction

        int nTotalSamples = (nLightSamples + nBRDFSamples);
        double alpha = (double)nLightSamples / nTotalSamples;

        // The direct illumination for chosen number of samples
        for (int i = 0; i < nLightSamples; i++) {

            vec3 outDir;
            float pdf;
            Light* lightSample = sampleLightSource(hit.position, outDir, pdf); // generate a light sample
            double distance2 = dot(outDir, outDir);
            double distance = sqrt(distance2);
            if (distance < epsilon) {
                continue;
            }
            outDir = outDir.normalize(); // normalize the direction
            double cosThetaLight = dot(lightSample->normal, -outDir);
            if (cosThetaLight > epsilon) {
                // visibility is not needed to handle, all lights are visible

                double pdfLightSourceSampling =
                    lightSample->pdf(hit.position, outDir, totalPower) * distance2 / cosThetaLight;
                double pdfBRDFSampling = hit.material->sampleProb(hit.normal, inDir, outDir);
                // the theta angle on the surface between normal and light direction
                double cosThetaSurface = dot(hit.normal, outDir);
                if (cosThetaSurface > 0) {
                    // yes, the light is visible and contributes to the output power
                    // The evaluation of rendering equation locally: (light power) * brdf * cos(theta)

                    vec3 brdf = hit.material->BRDF(hit.normal, inDir, outDir);

                    vec3 f = lightSample->get_illumination(outDir) * cosThetaSurface * brdf;
                    double pdfBRDFSampling = hit.material->sampleProb(hit.normal, inDir, outDir);
                    double p = pdfLightSourceSampling + pdfBRDFSampling;
                    // importance sample = 1/n . \sum (f/prob)                    
                    radianceLightSourceSampling +=  (f / p / nTotalSamples);
                } // if
            }
        } // for all the samples from light

        for (int i = 0; i < nBRDFSamples; i++) {
            // BRDF.cos(theta) sampling should be implemented first!
            vec3 outDir;
            // BRDF sampling with Russian roulette
            if (hit.material->sampleDirection(hit.normal, inDir, outDir)) {
                double pdfBRDFSampling = hit.material->sampleProb(hit.normal, inDir, outDir);
                double cosThetaSurface = dot(hit.normal, outDir);
                if (cosThetaSurface > 0) {
                    vec3 brdf = hit.material->BRDF(hit.normal, inDir, outDir);
                    // Trace a ray to the scene
                    Hit lightSource = firstIntersect(Ray(hit.position, outDir), hit.object);
                    // Do we hit a light source
                    if (lightSource.t > 0 && lightSource.material->Le.average() > 0) {
                        // squared distance between an illuminated point and light source
                        double distance2 = lightSource.t * lightSource.t;
                        double cosThetaLight = dot(lightSource.normal, -outDir);
                        if (cosThetaLight > epsilon) {
                            double pdfLightSourceSampling =
                                lightSource.object->pointSampleProb(totalPower) * distance2 / cosThetaLight;
                            // The evaluation of rendering equation locally: (light power) * brdf * cos(theta)
                            vec3 f = lightSource.material->Le * brdf * cosThetaSurface;
                            double p = pdfBRDFSampling + pdfLightSourceSampling;
                            radianceBRDFSampling += f / p / nTotalSamples;
                        }
                        else
                            printf("ERROR: Sphere hit from back\n");
                    }
                }
            }
        } // for i


        return radianceEmitted + radianceLightSourceSampling + radianceBRDFSampling;
    }

    // Only testing routine for debugging
    void testRay(int X, int Y) {
        nBRDFSamples = nLightSamples = 1000;
        vec3 current = trace(camera.getRay(X, Y));
        printf("Pixel %d, %d Value = %f, %f, %f\n", X, Y, current.x, current.y, current.z);
    }
};