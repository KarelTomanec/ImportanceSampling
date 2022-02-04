#pragma once

#define _USE_MATH_DEFINES
#include <cassert>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <vector>

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "intersectable.h"
#include "scene.h"


// Global variable
Scene scene;

void
onInitialization()
{
    for (int Y = 0; Y < screenHeight; Y++) {
        for (int X = 0; X < screenWidth; X++) {
            reference[Y * screenWidth + X] = image[Y * screenWidth + X] = vec3(0, 0, 0);
        }
    }
    // Read the reference image from binary file
    FILE* refImage = fopen("image.bin", "rb");
    if (!refImage) {
        printf("No reference file\n");
    }
    else {
        int sz = fread(reference, sizeof(vec3), screenWidth * screenHeight, refImage);
        fclose(refImage);
        for (int Y = 0; Y < screenHeight; Y++) {
            for (int X = 0; X < screenWidth; X++) {
                image[Y * screenWidth + X] = reference[Y * screenWidth + X];
            }
        }
    }
    glViewport(0, 0, screenWidth, screenHeight);
    scene.build(); // create scene objects
    method = LIGHT_SOURCE; // the method to compute an image


}

void
getPseudocolorRainbow(double val, double minVal, double maxVal,
    double& r, double& g, double& b)
{
    if (isnan(val) || isinf(val)) {
        r = g = b = 0; // black ... exception
        return;
    }
    if (val < minVal) val = minVal;
    if (val > maxVal) val = maxVal;
    double ratio = (val - minVal) / (maxVal - minVal);
    double value = 1.0f - ratio;
    float val4 = value * 4.0f;
    value = val4 - (int)val4;
    switch ((int)(val4)) {
    case 0: r = 1.0; g = value; b = 0.f; break;
    case 1: r = 1.0 - value; g = 1.0; b = 0.f; break;
    case 2: r = 0.f; g = 1.0; b = value; break;
    case 3: r = 0.f; g = 1.0 - value; b = 1.0; break;
    default: r = value * 1.0; g = 0.f; b = 1.0; break;
    }
    return;
}

void
getPseudocolorCoolWarm(double val, double minVal, double maxVal,
    double& r, double& g, double& b)
{
    if (isnan(val) || isinf(val)) {
        r = g = b = 0; // black ... exception
        return;
    }
    if (val < minVal) val = minVal;
    if (val > maxVal) val = maxVal;
    double ratio = (val - minVal) / (maxVal - minVal);
    int i = int(ratio * 31.999);
    assert(i < 33); assert(i >= 0);
    float alpha = i + 1.0 - (ratio * 31.999);
    r = pscols[4 * i + 1] * alpha + pscols[4 * (i + 1) + 1] * (1.0 - alpha);
    g = pscols[4 * i + 2] * alpha + pscols[4 * (i + 1) + 2] * (1.0 - alpha);
    b = pscols[4 * i + 3] * alpha + pscols[4 * (i + 1) + 3] * (1.0 - alpha);
    //printf("rgb=%f %f %f index=%d a=%g\n",r,g,b,i, alpha);
}

void
onDisplay()
{
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    static float displayImage[screenWidth * screenHeight * 3] = { 0 };
    // Image on the left
    for (int Y = 0; Y < screenHeight; Y++) {
        for (int X = 0; X < screenWidth; X++) {
            displayImage[3 * (Y * screenWidth + X)] = image[Y * screenWidth + X].x;
            displayImage[3 * (Y * screenWidth + X) + 1] = image[Y * screenWidth + X].y;
            displayImage[3 * (Y * screenWidth + X) + 2] = image[Y * screenWidth + X].z;
        }
    }
    glRasterPos2i(-1, -1);
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, displayImage);
    // Image on the right
    for (int Y = 0; Y < screenHeight; Y++) {
        for (int X = 0; X < screenWidth; X++) {
            if (showFlag == DIFF) {
                double diff = (image[Y * screenWidth + X] - reference[Y * screenWidth + X]).average();
                displayImage[3 * (Y * screenWidth + X)] = diff > 0 ? diff : 0;
                displayImage[3 * (Y * screenWidth + X) + 1] = diff < 0 ? -diff : 0;
                displayImage[3 * (Y * screenWidth + X) + 2] = 0;
            }
            if (showFlag == WEIGHT) { // black to white
                double w = weight[Y * screenWidth + X];
                displayImage[3 * (Y * screenWidth + X)] = w;
                displayImage[3 * (Y * screenWidth + X) + 1] = w;
                displayImage[3 * (Y * screenWidth + X) + 2] = w;
            }
            if (showFlag == WEIGHT_PSEUDOCOLOR) {
                double w = weight[Y * screenWidth + X];
                if (showBargraph && (X > 0.98 * screenWidth))
                    w = (double)Y / screenHeight; // thin bar on the right showing the mapping
                double r, g, b;
                if (rainbowPSC)
                    getPseudocolorRainbow(w, 0.0, 1.0, r, g, b); // is more common but wrong perceptually
                else
                    getPseudocolorCoolWarm(w, 0.0, 1.0, r, g, b); // is perceptually better
                displayImage[3 * (Y * screenWidth + X)] = r;
                displayImage[3 * (Y * screenWidth + X) + 1] = g;
                displayImage[3 * (Y * screenWidth + X) + 2] = b;
            }
        }
    }
    glRasterPos2i(0, -1);
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, displayImage);
    glutSwapBuffers();

    // Save TGA file for the image, simple format
    FILE* ofile = 0;
    switch (method) {
    case BRDF:	ofile = fopen("brdf.tga", "wb"); break;
    case LIGHT_SOURCE:	ofile = fopen("lightsource.tga", "wb"); break;
    case INF_AREA: ofile = fopen("infarea.tga", "wb"); break;
    case MULTIPLE: ofile = fopen("multiple.tga", "wb"); break;
    }
    if (!ofile) return;

    fputc(0, ofile); fputc(0, ofile); fputc(2, ofile);
    for (int i = 3; i < 12; i++) { fputc(0, ofile); }
    int width = screenWidth * 2, height = screenHeight;
    fputc(width % 256, ofile); fputc(width / 256, ofile);
    fputc(height % 256, ofile); fputc(height / 256, ofile);
    fputc(24, ofile);
    fputc(32, ofile);

    for (int Y = screenHeight - 1; Y >= 0; Y--) {
        for (int X = 0; X < width; X++) {
            double r, g, b;
            if (X < screenWidth) {
                r = image[Y * screenWidth + X].x;
                g = image[Y * screenWidth + X].y;
                b = image[Y * screenWidth + X].z;
            }
            else {
                int XX = X - screenWidth;
                double w = weight[Y * screenWidth + XX];
                if (showBargraph && (XX > 0.98 * screenWidth))
                    w = (double)Y / screenHeight; // thin bar on the right showing the mapping
                if (rainbowPSC)
                    getPseudocolorRainbow(w, 0.0, 1.0, r, g, b); // is more common but wrong perceptually
                else
                    getPseudocolorCoolWarm(w, 0.0, 1.0, r, g, b); // is perceptually better
            }
            int R = fmax(fmin(r * 255.5, 255), 0);
            int G = fmax(fmin(g * 255.5, 255), 0);
            int B = fmax(fmin(b * 255.5, 255), 0);
            fputc(B, ofile); fputc(G, ofile); fputc(R, ofile);
        }
    }
    fclose(ofile);
}

void Usage() {
    printf("Usage:\n");
    printf(" 'b': BRDF sampling \n");
    printf(" 'l': light source sampling \n");
    printf(" 'r': Show reference\n");
    printf(" 'w': Print current as a ground truth reference file for the future\n");
    printf(" 'o': Write output HDR file of rendered image\n\n");
    printf(" 'O': Write output HDR file render+pseudocolor visualization\n\n");
}

// Mouse click event
void
onMouse(int button, int state, int pX, int pY)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        // select GLUT_LEFT_BUTTON / GLUT_RIGHT_BUTTON and GLUT_DOWN / GLUT_UP
        scene.testRay(pX, screenHeight - pY);
    }
}

void
onKeyboard(unsigned char key, int x, int y)
{
    switch (key) {
    case 'l':
        method = LIGHT_SOURCE;
        printf("Light source sampling\n");
        scene.render();
        break;
    case 'b':
        method = BRDF;
        printf("BRDF sampling\n");
        scene.render();
        break;
    case 'w': {
        printf("Writing reference file\n");
        FILE* refImage = fopen("image.bin", "wb");
        if (refImage) {
            fwrite(image, sizeof(vec3), screenWidth * screenHeight, refImage);
            fclose(refImage);
        }
    }
    case 'i': {
        method = INF_AREA;
        printf("Infinite area light source sampling\n");
        scene.render();
        break;
    }
    case 'm': {
        method = MULTIPLE;
        printf("Infinite area light source sampling\n");
        scene.render();
        break;
    }
    case 'O':
    case 'o': {
        printf("Writing output HDR file (extension .hdr)\n");
        FILE* fp = nullptr;
        switch (method) {
        case LIGHT_SOURCE: fp = fopen("lightsource.hdr", "wb"); break;
        case BRDF: fp = fopen("brdf.hdr", "wb"); break;
        case INF_AREA: fp = fopen("infarea.hdr", "wb"); break;
        case MULTIPLE: fp = fopen("multiple.hdr", "wb"); break;
        }
        int width = screenWidth;
        bool psf = false;
        if (key == 'O') { width *= 2; psf = true; }
        if (fp) {
            size_t nmemb = width * screenHeight;
            typedef unsigned char RGBE[4];
            RGBE* data = new RGBE[nmemb];
            for (int ii = 0; ii < nmemb; ii++) {
                RGBE& rgbe = data[ii];
                //int x = (ii % width);
                //int y = screenHeight - (ii / width) + 1;
                int x = (ii % width);
                int y = screenHeight - (ii / width) - 1;
                vec3 vv;
                //std::cout << x << " " << y << std::endl;
                vv = image[y * screenWidth + x];
                if (psf) {
                    if (x < screenWidth) {
                        vv = image[y * screenWidth + x];
                    }
                    else {
                        x -= screenWidth;
                        double w = weight[y * screenWidth + x];
                        if (showBargraph && (x > 0.98 * screenWidth))
                            w = (double)y / screenHeight; // thin bar on the right showing the mapping
                        if (rainbowPSC)
                            getPseudocolorRainbow(w, 0.0, 1.0, vv.x, vv.y, vv.z); // is more common but wrong perceptually
                        else
                            getPseudocolorCoolWarm(w, 0.0, 1.0, vv.x, vv.y, vv.z); // is perceptually better
                    }
                }
                float v; int e;
                v = vv.x; if (vv.y > v) v = vv.y; if (vv.z > v) v = vv.z;
                if (v < 1e-32) {
                    rgbe[0] = rgbe[1] = rgbe[2] = rgbe[3] = 0x0;
                }
                else {
                    v = (float)(frexp(v, &e) * 256.0 / v);
                    rgbe[0] = (unsigned char)(vv.x * v);
                    rgbe[1] = (unsigned char)(vv.y * v);
                    rgbe[2] = (unsigned char)(vv.z * v);
                    rgbe[3] = (unsigned char)(e + 128);
                }
            }
            fflush(stdout);
            const char* programtype = "RADIANCE";
            if (fprintf(fp, "#?%s\n", programtype) < 0) { abort(); }
            float gamma = 2.2; float exposure = 1.0;
            if (fprintf(fp, "GAMMA=%g\n", gamma) < 0) { abort(); }
            if (fprintf(fp, "EXPOSURE=%g\n", exposure) < 0) { abort(); }
            if (fprintf(fp, "FORMAT=32-bit_rle_rgbe\n\n") < 0) { abort(); }
            if (fprintf(fp, "-Y %d +X %d\n", screenHeight, width) < 0) { abort(); }
            // Write data
            size_t kk = fwrite(data, (size_t)4, nmemb, fp);
            fclose(fp);
            if (kk != nmemb) {
                printf("ERROR - was not able to save all kk=%d entries to file, exiting\n",
                    (int)nmemb); fflush(stdout);
                abort(); // error
            }
        }
    }
    } // switch (key)
    Usage();
    glutPostRedisplay();
}

int
main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(2 * screenWidth, screenHeight);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("DEMO RSO 2018");

    Usage();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    onInitialization();

    glutDisplayFunc(onDisplay);
    glutKeyboardFunc(onKeyboard);
    glutMouseFunc(onMouse);
    glutMainLoop();
    return 0;
}
