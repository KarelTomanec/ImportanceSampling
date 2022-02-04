#pragma once
#define STB_IMAGE_IMPLEMENTATION  
#define STB_IMAGE_WRITE_IMPLEMENTATION  
#include "stb_image.h"
#include "stb_image_write.h"

#include <memory>

#include <algorithm>

#include <iostream>

class img {
public:
    //color* data;
    std::unique_ptr<color[]> data;
    float* raw_data;
    //unsigned char* raw_data
    int width, height;

	const static int bytes_per_pixel = 3;

    img(const char* filename) {
        auto components_per_pixel = bytes_per_pixel;
        stbi_ldr_to_hdr_gamma(1.0f);
        raw_data = stbi_loadf(
            filename, &width, &height, &components_per_pixel, 0);

        if (!raw_data) {
            std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
            width = height = 0;
            return;
        }

        //data = new color[width * height];
        data = std::make_unique<color[]>(width * height);

        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                data[i * width + j] = color(
                    raw_data[3 * i * width + j * 3], 
                    raw_data[3 * i * width + j * 3 + 1], 
                    raw_data[3 * i * width + j * 3 + 2]
                );
            }
        }

        //delete[] raw_data;
    }


    void write(float u, float v, color c) {
        //std::cout << u << " " << v << std::endl;
        assert(u <= 1);
        assert(v <= 1);
        data[this->width * int((float)(this->height - 1) * v) + int((float)(this->width - 1) * u)] = c;
    }

    void save(const char* filename) {
        std::copy(reinterpret_cast<double*>(data.get()), reinterpret_cast<double*>(data.get()) + width * height * 3, raw_data);
        stbi_write_hdr(filename, width, height, bytes_per_pixel, raw_data);
    }

    color lookup(int x, int y) {
        return data[this->width * y + x];
    }

    color lookupUV(float u, float v) {
        //std::cout << u << " " << v << std::endl;
        return data[this->width * int((this->height - 1) * v) + int((this->width - 1) * u)];
    }
};

