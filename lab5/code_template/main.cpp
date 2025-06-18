#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

struct Vec3 {
    double x, y, z;
};

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3& operator+=(Vec3& a, const Vec3& b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}
Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 operator*(double s, const Vec3& a) {
    return {s * a.x, s * a.y, s * a.z};
}
Vec3 operator*(const Vec3& a, double s) {
    return s * a;
}
double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 multiply(const double M[3][3], const Vec3& v) {
    return {
        M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z,
        M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z,
        M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z
    };
}

const double RGB_to_LMS_matrix[3][3] = {
    {0.3811, 0.5783, 0.0402},
    {0.1967, 0.7244, 0.0782},
    {0.0241, 0.1288, 0.8444}
};

const double LMS_to_lab_matrix[3][3] = {
    {1/sqrt(3.0), 0, 0},
    {0, 1/sqrt(6.0), 0},
    {0, 0, 1/sqrt(2.0)}
};

const double lab_to_LMS_matrix[3][3] = {
    {1, 1, 1},
    {1, 1, -2},
    {1, -1, 0}
};

const double LMS_to_RGB_matrix[3][3] = {
    {4.4679, -3.5873, 0.1193},
    {-1.2186, 2.3809, -0.1624},
    {0.0497, -0.2439, 1.2045}
};

const double LMS_to_lab_matrix_2[3][3] = {
    {1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)},
    {1/sqrt(6.0), 1/sqrt(6.0), -2/sqrt(6.0)},
    {1/sqrt(2.0), -1/sqrt(2.0), 0}
};


Vec3 rgb_to_lab(const Vec3& rgb) {
    Vec3 lms = multiply(RGB_to_LMS_matrix, rgb);
    lms.x = (lms.x > 0) ? log10(lms.x) : -10;
    lms.y = (lms.y > 0) ? log10(lms.y) : -10;
    lms.z = (lms.z > 0) ? log10(lms.z) : -10;
    
    Vec3 lab;
    lab.x = (1.0 / sqrt(3.0)) * (lms.x + lms.y + lms.z);
    lab.y = (1.0 / sqrt(6.0)) * (lms.x + lms.y - 2.0 * lms.z);
    lab.z = (1.0 / sqrt(2.0)) * (lms.x - lms.y);

    return lab;
}

Vec3 lab_to_rgb(const Vec3& lab) {
    Vec3 lms;
    lms.x = (1.0 / sqrt(3.0)) * lab.x + (1.0 / sqrt(6.0)) * lab.y + (1.0 / sqrt(2.0)) * lab.z;
    lms.y = (1.0 / sqrt(3.0)) * lab.x + (1.0 / sqrt(6.0)) * lab.y - (1.0 / sqrt(2.0)) * lab.z;
    lms.z = (1.0 / sqrt(3.0)) * lab.x - (2.0 / sqrt(6.0)) * lab.y;

    lms.x = pow(10, lms.x);
    lms.y = pow(10, lms.y);
    lms.z = pow(10, lms.z);

    return multiply(LMS_to_RGB_matrix, lms);
}


int main() {
	int W, H, C, W_model, H_model, C_model;
	
	unsigned char *input_img_data = stbi_load("imgA.jpg", &W, &H, &C, STBI_rgb);
	if (!input_img_data) {
        std::cerr << "Error loading input image." << std::endl;
        return 1;
    }
	unsigned char *model_img_data = stbi_load("redim.jpg", &W_model, &H_model, &C_model, STBI_rgb);
	if (!model_img_data) {
        std::cerr << "Error loading model image." << std::endl;
        stbi_image_free(input_img_data);
        return 1;
    }

    if (W != W_model || H != H_model) {
        std::cerr << "Images must have the same dimensions." << std::endl;
        stbi_image_free(input_img_data);
        stbi_image_free(model_img_data);
        return 1;
    }

	std::vector<Vec3> input_lab(W * H);
    for (int i = 0; i < W * H; ++i) {
        Vec3 rgb = {input_img_data[3 * i] / 255.0, input_img_data[3 * i + 1] / 255.0, input_img_data[3 * i + 2] / 255.0};
        input_lab[i] = rgb_to_lab(rgb);
    }

    std::vector<Vec3> model_lab(W * H);
    for (int i = 0; i < W * H; ++i) {
        Vec3 rgb = {model_img_data[3 * i] / 255.0, model_img_data[3 * i + 1] / 255.0, model_img_data[3 * i + 2] / 255.0};
        model_lab[i] = rgb_to_lab(rgb);
    }

    stbi_image_free(input_img_data);
    stbi_image_free(model_img_data);
	
    int nb_iter = 100;
    std::default_random_engine generator;
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    for (int iter = 0; iter < nb_iter; ++iter) {
        Vec3 v = {normal_dist(generator), normal_dist(generator), normal_dist(generator)};
        double norm = sqrt(dot(v, v));
        if (norm > 1e-9) {
            v.x /= norm;
            v.y /= norm;
            v.z /= norm;
        }

        std::vector<std::pair<double, int>> proj_input(W * H);
        std::vector<std::pair<double, int>> proj_model(W * H);

        for (int i = 0; i < W * H; ++i) {
            proj_input[i] = {dot(input_lab[i], v), i};
            proj_model[i] = {dot(model_lab[i], v), i};
        }

        std::sort(proj_input.begin(), proj_input.end());
        std::sort(proj_model.begin(), proj_model.end());

        for (int i = 0; i < W * H; ++i) {
            int original_idx = proj_input[i].second;
            double diff = proj_model[i].first - proj_input[i].first;
            input_lab[original_idx] += diff * v;
        }
    }
	
	std::vector<unsigned char> image_result(W * H * 3);
    for (int i = 0; i < W * H; ++i) {
        Vec3 rgb = lab_to_rgb(input_lab[i]);
        image_result[3 * i] = std::min(255, std::max(0, (int)(rgb.x * 255)));
        image_result[3 * i + 1] = std::min(255, std::max(0, (int)(rgb.y * 255)));
        image_result[3 * i + 2] = std::min(255, std::max(0, (int)(rgb.z * 255)));
    }
	
	stbi_write_png("color_transfer_result.png", W, H, 3, &image_result[0], W * 3);

	return 0;
}