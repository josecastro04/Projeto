#pragma once

#include <vector>
#include <variant>
#include <string>
#include <GL/glut.h>

struct Point {
    float x, y, z;
};

struct Translate {
    float x, y, z;
};

struct Rotate {
    float x, y, z, angle;
};

struct Scale {
    float x, y, z;
};

struct Anime_Rotate {
    float x, y, z;
    float time;
};

struct Anime_Translate {
    float time;
    bool align;
    std::vector<Point> points;
};

struct Transformation {
    std::variant<Translate, Anime_Translate, Rotate, Anime_Rotate, Scale> type;
    bool is_active = true;
    float start_time = 0.0f;
};

void normalize(float* p);
void cross(float* a, float* b, float* res);
void buildRotMatrix(float* x, float* y, float* z, float* m);
void multMatrixVector(float* m, float* v, float* res);

void getCatmullRomPoint(float t, const Point& p0, const Point& p1, const Point& p2, const Point& p3, float* pos, float* deriv);
void getGlobalCatmullRomPoint(float gt, const std::vector<Point>& points, float* pos, float* deriv);
void drawCatmullRomCurve(const std::vector<Point>& points);

void applyTransformation(Transformation& transformation);
