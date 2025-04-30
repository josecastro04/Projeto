#include "transformations.h"
#include <vector>
using namespace std;
#include <cmath>
#include <GL/glut.h>


void normalize(float *p) {
    float len = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    if (len > 0) {
        p[0] /= len;
        p[1] /= len;
        p[2] /= len;
    }
}
void cross(float *a, float *b, float *res) {
   res[0] = a[1] * b[2] - a[2] * b[1];
   res[1] = a[2] * b[0] - a[0] * b[2];
   res[2] = a[0] * b[1] - a[1] * b[0] ;
    
}
void buildRotMatrix(float *x, float *y, float *z, float *m) {
    m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
    m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
    m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
    m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}
void multMatrixVector(float *m, float *v, float *res) {
    for (int j = 0; j < 4; ++j) {
        res[j] = 0;
        for (int k = 0; k < 4; ++k) {
            res[j] += v[k] * m[j * 4 + k];
        }
    }
}
void getCatmullRomPoint(float t, const Point& p0, const Point& p1, const Point& p2, const Point& p3, float *pos, float *deriv){
   
    float m[16] = { -0.5f, 1.5f, -1.5f, 0.5f,
                    1.0f, -2.5f, 2.0f, -0.5f,
                    -0.5f, 0.0f, 0.5f, 0.0f,
                    0.0f, 1.0f, 0.0f, 0.0f };

    float a[4];
    float tt[4] = { t*t*t, t*t, t, 1 };
    float td[4] = { 3*t*t, 2*t, 1, 0 };

    float p[4][3] = {
        {p0.x, p0.y, p0.z},
        {p1.x, p1.y, p1.z},
        {p2.x, p2.y, p2.z},
        {p3.x, p3.y, p3.z}
    };

    for (int i = 0; i < 3; i++)
	{
		float pp[4] = {p[0][i], p[1][i], p[2][i], p[3][i]};
		multMatrixVector((float *)m, pp, a);
		pos[i] = tt[0] * a[0] +
				 tt[1] * a[1] +
				 tt[2] * a[2] +
				 tt[3] * a[3];
		deriv[i] = td[0] * a[0] +
				   td[1] * a[1] +
				   td[2] * a[2] +
				   td[3] * a[3];
	}
}

void getGlobalCatmullRomPoint(float gt, const vector<Point>& points, float *pos, float *deriv) {
    int point_count = points.size();
    float t = gt * point_count;
    int index = floor(t);
    t = t - index;

    int indices[4];
    indices[0] = (index + point_count - 1) % point_count;
    indices[1] = (indices[0] + 1) % point_count;
    indices[2] = (indices[1] + 1) % point_count;
    indices[3] = (indices[2] + 1) % point_count;

    getCatmullRomPoint(t, points[indices[0]], points[indices[1]], 
                       points[indices[2]], points[indices[3]], pos, deriv);
}

void drawCatmullRomCurve(const vector<Point>& points) {
   float pos[3], deriv[3];
    glBegin(GL_LINE_LOOP);

    for (float gt = 0; gt < 1; gt += 0.01) {
        
        getGlobalCatmullRomPoint(gt, points, pos, deriv);
        glVertex3f(pos[0], pos[1], pos[2]);
    }

    glEnd();
}

void applyTransformation(Transformation& transformation) {
    if (const auto* trans = get_if<Anime_Translate>(&transformation.type)) {
        float new_time = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
        float elapsed = new_time - transformation.start_time;
        float gt = fmod(elapsed, trans->time) / trans->time;
        
        float pos[3], deriv[3];
        getGlobalCatmullRomPoint(gt, trans->points, pos, deriv);
        
        glTranslatef(pos[0], pos[1], pos[2]);
        
        if (trans->align) {
            float X[3] = { deriv[0], deriv[1], deriv[2] };

            normalize(X);
            
            float Y[3] = { 0.0f, 1.0f, 0.0f };
            float Z[3];
            cross(X, Y,Z);
            
            normalize(Z);
            
            cross(Z, X, Y);
            normalize(Y);
            
            float m[16];
            buildRotMatrix(X, Y, Z, m);
            glMultMatrixf(m);

        
        }
        
     
    } 
    else if (const auto* trans = get_if<Translate>(&transformation.type)) {
        glTranslatef(trans->x, trans->y, trans->z);
    }
    else if (const auto* trans = get_if<Anime_Rotate>(&transformation.type)) {
        float elapsed = glutGet(GLUT_ELAPSED_TIME) / 1000.0f - transformation.start_time;
        float angle = 360.0f * fmod(elapsed, trans->time) / trans->time;
        glRotatef(angle, trans->x, trans->y, trans->z);
    }
    else if (const auto* trans = get_if<Rotate>(&transformation.type)) {
        glRotatef(trans->angle, trans->x, trans->y, trans->z);
    }
    else if (const auto* trans = get_if<Scale>(&transformation.type)) {
        glScalef(trans->x, trans->y, trans->z);
    }
}


