#include <cmath>
#include <fstream>
#include <iostream>
#include <GL/freeglut_std.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <list>
#include <sstream>
#include <string>
#include <variant>
#include <vector>
#include <math.h>
#include "TinyXML/tinyxml2.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

struct Point {
    float x, y, z;
};

struct Window {
    int width;
    int height;
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
    float x, y, z, time;
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

struct Models {
    std::vector<std::string> model;
    std::vector<Transformation> transformations;
    std::vector<Models> models;
};

struct Camera {
    Point position;
    Point lookAt;
    Point up;
    struct {
        float fov, near, far;
    } projection;
};

struct World {
    Window window;
    Camera camera;
    Models models;
};

World world;

float omega = 0.0f, alpha = 0.0f;
float k = 0.5;

Point calculateVector(Point L, Point P){
    return {L.x - P.x,
            L.y - P.y,
            L.z - P.z};
}

Point cross(const Point &a, const Point &b) {
    return { a.y * b.z - a.z * b.y, 
             a.z * b.x - a.x * b.z, 
             a.x * b.y - a.y * b.x };
}

void normalize(Point &p) {
    float len = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    if (len > 0) {
        p.x /= len;
        p.y /= len;
        p.z /= len;
    }
}

void SphericalToCartesianFPS() {
    Point direction;
    direction.x = cos(alpha) * cos(omega);
    direction.y = sin(-omega);
    direction.z = sin(alpha) * cos(omega);

    normalize(direction);

    world.camera.lookAt.x = world.camera.position.x + direction.x;
    world.camera.lookAt.y = world.camera.position.y + direction.y;
    world.camera.lookAt.z = world.camera.position.z + direction.z;
}



void key_press(unsigned char key, int x, int y) {
    Point d = calculateVector(world.camera.lookAt, world.camera.position);
    if(key == 'e'){
        alpha += 0.02f;
        SphericalToCartesianFPS();
    }
    if(key == 'q'){
        alpha -= 0.02f;
        SphericalToCartesianFPS();
    }
    if(key == 'x'){
        omega += 0.02f;
        SphericalToCartesianFPS();
    }
    if(key == 'z'){
        omega -= 0.02f;
        SphericalToCartesianFPS();
    }

    if(key == 'w'){
        world.camera.lookAt.x += k * d.x;
        world.camera.lookAt.y += k * d.y;
        world.camera.lookAt.z += k * d.z;
        world.camera.position.x += k * d.x;
        world.camera.position.y += k * d.y;
        world.camera.position.z += k * d.z;
    }
    if(key == 's'){ 
        world.camera.lookAt.x -= k * d.x;
        world.camera.lookAt.y -= k * d.y;
        world.camera.lookAt.z -= k * d.z;
        world.camera.position.x -= k * d.x;
        world.camera.position.y -= k * d.y;
        world.camera.position.z -= k * d.z;
    }
    if(key == 'd'){
        Point r = cross(d, world.camera.up);
        world.camera.lookAt.x += k * r.x;
        world.camera.lookAt.y += k * r.y;
        world.camera.lookAt.z += k * r.z;
        world.camera.position.x += k * r.x;
        world.camera.position.y += k * r.y;
        world.camera.position.z += k * r.z;
    }
    if(key == 'a'){ 
        Point r = cross(d, world.camera.up);
        world.camera.lookAt.x -= k * r.x;
        world.camera.lookAt.y -= k * r.y;
        world.camera.lookAt.z -= k * r.z;
        world.camera.position.x -= k * r.x;
        world.camera.position.y -= k * r.y;
        world.camera.position.z -= k * r.z;
    }
    if(key =='-'){
        k /=2;
        if(k < 0.1) k = 0.1;
    }
    if(key == '+'){
        k *=2;
        if(k > 25) k = 25;
    }

    glutPostRedisplay();
    
}

void drawAxis() {
    glBegin(GL_LINES);

    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, -50.0f, 0.0f);
    glVertex3f(0.0f, 150.0f, 0.0f);

    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-50.0f, 0.0f, 0.0f);
    glVertex3f(150.0f, 0.0f, 0.0f);

    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, -50.0f);
    glVertex3f(0.0f, 0.0f, 150.0f);

    glEnd();
    glColor3f(1.0f, 1.0f, 1.0f);
}

void drawFigure(string filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error when trying to open the file!" << endl;
        exit(0);
    }

    glBegin(GL_TRIANGLES);
    glColor3f(1.0f, 1.0f, 1.0f);
    string line;
    getline(file, line);
    while (getline(file, line)) {
        istringstream stream(line);
        float x, y, z;

        if (!(stream >> x >> y >> z)) {
            cerr << "Error when trying to read the values!" << endl;
            glEnd();
            file.close();
            return;
        }

        glVertex3f(x, y, z);
    }

    glEnd();
    file.close();
}

void multMatrixVector(float *m, float *v, float *res) {
    for (int j = 0; j < 4; ++j) {
        res[j] = 0;
        for (int k = 0; k < 4; ++k) {
            res[j] += v[k] * m[j * 4 + k];
        }
    }
}

void getCatmullRomPoint(float t, Point p0, Point p1, Point p2, Point p3, Point *pos, Point *deriv) {
   
    float m[16] = { -0.5f, 1.5f, -1.5f, 0.5f,
                    1.0f, -2.5f, 2.0f, -0.5f,
                    -0.5f, 0.0f, 0.5f, 0.0f,
                    0.0f, 1.0f, 0.0f, 0.0f };

    float T[4] = { t*t*t, t*t, t, 1 };
    float dT[4] = { 3*t*t, 2*t, 1, 0 };

  
    float px[4] = { p0.x, p1.x, p2.x, p3.x };
    float py[4] = { p0.y, p1.y, p2.y, p3.y };
    float pz[4] = { p0.z, p1.z, p2.z, p3.z };

    float a_x[4], a_y[4], a_z[4];
    multMatrixVector(m, px, a_x);
    multMatrixVector(m, py, a_y);
    multMatrixVector(m, pz, a_z);

    pos->x = T[0]*a_x[0] + T[1]*a_x[1] + T[2]*a_x[2] + T[3]*a_x[3];
    pos->y = T[0]*a_y[0] + T[1]*a_y[1] + T[2]*a_y[2] + T[3]*a_y[3];
    pos->z = T[0]*a_z[0] + T[1]*a_z[1] + T[2]*a_z[2] + T[3]*a_z[3];

   
    deriv->x = dT[0]*a_x[0] + dT[1]*a_x[1] + dT[2]*a_x[2] + dT[3]*a_x[3];
    deriv->y = dT[0]*a_y[0] + dT[1]*a_y[1] + dT[2]*a_y[2] + dT[3]*a_y[3];
    deriv->z = dT[0]*a_z[0] + dT[1]*a_z[1] + dT[2]*a_z[2] + dT[3]*a_z[3];
}

void getGlobalCatmullRomPoint(float gt, const vector<Point>& points, Point *pos, Point *deriv) {
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
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_LINE_LOOP);

    for (float gt = 0; gt < 1; gt += 0.01) {
        Point pos, deriv;
        getGlobalCatmullRomPoint(gt, points, &pos, &deriv);
        glVertex3f(pos.x, pos.y, pos.z);
    }

    glEnd();
}

void buildRotMatrix(Point x, Point y, Point z, float *m) {
    m[0] = x.x; m[1] = x.y; m[2] = x.z; m[3] = 0;
    m[4] = y.x; m[5] = y.y; m[6] = y.z; m[7] = 0;
    m[8] = z.x; m[9] = z.y; m[10] = z.z; m[11] = 0;
    m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}

void applyTransformation(const Transformation& transformation) {
    if (const auto* trans = get_if<Anime_Translate>(&transformation.type)) {
        float elapsed = glutGet(GLUT_ELAPSED_TIME) / 1000.0f - transformation.start_time;
        float gt = fmod(elapsed, trans->time) / trans->time;
        
        Point pos, deriv;
        getGlobalCatmullRomPoint(gt, trans->points, &pos, &deriv);
        
        glTranslatef(pos.x, pos.y, pos.z);
        
        if (trans->align) {
            Point X = deriv;
            normalize(X);
            
            Point Y = {0.0f, 1.0f, 0.0f};
            Point Z = cross(X, Y);
            normalize(Z);
            
            Y = cross(Z, X);
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

void drawModel(Models& models) {
    
    for (const auto& transformation : models.transformations) {
        if (holds_alternative<Anime_Translate>(transformation.type)) {
            const auto& trans = get<Anime_Translate>(transformation.type);
            drawCatmullRomCurve(trans.points);
        }
    }

    glPushMatrix();

    
    for (const auto& transformation : models.transformations) {
        applyTransformation(transformation);
    }

   
    for (const string& filename : models.model) {
        drawFigure(filename);
    }

   
    for (Models& child : models.models) {
        drawModel(child);
    }

    glPopMatrix();
}

void changeSize(int w, int h) {
    if (h == 0) h = 1;
    float ratio = w * 1.0 / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(world.camera.projection.fov, ratio, 
                  world.camera.projection.near, world.camera.projection.far);
    glMatrixMode(GL_MODELVIEW);
}

void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(world.camera.position.x, world.camera.position.y, world.camera.position.z,
              world.camera.lookAt.x, world.camera.lookAt.y, world.camera.lookAt.z,
              world.camera.up.x, world.camera.up.y, world.camera.up.z);

    glPolygonMode(GL_FRONT, GL_LINE);
    drawAxis();
    drawModel(world.models);
    
    glutSwapBuffers();
}

void parseGroup(tinyxml2::XMLElement *groupElement, Models &models) {
    using namespace tinyxml2;

    
    XMLElement *transformElement = groupElement->FirstChildElement("transform");
    if (transformElement) {
        for (XMLElement *child = transformElement->FirstChildElement(); child; child = child->NextSiblingElement()) {
            string name_trans = child->Name();

            if (name_trans == "translate") {
                if (child->Attribute("time")) {
                    
                    Anime_Translate at;
                    child->QueryFloatAttribute("time", &at.time);
                    const char* align = child->Attribute("align");
                    at.align = align && string(align) == "true";

                    for (XMLElement *point = child->FirstChildElement("point"); point; point = point->NextSiblingElement("point")) {
                        Point p;
                        point->QueryFloatAttribute("x", &p.x);
                        point->QueryFloatAttribute("y", &p.y);
                        point->QueryFloatAttribute("z", &p.z);
                        at.points.push_back(p);
                    }
                    models.transformations.push_back({at});
                } else {
                   
                    Translate t;
                    child->QueryFloatAttribute("x", &t.x);
                    child->QueryFloatAttribute("y", &t.y);
                    child->QueryFloatAttribute("z", &t.z);
                    models.transformations.push_back({t});
                }
            }
            else if (name_trans == "rotate") {
                if (child->Attribute("time")) {
                    
                    Anime_Rotate ar;
                    child->QueryFloatAttribute("time", &ar.time);
                    child->QueryFloatAttribute("x", &ar.x);
                    child->QueryFloatAttribute("y", &ar.y);
                    child->QueryFloatAttribute("z", &ar.z);
                    models.transformations.push_back({ar});
                } else {
                    
                    Rotate r;
                    child->QueryFloatAttribute("angle", &r.angle);
                    child->QueryFloatAttribute("x", &r.x);
                    child->QueryFloatAttribute("y", &r.y);
                    child->QueryFloatAttribute("z", &r.z);
                    models.transformations.push_back({r});
                }
            }
            else if (name_trans == "scale") {
                Scale s;
                child->QueryFloatAttribute("x", &s.x);
                child->QueryFloatAttribute("y", &s.y);
                child->QueryFloatAttribute("z", &s.z);
                models.transformations.push_back({s});
            }
        }
    }

    
    XMLElement *modelsElement = groupElement->FirstChildElement("models");
    if (modelsElement) {
        for (XMLElement *model = modelsElement->FirstChildElement("model"); model; model = model->NextSiblingElement("model")) {
            const char *file = model->Attribute("file");
            if (file) models.model.push_back(file);
        }
    }


    for (XMLElement *childGroup = groupElement->FirstChildElement("group"); childGroup; childGroup = childGroup->NextSiblingElement("group")) {
        Models childModels;
        parseGroup(childGroup, childModels);
        models.models.push_back(childModels);
    }
}

void parseInfo(char *filename) {
    using namespace tinyxml2;
    XMLDocument doc;
    if (doc.LoadFile(filename) != XML_SUCCESS) {
        cerr << "Error loading XML file" << endl;
        exit(1);
    }

    XMLNode *pRoot = doc.FirstChild();
    if (!pRoot) {
        cerr << "Empty XML file" << endl;
        exit(1);
    }

    
    XMLElement *window = pRoot->FirstChildElement("window");
    if (window) {
        window->QueryIntAttribute("width", &world.window.width);
        window->QueryIntAttribute("height", &world.window.height);
    }

    
    XMLElement *camera = pRoot->FirstChildElement("camera");
    if (camera) {
        XMLElement *pos = camera->FirstChildElement("position");
        if (pos) {
            pos->QueryFloatAttribute("x", &world.camera.position.x);
            pos->QueryFloatAttribute("y", &world.camera.position.y);
            pos->QueryFloatAttribute("z", &world.camera.position.z);
        }

        XMLElement *lookAt = camera->FirstChildElement("lookAt");
        if (lookAt) {
            lookAt->QueryFloatAttribute("x", &world.camera.lookAt.x);
            lookAt->QueryFloatAttribute("y", &world.camera.lookAt.y);
            lookAt->QueryFloatAttribute("z", &world.camera.lookAt.z);
        }

        XMLElement *up = camera->FirstChildElement("up");
        if (up) {
            up->QueryFloatAttribute("x", &world.camera.up.x);
            up->QueryFloatAttribute("y", &world.camera.up.y);
            up->QueryFloatAttribute("z", &world.camera.up.z);
        }

        XMLElement *proj = camera->FirstChildElement("projection");
        if (proj) {
            proj->QueryFloatAttribute("fov", &world.camera.projection.fov);
            proj->QueryFloatAttribute("near", &world.camera.projection.near);
            proj->QueryFloatAttribute("far", &world.camera.projection.far);
        }
    }

    
    float radius = sqrt(pow(world.camera.lookAt.x - world.camera.position.x, 2) + 
             pow(world.camera.lookAt.y - world.camera.position.y, 2) + 
             pow(world.camera.lookAt.z - world.camera.position.z, 2));
    omega = asin(world.camera.position.y / radius);
    alpha = atan2(world.camera.position.x, world.camera.position.z) + M_PI;

    XMLElement *group = pRoot->FirstChildElement("group");
    if (group) parseGroup(group, world.models);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <scene_file.xml>" << endl;
        return 1;
    }

    parseInfo(argv[1]);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(world.window.width, world.window.height);
    glutCreateWindow("CG@DI Phase 3");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutKeyboardFunc(key_press);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT, GL_LINE);

    glutMainLoop();

    return 0;
}
