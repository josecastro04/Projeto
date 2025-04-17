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

struct Point
{
    float x, y, z;
};

struct Window
{
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
struct Anime_Translate{
    float time;
    bool align;
    std::vector<Point> points;
};

using Transformation = std::variant<Translate, Anime_Translate, Rotate, Scale>;

struct Models {
    std::vector<std::string> model;
    std::vector<Transformation> transformations;
    std::vector<Models> models;
};


struct Camera
{
    Point position;
    Point lookAt;
    Point up;
    struct
    {
        float fov, near, far;
    } projection;
};

struct World
{
    Window window;
    Camera camera;
    Models models;
};

World world;

float omega = 0.0f, alpha = 0.0f, radius = 5.0f;
bool keyStates[256] = {false};

void key_press(unsigned char key, int x, int y)
{
    keyStates[key] = true;
}

void key_up(unsigned char key, int x, int y)
{
    keyStates[key] = false;
}

void SphericalToCartesian()
{
    world.camera.position.x = radius * cos(omega) * sin(alpha);
    world.camera.position.y = radius * sin(omega);
    world.camera.position.z = radius * cos(omega) * cos(alpha);
}

void drawAxis()
{
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


void drawFigure(string filename)
{
    ifstream file(filename);
    if (!file)
    {
        cerr << "Error when trying to open the file!" << endl;
        exit(0);
    }

    glBegin(GL_TRIANGLES);
    glColor3f(1.0f, 1.0f, 1.0f);
    string line;
    getline(file, line);
    while (getline(file, line))
    {
        istringstream stream(line);

        float x, y, z;

        if (!(stream >> x >> y >> z))
        {
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
void getCatmullRomPoint(float t, Point p0, Point p1, Point p2, Point p3, Point &res, Point &deriv) {
    float m[4][4] = {
        {-0.5f,  1.5f, -1.5f, 0.5f},
        {1.0f, -2.5f,  2.0f, -0.5f},
        {-0.5f,  0.0f,  0.5f, 0.0f},
        {0.0f,  1.0f,  0.0f, 0.0f}
    };

    float tVec[4] = {t * t * t, t * t, t, 1};
    float tDeriv[4] = {3 * t * t, 2 * t, 1, 0};

    for (int dim = 0; dim < 3; dim++) {
        float p[4] = {
            (dim == 0 ? p0.x : dim == 1 ? p0.y : p0.z),
            (dim == 0 ? p1.x : dim == 1 ? p1.y : p1.z),
            (dim == 0 ? p2.x : dim == 1 ? p2.y : p2.z),
            (dim == 0 ? p3.x : dim == 1 ? p3.y : p3.z)
        };

        float a[4];
        for (int i = 0; i < 4; i++) {
            a[i] = 0;
            for (int j = 0; j < 4; j++) {
                a[i] += m[i][j] * p[j];
            }
        }

        float val = 0, dVal = 0;
        for (int i = 0; i < 4; i++) {
            val += tVec[i] * a[i];
            dVal += tDeriv[i] * a[i];
        }

        if (dim == 0) { res.x = val; deriv.x = dVal; }
        else if (dim == 1) { res.y = val; deriv.y = dVal; }
        else { res.z = val; deriv.z = dVal; }
    }
}

void getGlobalCatmullRomPoint(float gt, const std::vector<Point>& points, Point &res, Point &deriv) {
    int pointCount = points.size();
    float t = gt * pointCount; // global t between 0 and pointCount
    int index = floor(t);
    t = t - index;

    int indices[4];
    indices[0] = (index + pointCount - 1) % pointCount;
    indices[1] = (indices[0] + 1) % pointCount;
    indices[2] = (indices[1] + 1) % pointCount;
    indices[3] = (indices[2] + 1) % pointCount;

    getCatmullRomPoint(t, points[indices[0]], points[indices[1]],
                          points[indices[2]], points[indices[3]], res, deriv);
}

void applyTransformation(const Transformation& transformation) {
    if (std::holds_alternative<Anime_Translate>(transformation)) {
        const auto& t = std::get<Anime_Translate>(transformation);
        float elapsed = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
        float timeRatio = fmod(elapsed, t.time) / t.time;

        Point pos, deriv;
        getGlobalCatmullRomPoint(timeRatio, t.points, pos, deriv);

        glTranslatef(pos.x, pos.y, pos.z);

        if (t.align) {
            // Normalize derivative
            float len = sqrt(deriv.x * deriv.x + deriv.y * deriv.y + deriv.z * deriv.z);
            Point X = { deriv.x / len, deriv.y / len, deriv.z / len };

            Point Y = {0, 1, 0};  // Up vector
            Point Z = {
                Y.y * X.z - Y.z * X.y,
                Y.z * X.x - Y.x * X.z,
                Y.x * X.y - Y.y * X.x
            };
        
            Y = {
                X.y * Z.z - X.z * Z.y,
                X.z * Z.x - X.x * Z.z,
                X.x * Z.y - X.y * Z.x
            };
        
            float m[16] = {
                X.x, X.y, X.z, 0,
                Y.x, Y.y, Y.z, 0,
                Z.x, Z.y, Z.z, 0,
                0,   0,   0,   1
            };
        
            glMultMatrixf(m);
            
        } else {
            for (const auto& point : t.points) {
                glTranslatef(point.x, point.y, point.z);
            }
        }
    } else if (std::holds_alternative<Translate>(transformation)) {
        const auto& t = std::get<Translate>(transformation);
        glTranslatef(t.x, t.y, t.z);

    
    } else if (std::holds_alternative<Rotate>(transformation)) {
        const auto& t = std::get<Rotate>(transformation);
        glRotatef(t.angle, t.x, t.y, t.z);
    } else if (std::holds_alternative<Scale>(transformation)) {
        const auto& t = std::get<Scale>(transformation);
        glScalef(t.x, t.y, t.z);
    }
}
void drawModel(Models& models) {
    glPushMatrix();
    
    for (const auto& transformation : models.transformations) {
        applyTransformation(transformation);
    }
    
    for (const std::string& filename : models.model) {
        drawFigure(filename);
    }
    
    for (Models& childModels : models.models) {
        drawModel(childModels);
    }
    
    glPopMatrix();
}


void changeSize(int w, int h)
{
    if (h == 0)
        h = 1;

    float ratio = w * 1.0 / h;

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    glViewport(0, 0, w, h);

    gluPerspective(world.camera.projection.fov, ratio, world.camera.projection.near, world.camera.projection.far);

    glMatrixMode(GL_MODELVIEW);
}

void renderScene(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // camera
    glLoadIdentity();
    gluLookAt(world.camera.position.x, world.camera.position.y, world.camera.position.z,
              world.camera.lookAt.x, world.camera.lookAt.y, world.camera.lookAt.z,
              world.camera.up.x, world.camera.up.y, world.camera.up.z);

    glPolygonMode(GL_FRONT,GL_LINE); 
    drawAxis(); 
    drawModel(world.models);
    
    
    glutSwapBuffers();
}

void processKeys()
{
    if (keyStates['w'])
    {
        omega += 0.02f;
        if (omega > M_PI / 2.0f)
        {
            omega = M_PI / 2.0f - 0.01f;
        }
    }
    if (keyStates['s'])
    {
        omega -= 0.02f;
        if (omega < -M_PI / 2.0f)
        {
            omega = -M_PI / 2.0f + 0.01f;
        }
    }
    if (keyStates['a'])
    {
        alpha -= 0.02f;
    }
    if (keyStates['d'])
    {
        alpha += 0.02f;
    }
    if (keyStates['-'])
    {
        radius -= 1.0f;

    }
    if (keyStates['+'])
    {
        radius += 1.0f;
    }

    SphericalToCartesian();

    glutPostRedisplay();
}
void parseGroup(tinyxml2::XMLElement *groupElement, Models &models) {
    using namespace tinyxml2;

    XMLElement *transformElement = groupElement->FirstChildElement("transform");
    if (transformElement) {
        for (XMLElement *child = transformElement->FirstChildElement(); child; child = child->NextSiblingElement()) {
            string name_trans = child->Value();

            if (name_trans == "translate") {
                if(child -> Attribute("time" ) != nullptr)
                {
                    float time;
                    bool align = false;
                    child->QueryFloatAttribute("time", &time);
                    const char* align_aux = child->Attribute("align");
                    if (align_aux && strcmp(align_aux, "true") == 0) {
                        align = true;
                    }
                    std::vector<Point> points;
                    for(XMLElement *child_time = child->FirstChildElement("point"); child_time; child_time = child_time->NextSiblingElement("point"))
                    {
                        float x , y , z ;
                        child_time->QueryFloatAttribute("x", &x);
                        child_time->QueryFloatAttribute("y", &y);
                        child_time->QueryFloatAttribute("z", &z);
                        points.emplace_back(Point{x, y, z});
                        
                    }
                    models.transformations.emplace_back(Anime_Translate{time, align, points});
                    
                }
                else
                {
                    float x , y , z ;
                    child->QueryFloatAttribute("x", &x);
                    child->QueryFloatAttribute("y", &y);
                    child->QueryFloatAttribute("z", &z);
                    models.transformations.emplace_back(Translate{x, y, z});
                }
            }
            else if (name_trans == "rotate") {
                float angle , x , y , z;
                child->QueryFloatAttribute("angle", &angle);
                child->QueryFloatAttribute("x", &x);
                child->QueryFloatAttribute("y", &y);
                child->QueryFloatAttribute("z", &z);
                models.transformations.emplace_back(Rotate{x, y, z, angle});

            }
            else if (name_trans == "scale") {
                float x , y, z ;
                child->QueryFloatAttribute("x", &x);
                child->QueryFloatAttribute("y", &y);
                child->QueryFloatAttribute("z", &z);
                models.transformations.emplace_back(Scale{x, y, z});
            }
        }
    }

    
    XMLElement *modelsElement = groupElement->FirstChildElement("models");
    if (modelsElement) {
        XMLElement *modelElement = modelsElement->FirstChildElement("model");
        while (modelElement) {
            const char *file = modelElement->Attribute("file");
            if (file) {
                models.model.push_back(std::string(file));
            }
            modelElement = modelElement->NextSiblingElement("model");
        }
    }

    
    XMLElement *childGroup = groupElement->FirstChildElement("group");
    while (childGroup) {
        Models childModels;
        parseGroup(childGroup, childModels);
        models.models.push_back(childModels);
        childGroup = childGroup->NextSiblingElement("group");
    }
}



void parseInfo(char *filename)
{
    using namespace tinyxml2;
    XMLDocument doc;
    XMLError eResult = doc.LoadFile(filename);

    if (eResult != XML_SUCCESS)
    {
        cout << "Error: " << eResult << endl;
        exit(0);
    }

    XMLNode *pRoot = doc.FirstChild();
    if (pRoot == nullptr)
    {
        cout << "Error: " << XML_ERROR_FILE_READ_ERROR << endl;
        exit(0);
    }

    XMLElement *pElement = pRoot->FirstChildElement("window");
    if (pElement)
    {
        pElement->QueryIntAttribute("width", &world.window.width);
        pElement->QueryIntAttribute("height", &world.window.height);
    }

    pElement = pRoot->FirstChildElement("camera");
    XMLElement *cameraElements = pElement->FirstChildElement("position");
    cameraElements->QueryFloatAttribute("x", &world.camera.position.x);
    cameraElements->QueryFloatAttribute("y", &world.camera.position.y);
    cameraElements->QueryFloatAttribute("z", &world.camera.position.z);

    radius = sqrt(pow(world.camera.position.x, 2) + pow(world.camera.position.y, 2) + pow(world.camera.position.z, 2));
    omega = asin(world.camera.position.y / radius);
    alpha = asin(world.camera.position.x / (radius * cos(omega)));
    XMLElement *cameraElement = pRoot->FirstChildElement("camera");
    if (cameraElement)
    {
        XMLElement *position = cameraElement->FirstChildElement("position");
        if (position)
        {
            position->QueryFloatAttribute("x", &world.camera.position.x);
            position->QueryFloatAttribute("y", &world.camera.position.y);
            position->QueryFloatAttribute("z", &world.camera.position.z);
        }

        XMLElement *lookAt = cameraElement->FirstChildElement("lookAt");
        if (lookAt)
        {
            lookAt->QueryFloatAttribute("x", &world.camera.lookAt.x);
            lookAt->QueryFloatAttribute("y", &world.camera.lookAt.y);
            lookAt->QueryFloatAttribute("z", &world.camera.lookAt.z);
        }

        XMLElement *up = cameraElement->FirstChildElement("up");
        if (up)
        {
            up->QueryFloatAttribute("x", &world.camera.up.x);
            up->QueryFloatAttribute("y", &world.camera.up.y);
            up->QueryFloatAttribute("z", &world.camera.up.z);
        }

        XMLElement *projection = cameraElement->FirstChildElement("projection");
        if (projection)
        {
            projection->QueryFloatAttribute("fov", &world.camera.projection.fov);
            projection->QueryFloatAttribute("near", &world.camera.projection.near);
            projection->QueryFloatAttribute("far", &world.camera.projection.far);
        }
    }
    XMLElement *modelsElement = pRoot->FirstChildElement("group");
    if(modelsElement)
    {
        parseGroup(modelsElement, world.models);
    }
    
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Someting went wrong!" << endl;
        return 0;
    }
    SphericalToCartesian();
    parseInfo(argv[1]);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(world.window.width, world.window.height);
    glutCreateWindow("Phase-1@CG");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);

    glutKeyboardFunc(key_press);
    glutKeyboardUpFunc(key_up);
    glutIdleFunc(processKeys);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT, GL_LINE);

    glutMainLoop();

    return 1;
}