#include <cmath>
#include <fstream>
#include <iostream>
#include <GL/freeglut_std.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <list>
#include <sstream>
#include <string>
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

struct Models
{
    list<string> model;
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

    glPolygonMode(GL_FRONT, GL_LINE);
    // desenhar models
    drawAxis();

    for (string filename : world.models.model)
    {
        drawFigure(filename);
    }
    glutSwapBuffers();
}

void processKeys()
{
    if (keyStates['w'])
    {
        omega += 0.02f;
        if (omega > M_PI / 2.0f)
        {
            omega = M_PI / 2.0f;
        }
    }
    if (keyStates['s'])
    {
        omega -= 0.02f;
        if (omega < -M_PI / 2.0f)
        {
            omega = -M_PI / 2.0f;
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
        radius -= 0.01f;
        if (radius < 3.0f)
        {
            radius = 3.0f;
        }
    }
    if (keyStates['+'])
    {
        radius += 0.01f;
        if (radius > 50.0f)
        {
            radius = 50.0f;
        }
    }

    SphericalToCartesian();

    glutPostRedisplay();
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
    if (modelsElement)
    {
        XMLElement *modelNode = modelsElement->FirstChildElement("models")->FirstChildElement("model");
        while (modelNode)
        {
            const char *modelFile = modelNode->Attribute("file");
            if (modelFile)
            {
                world.models.model.push_back(std::string(modelFile));
            }
            modelNode = modelNode->NextSiblingElement("model");
        }
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