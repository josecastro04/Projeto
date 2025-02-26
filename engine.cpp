#include <cmath>
#include <fstream>
#include <iostream>
#include <GL/freeglut_std.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <list>
#include <sstream>
#include <string>
#include "TinyXML/tinyxml2.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

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
    std::list<std::string> model;
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
float radius;
float alfa = 0.0f, beta = 0.0f;

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
}

void drawFigure(std::string filename)
{
    std::ifstream file(filename);
    if (!file)
    {
        std::cerr << "Error when trying to open the file!" << std::endl;
        return;
    }

    glBegin(GL_TRIANGLES);
    glColor3f(1.0f, 1.0f, 1.0f);
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line))
    {
        std::istringstream stream(line);

        float x, y, z;

        if (!(stream >> x >> y >> z))
        {
            std::cerr << "Error when trying to read the values!" << std::endl;
            glEnd();
            file.close();
            return;
        }

        glVertex3f(x, y, z);
    }

    glEnd();
    file.close();
}

void spherical2Cartesian(){
    world.camera.position.x = radius * cos(beta) * sin(alfa);
    world.camera.position.y = radius * sin(beta);
    world.camera.position.z = radius * cos(beta) * cos(alfa);
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

    for (std::string filename : world.models.model)
    {
        drawFigure(filename);
    }
    glutSwapBuffers();
}

void parseInfo(char *filename)
{
    using namespace tinyxml2;
    XMLDocument doc;
    XMLError eResult = doc.LoadFile(filename);

    if (eResult != XML_SUCCESS)
    {
        std::cout << "Error: " << eResult << std::endl;
        exit(0);
    }

    XMLNode *pRoot = doc.FirstChild();
    if (pRoot == nullptr)
    {
        std::cout << "Error: " << XML_ERROR_FILE_READ_ERROR << std::endl;
        exit(0);
    }

    XMLElement *pElement = pRoot->FirstChildElement("window");
    if (pElement == nullptr)
    {
        std::cout << "Error: " << XML_ERROR_FILE_READ_ERROR << std::endl;
        exit(0);
    }

    pElement->QueryIntAttribute("width", &world.window.width);
    pElement->QueryIntAttribute("height", &world.window.height);

    pElement = pRoot->FirstChildElement("camera");
    XMLElement *cameraElements = pElement->FirstChildElement("position");
    cameraElements->QueryFloatAttribute("x", &world.camera.position.x);
    cameraElements->QueryFloatAttribute("y", &world.camera.position.y);
    cameraElements->QueryFloatAttribute("z", &world.camera.position.z);

    cameraElements = pElement->FirstChildElement("lookAt");
    cameraElements->QueryFloatAttribute("x", &world.camera.lookAt.x);
    cameraElements->QueryFloatAttribute("y", &world.camera.lookAt.y);
    cameraElements->QueryFloatAttribute("z", &world.camera.lookAt.z);

    cameraElements = pElement->FirstChildElement("up");
    cameraElements->QueryFloatAttribute("x", &world.camera.up.x);
    cameraElements->QueryFloatAttribute("y", &world.camera.up.y);
    cameraElements->QueryFloatAttribute("z", &world.camera.up.z);

    cameraElements = pElement->FirstChildElement("projection");
    cameraElements->QueryFloatAttribute("fov", &world.camera.projection.fov);
    cameraElements->QueryFloatAttribute("near", &world.camera.projection.near);
    cameraElements->QueryFloatAttribute("far", &world.camera.projection.far);

    pElement = pRoot->FirstChildElement("group");
    XMLElement *modelsElements = pElement->FirstChildElement("models");
    XMLElement *listModelElements = modelsElements->FirstChildElement("model");

    while (listModelElements != nullptr)
    {
        const char *model_filename = nullptr;

        model_filename = listModelElements->Attribute("file");

        std::string model = model_filename;
        world.models.model.push_back(model);
        listModelElements = listModelElements->NextSiblingElement("model");
    }
}

void processSpecialKeys(int key, int xx, int yy){
    switch(key){
    case GLUT_KEY_RIGHT:
        alfa -= 0.1; break;

    case GLUT_KEY_LEFT:
        alfa += 0.1; break;

    case GLUT_KEY_UP:
        beta += 0.1; break;

    case GLUT_KEY_DOWN:
        beta -= 0.1f; break;

    case GLUT_KEY_PAGE_DOWN: radius -= 0.1f;
        if(radius < 0.1f)
            radius = 0.1f;
        break;

    case GLUT_KEY_PAGE_UP: radius += 0.1f; break;
    }

    spherical2Cartesian();
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Someting went wrong!" << std::endl;
        return 0;
    }
    parseInfo(argv[1]);

    radius = sqrt(pow(world.camera.position.x - world.camera.lookAt.x, 2)
            + pow(world.camera.position.y - world.camera.lookAt.y, 2)
            + pow(world.camera.position.z - world.camera.lookAt.z, 2));

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(world.window.width, world.window.height);
    glutCreateWindow("Phase-1@CG");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);

    glutSpecialFunc(processSpecialKeys);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT, GL_LINE);

    spherical2Cartesian();

    glutMainLoop();

    return 1;
}
