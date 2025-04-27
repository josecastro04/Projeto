#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <variant>
#include <tinyxml2.h>
#include "transformations.h"
#include <iostream>
#include <map>

using namespace std;
static float omega = 0.0f, alpha = 0.0f, radius = 5.0f;
bool keyStates[256] = {false};

struct Models
{
    std::vector<std::string> model;
    std::vector<Transformation> transformations;
    std::vector<Models> models;
};

struct ModelData
{
    GLuint vbo[2];
    int vertexCount;
};

std::map<std::string, ModelData> modelCache;

struct Window
{
    int width;
    int height;
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
void parseGroup(tinyxml2::XMLElement *groupElement, Models &models)
{
    using namespace tinyxml2;

    XMLElement *transformElement = groupElement->FirstChildElement("transform");
    if (transformElement)
    {
        for (XMLElement *child = transformElement->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            string name_trans = child->Name();

            if (name_trans == "translate")
            {
                if (child->Attribute("time"))
                {

                    Anime_Translate at;
                    child->QueryFloatAttribute("time", &at.time);
                    const char *align = child->Attribute("align");
                    at.align = align && string(align) == "true";

                    for (XMLElement *point = child->FirstChildElement("point"); point; point = point->NextSiblingElement("point"))
                    {
                        Point p;
                        point->QueryFloatAttribute("x", &p.x);
                        point->QueryFloatAttribute("y", &p.y);
                        point->QueryFloatAttribute("z", &p.z);
                        at.points.push_back(p);
                    }
                    models.transformations.push_back({at});
                }
                else
                {

                    Translate t;
                    child->QueryFloatAttribute("x", &t.x);
                    child->QueryFloatAttribute("y", &t.y);
                    child->QueryFloatAttribute("z", &t.z);
                    models.transformations.push_back({t});
                }
            }
            else if (name_trans == "rotate")
            {
                if (child->Attribute("time"))
                {

                    Anime_Rotate ar;
                    child->QueryFloatAttribute("time", &ar.time);
                    child->QueryFloatAttribute("x", &ar.x);
                    child->QueryFloatAttribute("y", &ar.y);
                    child->QueryFloatAttribute("z", &ar.z);
                    models.transformations.push_back({ar});
                }
                else
                {

                    Rotate r;
                    child->QueryFloatAttribute("angle", &r.angle);
                    child->QueryFloatAttribute("x", &r.x);
                    child->QueryFloatAttribute("y", &r.y);
                    child->QueryFloatAttribute("z", &r.z);
                    models.transformations.push_back({r});
                }
            }
            else if (name_trans == "scale")
            {
                Scale s;
                child->QueryFloatAttribute("x", &s.x);
                child->QueryFloatAttribute("y", &s.y);
                child->QueryFloatAttribute("z", &s.z);
                models.transformations.push_back({s});
            }
        }
    }

    XMLElement *modelsElement = groupElement->FirstChildElement("models");
    if (modelsElement)
    {
        for (XMLElement *model = modelsElement->FirstChildElement("model"); model; model = model->NextSiblingElement("model"))
        {
            const char *file = model->Attribute("file");
            if (file)
                models.model.push_back(file);
        }
    }

    for (XMLElement *childGroup = groupElement->FirstChildElement("group"); childGroup; childGroup = childGroup->NextSiblingElement("group"))
    {
        Models childModels;
        parseGroup(childGroup, childModels);
        models.models.push_back(childModels);
    }
}

void parseInfo(char *filename)
{
    using namespace tinyxml2;
    XMLDocument doc;
    if (doc.LoadFile(filename) != XML_SUCCESS)
    {
        cerr << "Error loading XML file" << endl;
        exit(1);
    }

    XMLNode *pRoot = doc.FirstChild();
    if (!pRoot)
    {
        cerr << "Empty XML file" << endl;
        exit(1);
    }

    XMLElement *window = pRoot->FirstChildElement("window");
    if (window)
    {
        window->QueryIntAttribute("width", &world.window.width);
        window->QueryIntAttribute("height", &world.window.height);
    }

    XMLElement *camera = pRoot->FirstChildElement("camera");
    if (camera)
    {
        XMLElement *pos = camera->FirstChildElement("position");
        if (pos)
        {
            pos->QueryFloatAttribute("x", &world.camera.position.x);
            pos->QueryFloatAttribute("y", &world.camera.position.y);
            pos->QueryFloatAttribute("z", &world.camera.position.z);
        }

        XMLElement *lookAt = camera->FirstChildElement("lookAt");
        if (lookAt)
        {
            lookAt->QueryFloatAttribute("x", &world.camera.lookAt.x);
            lookAt->QueryFloatAttribute("y", &world.camera.lookAt.y);
            lookAt->QueryFloatAttribute("z", &world.camera.lookAt.z);
        }

        XMLElement *up = camera->FirstChildElement("up");
        if (up)
        {
            up->QueryFloatAttribute("x", &world.camera.up.x);
            up->QueryFloatAttribute("y", &world.camera.up.y);
            up->QueryFloatAttribute("z", &world.camera.up.z);
        }

        XMLElement *proj = camera->FirstChildElement("projection");
        if (proj)
        {
            proj->QueryFloatAttribute("fov", &world.camera.projection.fov);
            proj->QueryFloatAttribute("near", &world.camera.projection.near);
            proj->QueryFloatAttribute("far", &world.camera.projection.far);
        }
    }

    radius = sqrt(pow(world.camera.position.x, 2) +
                  pow(world.camera.position.y, 2) +
                  pow(world.camera.position.z, 2));
    omega = asin(world.camera.position.y / radius);
    alpha = atan2(world.camera.position.x, world.camera.position.z);

    XMLElement *group = pRoot->FirstChildElement("group");
    if (group)
        parseGroup(group, world.models);
}

void drawFigureVBO(string filename) {

    if (modelCache.find(filename) == modelCache.end()) {
        ifstream file(filename);
        if (!file) {
            cerr << "Error when trying to open the file!" << endl;
            return;
        }

        int num_vertices = 0;
        string line;
        getline(file, line);
        istringstream stream(line);
        stream >> num_vertices;

        float *v = (float *)malloc(sizeof(float) * num_vertices * 3);
        float *n = (float *)malloc(sizeof(float) * num_vertices * 3);

        for (int i = 0; i < num_vertices; i++) {
            getline(file, line);
            istringstream stream(line);
            float x, y, z;

            if (!(stream >> x >> y >> z)) {
                cerr << "Error when trying to read the values!" << endl;
                free(v);
                free(n);
                file.close();
                return;
            }

            v[i*3] = x;
            v[i*3+1] = y;
            v[i*3+2] = z;
            n[i*3] = x;
            n[i*3+1] = y;
            n[i*3+2] = z;
        }
        file.close();

        // Cria VBOs
        GLuint buffers[2];
        glGenBuffers(2, buffers);
        
        // Buffer de vértices
        glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * num_vertices * 3, v, GL_STATIC_DRAW);
        
        // Buffer de normais
        glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * num_vertices * 3, n, GL_STATIC_DRAW);

        free(v);
        free(n);
        
        // Armazena os dados no cache
        ModelData data;
        data.vbo[0] = buffers[0];
        data.vbo[1] = buffers[1];
        data.vertexCount = num_vertices;
        modelCache[filename] = data;
    }

    // Obtém os dados do cache
    ModelData &data = modelCache[filename];

    // Desenha usando os VBOs
    glBindBuffer(GL_ARRAY_BUFFER, data.vbo[0]);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, data.vbo[1]);
    glNormalPointer(GL_FLOAT, 0, 0);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    
    glDrawArrays(GL_TRIANGLES, 0, data.vertexCount);
    
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
}

void drawModel(Models &models)
{

    for (const auto &transformation : models.transformations)
    {
        if (holds_alternative<Anime_Translate>(transformation.type))
        {
            const auto &trans = get<Anime_Translate>(transformation.type);
            drawCatmullRomCurve(trans.points);
        }
    }

    glPushMatrix();

    for (auto &transformation : models.transformations)
    {
        applyTransformation(transformation);
    }

    for (const string &filename : models.model)
    {
        glPushMatrix();
        drawFigureVBO(filename);
        glPopMatrix();
    }

    for (Models &child : models.models)
    {
        drawModel(child);
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
    gluPerspective(world.camera.projection.fov, ratio,
                   world.camera.projection.near, world.camera.projection.far);
    glMatrixMode(GL_MODELVIEW);
}

void renderScene()
{
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

void processKeys()
{
    if (keyStates['w'])
    {
        omega += 0.02f;
        if (omega > M_PI / 2.0f)
            omega = M_PI / 2.0f - 0.01f;
    }
    if (keyStates['s'])
    {
        omega -= 0.02f;
        if (omega < -M_PI / 2.0f)
            omega = -M_PI / 2.0f + 0.01f;
    }
    if (keyStates['a'])
        alpha -= 0.02f;
    if (keyStates['d'])
        alpha += 0.02f;
    if (keyStates['-'])
        radius -= 1.0f;
    if (keyStates['+'])
        radius += 1.0f;

    SphericalToCartesian();
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <scene_file.xml>" << endl;
        return 1;
    }

    parseInfo(argv[1]);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(world.window.width, world.window.height);
    glutCreateWindow("CG@DI Phase 3");

    glewInit();

    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);

    glutKeyboardFunc(key_press);
    glutKeyboardUpFunc(key_up);
    glutIdleFunc(processKeys);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT, GL_LINE);

    glutMainLoop();

    for (auto &entry : modelCache)
    {
        glDeleteBuffers(2, entry.second.vbo);
    }

    return 0;
}