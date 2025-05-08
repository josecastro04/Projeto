#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif
#include <GL/glut.h>
#include "transformations.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <variant>
#include <tinyxml2.h>
#include <map>

using namespace std;
bool solidMode = false;
static float omega, alpha, radius = 5.0f;
static float k = 0.5f;

struct Color
{
    float diffuse[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float ambient[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float specular[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float emissive[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float shininess = 0.0f;
};

struct Model
{
    std::string file;
    Color color;
};

struct Models
{
    std::vector<Model> model;
    std::vector<Transformation> transformations;
    std::vector<Models> models;
};

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

struct ModelData
{
    GLuint vbo[2];
    int vertexCount;
};

std::map<std::string, ModelData> modelCache;

Point calculateVector(Point L, Point P)
{
    return {L.x - P.x,
            L.y - P.y,
            L.z - P.z};
}

void SphericalToCartesianFPS()
{
    float direction[3];
    direction[0] = cos(alpha) * cos(omega);
    direction[1] = sin(omega);
    direction[2] = sin(alpha) * cos(omega);

    normalize(direction);

    world.camera.lookAt.x = world.camera.position.x + direction[0];
    world.camera.lookAt.y = world.camera.position.y + direction[1];
    world.camera.lookAt.z = world.camera.position.z + direction[2];
}
void key_press(unsigned char key, int x, int y)
{

    Point d1 = calculateVector(world.camera.lookAt, world.camera.position);
    float d[3];
    d[0] = d1.x;
    d[1] = d1.y;
    d[2] = d1.z;
    normalize(d);

    if (key == 'e')
    {
        alpha += 0.02f;
        SphericalToCartesianFPS();
    }
    if (key == 'q')
    {
        alpha -= 0.02f;
        SphericalToCartesianFPS();
    }
    if (key == 'x')
    {
        omega += 0.02f;
        SphericalToCartesianFPS();
    }
    if (key == 'z')
    {
        omega -= 0.02f;
        SphericalToCartesianFPS();
    }

    if (key == 'w')
    {
        world.camera.lookAt.x += k * d[0];
        world.camera.lookAt.y += k * d[1];
        world.camera.lookAt.z += k * d[2];
        world.camera.position.x += k * d[0];
        world.camera.position.y += k * d[1];
        world.camera.position.z += k * d[2];
    }
    if (key == 's')
    {
        world.camera.lookAt.x -= k * d[0];
        world.camera.lookAt.y -= k * d[1];
        world.camera.lookAt.z -= k * d[2];
        world.camera.position.x -= k * d[0];
        world.camera.position.y -= k * d[1];
        world.camera.position.z -= k * d[2];
    }
    if (key == 'a')
    {
        float right[3];
        float up[3] = {world.camera.up.x, world.camera.up.y, world.camera.up.z}; // (0,1,0) normalmente
        cross(up, d, right);
        normalize(right);

        world.camera.position.x += k * right[0];
        world.camera.position.y += k * right[1];
        world.camera.position.z += k * right[2];
        world.camera.lookAt.x += k * right[0];
        world.camera.lookAt.y += k * right[1];
        world.camera.lookAt.z += k * right[2];
    }
    if (key == 'd')
    {
        float right[3];
        float up[3] = {world.camera.up.x, world.camera.up.y, world.camera.up.z}; // (0,1,0) normalmente
        cross(up, d, right);
        normalize(right);

        world.camera.position.x -= k * right[0];
        world.camera.position.y -= k * right[1];
        world.camera.position.z -= k * right[2];
        world.camera.lookAt.x -= k * right[0];
        world.camera.lookAt.y -= k * right[1];
        world.camera.lookAt.z -= k * right[2];
    }
    if (key == '-')
    {
        k /= 2;
        if (k < 0.1)
            k = 0.1;
    }
    if (key == '+')
    {
        k *= 2;
        if (k > 25)
            k = 25;
    }

    if (key == 't')
    {
        
        solidMode = !solidMode;
        if (solidMode)
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }
        else
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
    }

    glutPostRedisplay();
}

void drawAxis()
{
    glDisable(GL_LIGHTING);
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
    glEnable(GL_LIGHTING);
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
            Model m;

            const char *file = model->Attribute("file");
            if (file)
                m.file = file;

            XMLElement *color = model->FirstChildElement("color");
            if (color)
            {
                XMLElement *diffuse = color->FirstChildElement("diffuse");
                if (diffuse)
                    diffuse->QueryFloatAttribute("R", &m.color.diffuse[0]),
                    diffuse->QueryFloatAttribute("G", &m.color.diffuse[1]),
                    diffuse->QueryFloatAttribute("B", &m.color.diffuse[2]);

                XMLElement *ambient = color->FirstChildElement("ambient");
                if (ambient)
                    ambient->QueryFloatAttribute("R", &m.color.ambient[0]),
                    ambient->QueryFloatAttribute("G", &m.color.ambient[1]),
                    ambient->QueryFloatAttribute("B", &m.color.ambient[2]);

                XMLElement *specular = color->FirstChildElement("specular");
                if (specular)
                    specular->QueryFloatAttribute("R", &m.color.specular[0]),
                    specular->QueryFloatAttribute("G", &m.color.specular[1]),
                    specular->QueryFloatAttribute("B", &m.color.specular[2]);

                XMLElement *emissive = color->FirstChildElement("emissive");
                if (emissive)
                    emissive->QueryFloatAttribute("R", &m.color.emissive[0]),
                    emissive->QueryFloatAttribute("G", &m.color.emissive[1]),
                    emissive->QueryFloatAttribute("B", &m.color.emissive[2]);

                XMLElement *shininess = color->FirstChildElement("shininess");
                if (shininess)
                    shininess->QueryFloatAttribute("value", &m.color.shininess);
            }

            models.model.push_back(m); 
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

    radius = sqrt(pow(world.camera.lookAt.x - world.camera.position.x, 2) +
                  pow(world.camera.lookAt.y - world.camera.position.y, 2) +
                  pow(world.camera.lookAt.z - world.camera.position.z, 2));
    omega = asin((world.camera.lookAt.y - world.camera.position.y) / radius);
    alpha = atan2(world.camera.lookAt.z - world.camera.position.z,
                  world.camera.lookAt.x - world.camera.position.x);

                  
    XMLElement *group = pRoot->FirstChildElement("group");
    if (group)
        parseGroup(group, world.models);
}

void drawFigureVBO(string filename)
{
    if (modelCache.find(filename) == modelCache.end())
    {

        ifstream file(filename);
        if (!file)
        {
            cerr << "Error when trying to open the file!" << endl;
            exit(0);
        }

        int num_vertices = 0;
        string line;
        getline(file, line);
        istringstream stream(line);
        stream >> num_vertices;

        float *v, *n;

        v = (float *)malloc(sizeof(float) * num_vertices * 3);
        n = (float *)malloc(sizeof(float) * num_vertices * 3);

        for (int i = 0; i < num_vertices; i++)
        {
            getline(file, line);
            istringstream stream(line);
            float x, y, z;

            if (!(stream >> x >> y >> z))
            {
                cerr << "Error when trying to read the values!" << endl;
                free(v);
                free(n);
                glEnd();
                file.close();
                return;
            }

            v[i * 3] = x;
            v[i * 3 + 1] = y;
            v[i * 3 + 2] = z;
            n[i * 3] = x;
            n[i * 3 + 1] = y;
            n[i * 3 + 2] = z;
        }

        GLuint buffers[2];

        glGenBuffers(2, buffers);

        glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * num_vertices * 3, v, GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * num_vertices * 3, n, GL_STATIC_DRAW);

        free(v);
        free(n);

        ModelData data;
        data.vbo[0] = buffers[0];
        data.vbo[1] = buffers[1];
        data.vertexCount = num_vertices;
        modelCache[filename] = data;
    }

    ModelData &data = modelCache[filename];

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
    glPushMatrix();
    for (const auto &transformation : models.transformations)
    {
        if (holds_alternative<Anime_Translate>(transformation.type))
        {
            const auto &trans = get<Anime_Translate>(transformation.type);
            drawCatmullRomCurve(trans.points);
        }
    }

    for (auto &transformation : models.transformations)
    {
        applyTransformation(transformation);
    }

    for (const Model &m : models.model)
    {

        glMaterialfv(GL_FRONT, GL_DIFFUSE, m.color.diffuse);
        glMaterialfv(GL_FRONT, GL_AMBIENT, m.color.ambient);
        glMaterialfv(GL_FRONT, GL_SPECULAR, m.color.specular);
        glMaterialfv(GL_FRONT, GL_EMISSION, m.color.emissive);
        glMaterialf(GL_FRONT, GL_SHININESS, m.color.shininess);

        drawFigureVBO(m.file);
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
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(world.camera.position.x, world.camera.position.y, world.camera.position.z,
              world.camera.lookAt.x, world.camera.lookAt.y, world.camera.lookAt.z,
              world.camera.up.x, world.camera.up.y, world.camera.up.z);
    drawAxis();
    if (solidMode)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    else
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
   

    drawModel(world.models);

    glutSwapBuffers();
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
    glutIdleFunc(renderScene);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glutMainLoop();

    for (auto &entry : modelCache)
    {
        glDeleteBuffers(2, entry.second.vbo);
    }

    return 0;
}