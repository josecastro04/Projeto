#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

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
#include <map>

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

struct Translate
{
    float x, y, z;
};

struct Rotate
{
    float x, y, z, angle;
};

struct Scale
{
    float x, y, z;
};

struct Anime_Rotate
{
    float x, y, z, time;
};

struct Anime_Translate
{
    float time;
    bool align;
    std::vector<Point> points;
};

struct Transformation
{
    std::variant<Translate, Anime_Translate, Rotate, Anime_Rotate, Scale> type;
    bool is_active = true;
    float start_time = 0.0f;
};

struct Models
{
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

struct ModelData
{
    GLuint vbo[2];
    int vertexCount;
};

std::map<std::string, ModelData> modelCache;

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

        // Liberar a memória alocada
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

void multMatrixVector(float *m, float *v, float *res)
{
    for (int j = 0; j < 4; ++j)
    {
        res[j] = 0;
        for (int k = 0; k < 4; ++k)
        {
            res[j] += v[k] * m[j * 4 + k];
        }
    }
}

void getCatmullRomPoint(float t, Point p0, Point p1, Point p2, Point p3, float *pos, float *deriv)
{

    float m[16] = {-0.5f, 1.5f, -1.5f, 0.5f,
                   1.0f, -2.5f, 2.0f, -0.5f,
                   -0.5f, 0.0f, 0.5f, 0.0f,
                   0.0f, 1.0f, 0.0f, 0.0f};
    float a[4];
    float T[4] = {t * t * t, t * t, t, 1};
    float dT[4] = {3 * t * t, 2 * t, 1, 0};
    float p[4][3] = {{p0.x, p0.y, p0.z},
                     {p1.x, p1.y, p1.z},
                     {p2.x, p2.y, p2.z},
                     {p3.x, p3.y, p3.z}};
    for ( int i = 0; i < 3; i++)
    {
        float pp[4] = {p[0][i], p[1][i], p[2][i], p[3][i]};
        multMatrixVector(m, pp, a);
        pos[i] = a[0] * T[0] + a[1] * T[1] + a[2] * T[2] + a[3] * T[3];
        deriv[i] = a[0] * dT[0] + a[1] * dT[1] + a[2] * dT[2] + a[3] * dT[3];

    }
}

void getGlobalCatmullRomPoint(float gt, const vector<Point> &points, float *pos, float *deriv)
{
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

void drawCatmullRomCurve(const vector<Point> &points)
{
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINE_LOOP);
    float pos[3], deriv[3];
    for (float gt = 0; gt < 1; gt += 0.01)
    {
        getGlobalCatmullRomPoint(gt, points, pos, deriv);
        glVertex3f(pos[0], pos[1], pos[2]);
    }

    glEnd();
}

void normalize(float *p)
{
    float length = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    if (length > 0)
    {
        p[0] /= length;
        p[1] /= length;
        p[2] /= length;
    }
}

void cross(float *a, float *b, float *result)
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

void buildRotMatrix(float *x, float *y, float *z, float *m)
{

        m[0] = x[0];
        m[1] = x[1];
        m[2] = x[2];
        m[3] = 0;
        m[4] = y[0];
        m[5] = y[1];
        m[6] = y[2];
        m[7] = 0;
        m[8] = z[0];
        m[9] = z[1];
        m[10] = z[2];
        m[11] = 0;
        m[12] = 0;
        m[13] = 0;
        m[14] = 0;
        m[15] = 1;
}



    void applyTransformation(const Transformation &transformation)
    {
        if (const auto *trans = get_if<Anime_Translate>(&transformation.type))
        {
            float elapsed = glutGet(GLUT_ELAPSED_TIME) / 1000.0f - transformation.start_time;
            float gt = fmod(elapsed, trans->time) / trans->time;

            float pos[3], deriv[3];
            getGlobalCatmullRomPoint(gt, trans->points, pos, deriv);

            glTranslatef(pos[0], pos[1], pos[2]);

            if (trans->align)
            {
                float X[3] = {deriv[0], deriv[1], deriv[2]};
                normalize(X);

                float Y[3] = {0.0f, 1.0f, 0.0f};
                float Z[3];
                cross(X, Y, Z);

                normalize(Z);

                cross(Z, X, Y);
                normalize(Y);

                float m[16];
                buildRotMatrix(X, Y, Z, m);
                glMultMatrixf(m);
            }
        }
        else if (const auto *trans = get_if<Translate>(&transformation.type))
        {
            glTranslatef(trans->x, trans->y, trans->z);
        }
        else if (const auto *trans = get_if<Anime_Rotate>(&transformation.type))
        {
            float elapsed = glutGet(GLUT_ELAPSED_TIME) / 1000.0f - transformation.start_time;
            float angle = 360.0f * fmod(elapsed, trans->time) / trans->time;
            glRotatef(angle, trans->x, trans->y, trans->z);
        }
        else if (const auto *trans = get_if<Rotate>(&transformation.type))
        {
            glRotatef(trans->angle, trans->x, trans->y, trans->z);
        }
        else if (const auto *trans = get_if<Scale>(&transformation.type))
        {
            glScalef(trans->x, trans->y, trans->z);
        }
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

    for (const auto &transformation : models.transformations)
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
