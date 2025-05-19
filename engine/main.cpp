#ifdef _APPLE_
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glext.h>

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
#include <IL/il.h>

using namespace std;
bool solidMode = true;
static float omega, alpha, radius = 5.0f;
static float k = 0.5f;
unsigned int figure = 0;
int i = 1;
static int mouseX = 0, mouseY = 0;
static float zoomFactor = 1.0f;
static bool mouseLeftDown = false, mouseRightDown = false;

struct Color
{
    float diffuse[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float ambient[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float specular[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float emissive[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float shininess = 0.0f;
};

struct Model
{
    std::string file;
    std::string filetextura;
    Color color;
    GLuint textureID = 0;

    GLuint vbo[3];
    int vertexCount;
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
struct Light
{
    enum Type
    {
        DIRECTIONAL,
        POINT,
        SPOT
    };
    Type type;
    float position[4];  
    float direction[4]; 
    float cutoff;       
};

struct World
{
    Window window;
    Camera camera;
    std::vector<Light> lights;
    Models models;
};
World world;

struct ModelData
{
    GLuint vbo[3];
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
        float up[3] = {world.camera.up.x, world.camera.up.y, world.camera.up.z}; 
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
        float up[3] = {world.camera.up.x, world.camera.up.y, world.camera.up.z}; 
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

int loadTexture(const char *filename)
{
    unsigned int t, tw, th;
    unsigned char *texData;
    unsigned int texID;

    ilInit();
    ilEnable(IL_ORIGIN_SET);
    ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
    ilGenImages(1, &t);
    ilBindImage(t);
    ilLoadImage((const char *)filename);
    tw = ilGetInteger(IL_IMAGE_WIDTH);
    th = ilGetInteger(IL_IMAGE_HEIGHT);
    ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);
    texData = ilGetData();

    glGenTextures(1, &texID);
    printf("Texture ID: %d, Width: %d, Height: %d\n", texID, tw, th);

    glBindTexture(GL_TEXTURE_2D, texID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tw, th, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
    glGenerateMipmap(GL_TEXTURE_2D);
    printf("TEXTid: %d\n", texID);

    glBindTexture(GL_TEXTURE_2D, 0);
    return texID;
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
            printf("loading model \n");
            Model m;

            const char *file = model->Attribute("file");
            if (file)
                m.file = file;

            XMLElement *color = model->FirstChildElement("color");
            if (color)
            {
                XMLElement *diffuse = color->FirstChildElement("diffuse");
                if (diffuse)
                {
                    diffuse->QueryFloatAttribute("R", &m.color.diffuse[0]);
                    diffuse->QueryFloatAttribute("G", &m.color.diffuse[1]);
                    diffuse->QueryFloatAttribute("B", &m.color.diffuse[2]);
                }
                m.color.diffuse[0] /= 255.0f;
                m.color.diffuse[1] /= 255.0f;
                m.color.diffuse[2] /= 255.0f;
                m.color.diffuse[3] = 1.0f;

                XMLElement *ambient = color->FirstChildElement("ambient");
                if (ambient)
                {
                    ambient->QueryFloatAttribute("R", &m.color.ambient[0]);
                    ambient->QueryFloatAttribute("G", &m.color.ambient[1]);
                    ambient->QueryFloatAttribute("B", &m.color.ambient[2]);
                }
                m.color.ambient[0] /= 255.0f;
                m.color.ambient[1] /= 255.0f;
                m.color.ambient[2] /= 255.0f;
                m.color.ambient[3] = 1.0f;
                XMLElement *specular = color->FirstChildElement("specular");
                if (specular)
                {
                    specular->QueryFloatAttribute("R", &m.color.specular[0]);
                    specular->QueryFloatAttribute("G", &m.color.specular[1]);
                    specular->QueryFloatAttribute("B", &m.color.specular[2]);
                }
                m.color.specular[0] /= 255.0f;
                m.color.specular[1] /= 255.0f;
                m.color.specular[2] /= 255.0f;
                m.color.specular[3] = 1.0f;

                XMLElement *emissive = color->FirstChildElement("emissive");
                if (emissive)
                {
                    emissive->QueryFloatAttribute("R", &m.color.emissive[0]);
                    emissive->QueryFloatAttribute("G", &m.color.emissive[1]);
                    emissive->QueryFloatAttribute("B", &m.color.emissive[2]);
                }
                m.color.emissive[0] /= 255.0f;
                m.color.emissive[1] /= 255.0f;
                m.color.emissive[2] /= 255.0f;
                m.color.emissive[3] = 1.0f;

                XMLElement *shininess = color->FirstChildElement("shininess");
                if (shininess)
                    shininess->QueryFloatAttribute("value", &m.color.shininess);
            }
            XMLElement *texture = model->FirstChildElement("texture");
            if (texture)
            {
                const char *textureFile = texture->Attribute("file");
                if (textureFile)
                {
                    m.filetextura = textureFile;
                   
                   
                }
            }

            models.model.push_back(m);
            printf("has now size %d\n", models.model.size());
        }
    }

    for (XMLElement *childGroup = groupElement->FirstChildElement("group"); childGroup; childGroup = childGroup->NextSiblingElement("group"))
    {
        Models childModels;
        parseGroup(childGroup, childModels);
        models.models.push_back(childModels);
        

    }
}

void parseInfo(char *filename){
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

    XMLElement *lights = pRoot->FirstChildElement("lights");
    if (lights)
    {
        for (XMLElement *lightElem = lights->FirstChildElement("light"); lightElem; lightElem = lightElem->NextSiblingElement("light"))
        {
            Light light; 

            const char *type = lightElem->Attribute("type");
            if (type)
            {
                if (strcmp(type, "directional") == 0)
                {
                    light.type = Light::DIRECTIONAL;
                    lightElem->QueryFloatAttribute("dirx", &light.direction[0]);
                    lightElem->QueryFloatAttribute("diry", &light.direction[1]);
                    lightElem->QueryFloatAttribute("dirz", &light.direction[2]);
                    float x = light.direction[0];
                    float y = light.direction[1];
                    float z = light.direction[2];
                    float length = sqrt(x * x + y * y + z * z);
                    if (length != 0.0f)
                    {
                        light.direction[0] = x / length;
                        light.direction[1] = y / length;
                        light.direction[2] = z / length;
                    }

                    light.direction[3] = 0.0f;
                }
                else if (strcmp(type, "point") == 0)
                {
                    light.type = Light::POINT;
                    lightElem->QueryFloatAttribute("posx", &light.position[0]);
                    lightElem->QueryFloatAttribute("posy", &light.position[1]);
                    lightElem->QueryFloatAttribute("posz", &light.position[2]);
                    light.position[3] = 1.0f; 

                    
                    light.direction[0] = light.direction[1] = light.direction[2] = 0.0f;
                    light.direction[3] = 0.0f;
                }
                else if (strcmp(type, "spot") == 0)
                {
                    light.type = Light::SPOT;
                    lightElem->QueryFloatAttribute("posx", &light.position[0]);
                    lightElem->QueryFloatAttribute("posy", &light.position[1]);
                    lightElem->QueryFloatAttribute("posz", &light.position[2]);
                    light.position[3] = 1.0f; // w = 1 (ponto)
                    lightElem->QueryFloatAttribute("dirx", &light.direction[0]);
                    lightElem->QueryFloatAttribute("diry", &light.direction[1]);
                    lightElem->QueryFloatAttribute("dirz", &light.direction[2]);
                    lightElem->QueryFloatAttribute("cutoff", &light.cutoff);
                    light.direction[3] = 0.0f; // w = 0 (direcional)
                }
            }
            world.lights.push_back(light); // Adiciona a luz ao vetor
        }
    }

    XMLElement *group = pRoot->FirstChildElement("group");
    if (group)
        parseGroup(group, world.models);
}


void drawFigureVBO(string filename, GLuint textureID)
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

        float *v, *n, *t;

        v = (float *)malloc(sizeof(float) * num_vertices * 3);
        n = (float *)malloc(sizeof(float) * num_vertices * 3);
        t = (float *)malloc(sizeof(float) * num_vertices * 2);

        for (int i = 0; i < num_vertices; i++)
        {
            
            getline(file, line);
            istringstream stream_vertex(line);
            float x, y, z;
            if (!(stream_vertex >> x >> y >> z))
            {
                cerr << "Error reading vertex coordinates!" << endl;
                free(v);
                free(n);
                file.close();
                return;
            }
            v[i * 3] = x;
            v[i * 3 + 1] = y;
            v[i * 3 + 2] = z;

           
            getline(file, line);
            istringstream stream_normal(line);
            float nx, ny, nz;
            if (!(stream_normal >> nx >> ny >> nz))
            {
                cerr << "Error reading normal coordinates!" << endl;
                free(v);
                free(n);
                file.close();
                return;
            }
            n[i * 3] = nx;
            n[i * 3 + 1] = ny;
            n[i * 3 + 2] = nz;

            
            getline(file, line);
            istringstream stream_texture(line);
            float tx, ty;
            if (!(stream_texture >> tx >> ty))
            {
                cerr << "Error reading texture coordinates!" << endl;
                free(v);
                free(n);
                file.close();
                return;
            }
            t[i * 2] = tx;
            t[i * 2 + 1] = ty;
        }

        GLuint buffers[3];
        glGenBuffers(3, buffers);

        glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * num_vertices * 3, v, GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * num_vertices * 3, n, GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * num_vertices * 2, t, GL_STATIC_DRAW);

        free(v);
        free(n);
        free(t);

        ModelData data;
        data.vbo[0] = buffers[0];
        data.vbo[1] = buffers[1];
        data.vbo[2] = buffers[2];
        data.vertexCount = num_vertices;
        modelCache[filename] = data;
    }

    ModelData &data = modelCache[filename];

    glBindTexture(GL_TEXTURE_2D, textureID);

    glBindBuffer(GL_ARRAY_BUFFER, data.vbo[0]);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, data.vbo[1]);
    glNormalPointer(GL_FLOAT, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, data.vbo[2]);
    glTexCoordPointer(2, GL_FLOAT, 0, 0);

    glDrawArrays(GL_TRIANGLES, 0, data.vertexCount);

    glBindTexture(GL_TEXTURE_2D, 0);

}

void drawModel(Models &models, bool colorPicking = false)
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

    for (Model &m : models.model)
    {
        if (m.filetextura.empty()) {
            
            glMaterialfv(GL_FRONT, GL_DIFFUSE, m.color.diffuse);
            glMaterialfv(GL_FRONT, GL_AMBIENT, m.color.ambient);
            glMaterialfv(GL_FRONT, GL_SPECULAR, m.color.specular);
            glMaterialfv(GL_FRONT, GL_EMISSION, m.color.emissive);
            glMaterialf(GL_FRONT, GL_SHININESS, m.color.shininess);
        }

        if (colorPicking)
        {
            float color = (float)i / 255.0f;
            glColor3f(color, color, color);
            i++;
        }

        if (!m.filetextura.empty() && m.textureID == 0)
        {
            m.textureID = loadTexture(m.filetextura.c_str());
            glBindTexture(GL_TEXTURE_2D, m.textureID);
        }
       
        drawFigureVBO(m.file, m.textureID);
    }

    for (Models &child : models.models)
    {
        drawModel(child, true);
    }

    glPopMatrix();
}

unsigned char picking(int x, int y)
{
    return 0;

    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    unsigned char res[4];
    GLint viewport[4];
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(world.camera.position.x, world.camera.position.y, world.camera.position.z,
              world.camera.lookAt.x, world.camera.lookAt.y, world.camera.lookAt.z,
              world.camera.up.x, world.camera.up.y, world.camera.up.z);
    glDepthFunc(GL_LEQUAL);

    drawModel(world.models, true);
    glDepthFunc(GL_LESS);

    glGetIntegerv(GL_VIEWPORT, viewport);
    glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, res);

    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);

    return res[0];
}

void processMouseButtons(int button, int state, int xx, int yy)
{
    if (button == GLUT_LEFT_BUTTON)
    {
        if (state == GLUT_DOWN)
        {
            mouseLeftDown = true;
            mouseX = xx;
            mouseY = yy;
        }
        else if (state == GLUT_UP)
        {
            mouseLeftDown = false;
        }
    }
    else if (button == GLUT_RIGHT_BUTTON)
    {
        if (state == GLUT_DOWN)
        {
            mouseRightDown = true;
            mouseX = xx;
            mouseY = yy;
        }
        else if (state == GLUT_UP)
        {
            mouseRightDown = false;
        }
    }
    if (state == GLUT_DOWN)
    {
        if (button == GLUT_LEFT_BUTTON)
        {
            i = 1;
            figure = picking(xx, yy);
            switch (figure)
            {
            case 1:
                cout << "Sol Clicado" << endl;
                break;
            case 2:
                cout << "Mercúrio Clicado" << endl;
                break;
            case 3:
                cout << "Vénus Clicado" << endl;
                break;
            case 4:
                cout << "Terra Clicada" << endl;
                break;
            case 5:
                cout << "Lua Clicada" << endl;
                break;
            case 6:
                cout << "Marte Clicado" << endl;
                break;
            case 7:
                cout << "Jupiter Clicado" << endl;
                break;
            case 8:
                cout << "Io Clicado" << endl;
                break;
            case 9:
                cout << "Europa Clicado" << endl;
                break;
            case 10:
                cout << "Ganimedes Clicado" << endl;
                break;
            case 11:
                cout << "Calisto Clicado" << endl;
                break;
            case 12:
                cout << "Saturno Clicado" << endl;
                break;
            case 13:
                cout << "Aneis de Saturno Clicados" << endl;
                break;
            case 14:
                cout << "Urano Clicado" << endl;
                break;
            case 15:
                cout << "Neptuno Clicado" << endl;
                break;
            case 16:
                cout << "Cometa Halley Clicado" << endl;
                break;
            default:
                cout << "Nada" << endl;
            }
            glutPostRedisplay();
        }
    }
}

void processMouseMotion(int xx, int yy)
{
    if (mouseLeftDown)
    {
        // Calcula a diferença de movimento
        int deltaX = xx - mouseX;
        int deltaY = yy - mouseY;

        // Atualiza os ângulos da câmera
        alpha += deltaX * 0.005f;
        omega += deltaY * 0.005f;

        // Limita o ângulo vertical para evitar inversões
        if (omega > M_PI / 2)
            omega = M_PI / 2;
        if (omega < -M_PI / 2)
            omega = -M_PI / 2;

        // Atualiza a posição da câmera
        SphericalToCartesianFPS();
    }

    // Atualiza a posição do rato para o próximo evento
    mouseX = xx;
    mouseY = yy;

    glutPostRedisplay();
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
void luz_ativa()
{
 

    for (size_t i = 0; i < world.lights.size() && i < GL_MAX_LIGHTS; ++i)
    {
        GLenum light_id = GL_LIGHT0 + i;

        
        Light &light = world.lights[i];
        GLfloat lightPos[4] = {light.position[0], light.position[1], light.position[2], light.position[3]};
        GLfloat lightDir[4] = {light.direction[0], light.direction[1], light.direction[2], light.direction[3]};

        float diffuse[4] = {1.0f, 1.0f, 1.0f, 1.0f};
        float ambient[4] = {0.5f, 0.5f, 0.5f, 1.0f};
        float specular[4] = {1.0f, 1.0f, 1.0f, 1.0f};

        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
        glLightfv(light_id, GL_DIFFUSE, diffuse);
        glLightfv(light_id, GL_AMBIENT, ambient);
        glLightfv(light_id, GL_SPECULAR, specular);

        switch (light.type)
        {
        case Light::DIRECTIONAL:
            glLightfv(light_id, GL_POSITION, light.direction); 
            break;
        case Light::POINT:
        {
            float position[4] = {light.position[0], light.position[1], light.position[2], 1.0f};
            glLightfv(light_id, GL_POSITION, position);
            glLightf(light_id, GL_CONSTANT_ATTENUATION, 1.0f);
            glLightf(light_id, GL_LINEAR_ATTENUATION, 0.05f);
            glLightf(light_id, GL_QUADRATIC_ATTENUATION, 0.001f);
            break;
        }
        case Light::SPOT:
        {
            float position[4] = {light.position[0], light.position[1], light.position[2], 1.0f};
            glLightfv(light_id, GL_POSITION, position);
            glLightfv(light_id, GL_SPOT_DIRECTION, light.direction);
            glLightf(light_id, GL_SPOT_CUTOFF, light.cutoff);
            glLightf(light_id, GL_SPOT_EXPONENT, 10.0f);
            break;
        }

            
            if (light.type == Light::SPOT)
            {
                glLightfv(GL_LIGHT0 + i, GL_SPOT_DIRECTION, lightDir);
                glLightf(GL_LIGHT0 + i, GL_SPOT_CUTOFF, light.cutoff);
            }
        }
    }
}

void renderScene()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(world.camera.position.x, world.camera.position.y, world.camera.position.z,
              world.camera.lookAt.x, world.camera.lookAt.y, world.camera.lookAt.z,
              world.camera.up.x, world.camera.up.y, world.camera.up.z);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glEnable(GL_TEXTURE_2D);

    luz_ativa();

    drawAxis();
    
    if (solidMode)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    else
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }


    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);


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
    glutMouseFunc(processMouseButtons);
    glutMotionFunc(processMouseMotion);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    for (int i = 0; i < world.lights.size(); i++)
    {
        glEnable(GL_LIGHT0 + i);
    }

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    

    glutMainLoop();

    for (auto &entry : modelCache)
    {
        glDeleteBuffers(2, entry.second.vbo);
    }

    return 0;
}