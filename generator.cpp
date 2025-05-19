#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstring>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
ofstream open_file(char *filename)
{
    ofstream file(filename);
    if (!file)
    {
        cerr << "Error when trying to open the file" << endl;
        exit(0);
    }

    return file;
}

void normalize(float *p)
{
    float len = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    if (len > 0)
    {
        p[0] /= len;
        p[1] /= len;
        p[2] /= len;
    }
}
void cross(float *a, float *b, float *res)
{
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}

void calculate_vectors(float p1[3], float p2[3], float p3[3], float *v1, float *v2)
{
    v1[0] = p2[0] - p1[0];
    v1[1] = p2[1] - p1[1];
    v1[2] = p2[2] - p1[2];

    v2[0] = p3[0] - p1[0];
    v2[1] = p3[1] - p1[1];
    v2[2] = p3[2] - p1[2];
}

void write_points_plane(ofstream &file, int points, int points_per_row, float coordinates[][3], int divisions, float normal[3], float texture[][2])
{
    for (int i = 0; i < divisions; i++)
    {
        int points = i * points_per_row;
        for (int j = 0; j < divisions; j++)
        {
            int index = points + j;
            file << coordinates[index][0] << " " << coordinates[index][1] << " " << coordinates[index][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << texture[index][0] << " " << texture[index][1] << "\n";
            file << coordinates[index + points_per_row][0] << " " << coordinates[index + points_per_row][1] << " " << coordinates[index + points_per_row][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << texture[index + points_per_row][0] << " " << texture[index +  points_per_row][1] << "\n";
            file << coordinates[index + points_per_row + 1][0] << " " << coordinates[index + points_per_row + 1][1] << " " << coordinates[index + points_per_row + 1][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << texture[index + points_per_row + 1][0] << " " << texture[index + points_per_row + 1][1] << "\n";

            file << coordinates[index][0] << " " << coordinates[index][1] << " " << coordinates[index][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << texture[index][0] << " " << texture[index][1] << "\n";
            file << coordinates[index + points_per_row + 1][0] << " " << coordinates[index + points_per_row + 1][1] << " " << coordinates[index + points_per_row + 1][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << texture[index + points_per_row + 1][0] << " " << texture[index + points_per_row + 1][1] << "\n";
            file << coordinates[index + 1][0] << " " << coordinates[index + 1][1] << " " << coordinates[index + 1][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << texture[index + 1][0] << " " << texture[index + 1][1] << "\n";
        }
    }
}

void generate_points_plane(float coordinates[][3], float start_x, float start_y, int points, int points_per_row, float distance)
{
    for (int row = 0; row < points_per_row; row++)
    {
        float space_z = start_x + (distance * row);
        int index = row * points_per_row;

        for (int column = 0; column < points_per_row; column++)
        {
            float space_x = start_x + (distance * column);
            coordinates[index + column][0] = space_x;
            coordinates[index + column][1] = start_y;
            coordinates[index + column][2] = space_z;
        }
    }
}

void generate_points_texture_plane(float textures_coordinates[][2], int points_per_row){
    float space = 1.0f / (points_per_row - 1);
    
    float x, y;
    for(int i = 0; i < points_per_row; i++){
        y = 1 - (i * space);
        for(int j = 0; j < points_per_row; j++){
            x = (j * space);
            textures_coordinates[i * points_per_row + j][0] = x;
            textures_coordinates[i * points_per_row + j][1] = y;
        }
    }
}

void generate_rotation_matrix_x(float rotation_matrix[4][4], float angle)
{
    rotation_matrix[0][0] = rotation_matrix[3][3] = 1;
    rotation_matrix[1][1] = rotation_matrix[2][2] = cos(angle);
    rotation_matrix[2][1] = sin(angle);
    rotation_matrix[1][2] = -1 * sin(angle);
}

void generate_rotation_matrix_z(float rotation_matrix[4][4], float angle)
{
    rotation_matrix[2][2] = rotation_matrix[3][3] = 1;
    rotation_matrix[0][0] = rotation_matrix[1][1] = cos(angle);
    rotation_matrix[0][1] = -1 * sin(angle);
    rotation_matrix[1][0] = sin(angle);
}

void generate_rotation_matrix_y(float rotation_matrix[4][4], float angle)
{
    rotation_matrix[1][1] = rotation_matrix[3][3] = 1;
    rotation_matrix[0][0] = rotation_matrix[2][2] = cos(angle);
    rotation_matrix[2][0] = -1 * sin(angle);
    rotation_matrix[0][2] = sin(angle);
}

void generate_plane(float length, int divisions, char *filename)
{
    int total_points = (divisions * divisions) * 6;
    int points_per_row = divisions + 1;
    float distance = length / divisions;
    float start = -length / 2;
    ofstream file = open_file(filename);

    file << total_points << "\n";
    int points = points_per_row * points_per_row;
    float coordinates[points][3];
    float textures_coordinates[points][2];
    float normal[3];
    float vector1[3];
    float vector2[3];
    generate_points_plane(coordinates, start, 0, points, points_per_row, distance);
    generate_points_texture_plane(textures_coordinates, points_per_row);
    calculate_vectors(coordinates[0], coordinates[points_per_row], coordinates[1], vector1, vector2);
    cross(vector1, vector2, normal);
    normalize(normal);
    write_points_plane(file, points, points_per_row, coordinates, divisions, normal, textures_coordinates);

    file.close();
}

void apply_rotation(float rotation_matrix[4][4], float coordinates[][3], float coordinates_rotation[][3], int points)
{
    for (int i = 0; i < points; i++)
    {
        float x = coordinates[i][0];
        float y = coordinates[i][1];
        float z = coordinates[i][2];

        coordinates_rotation[i][0] = x * rotation_matrix[0][0] + y * rotation_matrix[0][1] + z * rotation_matrix[0][2];
        coordinates_rotation[i][1] = x * rotation_matrix[1][0] + y * rotation_matrix[1][1] + z * rotation_matrix[1][2];
        coordinates_rotation[i][2] = x * rotation_matrix[2][0] + y * rotation_matrix[2][1] + z * rotation_matrix[2][2];
    }
}

void generate_box(float length, int divisions, char *filename)
{
    int total_points = (divisions * divisions) * 36;
    int points_per_row = divisions + 1;

    ofstream file = open_file(filename);

    file << total_points << "\n";

    float start = -length / 2;

    float distance = length / divisions;

    int points = points_per_row * points_per_row;

    float coordinates[points][3];
    float textures_coordinates[points][2];
    float normal[3];
    float vector1[3];
    float vector2[3];

    generate_points_plane(coordinates, start, -start, points, points_per_row, distance);
    generate_points_texture_plane(textures_coordinates, points_per_row);
    calculate_vectors(coordinates[0], coordinates[points_per_row], coordinates[1], vector1, vector2);
    cross(vector1, vector2, normal);
    normalize(normal);
    write_points_plane(file, points, points_per_row, coordinates, divisions, normal, textures_coordinates);

    float rotation_matrix_z[4][4] = {0};
    float temp_normal[1][3] = {{normal[0], normal[1], normal[2]}};
    for (int i = 1; i < 4; i++)
    {
        generate_rotation_matrix_z(rotation_matrix_z, i * M_PI / 2);

        float coordinates_rotation[points][3];
        apply_rotation(rotation_matrix_z, coordinates, coordinates_rotation, points);

        float temp_normal_rotation[1][3];

        apply_rotation(rotation_matrix_z, temp_normal, temp_normal_rotation, 1);
        write_points_plane(file, points, points_per_row, coordinates_rotation, divisions, temp_normal_rotation[0], textures_coordinates);
    }

    float rotation_matrix_x[4][4] = {0};
    for (int i = 0; i < 2; i++)
    {
        generate_rotation_matrix_x(rotation_matrix_x, (2 * i + 1) * M_PI / 2);

        float coordinates_rotation[points][3];

        apply_rotation(rotation_matrix_x, coordinates, coordinates_rotation, points);

        float temp_normal_rotation[1][3];

        apply_rotation(rotation_matrix_x, temp_normal, temp_normal_rotation, 1);
        normal[0] = temp_normal_rotation[0][0];
        normal[1] = temp_normal_rotation[0][1];
        normal[2] = temp_normal_rotation[0][2];
        write_points_plane(file, points, points_per_row, coordinates_rotation, divisions, temp_normal_rotation[0], textures_coordinates);
    }

    file.close();
}

void writeVertexWithNormalAndTexture(float x, float y, float z, ofstream &file, float tx, float ty)
{
    float len = sqrt(x * x + y * y + z * z);
    float nx = x / len;
    float ny = y / len;
    float nz = z / len;
    file << x << " " << y << " " << z << "\n";
    file << nx << " " << ny << " " << nz << "\n";
    file << tx << " " << ty << "\n";
}

void generate_sphere(float radius, int slices, int stacks, char *filename)
{
    float slice_angle = 2 * M_PI / slices;
    float stack_angle = M_PI / stacks;

    ofstream file = open_file(filename);

    file << slices * stacks * 6 << "\n"; // 2 triângulos por quadrado

    float space_x = 1.0f / slices;
    float space_y = 1.0f / stacks;

    for (int i = 0; i < slices; i++)
    {
        float theta1 = i * slice_angle;
        float theta2 = (i + 1) * slice_angle;

        for (int j = 0; j < stacks; j++)
        {
            float phi1 = j * stack_angle;
            float phi2 = (j + 1) * stack_angle;

            // Vértices dos 4 pontos do "quadrado" atual
            float x1 = radius * sin(phi1) * cos(theta1);
            float y1 = radius * cos(phi1);
            float z1 = radius * sin(phi1) * sin(theta1);

            float x2 = radius * sin(phi1) * cos(theta2);
            float y2 = radius * cos(phi1);
            float z2 = radius * sin(phi1) * sin(theta2);

            float x3 = radius * sin(phi2) * cos(theta1);
            float y3 = radius * cos(phi2);
            float z3 = radius * sin(phi2) * sin(theta1);

            float x4 = radius * sin(phi2) * cos(theta2);
            float y4 = radius * cos(phi2);
            float z4 = radius * sin(phi2) * sin(theta2);

            // Triângulo 1
            writeVertexWithNormalAndTexture(x1, y1, z1, file, space_x * i, 1 - (space_y * j));
            writeVertexWithNormalAndTexture(x3, y3, z3, file, space_x * i, 1 - (space_y * (j + 1)));
            writeVertexWithNormalAndTexture(x4, y4, z4, file, space_x * (i + 1), 1 - (space_y * (j + 1)));

            // Triângulo 2
            writeVertexWithNormalAndTexture(x1, y1, z1, file, space_x * i, 1 - (space_y * j));
            writeVertexWithNormalAndTexture(x4, y4, z4, file, space_x * (i + 1), 1 - (space_y * (j + 1)));
            writeVertexWithNormalAndTexture(x2, y2, z2, file, space_x * (i + 1), 1 - (space_y * j));
        }
    }

    file.close();
}

void generate_cone(float radius, float height, int slices, int stacks, char *filename)
{
    ofstream file = open_file(filename);
    float angle = 2 * M_PI / slices;
    float stack_height = height / stacks;
    float stack_radius_step = radius / stacks;

    int total_points = (slices * 3) + (slices * stacks * 6) + (slices * 3);
    file << total_points << "\n";

    // Base (normal para baixo)
    for (int i = 0; i < slices; i++)
    {
        float theta1 = i * angle;
        float theta2 = (i + 1) * angle;

        file << "0 0 0\n0 -1 0\n0.5 0.5\n";                                                     // Centro da base
        file << radius * cos(theta1) << " 0 " << radius * sin(theta1) << "\n0 -1 0\n" << 0.5 + (0.5 * cos(theta1)) << " " << 0.5 + (0.5 * sin(theta1)) << "\n"; // Ponto atual
        file << radius * cos(theta2) << " 0 " << radius * sin(theta2) << "\n0 -1 0\n" << 0.5 + (0.5 * cos(theta2)) << " " << 0.5 + (0.5 * sin(theta2)) << "\n"; // Próximo ponto
    }

    // Laterais (normais inclinadas)
    for (int j = 0; j < stacks; j++)
    {
        float current_radius = radius - (j * stack_radius_step);
        float next_radius = radius - ((j + 1) * stack_radius_step);
        float current_height = j * stack_height;
        float next_height = (j + 1) * stack_height;
        float slope = radius / height;
        float ny = 1 / sqrt(1 + slope * slope); // Componente Y da normal

        for (int i = 0; i < slices; i++)
        {
            float theta1 = i * angle;
            float theta2 = (i + 1) * angle;

            // Pontos e normais
            float x1 = current_radius * cos(theta1);
            float z1 = current_radius * sin(theta1);
            
            float nx1 = cos(theta1) * slope;
            float nz1 = sin(theta1) * slope;
            float normal1[3] = {nx1, ny, nz1};
            normalize(normal1);

            float x2 = current_radius * cos(theta2);
            float z2 = current_radius * sin(theta2);
            float nx2 = cos(theta2) * slope;
            float nz2 = sin(theta2) * slope;
            float normal2[3] = {nx2, ny, nz2};
            normalize(normal2);

            float u1 = (float) i / slices;
            float u2 = (float)(i + 1) / slices;
            float v1 = (float) j /stacks;
            float v2 = (float)(j + 1)/stacks;

            // Triângulos laterais
            file << x1 << " " << current_height << " " << z1 << "\n"
                 << normal1[0] << " " << normal1[1] << " " << normal1[2] << "\n"
                 << u1 << " " << v1 << "\n";
            file << x1 * (next_radius / current_radius) << " " << next_height << " " << z1 * (next_radius / current_radius) << "\n"
                 << normal1[0] << " " << normal1[1] << " " << normal1[2] << "\n"
                 << u1 << " " << v2 << "\n";
            file << x2 * (next_radius / current_radius) << " " << next_height << " " << z2 * (next_radius / current_radius) << "\n"
                 << normal2[0] << " " << normal2[1] << " " << normal2[2] << "\n"
                 << u2 << " " << v2 << "\n";   

            file << x1 << " " << current_height << " " << z1 << "\n"
                 << normal1[0] << " " << normal1[1] << " " << normal1[2] << "\n"
                 << u1 << " " << v1 << "\n";
            file << x2 * (next_radius / current_radius) << " " << next_height << " " << z2 * (next_radius / current_radius) << "\n"
                 << normal2[0] << " " << normal2[1] << " " << normal2[2] << "\n"
                 << u2 << " " << v2 << "\n";
            file << x2 << " " << current_height << " " << z2 << "\n"
                 << normal2[0] << " " << normal2[1] << " " << normal2[2] << "\n"
                 << u2 << " " << v1 << "\n";
        }
    }

    // Topo (normais horizontais)
    for (int i = 0; i < slices; i++)
    {
        float theta = i * angle;
        float nx = cos(theta);
        float nz = sin(theta);
        float normal[3] = {nx, 0, nz};

        file << radius * cos(theta) << " 0 " << radius * sin(theta) << "\n"
             << normal[0] << " " << normal[1] << " " << normal[2] << "\n"
             << (float)i/slices << " 0\n";
        file << "0 " << height << " 0\n0 1 0\n"
             << (float) i / slices << " 1\n"; // Ponto do topo
        file << radius * cos(theta + angle) << " 0 " << radius * sin(theta + angle) << "\n"
             << nx << " 0 " << nz << "\n"
             << (float)(i + 1) / slices << " 0\n";
    }

    file.close();
}

void generate_torus(float radius, float circle_radius, int slices, int divisions, char *filename) {
    float coordinates[divisions + 1][3];
    float normals[divisions + 1][3]; // Array para normais
    float angle = 2 * M_PI / divisions;

    ofstream file = open_file(filename);
    file << divisions * slices * 6 << endl;

    // Gera o círculo inicial com normais
    float rotation_matrix_z[4][4] = {0};
    for (int i = 0; i < divisions + 1; i++) {
        generate_rotation_matrix_z(rotation_matrix_z, angle * i);
        
        // Calcula coordenadas
        coordinates[i][0] = radius + (0 * rotation_matrix_z[0][0] + circle_radius * rotation_matrix_z[0][1] + 0 * rotation_matrix_z[0][2]);
        coordinates[i][1] = 0 * rotation_matrix_z[1][0] + circle_radius * rotation_matrix_z[1][1] + 0 * rotation_matrix_z[1][2];
        coordinates[i][2] = 0 * rotation_matrix_z[2][0] + circle_radius * rotation_matrix_z[2][1] + 0 * rotation_matrix_z[2][2];
        
        // Normais iniciais (direção radial do círculo menor)
        float theta = angle * i;
        normals[i][0] = cos(theta);
        normals[i][1] = sin(theta);
        normals[i][2] = 0;
    }

    angle = 2 * M_PI / slices;

    float previous_coordinates[divisions + 1][3];
    float previous_normals[divisions + 1][3];
    memcpy(previous_coordinates, coordinates, sizeof(coordinates));
    memcpy(previous_normals, normals, sizeof(normals));

    float rotation_matrix_y[4][4] = {0};
    float space_x = 1.0f / slices;
    float space_y = 1.0f / divisions;
    for (int i = 0; i < slices; i++) {
        generate_rotation_matrix_y(rotation_matrix_y, angle * (i + 1));

        // Rotaciona coordenadas e normais
        float next_coordinates[divisions + 1][3];
        float next_normals[divisions + 1][3];
        apply_rotation(rotation_matrix_y, coordinates, next_coordinates, divisions + 1);
        apply_rotation(rotation_matrix_y, normals, next_normals, divisions + 1);

        for (int j = 0; j < divisions; j++) {
            // Triângulo 1
            file << previous_coordinates[j][0] << " " << previous_coordinates[j][1] << " " << previous_coordinates[j][2] << "\n";
            file << previous_normals[j][0] << " " << previous_normals[j][1] << " " << previous_normals[j][2] << "\n";
            file << 1 - space_x * i << " " << 1 - space_y * j << "\n";
            file << next_coordinates[j + 1][0] << " " << next_coordinates[j + 1][1] << " " << next_coordinates[j + 1][2] << "\n";
            file << next_normals[j + 1][0] << " " << next_normals[j + 1][1] << " " << next_normals[j + 1][2] << "\n";
            file << 1 - space_x * (i + 1) << " " << 1 - space_y * (j + 1) << "\n";
            file << previous_coordinates[j + 1][0] << " " << previous_coordinates[j + 1][1] << " " << previous_coordinates[j + 1][2] << "\n";
            file << previous_normals[j + 1][0] << " " << previous_normals[j + 1][1] << " " << previous_normals[j + 1][2] << "\n";
            file << 1 - space_x * i << " " << 1 - space_y * (j + 1) << "\n";

            // Triângulo 2
            file << previous_coordinates[j][0] << " " << previous_coordinates[j][1] << " " << previous_coordinates[j][2] << "\n";
            file << previous_normals[j][0] << " " << previous_normals[j][1] << " " << previous_normals[j][2] << "\n";
            file << 1 - space_x * i << " " << 1 - space_y * j << "\n";
            file << next_coordinates[j][0] << " " << next_coordinates[j][1] << " " << next_coordinates[j][2] << "\n";
            file << next_normals[j][0] << " " << next_normals[j][1] << " " << next_normals[j][2] << "\n";
            file << 1 - space_x * (i + 1) << " " << 1 - space_y * j << "\n";
            file << next_coordinates[j + 1][0] << " " << next_coordinates[j + 1][1] << " " << next_coordinates[j + 1][2] << "\n";
            file << next_normals[j + 1][0] << " " << next_normals[j + 1][1] << " " << next_normals[j + 1][2] << "\n";
            file << 1 - space_x * (i + 1) << " " << 1 - space_y * (j + 1) << "\n";
        }

        memcpy(previous_coordinates, next_coordinates, sizeof(next_coordinates));
        memcpy(previous_normals, next_normals, sizeof(next_normals));
    }

    file.close();
}
void generate_cylinder(float radius, float height, int divisions, char *filename) {
    ofstream file = open_file(filename);
    float angle = 2 * M_PI / divisions;
    file << (divisions * 6 * 2) + (divisions * 6 * 2) << "\n"; 

    // Tampas (superior e inferior)
    for (int i = 0; i < divisions; i++) {
        float theta1 = i * angle;
        float theta2 = (i + 1) * angle;

        // Tampa superior (normal para cima)
        file << "0 " << height/2 << " 0\n0 1 0\n";
        file << radius * cos(theta1) << " " << height/2 << " " << radius * sin(theta1) << "\n0 1 0\n";
        file << radius * cos(theta2) << " " << height/2 << " " << radius * sin(theta2) << "\n0 1 0\n";

        // Tampa inferior (normal para baixo)
        file << "0 " << -height/2 << " 0\n0 -1 0\n";
        file << radius * cos(theta2) << " " << -height/2 << " " << radius * sin(theta2) << "\n0 -1 0\n";
        file << radius * cos(theta1) << " " << -height/2 << " " << radius * sin(theta1) << "\n0 -1 0\n";
    }

    // Laterais (normais radiais)
    for (int i = 0; i < divisions; i++) {
        float theta1 = i * angle;
        float theta2 = (i + 1) * angle;

        // Pontos e normais
        float x1 = radius * cos(theta1);
        float z1 = radius * sin(theta1);
        float nx1 = cos(theta1);
        float nz1 = sin(theta1);

        float x2 = radius * cos(theta2);
        float z2 = radius * sin(theta2);
        float nx2 = cos(theta2);
        float nz2 = sin(theta2);

        // Triângulos laterais
        file << x1 << " " << height/2 << " " << z1 << "\n" << nx1 << " 0 " << nz1 << "\n";
        file << x1 << " " << -height/2 << " " << z1 << "\n" << nx1 << " 0 " << nz1 << "\n";
        file << x2 << " " << -height/2 << " " << z2 << "\n" << nx2 << " 0 " << nz2 << "\n";

        file << x1 << " " << height/2 << " " << z1 << "\n" << nx1 << " 0 " << nz1 << "\n";
        file << x2 << " " << -height/2 << " " << z2 << "\n" << nx2 << " 0 " << nz2 << "\n";
        file << x2 << " " << height/2 << " " << z2 << "\n" << nx2 << " 0 " << nz2 << "\n";
    }

    file.close();
}

void multMatrixMatrix(float first[4][4], float second[4][4], float result[4][4])
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < 4; k++)
                result[i][j] += first[i][k] * second[k][j];
        }
}

void multVectorMatrix(float vector[4], float matrix[4][4], float result[4])
{
    for (int i = 0; i < 4; i++)
    {
        result[i] = 0;
        for (int j = 0; j < 4; j++)
        {
            result[i] += vector[j] * matrix[j][i];
        }
    }
}

void generate_bezier_patch(char *filename, int t, char *filename_destination) {
    ifstream file_patch(filename);
    if (!file_patch) {
        cerr << "Error when trying to open the file" << endl;
        exit(0);
    }

    string line;
    getline(file_patch, line);
    int n_patches = stoi(line);

    int patches[n_patches][16];
    for (int i = 0; i < n_patches; i++) {
        getline(file_patch, line);
        char line_buffer[line.size() + 1];
        strcpy(line_buffer, line.c_str());
        char *token = strtok(line_buffer, ", ");
        int j = 0;
        while (token != NULL) {
            patches[i][j++] = stoi(token);
            token = strtok(NULL, ", ");
        }
    }

    getline(file_patch, line);
    int n_control_points = stoi(line);

    float control_points[n_control_points][3];
    for (int i = 0; i < n_control_points; i++) {
        getline(file_patch, line);
        char line_buffer[line.size() + 1];
        strcpy(line_buffer, line.c_str());
        char *token = strtok(line_buffer, ", ");
        int j = 0;
        while (token != NULL) {
            control_points[i][j++] = stof(token);
            token = strtok(NULL, ", ");
        }
    }

    file_patch.close();

    float M[4][4] = {{-1, 3, -3, 1}, {3, -6, 3, 0}, {-3, 3, 0, 0}, {1, 0, 0, 0}};
    float delta = 1.0f / t;
    float ts[t + 1];
    for (int i = 0; i <= t; i++) ts[i] = delta * i;

    ofstream file = open_file(filename_destination);
    file << t * t * 6 * n_patches << endl;

    float textures_coordinates[(t + 1) * (t + 1)][2];
    float space = 1.0f / t;

    for(int i = 0; i <= t; i++){
        for(int j = 0; j <= t; j++){
                textures_coordinates[i * (t + 1) + j][0] = space * j;
                textures_coordinates[i * (t + 1) + j][1] = 1 - (space * i);
        }
    }

    for (int i = 0; i < n_patches; i++) {
        float px[4][4], py[4][4], pz[4][4];
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                int idx = patches[i][4 * j + k];
                px[k][j] = control_points[idx][0];
                py[k][j] = control_points[idx][1];
                pz[k][j] = control_points[idx][2];
            }
        }

        float aux[4][4], mx[4][4], my[4][4], mz[4][4];
        multMatrixMatrix(M, px, aux); multMatrixMatrix(aux, M, mx);
        multMatrixMatrix(M, py, aux); multMatrixMatrix(aux, M, my);
        multMatrixMatrix(M, pz, aux); multMatrixMatrix(aux, M, mz);

        // Arrays para pontos e normais
        float final_points[(t + 1) * (t + 1)][3];
        float final_normals[(t + 1) * (t + 1)][3];

        for (int j = 0; j <= t; j++) {
            for (int k = 0; k <= t; k++) {
                float u[4] = {powf(ts[j], 3), powf(ts[j], 2), ts[j], 1};
                float v[4] = {powf(ts[k], 3), powf(ts[k], 2), ts[k], 1};

                // Calcula pontos
                float aux_x[4], aux_y[4], aux_z[4];
                multVectorMatrix(u, mx, aux_x);
                multVectorMatrix(u, my, aux_y);
                multVectorMatrix(u, mz, aux_z);

                final_points[j*(t+1)+k][0] = aux_x[0]*v[0] + aux_x[1]*v[1] + aux_x[2]*v[2] + aux_x[3]*v[3];
                final_points[j*(t+1)+k][1] = aux_y[0]*v[0] + aux_y[1]*v[1] + aux_y[2]*v[2] + aux_y[3]*v[3];
                final_points[j*(t+1)+k][2] = aux_z[0]*v[0] + aux_z[1]*v[1] + aux_z[2]*v[2] + aux_z[3]*v[3];

                // Calcula derivadas parciais para normais
                float du[4] = {3*ts[j]*ts[j], 2*ts[j], 1, 0}; // Derivada de u
                float dv[4] = {3*ts[k]*ts[k], 2*ts[k], 1, 0}; // Derivada de v

                float du_x[4], du_y[4], du_z[4];
                multVectorMatrix(du, mx, du_x);
                multVectorMatrix(du, my, du_y);
                multVectorMatrix(du, mz, du_z);

                float dv_x[4], dv_y[4], dv_z[4];
                multVectorMatrix(dv, mx, dv_x);
                multVectorMatrix(dv, my, dv_y);
                multVectorMatrix(dv, mz, dv_z);

                // Vetores tangentes
                float tangent_u[3] = {
                    du_x[0]*v[0] + du_x[1]*v[1] + du_x[2]*v[2] + du_x[3]*v[3],
                    du_y[0]*v[0] + du_y[1]*v[1] + du_y[2]*v[2] + du_y[3]*v[3],
                    du_z[0]*v[0] + du_z[1]*v[1] + du_z[2]*v[2] + du_z[3]*v[3]
                };

                float tangent_v[3] = {
                    aux_x[0]*dv[0] + aux_x[1]*dv[1] + aux_x[2]*dv[2] + aux_x[3]*dv[3],
                    aux_y[0]*dv[0] + aux_y[1]*dv[1] + aux_y[2]*dv[2] + aux_y[3]*dv[3],
                    aux_z[0]*dv[0] + aux_z[1]*dv[1] + aux_z[2]*dv[2] + aux_z[3]*dv[3]
                };

                // Calcula e normaliza a normal
                float normal[3];
                cross(tangent_u, tangent_v, normal);
                normalize(normal);

                final_normals[j*(t+1)+k][0] = normal[0];
                final_normals[j*(t+1)+k][1] = normal[1];
                final_normals[j*(t+1)+k][2] = normal[2];
            }
        }

        // Escreve os triângulos com normais
        for (int j = 0; j < t; j++) {
            for (int k = 0; k < t; k++) {
                int idx = j*(t+1)+k;

                // Triângulo 1
                file << final_points[idx][0] << " " << final_points[idx][1] << " " << final_points[idx][2] << "\n";
                file << final_normals[idx][0] << " " << final_normals[idx][1] << " " << final_normals[idx][2] << "\n";
                file << textures_coordinates[idx][0] << " " << textures_coordinates[idx][1] << "\n";
                file << final_points[idx + t + 1][0] << " " << final_points[idx + t + 1][1] << " " << final_points[idx + t + 1][2] << "\n";
                file << final_normals[idx + t + 1][0] << " " << final_normals[idx + t + 1][1] << " " << final_normals[idx + t + 1][2] << "\n";
                file << textures_coordinates[idx + t + 1][0] << " " << textures_coordinates[idx + t + 1][1] << "\n";
                file << final_points[idx + t + 2][0] << " " << final_points[idx + t + 2][1] << " " << final_points[idx + t + 2][2] << "\n";
                file << final_normals[idx + t + 2][0] << " " << final_normals[idx + t + 2][1] << " " << final_normals[idx + t + 2][2] << "\n";
                file << textures_coordinates[idx + t + 2][0] << " " << textures_coordinates[idx + t + 2][1] << "\n";

                // Triângulo 2
                file << final_points[idx][0] << " " << final_points[idx][1] << " " << final_points[idx][2] << "\n";
                file << final_normals[idx][0] << " " << final_normals[idx][1] << " " << final_normals[idx][2] << "\n";
                file << textures_coordinates[idx][0] << " " << textures_coordinates[idx][1] << "\n";
                file << final_points[idx + t + 2][0] << " " << final_points[idx + t + 2][1] << " " << final_points[idx + t + 2][2] << "\n";
                file << final_normals[idx + t + 2][0] << " " << final_normals[idx + t + 2][1] << " " << final_normals[idx + t + 2][2] << "\n";
                file << textures_coordinates[idx + t + 2][0] << " " << textures_coordinates[idx + t + 2][1] << "\n";
                file << final_points[idx + 1][0] << " " << final_points[idx + 1][1] << " " << final_points[idx + 1][2] << "\n";
                file << final_normals[idx + 1][0] << " " << final_normals[idx + 1][1] << " " << final_normals[idx + 1][2] << "\n";
                file << textures_coordinates[idx + 1][0] << " " << textures_coordinates[idx + 1][1] << "\n";
            }
        }
    }

    file.close();
}

int main(int argc, char **argv)
{

    if (strcmp(argv[1], "box") == 0 && argc == 5)
    {
        float length = stof(argv[2]);
        int divisions = stoi(argv[3]);

        if (length <= 0 || divisions <= 0)
        {
            cout << "Length and divisions must be greater than 0!" << endl;
        }
        else
        {
            generate_box(length, divisions, argv[4]);
        }
    }
    else if (strcmp(argv[1], "sphere") == 0 && argc == 6)
    {
        float radius = stof(argv[2]);
        int slices = stoi(argv[3]);
        int stacks = stoi(argv[4]);

        if (radius <= 0 || slices <= 0 || stacks <= 0)
        {
            cout << "Radius, slices and stacks must be greater than 0!" << endl;
        }
        else
        {
            generate_sphere(radius, slices, stacks, argv[5]);
        }
    }
    else if (strcmp(argv[1], "cone") == 0 && argc == 7)
    {
        float radius = stof(argv[2]);
        float height = stof(argv[3]);
        int slices = stoi(argv[4]);
        int stacks = stoi(argv[5]);

        if (radius <= 0 || height <= 0 || slices <= 0 || stacks <= 0)
        {
            cout << "Radius, height, slices and stacks must be greater than 0!" << endl;
        }
        else
        {
            generate_cone(radius, height, slices, stacks, argv[6]);
        }
    }
    else if (strcmp(argv[1], "plane") == 0 && argc == 5)
    {
        float length = stof(argv[2]);
        int divisions = stoi(argv[3]);

        if (length <= 0 || divisions <= 0)
        {
            cout << "Length and divisions must be greater than 0!" << endl;
        }
        else
        {
            generate_plane(length, divisions, argv[4]);
        }
    }
    else if (strcmp(argv[1], "torus") == 0 && argc == 7)
    {
        float radius = stof(argv[2]);
        float circle_radius = stof(argv[3]);
        float slices = stoi(argv[4]);
        float divisions = stoi(argv[5]);

        if (radius < circle_radius || slices < 2 || divisions < 2)
        {
            cout << "Radius must be greater than circle_radius and slices and divisions must be greater than 2" << endl;
        }
        else
        {
            generate_torus(radius, circle_radius, slices, divisions, argv[6]);
        }
    }
    else if (strcmp(argv[1], "cylinder") == 0 && argc == 6)
    {
        float radius = stof(argv[2]);
        float height = stof(argv[3]);
        float divisions = stof(argv[4]);

        if (radius <= 0 || height <= 0 || divisions < 2)
        {
            cout << "Radius and height must be greater than 0 and divisions must be greater than 2" << endl;
        }
        else
        {
            generate_cylinder(radius, height, divisions, argv[5]);
        }
    }
    else if (strcmp(argv[1], "bezier_patch") == 0 && argc == 5)
    {
        int t = stoi(argv[3]);

        generate_bezier_patch(argv[2], t, argv[4]);
    }
    else
    {
        cout << "Something went wrong" << endl;
    }

    return 0;
}
