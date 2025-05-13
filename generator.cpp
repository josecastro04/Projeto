#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstring>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
ofstream open_file(char *filename){
    ofstream file(filename);
    if(!file){
        cerr << "Error when trying to open the file" << endl;
        exit(0);
    }

    return file;
}

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

void calculate_vectors(float p1[3], float p2[3], float p3[3], float *v1, float *v2){
    v1[0] = p2[0] - p1[0];
    v1[1] = p2[1] - p1[1];
    v1[2] = p2[2] - p1[2];

    v2[0] = p3[0] - p1[0];
    v2[1] = p3[1] - p1[1];
    v2[2] = p3[2] - p1[2];
}

void write_points_plane(ofstream &file, int points, int points_per_row, float coordinates[][3], int divisions, float normal[3])
{
    for (int i = 0; i < divisions; i++){
        int points = i * points_per_row;
        for(int j = 0; j < divisions; j++){
            int index = points + j;
            file << coordinates[index][0] << " " << coordinates[index][1] << " " << coordinates[index][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << coordinates[index + points_per_row][0] << " " << coordinates[index + points_per_row][1] << " " << coordinates[index + points_per_row][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << coordinates[index + points_per_row + 1][0] << " " << coordinates[index + points_per_row + 1][1] << " " << coordinates[index + points_per_row + 1][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";

            file << coordinates[index][0] << " " << coordinates[index][1] << " " << coordinates[index][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << coordinates[index + points_per_row + 1][0] << " " << coordinates[index + points_per_row + 1][1] << " " << coordinates[index + points_per_row + 1][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << coordinates[index + 1][0] << " " << coordinates[index + 1][1] << " " << coordinates[index + 1][2] << "\n";
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        }
    }
}

void generate_points_plane(float coordinates[][3], float start_x, float start_y, int points, int points_per_row , float distance)
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
    float start = - length / 2;
    ofstream file = open_file(filename);

    file << total_points << "\n";
    int points = points_per_row * points_per_row;
    float coordinates[points][3];
    float normal[3];
    float vector1[3];
    float vector2[3];
    generate_points_plane(coordinates, start, 0, points, points_per_row, distance);
    calculate_vectors(coordinates[0], coordinates[points_per_row], coordinates[1], vector1, vector2);
    cross(vector1, vector2, normal);
    normalize(normal);
    write_points_plane(file, points, points_per_row, coordinates, divisions, normal);

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
    float normal[3];
    float vector1[3];
    float vector2[3];

    generate_points_plane(coordinates, start, -start, points, points_per_row, distance);
    calculate_vectors(coordinates[0], coordinates[points_per_row], coordinates[1], vector1, vector2);
    cross(vector1, vector2, normal);
    normalize(normal);
    write_points_plane(file, points, points_per_row, coordinates, divisions, normal);

    float rotation_matrix_z[4][4] = {0};
    float temp_normal[1][3] =  {{normal[0], normal[1], normal[2]}};
    for (int i = 1; i < 4; i++){
        generate_rotation_matrix_z(rotation_matrix_z, i * M_PI / 2);

        float coordinates_rotation[points][3];
        apply_rotation(rotation_matrix_z, coordinates, coordinates_rotation, points);
        
        float temp_normal_rotation[1][3];

        apply_rotation(rotation_matrix_z, temp_normal, temp_normal_rotation, 1);
        write_points_plane(file, points, points_per_row, coordinates_rotation, divisions, temp_normal_rotation[0]);
    }

    float rotation_matrix_x[4][4] = {0};
    for (int i = 0; i < 2; i++){
        generate_rotation_matrix_x(rotation_matrix_x, (2 * i + 1) * M_PI / 2);

        float coordinates_rotation[points][3];

        apply_rotation(rotation_matrix_x, coordinates, coordinates_rotation, points);

        float temp_normal_rotation[1][3];

        apply_rotation(rotation_matrix_x, temp_normal, temp_normal_rotation, 1);
        normal[0] = temp_normal_rotation[0][0];
        normal[1] = temp_normal_rotation[0][1];
        normal[2] = temp_normal_rotation[0][2];
        write_points_plane(file, points, points_per_row, coordinates_rotation, divisions, temp_normal_rotation[0]);
    }

    file.close();
}

void generate_sphere(float radius, int slices, int stacks, char *filename)
{
    float x1, x2, x3, x4, y1, y2, z1, z2, z3, z4;
    float slices_angle = 2 * M_PI / slices;
    float stacks_angle = M_PI / stacks;

    ofstream file = open_file(filename);

    file << slices * stacks * 6 << "\n";

    for (int i = 0; i < slices; i++)
    {
        for (int j = 0; j < stacks; j++)
        {
            x1 = radius * sin(i * slices_angle) * cos(M_PI_2 - stacks_angle * j);
            x2 = radius * sin((i + 1) * slices_angle) * cos(M_PI_2 - stacks_angle * j);
            x3 = radius * sin(i * slices_angle) * cos(M_PI_2 - stacks_angle * (j + 1));
            x4 = radius * sin((i + 1) * slices_angle) * cos(M_PI_2 - stacks_angle * (j + 1));

            y1 = radius * sin(M_PI_2 - stacks_angle * j);
            y2 = radius * sin(M_PI_2 - stacks_angle * (j + 1));

            z1 = radius * cos(i * slices_angle) * cos(M_PI_2 - stacks_angle * j);
            z2 = radius * cos((i + 1) * slices_angle) * cos(M_PI_2 - stacks_angle * j);
            z3 = radius * cos(i * slices_angle) * cos(M_PI_2 - stacks_angle * (j + 1));
            z4 = radius * cos((i + 1) * slices_angle) * cos(M_PI_2 - stacks_angle * (j + 1));

            float vector1[3];
            float vector2[3];
            float normal1[3];
            float normal2[3];
            float normal3[3];
            float normal4[3];
            float point1[3];
            float point2[3];
            float point3[3];

            point1[0] = x1; point1[1] = y1; point1[2] = z1;
            point2[0] = x3; point2[1] = y2; point2[2] = z3;
            point3[0] = x2; point3[1] = y1; point3[2] = z2;

            calculate_vectors(point1, point2, point3, vector1, vector2);
            cross(vector1, vector2, normal1);
            normalize(normal1);

            point1[0] = x2; point1[1] = y1; point1[2] = z2;
            point2[0] = x1; point2[1] = y1; point2[2] = z1;
            point3[0] = x4; point3[1] = y2; point3[2] = z4;

            calculate_vectors(point1, point2, point3, vector1, vector2);
            cross(vector1, vector2, normal2);
            normalize(normal2);

            point1[0] = x3; point1[1] = y2; point1[2] = z3;
            point2[0] = x4; point2[1] = y2; point2[2] = z4;
            point3[0] = x1; point3[1] = y1; point3[2] = z1;

            calculate_vectors(point1, point2, point3, vector1, vector2);
            cross(vector1, vector2, normal3);
            normalize(normal3);

            point1[0] = x4; point1[1] = y2; point1[2] = z4;
            point2[0] = x2; point2[1] = y1; point2[2] = z2;
            point3[0] = x3; point3[1] = y2; point3[2] = z3;

            calculate_vectors(point1, point2, point3, vector1, vector2);
            cross(vector1, vector2, normal4);
            normalize(normal4);

            file << x1 << " " << y1 << " " << z1 << "\n";
            file << normal1[0] << " " << normal1[1] << " " << normal1[2] << "\n";
            file << x3 << " " << y2 << " " << z3 << "\n";
            file << normal3[0] << " " << normal3[1] << " " << normal3[2] << "\n";
            file << x4 << " " << y2 << " " << z4 << "\n";
            file << normal4[0] << " " << normal4[1] << " " << normal4[2] << "\n";

            file << x1 << " " << y1 << " " << z1 << "\n";
            file << normal1[0] << " " << normal1[1] << " " << normal1[2] << "\n";
            file << x4 << " " << y2 << " " << z4 << "\n";
            file << normal4[0] << " " << normal4[1] << " " << normal4[2] << "\n";
            file << x2 << " " << y1 << " " << z2 << "\n";
            file << normal2[0] << " " << normal2[1] << " " << normal4[2] << "\n";
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

    int total_points = (slices * 6) + (slices * stacks * 6) + (slices * 3);
    file << total_points << "\n";

    
    float base_coordinates[slices][3];
    float rotation_matrix_y[4][4] = {0};
    for (int i = 0; i < slices; i++){
        generate_rotation_matrix_y(rotation_matrix_y, i * angle);

        base_coordinates[i][0] = radius * rotation_matrix_y[2][0];
        base_coordinates[i][1] = 0;
        base_coordinates[i][2] = radius * rotation_matrix_y[0][0];
    }

    
    for (int i = 0; i < slices; i++)
    {
        int next = (i + 1) % slices;
        file << "0 0 0\n";
        file << base_coordinates[i][0] << " 0 " << base_coordinates[i][2] << "\n";
        file << base_coordinates[next][0] << " 0 " << base_coordinates[next][2] << "\n";
    }

    float next_coordinates[slices][3];
    
    for (int j = 0; j < stacks; j++)
    {
        float current_radius = radius - (stack_radius_step * j);
        float next_radius = radius - (stack_radius_step * (j + 1));
        float current_height = stack_height * j;
        float next_height = stack_height * (j + 1);

        float current_coordinates[slices][3];
        

        for (int i = 0; i < slices; i++)
        {
            generate_rotation_matrix_y(rotation_matrix_y, i * angle);

            current_coordinates[i][0] = current_radius * rotation_matrix_y[2][0];
            current_coordinates[i][1] = current_height;
            current_coordinates[i][2] = current_radius * rotation_matrix_y[0][0];

            next_coordinates[i][0] = next_radius * rotation_matrix_y[2][0];
            next_coordinates[i][1] = next_height;
            next_coordinates[i][2] = next_radius * rotation_matrix_y[0][0];
        }

        for (int i = 0; i < slices; i++)
        {
            int next = (i + 1) % slices;
            
            file << current_coordinates[i][0] << " " << current_coordinates[i][1] << " " << current_coordinates[i][2] << "\n";
            file << next_coordinates[i][0] << " " << next_coordinates[i][1] << " " << next_coordinates[i][2] << "\n";
            file << next_coordinates[next][0] << " " << next_coordinates[next][1] << " " << next_coordinates[next][2] << "\n";
 
            file << current_coordinates[i][0] << " " << current_coordinates[i][1] << " " << current_coordinates[i][2] << "\n";
            file << next_coordinates[next][0] << " " << next_coordinates[next][1] << " " << next_coordinates[next][2] << "\n";
            file << current_coordinates[next][0] << " " << current_coordinates[next][1] << " " << current_coordinates[next][2] << "\n";
        }
    }


    for (int i = 0; i < slices; i++)
    {
        int next = (i + 1) % slices;
        
        
        file << next_coordinates[i][0] << " " << next_coordinates[i][1] << " " << next_coordinates[i][2] << "\n";
        file << "0 " << height << " 0\n";
        file << next_coordinates[next][0] << " " << next_coordinates[next][1] << " " << next_coordinates[next][2] << "\n";
    }
    

    file.close();
}

void generate_torus(float radius, float circle_radius, int slices, int divisions, char *filename){
    float coordinates[divisions + 1][3];
    float angle = 2 * M_PI / divisions;

    ofstream file = open_file(filename);
    file << divisions * slices * 6 << endl;

    float rotation_matrix_z[4][4] = {0};
    for(int i = 0; i < divisions + 1; i++){
        generate_rotation_matrix_z(rotation_matrix_z, angle * i);

        coordinates[i][0] = radius  + (0 * rotation_matrix_z[0][0]) + circle_radius * rotation_matrix_z[0][1] + 0 * rotation_matrix_z[0][2];
        coordinates[i][1] = 0 * rotation_matrix_z[1][0] + circle_radius * rotation_matrix_z[1][1] + 0 * rotation_matrix_z[1][2];
        coordinates[i][2] = 0 * rotation_matrix_z[2][0] + circle_radius * rotation_matrix_z[2][1] + 0 * rotation_matrix_z[2][2];

    }
    
    angle = 2 * M_PI / slices;

    float previous_coordinates[divisions + 1][3];
    memcpy(previous_coordinates, coordinates, sizeof(coordinates));
    float rotation_matrix_y[4][4] = {0};
    float vector1[3];
    float vector2[3];
    float normal[3];
    for(int i = 0; i < slices; i++){
        generate_rotation_matrix_y(rotation_matrix_y, angle * (i + 1));

        float next_coordinates[divisions + 1][3];
        apply_rotation(rotation_matrix_y, coordinates, next_coordinates, divisions + 1);  

        for(int j = 0; j < divisions; j++){
            file << previous_coordinates[j][0] << " " << previous_coordinates[j][1] << " " << previous_coordinates[j][2] << "\n";
            calculate_vectors(previous_coordinates[j], previous_coordinates[j + 1], next_coordinates[j]);
            croos(vector1, vector2, normal);
            normalize(normal);
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << next_coordinates[j + 1][0] << " " << next_coordinates[j + 1][1] << " " << next_coordinates[j + 1][2] << "\n";
            calculate_vectors(next_coordinates[j + 1], previous_coordinates[j + 1], next_coordinates[j]);
            croos(vector1, vector2, normal);
            normalize(normal);
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << previous_coordinates[j + 1][0] << " " << previous_coordinates[j + 1][1] << " " << previous_coordinates[j + 1][2] << "\n";
            calculate_vectors(previous_coordinates[j + 1], previous_coordinates[j + 1], next_coordinates[j + 1]);
            croos(vector1, vector2, normal);
            normalize(normal);
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";

            file << previous_coordinates[j][0] << " " << previous_coordinates[j][1] << " " << previous_coordinates[j][2] << "\n";
            calculate_vectors(previous_coordinates[j], previous_coordinates[j + 1], next_coordinates[j]);
            croos(vector1, vector2, normal);
            normalize(normal);
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << next_coordinates[j][0] << " " << next_coordinates[j][1] << " " << next_coordinates[j][2] << "\n";
            calculate_vectors(new_coordinates[j], new_coordinates[j + 1], previous_coordinates[j]);
            croos(vector1, vector2, normal);
            normalize(normal);
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            file << next_coordinates[j + 1][0] << " " << next_coordinates[j + 1][1] << " " << next_coordinates[j + 1][2] << "\n";
            calculate_vectors(next_coordinates[j], previous_coordinates[j + 1], next_coordinates[j]);
            croos(vector1, vector2, normal);
            normalize(normal);
            file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        }
        memcpy(previous_coordinates, next_coordinates, sizeof(next_coordinates));

    }

    file.close();
}   

void generate_cylinder(float radius, float height, int divisions, char *filename){
    ofstream file = open_file(filename);
    
    float angle = 2 * M_PI / divisions;
    float coordinates[2][3] = {{radius, height / 2, 0}, {radius, -height / 2, 0}};
    file << (3 * divisions) * 2 + 6 * divisions << "\n";
    
    float rotation_matrix_y[4][4] = {0};
    float previous_coordinates[2][3];

    memcpy(previous_coordinates, coordinates, sizeof(coordinates));
    for(int i = 1; i <= divisions; i++){
        float new_coordinates[2][3];
        generate_rotation_matrix_y(rotation_matrix_y, angle * i);
        apply_rotation(rotation_matrix_y, coordinates, new_coordinates, 2);

        file << previous_coordinates[0][0] << " " << previous_coordinates[0][1] << " " << previous_coordinates[0][2] << "\n";
        file << new_coordinates[0][0] << " " << new_coordinates[0][1] << " " << new_coordinates[0][2] << "\n";
        file << 0 << " " << previous_coordinates[0][1] << " " << 0 << "\n";

        file << new_coordinates[1][0] << " " << new_coordinates[1][1] << " " << new_coordinates[1][2] << "\n";
        file << previous_coordinates[1][0] << " " << previous_coordinates[1][1] << " " << previous_coordinates[1][2] << "\n";
        file << 0 << " " << new_coordinates[1][1] << " " << 0 << "\n";

        file << previous_coordinates[0][0] << " " << previous_coordinates[0][1] << " " << previous_coordinates[0][2] << "\n";
        file << previous_coordinates[1][0] << " " << previous_coordinates[1][1] << " " << previous_coordinates[1][2] << "\n";
        file << new_coordinates[1][0] << " " << new_coordinates[1][1] << " " << new_coordinates[1][2] << "\n";

        file << previous_coordinates[0][0] << " " << previous_coordinates[0][1] << " " << previous_coordinates[0][2] << "\n";
        file << new_coordinates[1][0] << " " << new_coordinates[1][1] << " " << new_coordinates[1][2] << "\n";
        file << new_coordinates[0][0] << " " << new_coordinates[0][1] << " " << new_coordinates[0][2] << "\n";

        memcpy(previous_coordinates, new_coordinates, sizeof(new_coordinates));
    }  
    
    file.close();
}

void multMatrixMatrix(float first[4][4], float second[4][4], float result[4][4]){
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++){
            result[i][j] = 0;
            for(int k = 0; k < 4; k++)
                result[i][j] += first[i][k] * second[k][j];
        }
}

void multVectorMatrix(float vector[4], float matrix[4][4], float result[4]){
    for(int i = 0; i < 4; i++){
        result[i] = 0;
        for(int j = 0; j < 4; j++){
            result[i] += vector[j] * matrix[j][i];
        }
    }
}

void generate_bezier_patch(char* filename, int t, char* filename_destination){

    ifstream file_patch(filename);
    if(!file_patch){
        cerr << "Error when trying to open the file" << endl;
        exit(0);
    }
    string line;
    getline(file_patch, line);

    int n_patches = stoi(line);

    int patches[n_patches][16];

    for(int i = 0; i < n_patches; i++){
        getline(file_patch, line);
        char line_buffer[line.size() + 1];
        strcpy(line_buffer, line.c_str());
        char* token = strtok(line_buffer, ", ");
        int j = 0;
        while(token != NULL){
            patches[i][j++] = stoi(token); 
            token = strtok(NULL, ", ");
        }
    }

    getline(file_patch, line);
    int n_control_points = stoi(line);

    float control_points[n_control_points][3];

    for(int i = 0; i < n_control_points; i++){
        getline(file_patch, line);
        char line_buffer[line.size() + 1];
        strcpy(line_buffer, line.c_str());
        char* token = strtok(line_buffer, ", ");
        int j = 0;
        while(token != NULL){
            control_points[i][j++] = stof(token); 
            token = strtok(NULL, ", ");
        }
    }

    file_patch.close();

    float M[4][4] = {{-1, 3, -3, 1},
                  {3, -6, 3, 0},
                  {-3, 3, 0, 0},
                  {1, 0, 0, 0,}
                };

    float delta = 1.0f / t;
    float ts[t + 1];
    for(int i = 0; i <= t; i++){
        ts[i] = delta * i;
    }

    ofstream file = open_file(filename_destination);
    file << t * t * 6 * n_patches << endl;
    int index;
    for(int i = 0; i < n_patches; i++){
        float px[4][4];
        float py[4][4];
        float pz[4][4];

        for(int j = 0; j < 4; j++){
            for(int k = 0; k < 4; k++){
                px[k][j] = control_points[patches[i][4 * j + k]][0];
                py[k][j] = control_points[patches[i][4 * j + k]][1];
                pz[k][j] = control_points[patches[i][4 * j + k]][2];
            }
        }

        float aux[4][4];
        float mx[4][4];
        float my[4][4];
        float mz[4][4];
        multMatrixMatrix(M, px, aux);
        multMatrixMatrix(aux, M, mx);

        multMatrixMatrix(M, py, aux);
        multMatrixMatrix(aux, M, my);

        multMatrixMatrix(M, pz, aux);
        multMatrixMatrix(aux, M, mz);

        float final_points[(t + 1) * (t + 1)][3] = {0};

        for(int j = 0; j <= t; j++){
            for(int k = 0; k <= t; k++){
                index = j * (t + 1) + k;
                float u[4] = {powf(ts[j],3), powf(ts[j],2), ts[j], 1};
                float v[4] = {powf(ts[k],3), powf(ts[k],2), ts[k], 1};

                float aux_x[4];
                float aux_y[4];
                float aux_z[4];
                multVectorMatrix(u, mx, aux_x);
                multVectorMatrix(u, my, aux_y);
                multVectorMatrix(u, mz, aux_z);

                final_points[index][0] = aux_x[0] * v[0] + aux_x[1] * v[1] + aux_x[2] * v[2] + aux_x[3] * v[3];
                final_points[index][1] = aux_y[0] * v[0] + aux_y[1] * v[1] + aux_y[2] * v[2] + aux_y[3] * v[3];
                final_points[index][2] = aux_z[0] * v[0] + aux_z[1] * v[1] + aux_z[2] * v[2] + aux_z[3] * v[3];
            }
        }

        for(int j = 0; j < t; j++){
            for(int k = 0; k < t; k++){
                file << final_points[j * (t + 1) + k][0] << " " << final_points[j * (t + 1) + k][1] << " " << final_points[j * (t + 1) + k][2] << endl;
                file << final_points[(j + 1) * (t + 1) + k][0] << " " << final_points[(j + 1) * (t + 1) + k][1] << " " << final_points[(j + 1) * (t + 1) + k][2] << endl;
                file << final_points[(j + 1) * (t + 1) + k + 1][0] << " " << final_points[(j + 1) * (t + 1) + k + 1][1] << " " << final_points[(j + 1) * (t + 1) + k + 1][2] << endl;

                file << final_points[j * (t + 1) + k][0] << " " << final_points[j * (t + 1) + k][1] << " " << final_points[j * (t + 1) + k][2] << endl;
                file << final_points[(j + 1) * (t + 1) + k + 1][0] << " " << final_points[(j + 1) * (t + 1) + k + 1][1] << " " << final_points[(j + 1) * (t + 1) + k + 1][2] << endl;
                file << final_points[j * (t + 1) + k + 1][0] << " " << final_points[j * (t + 1) + k + 1][1] << " " << final_points[j * (t + 1) + k + 1][2] << endl;
            }
        }
    }
    
    file.close();
}

int main(int argc, char **argv){

    if (strcmp(argv[1], "box") == 0 && argc == 5)
    {
        float length = stof(argv[2]);
        int divisions = stoi(argv[3]);
        
        if(length <= 0 || divisions <= 0){
            cout << "Length and divisions must be greater than 0!" << endl;
        }else{
            generate_box(length, divisions, argv[4]);
        }
    }
    else if (strcmp(argv[1], "sphere") == 0 && argc == 6)
    {
        float radius = stof(argv[2]);
        int slices = stoi(argv[3]);
        int stacks = stoi(argv[4]);

        if(radius <= 0 || slices <= 0 || stacks <= 0) {
            cout << "Radius, slices and stacks must be greater than 0!" << endl;
        }else{
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
        
        if(length <= 0 || divisions <= 0){
            cout << "Length and divisions must be greater than 0!" << endl;
        }else{
            generate_plane(length, divisions, argv[4]);
        }
    }else if(strcmp(argv[1], "torus") == 0 && argc == 7){
        float radius = stof(argv[2]);
        float circle_radius = stof(argv[3]);
        float slices = stoi(argv[4]);
        float divisions = stoi(argv[5]);
        
        if(radius < circle_radius || slices < 2 || divisions < 2){
            cout << "Radius must be greater than circle_radius and slices and divisions must be greater than 2" << endl;
        }else{
            generate_torus(radius, circle_radius, slices, divisions, argv[6]);
        }
    }else if(strcmp(argv[1], "cylinder") == 0  && argc == 6){
        float radius = stof(argv[2]);
        float height = stof(argv[3]);
        float divisions = stof(argv[4]);

        if(radius <= 0 || height <= 0 || divisions < 2){
            cout << "Radius and height must be greater than 0 and divisions must be greater than 2" << endl;
        }else{
            generate_cylinder(radius, height, divisions, argv[5]);
        }
    }else if(strcmp(argv[1], "bezier_patch") == 0 && argc == 5){
        int t = stoi(argv[3]);

        generate_bezier_patch(argv[2], t, argv[4]);
    }else{
        cout << "Something went wrong" << endl;
    }

    return 0;
}
