#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstring>

#define _USE_MATH_DEFINES
#include <math.h>
std::ofstream open_file(char *filename){
    std::ofstream file(filename);
    if(!file){
         std::cerr << "Error when trying to open the file" << std::endl;
        exit(0);
    }

    return file;
}
void write_points_plane(std::ofstream &file, int points, int points_per_row, float coordinates[][3], int divisions)
{
    for (int i = 0; i < points - points_per_row; i++)
    {
        if (i % points_per_row == divisions)
            continue;

        file << coordinates[i][0] << " " << coordinates[i][1] << " " << coordinates[i][2] << "\n";
        file << coordinates[i + points_per_row][0] << " " << coordinates[i + points_per_row][1] << " " << coordinates[i + points_per_row][2] << "\n";
        file << coordinates[i + points_per_row + 1][0] << " " << coordinates[i + points_per_row + 1][1] << " " << coordinates[i + points_per_row + 1][2] << "\n";

        file << coordinates[i][0] << " " << coordinates[i][1] << " " << coordinates[i][2] << "\n";
        file << coordinates[i + points_per_row + 1][0] << " " << coordinates[i + points_per_row + 1][1] << " " << coordinates[i + points_per_row + 1][2] << "\n";
        file << coordinates[i + 1][0] << " " << coordinates[i + 1][1] << " " << coordinates[i + 1][2] << "\n";
    }
}

void generate_points_plane(float coordinates[][3], float start_x, float start_y, int points, int divisions, float distance)
{
    for (int row = 0; row <= divisions; row++)
    {
        float space_z = start_x + (distance * row);
        int index = row * (divisions + 1);

        for (int column = 0; column <= divisions; column++)
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
    std::ofstream file = open_file(filename);

    file << total_points << "\n";
    int points = points_per_row * points_per_row;
    float coordinates[points][3];

    generate_points_plane(coordinates, start, 0, points, divisions, distance);
    write_points_plane(file, points, points_per_row, coordinates, divisions);

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

    std::ofstream file = open_file(filename);

    file << total_points << "\n";

    float start = -length / 2;

    float distance = length / divisions;

    int points = points_per_row * points_per_row;

    float coordinates[points][3];

    generate_points_plane(coordinates, start, -start, points, divisions, distance);
    write_points_plane(file, points, points_per_row, coordinates, divisions);

    float rotation_matrix_z[4][4] = {0};

    for (int i = 1; i < 4; i++){
        generate_rotation_matrix_z(rotation_matrix_z, i * M_PI / 2);

        float coordinates_rotation[points][3];

        apply_rotation(rotation_matrix_z, coordinates, coordinates_rotation, points);

        write_points_plane(file, points, points_per_row, coordinates_rotation, divisions);
    }

    float rotation_matrix_x[4][4] = {0};
    for (int i = 0; i < 2; i++){
        generate_rotation_matrix_x(rotation_matrix_x, (2 * i + 1) * M_PI / 2);

        float coordinates_rotation[points][3];

        apply_rotation(rotation_matrix_x, coordinates, coordinates_rotation, points);
        write_points_plane(file, points, points_per_row, coordinates_rotation, divisions);
    }

    file.close();
}

void generate_sphere(float radius, int slices, int stacks, char *filename)
{
    float x1, x2, x3, x4, y1, y2, z1, z2, z3, z4;
    float slices_angle = 2 * M_PI / slices;
    float stacks_angle = M_PI / stacks;

    std::ofstream file = open_file(filename);

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

            file << x1 << " " << y1 << " " << z1 << "\n";
            file << x3 << " " << y2 << " " << z3 << "\n";
            file << x4 << " " << y2 << " " << z4 << "\n";

            file << x1 << " " << y1 << " " << z1 << "\n";
            file << x4 << " " << y2 << " " << z4 << "\n";
            file << x2 << " " << y1 << " " << z2 << "\n";
            
        }
    }
    file.close();
}
void generate_cone(float radius, float height, int slices, int stacks, char *filename)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Error when trying to open the file!" << std::endl;
        exit(0);
    }

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

    std::ofstream file = open_file(filename);
    file << divisions * slices * 6 << std::endl;

    float rotation_matrix_z[4][4] = {0};
    for(int i = 0; i < divisions + 1; i++){
        generate_rotation_matrix_z(rotation_matrix_z, angle * i);

        coordinates[i][0] = radius  + (0 * rotation_matrix_z[0][0]) + circle_radius * rotation_matrix_z[0][1] + 0 * rotation_matrix_z[0][2];
        coordinates[i][1] = 0 * rotation_matrix_z[1][0] + circle_radius * rotation_matrix_z[1][1] + 0 * rotation_matrix_z[1][2];
        coordinates[i][2] = 0 * rotation_matrix_z[2][0] + circle_radius * rotation_matrix_z[2][1] + 0 * rotation_matrix_z[2][2];

    }
    
    angle = 2 * M_PI / slices;

    float previous_coordinates[divisions + 1][3];
    std::memcpy(previous_coordinates, coordinates, sizeof(coordinates));
    float rotation_matrix_y[4][4] = {0};
    for(int i = 0; i < slices; i++){
        generate_rotation_matrix_y(rotation_matrix_y, angle * (i + 1));

        float next_coordinates[divisions + 1][3];
        apply_rotation(rotation_matrix_y, coordinates, next_coordinates, divisions + 1);  

        for(int j = 0; j < divisions; j++){
            file << previous_coordinates[j][0] << " " << previous_coordinates[j][1] << " " << previous_coordinates[j][2] << "\n";
            file << next_coordinates[j + 1][0] << " " << next_coordinates[j + 1][1] << " " << next_coordinates[j + 1][2] << "\n";
            file << previous_coordinates[j + 1][0] << " " << previous_coordinates[j + 1][1] << " " << previous_coordinates[j + 1][2] << "\n";

            file << previous_coordinates[j][0] << " " << previous_coordinates[j][1] << " " << previous_coordinates[j][2] << "\n";
            file << next_coordinates[j][0] << " " << next_coordinates[j][1] << " " << next_coordinates[j][2] << "\n";
            file << next_coordinates[j + 1][0] << " " << next_coordinates[j + 1][1] << " " << next_coordinates[j + 1][2] << "\n";
        }
        std::memcpy(previous_coordinates, next_coordinates, sizeof(next_coordinates));

    }

    file.close();
}   

void generate_cylinder(float radius, float height, int divisions, char *filename){
    std::ofstream file = open_file(filename);
    
    float angle = 2 * M_PI / divisions;
    float coordinates[2][3] = {{radius, height / 2, 0}, {radius, -height / 2, 0}};
    file << (3 * divisions) * 2 + 6 * divisions << "\n";
    
    float rotation_matrix_y[4][4] = {0};
    float previous_coordinates[2][3];

    std::memcpy(previous_coordinates, coordinates, sizeof(coordinates));
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

        std::memcpy(previous_coordinates, new_coordinates, sizeof(new_coordinates));
    }  
    
    file.close();
}

int main(int argc, char **argv){

    if (std::strcmp(argv[1], "box") == 0 && argc == 5)
    {
        float length = std::stof(argv[2]);
        int divisions = std::stoi(argv[3]);
        
        if(length <= 0 || divisions <= 0){
            std::cout << "Length and divisions must be greater than 0!" << std::endl;
        }else{
            generate_box(length, divisions, argv[4]);
        }
    }
    else if (strcmp(argv[1], "sphere") == 0 && argc == 6)
    {
        float radius = std::stof(argv[2]);
        int slices = std::stoi(argv[3]);
        int stacks = std::stoi(argv[4]);

        if(radius <= 0 || slices <= 0 || stacks <= 0) {
            std::cout << "Radius, slices and stacks must be greater than 0!" << std::endl;
        }else{
            generate_sphere(radius, slices, stacks, argv[5]);
        }
    }
    else if (strcmp(argv[1], "cone") == 0 && argc == 7)
    {
        float radius = std::stof(argv[2]);
        float height = std::stof(argv[3]);
        int slices = std::stoi(argv[4]);
        int stacks = std::stoi(argv[5]);

        if (radius <= 0 || height <= 0 || slices <= 0 || stacks <= 0)
        {
            std::cout << "Radius, height, slices and stacks must be greater than 0!" << std::endl;
        }
        else
        {
            generate_cone(radius, height, slices, stacks, argv[6]);
        }
    }
    else if (strcmp(argv[1], "plane") == 0 && argc == 5)
    {
        float length = std::stof(argv[2]);
        int divisions = std::stoi(argv[3]);
        
        if(length <= 0 || divisions <= 0){
            std::cout << "Length and divisions must be greater than 0!" << std::endl;
        }else{
            generate_plane(length, divisions, argv[4]);
        }
    }else if(strcmp(argv[1], "torus") == 0 && argc == 7){
        float radius = std::stof(argv[2]);
        float circle_radius = std::stof(argv[3]);
        float slices = std::stoi(argv[4]);
        float divisions = std::stoi(argv[5]);
        
        if(radius < circle_radius || slices < 2 || divisions < 2){
            std::cout << "Radius must be greater than circle_radius and slices and divisions must be greater than 2" << std::endl;
        }else{
            generate_torus(radius, circle_radius, slices, divisions, argv[6]);
        }
    }else if(strcmp(argv[1], "cylinder") == 0  && argc == 6){
        float radius = std::stof(argv[2]);
        float height = std::stof(argv[3]);
        float divisions = std::stof(argv[4]);

        if(radius <= 0 || height <= 0 || divisions < 2){
            std::cout << "Radius and height must be greater than 0 and divisions must be greater than 2" << std::endl;
        }else{
            generate_cylinder(radius, height, divisions, argv[5]);
        }
    }else{
        std::cout << "Something went wrong" << std::endl;
    }

    return 0;
}
