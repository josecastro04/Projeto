#include <cstring>
#include <iostream>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

void write_points(std::ofstream& file, int points, int points_per_row, float coordinates[][3], int divisions){
    for(int i = 0; i < points - points_per_row ; i++){
        if(i % points_per_row == divisions) continue;

        file << coordinates[i][0] << " " << coordinates[i][1] << " " << coordinates[i][2] << "\n";
        file << coordinates[i + points_per_row][0] << " " << coordinates[i + points_per_row][1] << " " << coordinates[i + points_per_row][2] << "\n";
        file << coordinates[i + points_per_row + 1][0] << " " << coordinates[i + points_per_row + 1][1] << " " << coordinates[i + points_per_row + 1][2] << "\n";

        file << coordinates[i][0] << " " << coordinates[i][1] << " " << coordinates[i][2] << "\n";
        file << coordinates[i + points_per_row + 1][0] << " " << coordinates[i + points_per_row + 1][1] << " " << coordinates[i + points_per_row + 1][2] << "\n";
        file << coordinates[i + 1][0] << " " << coordinates[i + 1][1] << " " << coordinates[i + 1][2] << "\n"; 
    }
}

void generate_points_plane(float coordinates[][3], float start_x, float start_y, int points, int divisions, float distance){
    for(int row = 0; row <= divisions; row++){
        float space_z = start_x + (distance * row);
        int index = row * (divisions + 1);

        for(int column = 0; column <= divisions; column++){
            float space_x = start_x + (distance * column);
            coordinates[index + column][0] = space_x;
            coordinates[index + column][1] = start_y;
            coordinates[index + column][2] = space_z;
        }
    }
}

void generate_rotation_matrix_x(float rotation_matrix[4][4], float angle){
    rotation_matrix[0][0] = rotation_matrix[3][3] = 1;
    rotation_matrix[1][1] = rotation_matrix[2][2] = cos(angle); 
    rotation_matrix[2][1] = sin(angle);
    rotation_matrix[1][2] = -1 * sin(angle);
}

void generate_rotation_matrix_z(float rotation_matrix[4][4], float angle){
    rotation_matrix[2][2] = rotation_matrix[3][3] = 1;
    rotation_matrix[0][0] = rotation_matrix[1][1] = cos(angle);
    rotation_matrix[0][1] = -1 * sin(angle);
    rotation_matrix[1][0] = sin(angle);
}

void generate_rotation_matrix_y(float rotation_matrix[4][4], float angle){
    rotation_matrix[1][1] = rotation_matrix[3][3] = 1;
    rotation_matrix[0][0] = rotation_matrix[2][2] = cos(angle);
    rotation_matrix[2][1] = -1 * sin(angle);
    rotation_matrix[0][2] = sin(angle);
}

void generate_plane(float length, int divisions, char* filename){
    int total_points = (divisions * divisions) * 6;
    int points_per_row = divisions + 1;
    std::ofstream file(filename);
    if(!file){
        std::cerr << "Error when trying to open the file" << std::endl;
        exit(0);
    }

    float start = -length / 2;
    float distance = length / divisions;

    file << total_points << "\n";
    int points = points_per_row * points_per_row; 
    float coordinates[points][3];


    generate_points_plane(coordinates, start, 0, points, divisions, distance);

    write_points(file, points, points_per_row, coordinates, divisions);

    file.close();

}

void apply_rotation(float rotation_matrix[4][4], float coordinates[][3], float coordinates_rotation[][3], int points){
    for(int i = 0; i < points; i++){
        float x = coordinates[i][0];
        float y = coordinates[i][1];
        float z = coordinates[i][2];


        coordinates_rotation[i][0] = x * rotation_matrix[0][0] + y * rotation_matrix[0][1] + z * rotation_matrix[0][2];
        coordinates_rotation[i][1] = x * rotation_matrix[1][0] + y * rotation_matrix[1][1] + z * rotation_matrix[1][2]; 
        coordinates_rotation[i][2] = x * rotation_matrix[2][0] + y * rotation_matrix[2][1] + z * rotation_matrix[2][2];
    }
}


void generate_box(float length, int divisions, char* filename){
    int total_points = (divisions * divisions) * 36;
    int points_per_row = divisions + 1;

    std::ofstream file(filename);
    if(!file){
        std::cerr << "Error when trying to open the file!" << std::endl;
        exit(0);
    }

    file << total_points << "\n";

    float start = -length / 2;

    float distance = length / divisions;

    int points = points_per_row * points_per_row;

    float coordinates[points][3];

    generate_points_plane(coordinates, start, -start, points, divisions, distance);
    write_points(file, points, points_per_row, coordinates, divisions); 

    for(int i = 1; i < 4; i++){
        float rotation_matrix[4][4] = {0};
        generate_rotation_matrix_z(rotation_matrix, i * M_PI/2);

        float coordinates_rotation[points][3];

        apply_rotation(rotation_matrix, coordinates, coordinates_rotation, points);
        
        write_points(file, points, points_per_row, coordinates_rotation, divisions);
    }

    for(int i = 0; i < 2; i++){
        float rotation_matrix[4][4] = {0};
        generate_rotation_matrix_x(rotation_matrix, (2 * i + 1) * M_PI/2);

        float coordinates_rotation[points][3];

        apply_rotation(rotation_matrix, coordinates, coordinates_rotation, points);

        write_points(file, points, points_per_row, coordinates_rotation, divisions);
    }

    file.close();

}

int main(int argc, char **argv){

    if(strcmp(argv[1], "box") == 0 && argc == 5){
        float length = std::stof(argv[2]);
        int divisions = std::stoi(argv[3]);
        generate_box(length, divisions, argv[4]);
    }else if(strcmp(argv[1], "sphere") == 0 && argc == 6){
    }else if(strcmp(argv[1], "cone") == 0 && argc == 7){
    }else if(strcmp(argv[1], "plane") == 0 && argc == 5){
        float length = std::stof(argv[2]);
        int divisions = std::stoi(argv[3]);
        generate_plane(length, divisions, argv[4]);
    }

    return 0;
}
