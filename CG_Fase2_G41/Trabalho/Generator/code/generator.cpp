#include <iostream>
#include <fstream>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <vector>

// Tipo da figura passada como argumento
char* figure;


// Função que calcula os pontos do plane de acordo com os parâmetros introduzidos
void plane(float dimension, float divisions, string file) {
	// Abre e escreve no ficheiro
	ofstream write(file);

	float vertices = 2 * 3 * divisions * divisions;
	// Escreve o número de vértices da figura
	write << vertices << endl;

	// Variáveis
	float lado = dimension / divisions;
	float d = dimension / static_cast<float>(2);

	// Triângulos de cima
	for (float z = -d; z < d; z += lado) {
		for (float x = -d; x < d; x += lado) {
			// Primeiro ponto do triângulo
			write << x + lado << " " << 0.0f << " " << z << endl;

			// Segundo ponto do triângulo
			write << x << " " << 0.0f << " " << z << endl;

			// Terceiro ponto do triângulo
			write << x << " " << 0.0f << " " << z + lado << endl;
		}
	}
	// Triângulos de baixo
	for (float z = d; z > -d; z -= lado) {
		for (float x = d; x > -d; x -= lado) {
			// Primeiro ponto do triângulo
			write << x - lado << " " << 0.0f << " " << z << endl;

			// Segundo ponto do triângulo
			write << x << " " << 0.0f << " " << z << endl;

			// Terceiro ponto do triângulo
			write << x << " " << 0.0f << " " << z - lado << endl;
		}
	}
	write.close();
}


// Função que calcula os pontos da caixa de acordo com os parâmetros introduzidos
void box(float dimension, float divisions, string file) {
	// Abre e escreve no ficheiro
	ofstream write(file);

	// Escreve o número de vértices da figura
	float vertices = 2 * 3 * divisions * divisions * 6;
	write << vertices << endl;

	// Variáveis
	float lado = dimension / divisions;
	float d = dimension / static_cast<float>(2);

	// Face de cima
	for (float z = -d; z < d; z += lado) {
		for (float x = -d; x < d; x += lado) {
			write << x + lado << " " << d << " " << z << endl;
			write << x << " " << d << " " << z << endl;
			write << x << " " << d << " " << z + lado << endl;
		}
	}
	for (float z = d; z > -d; z -= lado) {
		for (float x = d; x > -d; x -= lado) {
			write << x - lado << " " << d << " " << z << endl;
			write << x << " " << d << " " << z << endl;
			write << x << " " << d << " " << z - lado << endl;
		}
	}


	// Face de baixo
	for (float z = -d; z < d; z += lado) {
		for (float x = -d; x < d; x += lado) {
			write << x << " " << -d << " " << z << endl;
			write << x + lado << " " << -d << " " << z << endl;
			write << x << " " << -d << " " << z + lado << endl;
		}
	}
	for (float z = d; z > -d; z -= lado) {
		for (float x = d; x > -d; x -= lado) {
			write << x << " " << -d << " " << z << endl;
			write << x - lado << " " << -d << " " << z << endl;
			write << x << " " << -d << " " << z - lado << endl;
		}
	}


	// Face de trás
	for (float y = -d; y < d; y += lado) {
		for (float x = -d; x < d; x += lado) {
			write << x + lado << " " << y << " " << -d << endl;
			write << x << " " << y << " " << -d << endl;
			write << x << " " << y + lado << " " << -d << endl;
		}
	}
	for (float y = d; y > -d; y -= lado) {
		for (float x = d; x > -d; x -= lado) {
			write << x - lado << " " << y << " " << -d << endl;
			write << x << " " << y << " " << -d << endl;
			write << x << " " << y - lado << " " << -d << endl;
		}
	}


	// Face da frente
	for (float y = -d; y < d; y += lado) {
		for (float x = -d; x < d; x += lado) {
			write << x << " " << y << " " << d << endl;
			write << x + lado << " " << y << " " << d << endl;
			write << x << " " << y + lado << " " << d << endl;
		}
	}
	for (float y = d; y > -d; y -= lado) {
		for (float x = d; x > -d; x -= lado) {
			write << x << " " << y << " " << d << endl;
			write << x - lado << " " << y << " " << d << endl;
			write << x << " " << y - lado << " " << d << endl;
		}
	}


	// Face da esquerda
	for (float y = -d; y < d; y += lado) {
		for (float z = -d; z < d; z += lado) {
			write << -d << " " << y << " " << z << endl;
			write << -d << " " << y << " " << z + lado << endl;
			write << -d << " " << y + lado << " " << z << endl;
		}
	}
	for (float y = d; y > -d; y -= lado) {
		for (float z = d; z > -d; z -= lado) {
			write << -d << " " << y << " " << z << endl;
			write << -d << " " << y << " " << z - lado << endl;
			write << -d << " " << y - lado << " " << z << endl;
		}
	}


	// Face da direita
	for (float y = -d; y < d; y += lado) {
		for (float z = -d; z < d; z += lado) {
			write << d << " " << y << " " << z + lado << endl;
			write << d << " " << y << " " << z << endl;
			write << d << " " << y + lado << " " << z << endl;
		}
	}
	for (float y = d; y > -d; y -= lado) {
		for (float z = d; z > -d; z -= lado) {
			write << d << " " << y << " " << z - lado << endl;
			write << d << " " << y << " " << z << endl;
			write << d << " " << y - lado << " " << z << endl;
		}
	}

	write.close();
}


// Função que calcula os pontos da esfera de acordo com os parâmetros introduzidos
void sphere(float radius, float slices, float stacks, string file) {
	// Abre e escreve no ficheiro
	ofstream write(file);

	// Escreve o número de vértices da figura
	float vertices =  2 * 6 * (stacks / 2) * slices;
	write << vertices << endl;

	// Variáveis destinadas às coordenas e aos ângulos
	float next_theta, next_alpha, theta = 0, alpha = 0, x, y, z, next_y;

	for (int i = 0; i < stacks / 2; i++) {
		theta = 0;
		next_alpha = alpha + (M_PI / stacks);

		for (int j = 0; j < slices; j++) {
			next_theta = theta + ((2 * M_PI) / slices);

			// Metade de cima
			x = radius * cos(alpha) * sin(theta);
			y = radius * sin(alpha);
			z = radius * cos(theta) * cos(alpha);
			write << x << " " << y << " " << z << endl;

			x = radius * cos(alpha) * sin(next_theta);
			y = radius * sin(alpha);
			z = radius * cos(next_theta) * cos(alpha);
			write << x << " " << y << " " << z << endl;

			x = radius * cos(next_alpha) * sin(theta);
			y = radius * sin(next_alpha);
			z = radius * cos(theta) * cos(next_alpha);
			write << x << " " << y << " " << z << endl;

			x = radius * cos(next_alpha) * sin(next_theta);
			y = radius * sin(next_alpha);
			z = radius * cos(next_theta) * cos(next_alpha);
			write << x << " " << y << " " << z << endl;

			x = radius * cos(next_alpha) * sin(theta);
			y = radius * sin(next_alpha);
			z = radius * cos(theta) * cos(next_alpha);
			write << x << " " << y << " " << z << endl;

			x = radius * cos(alpha) * sin(next_theta);
			y = radius * sin(alpha);
			z = radius * cos(next_theta) * cos(alpha);
			write << x << " " << y << " " << z << endl;


			// Metade de baixo
			x = radius * cos(alpha) * sin(theta);
			y = radius * sin(alpha);
			z = radius * cos(theta) * cos(alpha);
			write << x << " " << -y << " " << z << endl;

			x = radius * cos(next_alpha) * sin(theta);
			y = radius * sin(next_alpha);
			z = radius * cos(theta) * cos(next_alpha);
			write << x << " " << -y << " " << z << endl;

			x = radius * cos(next_alpha) * sin(next_theta);
			y = radius * sin(next_alpha);
			z = radius * cos(next_theta) * cos(next_alpha);
			write << x << " " << -y << " " << z << endl;

			x = radius * cos(next_alpha) * sin(next_theta);
			y = radius * sin(next_alpha);
			z = radius * cos(next_theta) * cos(next_alpha);
			write << x << " " << -y << " " << z << endl;

			x = radius * cos(alpha) * sin(next_theta);
			y = radius * sin(alpha);
			z = radius * cos(next_theta) * cos(alpha);
			write << x << " " << -y << " " << z << endl;

			x = radius * cos(alpha) * sin(theta);
			y = radius * sin(alpha);
			z = radius * cos(theta) * cos(alpha);
			write << x << " " << -y << " " << z << endl;

			theta = next_theta;
		}
		alpha = next_alpha;
	}
}


// Função que calcula os pontos do cone de acordo com os parâmetros introduzidos
void cone(float radius, float height, float slices, float stacks, string file) {
	// Abre e escreve no ficheiro
	ofstream write(file);

	// Escreve o número de vértices da figura
	float vertices = 3 * slices + 6 * stacks * slices;
	write << vertices << endl;

	// Variáveis destinadas às coordenas e aos ângulos
	float step = 360.0 / slices;
	float x, y = 0, z, next_y, next_r, r = radius;

	// Base
	for (int i = 0; i < slices; i++) {
		write << cos((i + 1) * step * M_PI / 180.0) * radius << " " << 0.0f << " " << -sin((i + 1) * step * M_PI / 180.0) * radius << endl;
		write << cos(i * step * M_PI / 180.0) * radius << " " << 0.0f << " " << -sin(i * step * M_PI / 180.0) * radius << endl;
		write << 0.0f << " " << 0.0f << " " << 0.0f << endl;
	}

	// Lateral
	for (int i = 0; i < stacks; i++) {
		next_y = y + height / stacks;
		next_r = r - (radius / stacks);
		for (int j = 0; j < slices; j++) {
			write << r * sin(j * ((2 * M_PI) / slices)) << " " << y << " " << r * cos(j * ((2 * M_PI) / slices)) << endl;
			write << r * sin((j + 1) * ((2 * M_PI) / slices)) << " " << y << " " << r * cos((j + 1) * ((2 * M_PI) / slices)) << endl;
			write << next_r * sin(j * ((2 * M_PI) / slices)) << " " << next_y << " " << next_r * cos(j * ((2 * M_PI) / slices)) << endl;

			write << next_r * sin(j * ((2 * M_PI) / slices)) << " " << next_y << " " << next_r * cos(j * ((2 * M_PI) / slices)) << endl;
			write << r * sin((j + 1) * ((2 * M_PI) / slices)) << " " << y << " " << r * cos((j + 1) * ((2 * M_PI) / slices)) << endl;
			write << next_r * sin((j + 1) * ((2 * M_PI) / slices)) << " " << next_y << " " << next_r * cos((j + 1) * ((2 * M_PI) / slices)) << endl;
		}
		r = next_r;
		y = next_y;
	}
	write.close();
}


// Função que calcula os pontos do torus de acordo com os parâmetros introduzidos
void torus(float inner_radius, float outer_radius, int slices, int stacks, string file) {
	// Abre e escreve no ficheiro
	ofstream write(file);

	// Escreve o número de vértices da figura
	float vertices = 6 * stacks * slices;
	write << vertices << endl;

	// Variação do ângulo entre slices
	float delta_slices = 2 * M_PI / slices;

	// Variação do ângulo entre stacks
	float delta_stacks = 2 * M_PI / stacks;

	float phi = 0;
	float theta = 0;

	for (int i = 0; i < slices; i++) {
		for (int j = 0; j < stacks; j++) {
			write << (outer_radius + inner_radius * cos(phi)) * cos(theta) << " "
				  << (outer_radius + inner_radius * cos(phi)) * sin(theta) << " "
				  << inner_radius * sin(phi) << endl;
			write << (outer_radius + inner_radius * cos(phi)) * cos(theta + delta_slices) << " "
				  << (outer_radius + inner_radius * cos(phi)) * sin(theta + delta_slices) << " "
				  << inner_radius * sin(phi) << endl;
			write << (outer_radius + inner_radius * cos(phi + delta_stacks)) * cos(theta + delta_slices) << " "
				  << (outer_radius + inner_radius * cos(phi + delta_stacks)) * sin(theta + delta_slices) << " "
				  << inner_radius * sin(phi + delta_stacks) << endl;

			write << (outer_radius + inner_radius * cos(phi + delta_stacks)) * cos(theta + delta_slices) << " "
				  << (outer_radius + inner_radius * cos(phi + delta_stacks)) * sin(theta + delta_slices) << " "
				  << inner_radius * sin(phi + delta_stacks) << endl;
			write << (outer_radius + inner_radius * cos(phi + delta_stacks)) * cos(theta) << " "
				  << (outer_radius + inner_radius * cos(phi + delta_stacks)) * sin(theta) << " "
				  << inner_radius * sin(phi + delta_stacks) << endl;
			write << (outer_radius + inner_radius * cos(phi)) * cos(theta) << " "
				  << (outer_radius + inner_radius * cos(phi)) * sin(theta) << " "
				  << inner_radius * sin(phi) << endl;

			phi = delta_stacks * (j + 1);
		}
		theta = delta_slices * (i + 1);
	}
	write.close();
}


// Função principal do programa
int main(int argc, char** argv) {
	// Armazena a figura e os seus parâmetros em variaveis globais
	if (argc > 1) {
		figure = argv[1];

		if (!strcmp(figure, "plane") && argc == 5) {
			float dimension = atof(argv[2]);
			float divisions = atof(argv[3]);
			string file = "../models/" + string(argv[4]);
			plane(dimension, divisions, file);
		}
		else if (!strcmp(figure, "box") && argc == 5) {
			float dimension = atof(argv[2]);
			float divisions = atof(argv[3]);
			string file = "../models/" + string(argv[4]);
			box(dimension, divisions, file);
		}
		else if (!strcmp(figure, "sphere") && argc == 6) {
			float radius = atof(argv[2]);
			float stacks = atof(argv[3]);
			float slices = atof(argv[4]);
			string file = "../models/" + string(argv[5]);
			sphere(radius, slices, stacks, file);
		}
		else if (!strcmp(figure, "cone") && argc == 7) {
			float radius = atof(argv[2]);
			float height = atof(argv[3]);
			float slices = atof(argv[4]);
			float stacks = atof(argv[5]);
			string file = "../models/" + string(argv[6]);
			cone(radius, height, slices, stacks, file);
		}
		else if (!strcmp(figure, "rings") && argc == 7) {
			float inner_radius = atof(argv[2]);
			float outer_radius = atof(argv[3]);
			int slices = atoi(argv[4]);
			int stacks = atoi(argv[5]);
			string file = "../models/" + string(argv[6]);
			torus(inner_radius, outer_radius, slices, stacks, file);
		}
		else {
			std::cout << "Erro: Argumentos invalidos";
			return 1;
		}
	}
	else {
		std::cout << "Erro: Figura nao suportada";
		return 1;
	}

	return 0;
}
