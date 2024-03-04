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

// 
void cross(float* a, float* b, float* res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

// Função que multiplica uma matriz por um vetor
void multMatrixVector(float* m, float* v, float* res) {
	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}
}

void getBezierPoint(float u, float v, float(*points)[3], float* pos, float* du, float* dv) {
	float m[4][4] = { {-1.0f, 3.0f, -3.0f, 1.0f},
					  {3.0f, -6.0f, 3.0f, 0.0f},
					  {-3.0f, 3.0f, 0.0f, 0.0f},
					  {1.0f,  0.0f, 0.0f, 0.0f} };

	// Vetores U e V
	float U[4] = { u * u * u, u * u, u, 1 };
	float V[4] = { v * v * v, v * v, v, 1 };

	// Vetores U' e V'
	float U_linha[4] = { 3 * u * u, 2 * u, 1, 0 };
	float V_linha[4] = { 3 * v * v, 2 * v, 1, 0 };

	// Armazena cada componente dos 16 pontos
	float p[3][16];
	for (int i = 0; i < 16; i++) {
		p[0][i] = points[i][0];
		p[1][i] = points[i][1];
		p[2][i] = points[i][2];
	}

	// Cálculo para cada componente x, y e z do ponto
	for (int i = 0; i < 3; i++) { // x, y, z
		float A1[4];
		float A2[4];
		float A1_linha[4];
		float A2_linha[4];

		// Compute A1 = U * M
		multMatrixVector((float*)m, U, A1);

		// Compute A2 = A1 * P 
		multMatrixVector(p[i], A1, A2);


		// Compute A1 = A2 * MT (where MT = M)
		multMatrixVector((float*)m, A2, A1);

		// Mesma coisa mas com as derivadas para calcular as normais
		multMatrixVector((float*)m, U_linha, A1_linha);
		multMatrixVector(p[i], A1_linha, A2_linha);
		multMatrixVector((float*)m, A2_linha, A1_linha);

		// Compute pos[i] = A1 * V  (i.e. U * M * P * MT * V)
		pos[i] = V[0] * A1[0] + V[1] * A1[1] + V[2] * A1[2] + V[3] * A1[3];
		du[i] = V[0] * A1_linha[0] + V[1] * A1_linha[1] + V[2] * A1_linha[2] + V[3] * A1_linha[3];
		dv[i] = V_linha[0] * A1[0] + V_linha[1] * A1[1] + V_linha[2] * A1[2] + V_linha[3] * A1[3];
	}
}

void patch(string input_file, int tesselation, string file) {
	int n_patches, n_cpoints;

	// Abre o ficheiro
	ifstream patch_file(input_file);
	if (!patch_file) {
		cerr << "Erro ao abrir o ficheiro: " << input_file << endl;
		return;
	}

	// Armazena o número de patches
	patch_file >> n_patches;

	// Armazena os 16 indices dos pontos de controle de cada Patch -> cada 4 pontos formam uma curva de Bezier
	vector<vector<int>> patches(n_patches);
	for (int i = 0; i < n_patches; i++) {
		for (int j = 0; j < 16; j++) {
			string line;
			patch_file >> line;
			patches[i].push_back(stoi(line));
		}
	}

	// Armazena o número de pontos de controlo
	patch_file >> n_cpoints;

	// Armazena as coordenadas dos pontos de controlo
	vector<vector<float>> cpoints(n_cpoints);
	for (int i = 0; i < n_cpoints; i++) {
		for (int j = 0; j < 3; j++) {
			string line;
			patch_file >> line;
			cpoints[i].push_back(stof(line));
		}
	}

	// Fecha o ficheiro de input
	patch_file.close();

	// Abre e escreve no ficheiro de output
	ofstream write(file);

	// Escreve o número de vértices da figura
	int vertices = n_patches * 6 * tesselation * tesselation;
	write << vertices << endl;

	//comet(n_patches,patches,n_cpoints,cpoints);

	// Definir a variação de u e de v
	float delta_u, delta_v;
	delta_u = delta_v = (float)1 / (float)tesselation;

	// Guardar as coordenadas de cada ponto calculado
	float pos[3], du[3], dv[3], normal[3];
	for (int i = 0; i < n_patches; i++) {
		// Armazena os pontos correspondentes ao índices de cada patch
		float points[16][3];
		for (int k = 0; k < 16; k++) {
			int index = patches[i][k];
			for (int l = 0; l < 3; l++) {
				points[k][l] = cpoints[index][l];
			}
		}
		// Calcula os pontos dos quadrados
		for (float u = 0; u < 1; u += delta_u) {
			for (float v = 0; v < 1; v += delta_v) {
				// Triangulo 1
				getBezierPoint(u, v, points, pos, du, dv);
				cross(du, dv, normal);
				write << pos[0] << " " << pos[1] << " " << pos[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << " " << u << " "<< v << endl;

				getBezierPoint(u + delta_u, v, points, pos, du, dv);
				cross(du, dv, normal);
				write << pos[0] << " " << pos[1] << " " << pos[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << " " << u + delta_u << " " << v << endl;

				getBezierPoint(u, v + delta_v, points, pos, du, dv);
				cross(du, dv, normal);
				write << pos[0] << " " << pos[1] << " " << pos[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << " " << u << " " << v + delta_v << endl;

				// Triangulo 2
				getBezierPoint(u + delta_u, v, points, pos, du, dv);
				cross(du, dv, normal);
				write << pos[0] << " " << pos[1] << " " << pos[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << " " << u + delta_u << " " << v << endl;

				getBezierPoint(u + delta_u, v + delta_v, points, pos, du, dv);
				cross(du, dv, normal);
				write << pos[0] << " " << pos[1] << " " << pos[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << " " << u + delta_u << " " << v + delta_v << endl;

				getBezierPoint(u, v + delta_v, points, pos, du, dv);
				cross(du, dv, normal);
				write << pos[0] << " " << pos[1] << " " << pos[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << " " << u << " " << v + delta_v << endl;
			}
		}
	}
	write.close();
}


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
	for (float z = 0; z < dimension; z += lado) {
		for (float x = 0; x < dimension; x += lado) {
			// Primeiro ponto do triângulo
			write << x + lado - d << " " << 0.0f << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << (x+lado)/dimension << " " << z/dimension << endl;

			// Segundo ponto do triângulo
			write << x - d << " " << 0.0f << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << z / dimension << endl;

			// Terceiro ponto do triângulo
			write << x - d << " " << 0.0f << " " << z + lado - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << (z + lado) / dimension << endl;
		}
	}
	// Triângulos de baixo
	for (float z = dimension; z > 0; z -= lado) {
		for (float x = dimension; x > 0; x -= lado) {
			// Primeiro ponto do triângulo
			write << x - lado - d << " " << 0.0f << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << (x - lado) / dimension << " " << z / dimension << endl;

			// Segundo ponto do triângulo
			write << x - d << " " << 0.0f << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << z / dimension << endl;

			// Terceiro ponto do triângulo
			write << x - d << " " << 0.0f << " " << z - lado - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << (z - lado) / dimension << endl;
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
	for (float z = 0; z < dimension; z += lado) {
		for (float x = 0; x < dimension; x += lado) {
			write << x + lado - d << " " << d << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << (x+lado)/dimension << " " << z/dimension << endl;
			write << x - d << " " << d << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << z / dimension << endl;
			write << x - d << " " << d << " " << z + lado - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << (z + lado) / dimension << endl;
		}
	}
	for (float z = dimension; z > 0; z -= lado) {
		for (float x = dimension; x > 0; x -= lado) {
			write << x - lado - d << " " << d << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << (x - lado) / dimension << " " << z / dimension << endl;
			write << x - d << " " << d << " " << z - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << z / dimension << endl;
			write << x - d << " " << d << " " << z - lado - d << " " << 0.0f << " " << 1.0f << " " << 0.0f << " " << x / dimension << " " << (z - lado) / dimension << endl;
		}
	}

	// Face de baixo
	for (float z = 0; z < dimension; z += lado) {
		for (float x = 0; x < dimension; x += lado) {
			write << x - d << " " << -d << " " << z - d << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << x / dimension << " " << z / dimension << endl;
			write << x + lado -d << " " << -d << " " << z - d << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << (x + lado) / dimension << " " << z / dimension << endl;
			write << x - d << " " << -d << " " << z + lado - d << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << x / dimension << " " << (z + lado) / dimension << endl;
		}
	}
	for (float z = dimension; z > 0; z -= lado) {
		for (float x = dimension; x > 0; x -= lado) {
			write << x - d << " " << -d << " " << z - d << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << x / dimension << " " << z / dimension << endl;
			write << x - lado - d << " " << -d << " " << z - d << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << (x - lado) / dimension << " " << z / dimension << endl;
			write << x - d << " " << -d << " " << z - lado - d << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << x / dimension << " " << (z - lado) / dimension << endl;
		}
	}


	// Face de trás
	for (float y = 0; y < dimension; y += lado) {
		for (float x = 0; x < dimension; x += lado) {
			write << x + lado - d << " " << y - d << " " << -d << " " << 0.0f << " " << 0.0f << " " << -1.0f << " " << (x + lado) / dimension << " " << y / dimension << endl;
			write << x - d << " " << y - d << " " << -d << " " << 0.0f << " " << 0.0f << " " << -1.0f << " " << x / dimension << " " << y / dimension << endl;
			write << x - d << " " << y + lado - d << " " << -d << " " << 0.0f << " " << 0.0f << " " << -1.0f << " " << x / dimension << " " << (y + lado) / dimension << endl;
		}
	}
	for (float y = dimension; y > 0; y -= lado) {
		for (float x = dimension; x > 0; x -= lado) {
			write << x - lado - d << " " << y - d << " " << -d << " " << 0.0f << " " << 0.0f << " " << -1.0f << " " << (x - lado) / dimension << " " << y / dimension << endl;
			write << x - d << " " << y - d << " " << -d << " " << 0.0f << " " << 0.0f << " " << -1.0f << " " << x / dimension << " " << y / dimension << endl;
			write << x - d << " " << y - lado - d << " " << -d << " " << 0.0f << " " << 0.0f << " " << -1.0f << " " << x / dimension << " " << (y - lado) / dimension << endl;
		}
	}


	// Face da frente
	for (float y = 0; y < dimension; y += lado) {
		for (float x = 0; x < dimension; x += lado) {
			write << x - d << " " << y - d << " " << d << " " << 0.0f << " " << 0.0f << " " << 1.0f << " " << x / dimension << " " << y / dimension << endl;
			write << x + lado - d << " " << y - d << " " << d << " " << 0.0f << " " << 0.0f << " " << 1.0f << " " << (x + lado) / dimension << " " << y / dimension << endl;
			write << x - d << " " << y + lado - d << " " << d << " " << 0.0f << " " << 0.0f << " " << 1.0f << " " << x / dimension << " " << (y + lado) / dimension << endl;
		}
	}
	for (float y = dimension; y > 0; y -= lado) {
		for (float x = dimension; x > 0; x -= lado) {
			write << x - d << " " << y - d << " " << d << " " << 0.0f << " " << 0.0f << " " << 1.0f << " " << x / dimension << " " << y / dimension << endl;
			write << x - lado - d << " " << y - d << " " << d << " " << 0.0f << " " << 0.0f << " " << 1.0f << " " << (x - lado) / dimension << " " << y / dimension << endl;
			write << x - d << " " << y - lado - d << " " << d << " " << 0.0f << " " << 0.0f << " " << 1.0f << " " << x / dimension << " " << (y - lado) / dimension << endl;
		}
	}


	// Face da esquerda
	for (float y = 0; y < dimension; y += lado) {
		for (float z = 0; z < dimension; z += lado) {
			write << -d << " " << y - d << " " << z - d << " " << -1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << y / dimension << endl;
			write << -d << " " << y - d << " " << z + lado - d << " " << -1.0f << " " << 0.0f << " " << 0.0f << " " << (z + lado) / dimension << " " << y / dimension << endl;
			write << -d << " " << y + lado - d << " " << z - d << " " << -1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << (y + lado) / dimension << endl;
		}
	}
	for (float y = dimension; y > 0; y -= lado) {
		for (float z = dimension; z > 0; z -= lado) {
			write << -d << " " << y - d << " " << z - d << " " << -1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << y / dimension << endl;
			write << -d << " " << y - d << " " << z - lado - d << " " << -1.0f << " " << 0.0f << " " << 0.0f << " " << (z - lado) / dimension << " " << y / dimension << endl;
			write << -d << " " << y - lado - d << " " << z - d << " " << -1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << (y - lado) / dimension << endl;
		}
	}


	// Face da direita
	for (float y = 0; y < dimension; y += lado) {
		for (float z = 0; z < dimension; z += lado) {
			write << d << " " << y - d << " " << z + lado - d << " " << 1.0f << " " << 0.0f << " " << 0.0f << " " << (z+lado) / dimension << " " << y / dimension << endl;
			write << d << " " << y - d << " " << z - d << " " << 1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << y / dimension << endl;
			write << d << " " << y + lado - d << " " << z - d << " " << 1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << (y + lado) / dimension << endl;
		}
	}
	for (float y = dimension; y > 0; y -= lado) {
		for (float z = dimension; z > 0; z -= lado) {
			write << d << " " << y - d << " " << z - lado - d << " " << 1.0f << " " << 0.0f << " " << 0.0f << " " << (z - lado) / dimension << " " << y / dimension << endl;
			write << d << " " << y - d << " " << z - d << " " << 1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << y / dimension << endl;
			write << d << " " << y - lado - d << " " << z - d << " " << 1.0f << " " << 0.0f << " " << 0.0f << " " << z / dimension << " " << (y-lado) / dimension << endl;
		}
	}
	write.close();
}



// Função que calcula os pontos da esfera de acordo com os parâmetros introduzidos
void sphere(float radius, float slices, float stacks, string file) {
	// Abre e escreve no ficheiro
	ofstream write(file);

	// Escreve o número de vértices da figura
	float vertices = 6 * stacks * slices;
	write << vertices << endl;

	// Variáveis destinadas às coordenas e aos ângulos
	float next_theta, next_alpha, theta = 0, alpha = -M_PI/(float)2, x, y, z, next_y, nx, ny, nz, t1, t2;

	for (int i = 0; i < stacks; i++) {
		theta = 0;
		next_alpha = alpha + (M_PI / stacks);

		for (int j = 0; j < slices; j++) {
			next_theta = theta + ((2 * M_PI) / slices);

			x = radius * cos(alpha) * sin(theta);
			y = radius * sin(alpha);
			z = radius * cos(theta) * cos(alpha);
			nx = cos(alpha) * sin(theta);
			ny = sin(alpha);
			nz = cos(theta) * cos(alpha);
			t1 = (float)j/slices;
			t2 = (float)i/stacks;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t1 << " " << t2 << endl;

			x = radius * cos(alpha) * sin(next_theta);
			y = radius * sin(alpha);
			z = radius * cos(next_theta) * cos(alpha);
			nx = cos(alpha) * sin(next_theta);
			ny = sin(alpha);
			nz = cos(next_theta) * cos(alpha);
			t1 = (float)(j+1) / slices;
			t2 = (float)i / stacks;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t1 << " " << t2 << endl;

			x = radius * cos(next_alpha) * sin(theta);
			y = radius * sin(next_alpha);
			z = radius * cos(theta) * cos(next_alpha);
			nx = cos(next_alpha) * sin(theta);
			ny = sin(next_alpha);
			nz = cos(theta) * cos(next_alpha);
			t1 = (float)j / slices;
			t2 = (float)(i+1) / stacks;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t1 << " " << t2 << endl;

			x = radius * cos(next_alpha) * sin(next_theta);
			y = radius * sin(next_alpha);
			z = radius * cos(next_theta) * cos(next_alpha);
			nx = cos(next_alpha) * sin(next_theta);
			ny = sin(next_alpha);
			nz = cos(next_theta) * cos(next_alpha);
			t1 = (float)(j+1) / slices;
			t2 = (float)(i+1) / stacks;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t1 << " " << t2 << endl;

			x = radius * cos(next_alpha) * sin(theta);
			y = radius * sin(next_alpha);
			z = radius * cos(theta) * cos(next_alpha);
			nx = cos(next_alpha) * sin(theta);
			ny = sin(next_alpha);
			nz = cos(theta) * cos(next_alpha);
			t1 = (float)j / slices;
			t2 = (float)(i+1) / stacks;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t1 << " " << t2 << endl;

			x = radius * cos(alpha) * sin(next_theta);
			y = radius * sin(alpha);
			z = radius * cos(next_theta) * cos(alpha);
			nx = cos(alpha) * sin(next_theta);
			ny = sin(alpha);
			nz = cos(next_theta) * cos(alpha);
			t1 = (float)(j+1) / slices;
			t2 = (float)i / stacks;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t1 << " " << t2 << endl;

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
	float x, y = 0, z, t1,t2, next_y, next_r, r = radius;

	// Base
	for (int i = 0; i < slices; i++) {
		x = cos((i + 1) * step * M_PI / 180.0) * radius;
		z = -sin((i + 1) * step * M_PI / 180.0) * radius;
		t1 = (float) (i+1) /slices;
		t2 = 0.0f;
		write << x << " " << 0.0f << " " << z << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << t1 << " " << t2 << endl;

		x = cos(i * step * M_PI / 180.0) * radius;
		z = -sin(i * step * M_PI / 180.0) * radius;
		t1 = (float)i / slices ;
		t2 = 0.0f;
		write << x << " " << 0.0f << " " << z << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << t1 << " " << t2 << endl;

		t1 = 0.5f;
		t2 = 0.5f;
		write << 0.0f << " " << 0.0f << " " << 0.0f << " " << 0.0f << " " << -1.0f << " " << 0.0f << " " << t1 << " " << t2 << endl;
	}

	// Lateral
	for (int i = 0; i < stacks; i++) {
		next_y = y + height / stacks;
		next_r = r - (radius / stacks);
		for (int j = 0; j < slices; j++) {
			x = r * sin(j * ((2 * M_PI) / slices));
			z = r * cos(j * ((2 * M_PI) / slices));
			t1 = (float)j/slices;
			t2 = (float)i/stacks;
			write << x << " " << y << " " << z << " " << sin(j * ((2 * M_PI) / slices)) << " " << height/stacks << " " << cos(j * ((2 * M_PI) / slices)) << " " << t1 << " " << t2 << endl;
			
			x = r * sin((j + 1) * ((2 * M_PI) / slices));
			z = r * cos((j + 1) * ((2 * M_PI) / slices));
			t1 = (float)(j+1) / slices;
			t2 = (float)i / stacks;
			write << x << " " << y << " " << z << " " << sin((j + 1) * ((2 * M_PI) / slices)) << " " << height / stacks << " " << cos((j + 1) * ((2 * M_PI) / slices)) << " " << t1 << " " << t2 << endl;
			
			x = next_r * sin(j * ((2 * M_PI) / slices)); 
			z = next_r * cos(j * ((2 * M_PI) / slices)); 
			t1 = (float)j / slices;
			t2 = (float)(i+1) / stacks;
			write << x << " " << next_y << " " << z << " " << sin(j * ((2 * M_PI) / slices)) << " " << height / stacks << " " << cos(j * ((2 * M_PI) / slices)) << " " << t1 << " " << t2 << endl;

			x = next_r * sin(j * ((2 * M_PI) / slices));
			z = next_r * cos(j * ((2 * M_PI) / slices));
			t1 = (float)j / slices;
			t2 = (float)(i+1) / stacks;
			write << x << " " << next_y << " " << z << " " << sin(j * ((2 * M_PI) / slices)) << " " << height / stacks << " " << cos(j * ((2 * M_PI) / slices)) << " " << t1 << " " << t2 << endl;
			
			x = r * sin((j + 1) * ((2 * M_PI) / slices)); 
			z = r * cos((j + 1) * ((2 * M_PI) / slices)); 
			t1 = (float)(j+1) / slices;
			t2 = (float)i / stacks;
			write << x << " " << y << " " << z << " " << sin((j + 1) * ((2 * M_PI) / slices)) << " " << height / stacks << " " << cos((j + 1) * ((2 * M_PI) / slices)) << " " << t1 << " " << t2 << endl;
			
			x = next_r * sin((j + 1) * ((2 * M_PI) / slices));
			z = next_r * cos((j + 1) * ((2 * M_PI) / slices));
			t1 = (float)(j+1) / slices;
			t2 = (float)(i+1) / stacks;
			write << x << " " << next_y << " " << z << " " << sin((j + 1) * ((2 * M_PI) / slices)) << " " << height / stacks << " " << cos((j + 1) * ((2 * M_PI) / slices)) << " " << t1 << " " << t2 << endl;
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
	float x, y, z, nx, ny, nz, t1, t2;

	for (int i = 0; i < slices; i++) {
		for (int j = 0; j < stacks; j++) {
			nx = 0.0f;
			ny = 0.1f;
			nz = 0.0f;
			t2 = (float)j / stacks;

			x = (outer_radius + inner_radius * cos(phi)) * cos(theta);
			y = (outer_radius + inner_radius * cos(phi)) * sin(theta);
			z = inner_radius * sin(phi);
			t1 = (float)i / slices;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t2 << " " << t1 << endl;
			
			x = (outer_radius + inner_radius * cos(phi)) * cos(theta + delta_slices); 
			y = (outer_radius + inner_radius * cos(phi)) * sin(theta + delta_slices);
			z = inner_radius * sin(phi);
			t1 = (float)(i+1) / slices;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t2 << " " << t1 << endl;

			x = (outer_radius + inner_radius * cos(phi + delta_stacks)) * cos(theta + delta_slices);
			y = (outer_radius + inner_radius * cos(phi + delta_stacks)) * sin(theta + delta_slices);
			z = inner_radius * sin(phi + delta_stacks);
			t1 = (float)(i+1) / slices;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t2 << " " << t1 << endl;

			x = (outer_radius + inner_radius * cos(phi + delta_stacks)) * cos(theta + delta_slices);
			y = (outer_radius + inner_radius * cos(phi + delta_stacks)) * sin(theta + delta_slices);
			z = inner_radius * sin(phi + delta_stacks);
			t1 = (float)(i+1) / slices;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t2 << " " << t1 << endl;

			x = (outer_radius + inner_radius * cos(phi + delta_stacks)) * cos(theta);
			y = (outer_radius + inner_radius * cos(phi + delta_stacks)) * sin(theta);
			z = inner_radius * sin(phi + delta_stacks);
			t1 = (float)i / slices;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t2 << " " << t1 << endl;
			
			x = (outer_radius + inner_radius * cos(phi)) * cos(theta);
			y = (outer_radius + inner_radius * cos(phi)) * sin(theta);
			z = inner_radius * sin(phi);
			t1 = (float)i / slices;
			write << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << t2 << " " << t1 << endl;

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
		else if (!strcmp(figure, "patch") && argc == 5) {
			string input_file = "../" + string(argv[2]);
			float tesselation = atoi(argv[3]);
			string file = "../models/" + string(argv[4]);
			patch(input_file, tesselation, file);
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