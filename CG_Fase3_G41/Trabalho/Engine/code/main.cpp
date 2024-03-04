#include <iostream>
#include <string>
#include <vector>
#include "group.h"
#include "point.h"
#include "figure.h"
#include "transformation.h"
#include "../RapidXML/rapidxml.hpp"
#include "../RapidXML/rapidxml_utils.hpp"
#include "../RapidXML/rapidxml_print.hpp"
#include "../RapidXML/rapidxml_iterators.hpp"


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif


#define _USE_MATH_DEFINES
#include <math.h>


using namespace std;
using namespace rapidxml;

// Variável global que irá armazenar a informação relevante do ficheiro XML
Group group;

// Variáveis para definir a janela no modo default
int width = 512, height = 512;

// Variáveis para definir a camera no modo default
float x_position = 20.0f, y_position = 10.0f, z_position = 20.0f;
float x_lookAt = 0.0f, y_lookAt = 0.0f, z_lookAt = 0.0f;
float x_up = 0.0f, y_up = 1.0f, z_up = 0.0f;
float fov = 60.0f, nearPlane = 1.0f, farPlane = 1000.0f;

// Variáveis para movimentar a camera
float alfa = 0.0f, beta = 0.0f, radius = 25.0f;

float rot_x = 0.0, rot_y = 0.0;
float camX = 0.0, camY = 0.0, camZ = 0.0;

float last_x, last_y;
bool valid = false;

// Variáveis para modificar outras configurações
int modoDesenho = GL_LINE;
bool axis = true;
bool orbits = true;

// Variável que define o nível de tesselação da curva
float tesselation = 0.01;
// Variável que define a posição do objeto
float Pos[3] = { 0,0,0 };
// Variável que define a derivada do objeto
float Deriv[3] = { 0,0,0 };

// VBO para armazenar os pontos das figuras
vector<float> vertexB;
GLuint buffers[1];
GLuint ptr = 0;


//
void changeSize(int w, int h) {
	// Prevent a divide by zero, when window is too short
	// (you cant make a window with zero width).
	if(h == 0)
		h = 1;

	// compute window's aspect ratio 
	float ratio = w * 1.0 / h;

	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);

	// Load Identity Matrix
	glLoadIdentity();
	
	// Set the viewport to be the entire window
    glViewport(0, 0, w, h);

	// Set perspective
	gluPerspective(45.0f ,ratio, 1.0f ,1000.0f);

	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

//
void buildRotMatrix(float* x, float* y, float* z, float* m) {
	m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
	m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
	m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}

// Função que calcula o produto de 2 vetores
void cross(float* a, float* b, float* res) {
	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

// Função que normaliza um vetor
void normalize(float* a) {
	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
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


void getCatmullRomPoint(float t, float* p0, float* p1, float* p2, float* p3, float* pos, float* deriv) {
	// Matriz Catmull-Rom
	float m[4][4] = {{-0.5f,  1.5f, -1.5f,  0.5f},
					 { 1.0f, -2.5f,  2.0f, -0.5f},
					 {-0.5f,  0.0f,  0.5f,  0.0f},
					 { 0.0f,  1.0f,  0.0f,  0.0f}};

	// 
	float T[4] = { t * t * t, t * t, t, 1 };
	float T_deriv[4] = { 3 * t * t,2 * t, 1, 0 };


	for (int i = 0; i < 3; i++) { // x, y, z
		float p[4] = { p0[i], p1[i], p2[i], p3[i] };
		float A[4];

		//Compute vector A = M * P 
		multMatrixVector((float*)m, p, A);

		//Compute pos[i] = T * A
		pos[i] = T[0] * A[0] + T[1] * A[1] + T[2] * A[2] + T[3] * A[3];

		//compute deriv[i] = T' * A
		deriv[i] = T_deriv[0] * A[0] + T_deriv[1] * A[1] + T_deriv[2] * A[2] + T_deriv[3] * A[3];
	}
}


// Given  a time, returns the point in the curve
void getGlobalCatmullRomPoint(float gt, vector<float> px, vector<float> py, vector<float> pz, float* pos, float* deriv) {
	int POINT_COUNT = (int)px.size();

	// Tempo local (real)
	float t = gt * POINT_COUNT;
	int index = floor(t);		
	t = t - index;			

	// Determina os indices dos pontos de controlo
	int indices[4];
	indices[0] = (index + POINT_COUNT - 1) % POINT_COUNT;
	indices[1] = (indices[0] + 1) % POINT_COUNT;
	indices[2] = (indices[1] + 1) % POINT_COUNT;
	indices[3] = (indices[2] + 1) % POINT_COUNT;

	// Determina os pontos de controlo a serem utilizados
	float p[4][3];
	p[0][0] = px[indices[0]];
	p[0][1] = py[indices[0]];
	p[0][2] = pz[indices[0]];

	p[1][0] = px[indices[1]];
	p[1][1] = py[indices[1]];
	p[1][2] = pz[indices[1]];

	p[2][0] = px[indices[2]];
	p[2][1] = py[indices[2]];
	p[2][2] = pz[indices[2]];

	p[3][0] = px[indices[3]];
	p[3][1] = py[indices[3]];
	p[3][2] = pz[indices[3]];

	getCatmullRomPoint(t, p[0], p[1], p[2], p[3], pos, deriv);
}

void renderCatmullRomCurve(vector<float> px, vector<float> py, vector<float> pz) {
	// draw curve using line segments with GL_LINE_LOOP
	glBegin(GL_LINE_LOOP);
	for (float i = 0; i < 1; i += tesselation) {
		getGlobalCatmullRomPoint(i, px, py, pz, Pos, Deriv);
		glVertex3f(Pos[0], Pos[1], Pos[2]);
	}
	glEnd();
}

// Função que permite a rotação da camera
void camera() {
	glRotatef(rot_x, 1.0, 0.0, 0.0);	// Along X axis
	glRotatef(rot_y, 0.0, 1.0, 0.0);    // Along Y axis
	glTranslatef(-camX, -camY, -camZ);  // Move the camera
}


// Função que converte coordenadas
void spherical2Cartesian() {
	x_position = radius * cos(beta) * sin(alfa);
	y_position = radius * sin(beta);
	z_position = radius * cos(beta) * cos(alfa);
}


// Função que desenha os eixos cartesianos
void drawAxis() {
	glBegin(GL_LINES);
		// X in red
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(-100.0f, 0.0f, 0.0f);
		glVertex3f(100.0f, 0.0f, 0.0f);

		// Y in green
		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(0.0f, -100.0f, 0.0f);
		glVertex3f(0.0f, 100.0f, 0.0f);

		// Z in blue
		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, 0.0f, -100.0f);
		glVertex3f(0.0f, 0.0f, 100.0f);
	glEnd();
}

//Desenha a orbita dos planetas
void drawOrbit(float raio) {
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i <= 360; i++) {
		float angle = (i * M_PI) / 180; // radianos
		glVertex3f(cos(angle) * raio, 0, sin(angle) * raio);
	}
	glEnd();
}

// Função que desenha o anel de asteroides
void drawAsteroids() {
	glColor3f(0.478f, 0.458f, 0.458f);

	srand(100000);
	for (int i = 0; i < 500; i++) {
		float x = 30 * (rand() / (float)RAND_MAX) - 15;  // x coordinate between -15 and 15
		float y = 4 * (rand() / (float)RAND_MAX) - 2;    // y coordinate between -2 and 2
		float z = 30 * (rand() / (float)RAND_MAX) - 15;  // z coordinate between -15 and 15

		if (169 < x*x + z*z && x*x + z*z <= 225) {
			glPushMatrix();
			glTranslatef(x, y, z);
			glScalef(0.01, 0.01, 0.01);
			glutSolidSphere(4, 80, 80);
			glPopMatrix();
		}
	}
	glClearColor;
}

// Função que desenha as figuras presentes em cada grupo
void drawFigures(Group g) {
	// Guarda o estado atual da matriz
	glPushMatrix();

	// Percorremos e aplicamos as transformações
	for (int k = 0; k < g.transformations.size(); k++) {
		Transformation transf = g.transformations[k];
		vector<float> params = transf.getParams();

		if (transf.getTipo() == "translate") {
			float time = params[0];

			// Objeto estático
			if (time == 0) {
				// Desenha as órbitas dos planetas e das luas tendo em conta o raio(x)
				if(orbits)
					drawOrbit(params[1]);
				glTranslatef(params[1], params[2], params[3]);
			}
			else { // Objeto com movimento
				vector<float> px = transf.getPx();
				vector<float> py = transf.getPy();
				vector<float> pz = transf.getPz();
				int align = transf.getAlign();

				if (orbits)
					renderCatmullRomCurve(px, py, pz);
				
				// Tempo desde o ínicio do programa / tempo de duração da animação (em milissegundos)
				float gt = (float)glutGet(GLUT_ELAPSED_TIME) / (time * 1000);
				getGlobalCatmullRomPoint(gt, px, py, pz, Pos, Deriv);
				glTranslatef(Pos[0], Pos[1], Pos[2]);

				if (align) {
					float rot[16];
					float* X = Deriv;
					float Y[3] = { 0,1,0 }; // Variável que define a normal do objeto
					float Z[3];

					// Z = X x Y
					cross(X, Y, Z);
					// Y = Z x X
					cross(Z, X, Y);

					// Normaliza os vetores
					normalize(X);
					normalize(Y);
					normalize(Z);

					buildRotMatrix(X, Y, Z, rot);

					glMultMatrixf(rot);
				}
			}
		}
		else if (transf.getTipo() == "scale") {
			float x = params[0];
			float y = params[1];
			float z = params[2];
			glScalef(x, y, z);
		}
		else if (transf.getTipo() == "rotate") {
			float angle = params[0];
			float time = params[1];
			float x = params[2];
			float y = params[3];
			float z = params[4];

			// Objeto estático
			if (time == 0) {
				glRotatef(angle, x, y, z);
			}
			else{ // Objeto com movimento
				float t = (float)glutGet(GLUT_ELAPSED_TIME) / (time * 1000);
				glRotatef(t * 360, x, y, z);
			}
		}
		else if (transf.getTipo() == "color") {
			float red = params[0];
			float green = params[1];
			float blue = params[2];
			glColor3f(red, green, blue);
		}
	}

	// Percorremos e desenhamos as figuras com VBOs
	for (int i = 0; i < g.figures.size(); i++) {
		Figure fig = g.figures[i];
		int n_vertices = fig.getNrVertices();

		// Ativa e define o tipo de buffer (array)
		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		// Forma como organizamos os vértices
		glVertexPointer(3, GL_FLOAT, 0, 0);
		// Desenha os triângulos (conjuntos de 3 vértices)
		glDrawArrays(GL_TRIANGLES, ptr, n_vertices);
		
		// O pointer passa para o próximo conjunto de vértices
		ptr = ptr + n_vertices;
	}

	// Percorremos e analisamos os subgrupos
	if (!g.subGroups.empty()) {
		for (int i = 0; i < g.subGroups.size(); i++) {
			drawFigures(g.subGroups[i]);
		}
	}

	// Voltamos ao estado anterior
	glPopMatrix();
}


// Função que renderiza uma cena
void renderScene(void) {
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set the camera
	glLoadIdentity();

	spherical2Cartesian();

	gluLookAt(x_position, y_position, z_position,
			  x_lookAt, y_lookAt, z_lookAt,
			  x_up, y_up, z_up);

	// Draw axis
	if (axis)
		drawAxis();

	// Inicial settings
	glPolygonMode(GL_FRONT, modoDesenho);

	// Camera rotation
	camera();

	// Draw figures
	drawFigures(group);
	
	// Reinicializa o pointer
	ptr = 0;

	// Draw Asteroids
	drawAsteroids();
	
	// End of frame
	glutSwapBuffers();
}

// Função que faz o parser dos ficheiros .3d
vector<Figure> readFile(string file_3d) {
		vector<Point> points;
		vector<Figure> figures;
		int numVertices;
		float x, y, z;

		// Abre o ficheiro
		ifstream ficheiro(file_3d);
		if (!ficheiro.is_open()) {
			cout << "Erro ao abrir o ficheiro " << file_3d << endl;
			return figures;
		}

		// Armazena o número de vértices que a figura irá ter
		ficheiro >> numVertices;
		for (int i = 0; i < numVertices; i++) {
			// Lê do ficheiro as coordenadas
			ficheiro >> x >> y >> z;
			
			// Cria um objeto ponto com as coordenadas lidas
			Point p(x,y,z);
			
			// Armazena cada ponto no vetor points
			points.push_back(p);

			// Armazena os vértices num VBO
			vertexB.push_back(x);
			vertexB.push_back(y);
			vertexB.push_back(z);
		}

		// Cria uma figura através de um vetor de pontos
		Figure fig(points);

		// Armazena a figura criada num vetor de figuras
		figures.push_back(fig);
		
		ficheiro.close();

		// Retorna um vetor de figuras
		return figures;
}


// Função que trata um nodo grupo
Group loadGroup(xml_node<>* group_node) {
	// Vetor para armazenar as transformações do grupo
	vector<Transformation> transformations;

	xml_node<>* transform_node = group_node->first_node("transform");
	while (transform_node) {
		// Analisa cada nodo dentro do bloco transform
		xml_node<>* node = transform_node->first_node();
		while (node) {
			// Transformação
			Transformation transform;

			// Vetor para armazenar os parâmetros da transformação
			vector <float> params;
			vector <float> px;
			vector <float> py;
			vector <float> pz;

			// Parse do elemneto translate
			if (strcmp(node->name(), "translate") == 0) {
				xml_attribute<>* x = node->first_attribute("x");
				if (x) { // Elemento estático
					params.push_back(0.0f); // time = 0.0f
					params.push_back(atof(x->value()));
					params.push_back(atof(node->first_attribute("y")->value()));
					params.push_back(atof(node->first_attribute("z")->value()));
				}
				else { // Elemento com movimento
					params.push_back(atof(node->first_attribute("time")->value()));
					char* align_prev = node->first_attribute("align")->value();
					int align = 1;

					if (!strcmp(align_prev, "false")) {
						align = 0;
					}

					xml_node<>* point = node->first_node("point");
					while (point) {
						px.push_back(atof(point->first_attribute("x")->value()));
						py.push_back(atof(point->first_attribute("y")->value()));
						pz.push_back(atof(point->first_attribute("z")->value()));

						point = point->next_sibling("point");
					}

					transform.setPx(px);
					transform.setPy(py);
					transform.setPz(pz);
					transform.setAlign(align);
				}

				transform.setTipo("translate");
				transform.setParams(params);
				transformations.push_back(transform);
			}

			// Parse do elemento rotate
			if (strcmp(node->name(), "rotate") == 0) {
				xml_attribute<>* rotate = node->first_attribute("angle");
				if (rotate){ // Elemento estático
					params.push_back(atof(rotate->value()));
				}
				else { // Elemento com movimento
					params.push_back(0.0f);
				}

				xml_attribute<>* time = node->first_attribute("time");
				if (time) {
					params.push_back(atof(time->value()));
				}
				else {
					params.push_back(0.0f);
				}

				params.push_back(atof(node->first_attribute("x")->value()));
				params.push_back(atof(node->first_attribute("y")->value()));
				params.push_back(atof(node->first_attribute("z")->value()));
				transform.setTipo("rotate");
				transform.setParams(params);
				transformations.push_back(transform);
			}

			// Parse do elemento scale
			if (strcmp(node->name(), "scale") == 0) {
				params.push_back(atof(node->first_attribute("x")->value()));
				params.push_back(atof(node->first_attribute("y")->value()));
				params.push_back(atof(node->first_attribute("z")->value()));
				transform.setTipo("scale");
				transform.setParams(params);
				transformations.push_back(transform);
			}

			// Parse do elemento color
			if (strcmp(node->name(), "color") == 0) {
				params.push_back(atof(node->first_attribute("x")->value()));
				params.push_back(atof(node->first_attribute("y")->value()));
				params.push_back(atof(node->first_attribute("z")->value()));
				transform.setTipo("color");
				transform.setParams(params);
				transformations.push_back(transform);
			}
			node = node->next_sibling();
		}
		transform_node = transform_node->next_sibling("transform");
	}

	// Vetor para armazenar as figuras do grupo
	vector<Figure> figures;

	xml_node<>* models = group_node->first_node("models");
	if (models) {
		xml_node<>* model = models->first_node("model");

		while (model) {
			string figure = model->first_attribute("file")->value();

			// Abrir e guardar a informação do ficheiro
			string file_3d = "../../Generator/models/" + figure;

			figures = readFile(file_3d);
			if (figure.empty()) {
				cout << "Erro na leitura do ficheiro: " + figure << endl;
			}

			model = model->next_sibling();
		}
	}

	// Vetor para armazenar os subgrupos do grupo
	vector<Group> subGroups;

	xml_node<>* subGroup = group_node->first_node("group");
	while (subGroup) {
		Group g = loadGroup(subGroup);
		subGroups.push_back(g);
		subGroup = subGroup->next_sibling("group");
	}

	// Criação do grupo com os seus elementos
	Group grupo(figures, subGroups, transformations);
	
	return grupo;
}


// Função que faz o parser de um ficheiro xml através do rapidXML
bool parseDocument(string input_file) {
	// Carrega o conteúdo do arquivo XML para a memória
	string path = "../test_files/" + input_file;
	ifstream file(path);
	if (!file) {
		cerr << "Erro ao abrir o arquivo XML: " << input_file << endl;
		return false;
	}

	string xml;
	file.seekg(0, ios::end);
	xml.reserve(file.tellg());
	file.seekg(0, ios::beg);
	xml.assign((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());

	// Faz o parse do documento XML
	xml_document<> doc;
	doc.parse<0>(&xml[0]);

	// Obtém o elemento raiz do documento
	xml_node<>* root = doc.first_node("world");
	if (root == nullptr) {
		cerr << "Erro na estrutura do arquivo XML" << endl;
		return false;
	}

	// Obtém os valores dos atributos do elemento <window>
	xml_node<>* window = root->first_node("window");
	if (window) {
		width = atoi(window->first_attribute("width")->value());
		height = atoi(window->first_attribute("height")->value());
	}

	// Obtém os valores dos atributos do elemento <camera>
	xml_node<>* camera = root->first_node("camera");
	if (camera) {
		// Obtém os valores dos atributos do elementos <position>
		xml_node<>* position = camera->first_node("position");
		if (position) {
			x_position = atof(position->first_attribute("x")->value());
			y_position = atof(position->first_attribute("y")->value());
			z_position = atof(position->first_attribute("z")->value());
		}

		// Obtém os valores dos atributos dos elementos <lookAt>
		xml_node<>* lookAt = camera->first_node("lookAt");
		if (lookAt) {
			x_lookAt = atof(lookAt->first_attribute("x")->value());
			y_lookAt = atof(lookAt->first_attribute("y")->value());
			z_lookAt = atof(lookAt->first_attribute("z")->value());
		}

		// Obtém os valores dos atributos dos elementos <up>
		xml_node<>* up = camera->first_node("up");
		if (up) {
			x_up = atof(up->first_attribute("x")->value());
			y_up = atof(up->first_attribute("y")->value());
			z_up = atof(up->first_attribute("z")->value());
		}

		// Obtém os valores dos atributos do elemento <projection>
		xml_node<>* projection = camera->first_node("projection");
		if (projection) {
			fov = atof(projection->first_attribute("fov")->value());
			nearPlane = atof(projection->first_attribute("near")->value());
			farPlane = atof(projection->first_attribute("far")->value());
		}
	}

	// Obtém os de todos os elementos de todos os grupos
	xml_node<>* group_node = root->first_node("group");
	if(group_node == nullptr){
		cerr << "Erro: leitura do grupo no arquivo XML" << endl;
		return false;
	}

	// Trata os elementos de um grupo
	group = loadGroup(group_node);

	return true;
}


// Função que age quando determinada tecla é pressionada
void processKeys(unsigned char key, int x, int y) {
	float xrotrad, yrotrad;

	switch (key) {
		// Sai da aplicação de clicar no Esc
		case 27:
			exit(0);

		// Muda o modo de desenho
		case '1':
			modoDesenho = GL_LINE;
			break;
		case '2':
			modoDesenho = GL_POINT;
			break;
		case '3':
			modoDesenho = GL_FILL;
			break;
		
		// Mostra ou esconde os eixos:
		case 'e':
			axis = !axis;
			break;

		// Mostra ou esconde as órbitas dos astros:
		case 'o':
			orbits = !orbits;
			break;

		// Movimento dos objetos
		case 'w':
			xrotrad = (rot_x / 180 * M_PI);
			yrotrad = (rot_y / 180 * M_PI);
			camX += float(sin(yrotrad));
			camZ -= float(cos(yrotrad));
			camY -= float(sin(xrotrad));
			break;
		case 's':
			yrotrad = (rot_y / 180 * M_PI);
			xrotrad = (rot_x / 180 * M_PI);
			camX -= float(sin(yrotrad));
			camZ += float(cos(yrotrad));
			camY += float(sin(xrotrad));
			break;
		case 'd':
			yrotrad = (rot_y / 180 * M_PI);
			camX += float(cos(yrotrad)) * 0.2;
			camZ += float(sin(yrotrad)) * 0.2;
			break;
		case 'a':
			yrotrad = (rot_y / 180 * M_PI);
			camX -= float(cos(yrotrad)) * 0.2;
			camZ -= float(sin(yrotrad)) * 0.2;
			break;
		default:
			break;
	}
	glutPostRedisplay();
}


// Função que age quando determinada tecla especial é pressionada
void processSpecialKeys(int key, int xx, int yy) {
	// Movimento da camera
	switch (key) {
		case GLUT_KEY_RIGHT:
			alfa -= 0.1;
			break;
		case GLUT_KEY_LEFT:
			alfa += 0.1;
			break;
		case GLUT_KEY_UP:
			beta += 0.1f;
			if (beta > 1.5f)
				beta = 1.5f;
			break;
		case GLUT_KEY_DOWN:
			beta -= 0.1f;
			if (beta < -1.5f)
				beta = -1.5f;
			break;
	}

	glutPostRedisplay();
}


// Função que o utilizador iterage com o rato
void processMouseButtons(int button, int state, int x, int y) {
	last_x = x;
	last_y = y;
	valid = state == GLUT_DOWN;
}


// Função que o utilizador iterage com o rato
void processMouseMotion(int x, int y) {
	if (valid) {
		int diffx = x - last_x;				//check the difference between the current x and the last x position
		int diffy = y - last_y;				//check the difference between the current y and the last y position
		last_x = x;							//set last_x to the current x position
		last_y = y;							//set last_y to the current y position
		if (rot_x < -87) rot_x = -87;
		else {
			if (rot_x > 87) rot_x = 87;
			else {
				rot_x -= (float)diffy / 10;  //set the rot_x to rot_x with the addition of the difference in the y position
			}
		}
		rot_y -= (float)diffx / 10;			//set the rot_x to rot_y with the addition of the difference in the x position
		glutPostRedisplay();
	}
}


// Função que imprime no ecrã um menu de ajuda com as possiveis ações
void printHelp() {
	cout << "------------ HELP ------------" << endl;
	cout << "Press the key:" << endl << endl;
	cout << "[w] - move forward" << endl;
	cout << "[s] - move backward" << endl;
	cout << "[d] - move to the right" << endl;
	cout << "[s] - move to the left" << endl;
	cout << "[1] - change drawing mode to line" << endl;
	cout << "[2] - change drawing mode to point" << endl;
	cout << "[3] - change drawing mode to fill" << endl;
	cout << "[e] - show/hide axis" << endl;
	cout << "[o] - show/hide orbits" << endl;
	cout << "[key_up] - move the camera up" << endl;
	cout << "[key_down] - move the camera down" << endl;
	cout << "[key_right] - move the camera to the right" << endl;
	cout << "[key_left] - move the camera to the left" << endl;
	cout << "mouse - move camera in all directions" << endl;
}	

void init(int argc, char** argv) {
	// Inicializa o GLUT e a janela
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow("TP-CG");

	// Print help menu
	printHelp();

	// Registo das callbacks
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);

	// Registo das callbacks quando uma tecla é pressionada
	glutKeyboardFunc(processKeys);
	glutSpecialFunc(processSpecialKeys);

	// Registo das callbacks quando o utilizador interage com o rato
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);

	glewInit();

	// OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT, GL_LINE);
	glEnableClientState(GL_VERTEX_ARRAY);

	// Cria o VBO das figuras
	glGenBuffers(1, buffers);
	glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertexB.size(), vertexB.data(), GL_STATIC_DRAW);
}

// Função principal do programa
int main(int argc, char **argv) {
	// Verifica se algum ficheiro xml foi passado como argumento
	if (argc < 2) {
		cout << "Nenhum ficheiro XML indicado" << endl;
		return 1;
	}

	// Faz o parse do ficheiro XML
	if (!parseDocument(argv[1])) {
		cout << "Erro na leitura do ficheiro XML" << endl;
		return 1;
	}

	// Inicializa o Glut, o Glew, o VBO e a janela
	init(argc,argv);

	// Ciclo principal do GLUT
	glutMainLoop();

	return 0;
}
