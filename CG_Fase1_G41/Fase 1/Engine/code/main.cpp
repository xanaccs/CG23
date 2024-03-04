#include <iostream>
#include <string>
#include <vector>
#include "point.h"
#include "figure.h"
#include "../RapidXML/rapidxml.hpp"
#include "../RapidXML/rapidxml_utils.hpp"
#include "../RapidXML/rapidxml_print.hpp"
#include "../RapidXML/rapidxml_iterators.hpp"


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
using namespace rapidxml;

// Vector that stores .3d filenames.
vector<string> files;

// Vector that stores all figures.
vector<Figure> figures;

// Vari�veis para manipular a janela
int width, height;

// Vari�veis para manipular a camara
float x_position, y_position, z_position, x_lookAt, y_lookAt, z_lookAt, x_up, y_up, z_up;

// Vari�veis para mover a camara
float alfa = 0.0f, beta = 0.0f, radius = 5.0f;


// Fun��o que permite a rota��o da camara
void spherical2Cartesian() {
	x_position = radius * cos(beta) * sin(alfa);
	y_position = radius * sin(beta);
	z_position = radius * cos(beta) * cos(alfa);
}


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


// Fun��o que desenha os eixos cartesianos
void drawAxis() {
	glBegin(GL_LINES);
		// put axis drawing in here
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


// Fun��o que desenha todas a figuras armazenadas previamente num vetor
void drawFigures() {
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < figures.size(); i++) {
		for (int j = 0; j < figures[i].getNrVertices(); j++) {
			glColor3f(1.0f, 1.0f, 1.0f);
			Point point = figures[i].getPoint(j);
			glVertex3f(point.getX(), point.getY(), point.getZ());
		}
	}
	glEnd();
}


// Fun��o que renderiza uma cena
void renderScene(void) {
	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();
	gluLookAt(x_position, y_position, z_position,
			  x_lookAt, y_lookAt, z_lookAt,
			  x_up,y_up,z_up);

	drawAxis();

	// Inicial settings
	glPolygonMode(GL_FRONT, GL_LINE);
	glColor3f(1.0f, 1.0f, 1.0f);

	// Draw figures
	drawFigures();

	// End of frame
	glutSwapBuffers();
}


// Fun��o que age quando determinada tecla � pressionada
void processKeys(unsigned char key, int xx, int yy) {
	// Sai da aplica��o de clicar no Esc
	if (key == 27)
		exit(0);
}


// Fun��o que age quando determinada tecla especial � pressionada
void processSpecialKeys(int key, int xx, int yy) {
	switch (key) {
		case GLUT_KEY_RIGHT:
			alfa -= 0.1; break;
		case GLUT_KEY_LEFT:
			alfa += 0.1; break;
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
		case GLUT_KEY_F2:
			radius -= 0.1f;
			if (radius < 0.1f)
				radius = 0.1f;
			break;

		case GLUT_KEY_F1:
			radius += 0.1f;
			break;
	}
	spherical2Cartesian();
	glutPostRedisplay();
}


// Fun��o que faz o parser dos ficheiros .3d
bool readFile(string file_3d) {
		vector<Point> points;
		int numVertices;
		float x, y, z;

		// Abre o ficheiro
		ifstream ficheiro(file_3d);
		if (!ficheiro.is_open()) {
			cout << "Erro ao abrir o ficheiro " << file_3d << endl;
			return false;
		}

		// Armazena o n�mero de v�rtices que a figura ir� ter
		ficheiro >> numVertices;
		for (int i = 0; i < numVertices; i++) {
			// L� do ficheiro as coordenadas
			ficheiro >> x >> y >> z;
			
			// Cria um objeto ponto com as coordenadas lidas
			Point p(x,y,z);

			// Armazena cada ponto no vetor points
			points.push_back(p);
		}
		// Cria uma figura atrav�s de um vetor de pontos
		Figure fig(points);

		// Armazena a figura criada num vetor de figuras
		figures.push_back(fig);
		
		ficheiro.close();
		return true;
}


// Fun��o que faz o parser de um ficheiro xml atrav�s do rapidXML
bool parseDocument(string input_file) {
	// Carrega o conte�do do arquivo XML para a mem�ria
	string path = "../test_files/" + input_file;
	ifstream file(path);
	if (!file) {
		cerr << "Erro ao abrir o arquivo XML" << endl;
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

	// Obt�m o elemento raiz do documento
	xml_node<>* root = doc.first_node("world");

	// Obt�m os valores dos atributos do elemento <window>
	xml_node<>* window = root->first_node("window");
	width = atoi(window->first_attribute("width")->value());
	height = atoi(window->first_attribute("height")->value());

	// Obt�m os valores dos atributos dos elementos <position>, <lookAt> e <up> dentro do elemento <camera>
	xml_node<>* camera = root->first_node("camera");
	x_position = atof(camera->first_node("position")->first_attribute("x")->value());
	y_position = atof(camera->first_node("position")->first_attribute("y")->value());
	z_position = atof(camera->first_node("position")->first_attribute("z")->value());

	x_lookAt = atof(camera->first_node("lookAt")->first_attribute("x")->value());
	y_lookAt = atof(camera->first_node("lookAt")->first_attribute("y")->value());
	z_lookAt = atof(camera->first_node("lookAt")->first_attribute("z")->value());

	x_up = atof(camera->first_node("up")->first_attribute("x")->value());
	y_up = atof(camera->first_node("up")->first_attribute("y")->value());
	z_up = atof(camera->first_node("up")->first_attribute("z")->value());

	// Obt�m os valores dos atributos do elemento <projection>
	xml_node<>* projection = camera->first_node("projection");
	float fov = atof(projection->first_attribute("fov")->value());
	float nearPlane = atof(projection->first_attribute("near")->value());
	float farPlane = atof(projection->first_attribute("far")->value());

	// Obt�m os nomes dos ficheiros .3d e armazena-os num vetor
	xml_node<>* group = root->first_node("group");
	xml_node<>* models = group->first_node("models");
	xml_node<>* model = models->first_node("model");
	while(model) {
		string figure = model->first_attribute("file")->value();
		files.push_back(figure);
		model = model->next_sibling();
	}
	
	// Libera o documento XML da mem�ria
	//delete doc;

	return true;
}


// Fun��o principal do programa
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

	// Faz o parse dos ficheiros .3d
	for (int i = 0; i < files.size(); i++) {
		string file_3d = "../../Generator/models/" + files[i];                
		bool res = readFile(file_3d);
		if (!res) {
			cout << "Erro na leitura de um ficheiro .3d" << endl;
		}
	}

	// Inicializa o GLUT e a janela
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow("TP-CG");

	// Registo das callbacks
	glutIdleFunc(renderScene);   
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);

	// Registo das callbacks quando uma tecla � pressionada
	glutKeyboardFunc(processKeys);
	glutSpecialFunc(processSpecialKeys);

	// OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	// Ciclo principal do GLUT
	glutMainLoop();
	
	return 0;
}
