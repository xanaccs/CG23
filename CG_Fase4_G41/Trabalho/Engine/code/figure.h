#ifndef figure_h
#define figure_h

#include <vector>
#include "point.h"
#include "color.h"

using namespace std;

class Figure {

private:
	vector <Point> vertices;
	Color color;
	string texture_file;

public:
	Figure() = default;

	Figure(vector<Point> vertices, Color color, string texture_file) {
		this->vertices = vertices; 
		this->color = color;
		this->texture_file = texture_file;
	}

	std::vector<Point> getVertices() {
		return vertices;
	}

	int getNrVertices() {
		return vertices.size();
	}

	Color getColor() {
		return color;
	}

	string getTextureID() {
		return texture_file;
	}
};

#endif