#ifndef figure_h
#define figure_h

#include <vector>
#include "point.h"

using namespace std;

class Figure {

private:
	vector <Point> vertices;

public:
	Figure() = default;

	Figure(vector<Point> vertices) {    
		this->vertices = vertices; 
	}

	std::vector<Point> getVertices() {
		return vertices;
	}

	int getNrVertices() {
		return vertices.size();
	}

	Point getPoint(int index) {
		return vertices[index];
	}

	void addPoint(Point p) {
		vertices.push_back(p);
	}
};

#endif