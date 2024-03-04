#ifndef Group_h
#define Group_h


#include <stdio.h>
#include <vector>
#include "transformation.h"
#include "figure.h"

using namespace std;


class Group {

public:
	vector<Transformation> transformations;
	vector<Figure> figures;
	vector<Group> subGroups;

	Group() = default;

	Group(vector<Figure> figures, vector<Group> subGroups, vector<Transformation> transformations) {
		this->transformations = transformations;
		this->figures = figures;
		this->subGroups = subGroups;
	}
};

#endif 
