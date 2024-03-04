#ifndef light_h
#define light_h

#include "point.h"
#include <string>

using namespace std;

class Light {

private:
	string type;
	Point coordsPos;
	Point coordsDir;
	float cutoff;

public:
	Light() = default;

	Light(string type, Point coordsPos,Point coordsDir, float cutoff) {
		this->type = type;
		this->coordsPos = coordsPos;
		this->coordsDir = coordsDir;
		this->cutoff = cutoff;
	}

	string getType() {
		return type;
	}

	Point getCoordsPos() {
		return coordsPos;
	}

	Point getCoordsDir() {
		return coordsDir;
	}

	float getCutoff() {
		return cutoff;
	}
};

#endif