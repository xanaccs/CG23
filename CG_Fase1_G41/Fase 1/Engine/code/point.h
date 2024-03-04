#ifndef Point_h
#define Point_h

class Point {
	float x;
	float y;
	float z;

public:
	Point() {
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}

	Point(float x, float y, float z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	float getX() {
		return x;
	}

	float getY() {
		return y;
	}

	float getZ() {
		return z;
	}

	void setX(float x) {
		this->x = x;
	}

	void setY(float y) {
		this->y = y;
	}

	void setZ(float z) {
		this->z = z;
	}
};

#endif