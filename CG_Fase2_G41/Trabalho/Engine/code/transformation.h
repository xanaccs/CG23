#ifndef Transformation_h
#define Transformation_h
#include <string>
#include <vector>
using namespace std;

class Transformation {

private:
	string tipo;
	vector <float> params;

public:
	Transformation() = default;

	Transformation(string tipo, vector<float> params) {
		this->tipo = tipo;
		this->params = params;
	}

	string getTipo() {
		return tipo;
	}

	vector<float> getParams() {
		return params;
	}

	void setTipo(string tipo) {
		this->tipo = tipo;
	}

	void setParams(vector <float> params) {
		this->params = params;
	}
};

#endif
