#ifndef Transformation_h
#define Transformation_h
#include <string>
#include <vector>
using namespace std;

class Transformation {

private:
	string tipo;
	vector <float> params;

	int align;

	// Vectores para armazenar os pontos de controle
	vector<float> px;
	vector<float> py;
	vector<float> pz;

	//float normal[3] = { 0,1,0 };

public:
	Transformation() = default;

	Transformation(string tipo, vector<float> params, vector<float> px, vector<float> py, vector<float> pz, int align) {
		this->tipo = tipo;
		this->params = params;
		this->px = px;
		this->py = py;
		this->pz = pz;
		this->align = align;
	}

	string getTipo() {
		return tipo;
	}

	vector<float> getParams() {
		return params;
	}

	vector<float> getPx() {
		return px;
	}

	vector<float> getPy() {
		return py;
	}

	vector<float> getPz() {
		return pz;
	}

	int getAlign() {
		return align;
	}

	void setTipo(string tipo) {
		this->tipo = tipo;
	}

	void setParams(vector <float> params) {
		this->params = params;
	}

	void setPx(vector<float> px) {
		this->px = px;
	}

	void setPy(vector<float> py) {
		this->py = py;
	}

	void setPz(vector<float> pz) {
		this->pz = pz;
	}

	void setAlign(int align) {
		this->align = align;
	}
};

#endif
