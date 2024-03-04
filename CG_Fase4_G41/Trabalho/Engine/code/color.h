#include "point.h"
#ifndef Color_h
#define Color_h

class Color {

private:
	Point diff;
	bool hasDiff;
	Point amb;
	bool hasAmb;
	Point spec;
	bool hasSpec;
	Point emi;
	bool hasEmi;
	float shininess;
	bool hasShininess;

public:
	Color() = default;

	Color(Point diff, bool hasDiff, Point amb, bool hasAmb, Point spec, bool hasSpec, Point emi, bool hasEmi, float shininess, bool hasShininess) {
		this->diff = diff;
		this->hasDiff = hasDiff;
		this->amb = amb;
		this->hasAmb = hasAmb;
		this->spec = spec;
		this->hasSpec = hasSpec;
		this->emi = emi;
		this->hasEmi = hasEmi;
		this->shininess = shininess;
		this->hasShininess = hasShininess;
	}

	Point getDiff(){
		return diff;
	}

	Point getAmb(){
		return amb;
	}

	Point getSpec(){
		return spec;
	}

	Point getEmi(){
		return emi;
	}

	float getShininess(){
		return shininess;
	}

	bool getHasDiff() {
		return hasDiff;
	}

	bool getHasAmb() {
		return hasAmb;
	}

	bool getHasSpec() {
		return hasSpec;
	}

	bool getHasEmi() {
		return hasEmi;
	}

	bool getHasShininess() {
		return hasShininess;
	}

	void setDiff(Point diff) {
		this->diff = diff;
	}

	void setAmb(Point amb) {
		this->amb = amb;
	}

	void setSpec(Point spec) {
		this->spec = spec;
	}

	void setEmi(Point emi) {
		this->emi = emi;
	}

	void setShininess(float shi) {
		this->shininess = shi;
	}

	void setHasDiff(bool hasDiff) {
		this->hasDiff = hasDiff;
	}

	void setHasAmb(bool hasAmb) {
		this->hasAmb = hasAmb;
	}

	void setHasSpec(bool hasSpec) {
		this->hasSpec = hasSpec;
	}

	void setHasEmi(bool hasEmi) {
		this->hasEmi = hasEmi;
	}

	void setHasShininess(float hasShi) {
		this->hasShininess = hasShi;
	}

};

#endif
