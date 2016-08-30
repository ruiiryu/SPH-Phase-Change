#include <cmath>
using namespace std;


#include "vec3.h"


Vector3 Point3::operator-(const Point3& p) const{
    return Vector3(v[0]-p[0], v[1]-p[1], v[2]-p[2]);
}

Point3 Point3::operator+(const Vector3& w) const{
    return Point3(v[0]+w[0], v[1]+w[1], v[2]+w[2]);
}

void Point3::operator+=(const Vector3& w){
    v[0]+=w[0]; v[1]+=w[1]; v[2]+=w[2];
}

Vector3 Vector3::operator/(double r){
    return Vector3(v[0]/r, v[1]/r, v[2]/r);
}

Vector3 Vector3::operator+(const Vector3& w){
    return Vector3(v[0]+w[0], v[1]+w[1], v[2]+w[2]);
}

void Vector3::operator+=(const Vector3& w){
    v[0]+=w[0]; v[1]+=w[1]; v[2]+=w[2];
}

void Vector3::operator-=(const Vector3& w){
    v[0]-=w[0]; v[1]-=w[1]; v[2]-=w[2];
}

double Vector3::squared_length(){
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

double Vector3::length(){
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

Vector3 Vector3::normalize(){
    return (*this)/length();
}

Vector3 operator*(const Vector3& v, double r){
    return Vector3(v[0]*r, v[1]*r, v[2]*r);
}

Vector3 operator*(double r, const Vector3& v){
    return Vector3(v[0]*r, v[1]*r, v[2]*r);
}

