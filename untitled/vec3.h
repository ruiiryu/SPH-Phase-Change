#ifndef VEC3_H
#define VEC3_H

#ifndef NULL
#define NULL ((void *)0)
#endif

#ifdef _DEBUG
# include "assert.h"
#else
# define assert(x) { }
#endif

#include <GL/gl.h>
#include <math.h>

class Vector3;
class Point3{
public:
    Point3(float x=0, float y=0, float z=0){
        v[0]=x; v[1]=y; v[2]=z;
    }

    Point3(float r[3]){
        v[0]=r[0]; v[1]=r[1]; v[2]=r[2];
    }

    Point3(const Point3& p){
        v[0]=p[0]; v[1]=p[1]; v[2]=p[2];
    }


    float x(){ return v[0]; }
    float y(){ return v[1]; }
    float z(){ return v[2]; }
    float& operator[] (unsigned i){ return v[i]; }
    const float& operator[] (unsigned i) const { return v[i]; }

    Vector3 operator- (const Point3& p) const;
    Point3 operator+ (const Vector3& w) const;

    void operator+=(const Vector3& v);

    float v[3];
};

class Vector3 : public Point3{
public:
    Vector3(double x=0, double y=0, double z=0):Point3(x, y, z){
    }

    Vector3 operator/(double r);
    Vector3 operator+(const Vector3& w);

    void operator+=(const Vector3& v);
    void operator-=(const Vector3& v);

    double squared_length();
    double length();
    Vector3 normalize();

    friend Vector3 operator*(const Vector3& v, double r);
    friend Vector3 operator*(double r, const Vector3& v);
};

inline Vector3 cross(const Vector3& v, const Vector3& w){
    return Vector3(v.v[1]*w.v[2]-v.v[2]*w.v[1], v.v[2]*w.v[0]-v.v[0]*w.v[2], v.v[0]*w.v[1]-v.v[1]*w.v[0]);
}

inline double dot(const Vector3& v, const Vector3& w){
    return v.v[0]*w.v[0]+v.v[1]*w.v[1]+v.v[2]*w.v[2];
}



#endif // VEC3_H
