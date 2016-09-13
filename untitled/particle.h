#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <algorithm>    // std::find
#include <cmath>
#include "vec3.h"
using namespace std;

struct Particle{
    Particle(double x, double y, double z, double r){
        pos = Point3(x,y,z);
        pre_pos = Point3(x,y,z);
        density = 0.9983f;
        mass = MASS;//unit: kg
        pressure = 0.0f;
        color = Point3(0.0f, 0.0f, 1.0f);
        velocity = velocity_delta = Vector3(0, 0, 0);
        radius = r;
        darkmatter = false;
        //type = Particle::ICE;
        type = Particle::WATER;
        heatvalue = 0;
        T[0] = 273.0f;  T[1] = 0;
        nbs_notAir[0] = 0;  nbs_notAir[1] = 0;  nbs_notAir[2] = 0;
        nbs_notAir[3] = 0;  nbs_notAir[4] = 0;  nbs_notAir[5] = 0;
    }
    unsigned id;
    Point3 pos;
    Point3 pre_pos;
    Point3 color;
    Vector3 acc;
    Vector3 pre_velocity;
    Vector3 velocity;
    Vector3 velocity_delta;

    Vector3 angular_velocity;

    Vector3 inter_tension_ice;
    double density;//actual density
    double mass;
    double pressure;
    double radius;
    Vector3 force_pressure;
    Vector3 force_vis;

    double ci;//color for interface tension
    double cs;//color for surface tension

    double T[2];
    unsigned nbs_notAir[6];

    double dist_center;
    double heatvalue;

    Vector3 relative_pos;
    double area;
    Vector3 Fn;
    Vector3 Fs;
    Vector3 Fd;

    Vector3 vT;

    enum {ICE,WATER,AIR,HEAT};
    unsigned type;
    bool darkmatter;
    vector<pair<unsigned, double> > nbs;
    bool high;


    static const double MASS;
    static double max_v;
    static double min_v;
    static double max_pressure;
    static double min_pressure;
    static double max_density;
    static double min_density;
    static double max_temperature;
    static double min_temperature;

};


#endif // PARTICLE_H
