#ifndef SPH_H
#define SPH_H

#include <vector>
#include <algorithm>    // std::find
using namespace std;


#include "particle.h"
#include "grid.h"


const double PI=3.1415926535;


struct Rigid_Body{
    vector<unsigned> listICEparticle;
    double mass;
    double radii;
    Point3 center;
    vector<unsigned> side1;
    vector<unsigned> side2;
    vector<unsigned> side3;
    vector<unsigned> side4;
    vector<unsigned> side5;
    vector<unsigned> side6;

    void remove(vector<unsigned>& vec, unsigned item){
        vector<unsigned>::iterator it;
        it = find(vec.begin(), vec.end(), item);
        if (it != vec.end()) vec.erase(vec.begin(),it);
        else std::cout << "Element not found in myvector\n";
    }

    bool finditem(vector<unsigned>& vec, unsigned item){
        std::vector<unsigned>::iterator it;
        it = std::find(vec.begin(), vec.end(), item);
        if (it != vec.end()) return true;
        else return false;
    }

};

struct Liquid{
    vector<unsigned> listWaterparticle;
    vector<unsigned> side1;
    vector<unsigned> side2;
    vector<unsigned> side3;
    vector<unsigned> side4;
    vector<unsigned> side5;
    vector<unsigned> side6;

    bool finditem(vector<unsigned>& vec, unsigned item){
        std::vector<unsigned>::iterator it;
        it = std::find(vec.begin(), vec.end(), item);
        if (it != vec.end()) return true;
        else return false;
    }
};

class SPHParticleSystem{
public:
    SPHParticleSystem(vector<Particle> &particles , double xmin, double xmax,
                      double ymin, double ymax,
                      double zmin, double zmax, Point3 center,double ixmin, double ixmax,
                      double iymin, double iymax,
                      double izmin, double izmax);
    void InitSPH();
    vector<Particle> GetParticles(){return list_particles;}
    unsigned GetCurrentFrame(){return frame;}
    bool savefile();

    double SPHStep();
    void SearchParticleNeighbors();
    //Heat Transfer
    void UpdateParticleTemperature(unsigned index, unsigned c);//c is 0 or 1 : 0=not increse T, 1=increse T
    void UpdateParticleEnergy();
    //Fluid Simulation
    void UpdateParticlePressure();
    void UpdateParticleAcceleration();
    //Rigid Body Simulation

    //Movtion
    void DetectingCollision(unsigned index);
    void ResolvingCollisions(Particle& p1,Particle& p2);
    void UpdateParticlePosition();
    void EnforceVelocityConstraints(unsigned index);
    void EnforceForceConstraints(unsigned index);

    unsigned GetIndex(double p1, double p2){
        if(p1 > p2){//-1
            return 0;
        }else if(p1 < p2){//1
            return 2;
        }else{
            return 1;//0
        }
    }

    double time_current;
    double bbox[6];
    vector<unsigned> listWaterparticle;
    vector<Particle> list_particles;
    vector<pair<Point3,Point3> > layer;//index=edge no.
    vector<unsigned> energy;
    Rigid_Body ice_obj;
    Liquid water;
protected:
    double time_delta;
    double time_end;
    bool change_const;
    Grid* p_grid;
    double h;
    unsigned frame;
    unsigned num_w;
    unsigned num_ice;
    //for SPH
    double reden;//rest density
    double vis;//viscosity
    double vis_ice;
    double k;//gas constant (stiffness)
    double hh;
    double gamma;
    double sigma;
    double zeta;
    double max_velocity;
    double max_velocity_sqrt;
    double max_acceleration;
    double max_acceleration_sqrt;
    double damp;
    Vector3 ambient_gravity;

    double T_room;

    //for Heat Tranfer
    double d;//voxel width
    double c_ice;//ice thermal diffusion
    double c_w;//water thermal diffusion
    double hc;//thermal conductivity
    double sb_cons;//Stefan Boltzmann constant
    double k_w;//interfacial tension coefficient for w-w particle
    double k_ice;//interfacial tension coefficient for w-ice particle
    double ch;//specific heat capacity of ice
    double LH;//Latent heat

    double uk;//friction coefficient
    double kn;//spring coefficient



    Vector3 Fice;

    double kernel_poly6;
    double kernel_poly6_diff;
    double kernel_poly6_2diff;
    double kernel_spiky_diff;
    double kernel_laplacian;
};

#endif // SPH_H
