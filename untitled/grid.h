#ifndef GRID_H
#define GRID_H

#include <vector>
using namespace std;

#include "particle.h"

class Grid{
public:
    Grid(double* boundingbox, double h, unsigned particle_num, Particle* p_particles);
    Grid(double* boundingbox, double h, const vector<Particle>& particles);
    ~Grid();

    //flag is true: point the same as the input point well be included in nbs, otherwise not.
    void Search(const Point3& point, vector<pair<unsigned, double> >& nbs, bool selfcounting = false);
    void Search(unsigned index, vector<pair<unsigned, double> >& nbs, bool selfcounting = false);
    void Search(unsigned index, vector<pair<unsigned,double> >& nbs, vector<unsigned> mesh, bool selfcounting = false);
private:
    void Locate(const Point3& pos, unsigned cell_index[3]);
    unsigned IndexToId(unsigned x, unsigned y, unsigned z);
    void GetNeighborCells(unsigned cell_index[3], vector<unsigned>& cells);


    vector<unsigned>* p_indices;
    Point3* p_points;
    unsigned points_num;
    double bbox[6];
    double cell_size, squared_cell_size;
    unsigned xdim, ydim, zdim, yzdim;
};

#endif // GRID_H
