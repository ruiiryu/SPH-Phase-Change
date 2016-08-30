#include <cmath>
#include <cassert>
using namespace std;

#include "grid.h"

Grid::Grid(double* boundingbox, double h, unsigned particle_num, Particle* p_particles){
    for(unsigned i=0; i<6; i++){
        if(i%2==0) bbox[i]=boundingbox[i]-h;
        else bbox[i]=boundingbox[i]+h;
    }

    cell_size = h;
    squared_cell_size = h*h;

    xdim = (unsigned)ceilf((bbox[1]-bbox[0])/cell_size);
    ydim = (unsigned)ceilf((bbox[3]-bbox[2])/cell_size);
    zdim = (unsigned)ceilf((bbox[5]-bbox[4])/cell_size);
    yzdim = ydim*zdim;

    p_indices = new vector<unsigned>[xdim*ydim*zdim];
    p_points = new Point3[particle_num];
    unsigned cell_index[3];
    for(unsigned i=0; i<particle_num; i++){
        p_points[i] = p_particles[i].pos;
        Locate(p_points[i], cell_index);
        unsigned id = IndexToId(cell_index[0], cell_index[1], cell_index[2]);
        p_indices[id].push_back(i);
    }
    points_num = particle_num;
}

Grid::Grid(double* boundingbox, double h, const vector<Particle>& particles){
    for(unsigned i=0; i<6; i++){
        if(i%2==0) bbox[i]=boundingbox[i]-h;
        else bbox[i]=boundingbox[i]+h;
    }
    cell_size = h;
    squared_cell_size = h*h;

    xdim = (unsigned)ceilf((bbox[1]-bbox[0])/cell_size);
    ydim = (unsigned)ceilf((bbox[3]-bbox[2])/cell_size);
    zdim = (unsigned)ceilf((bbox[5]-bbox[4])/cell_size);
    yzdim = ydim*zdim;

    p_indices = new vector<unsigned>[xdim*ydim*zdim];
    p_points = new Point3[particles.size()];
    unsigned cell_index[3];
    for(unsigned i=0; i<particles.size(); i++){
        p_points[i] = particles[i].pos;
        Locate(p_points[i], cell_index);
        unsigned id = IndexToId(cell_index[0], cell_index[1], cell_index[2]);
        p_indices[id].push_back(i);
    }
    points_num = particles.size();
}

Grid::~Grid(){
    delete [] p_indices;
    delete [] p_points;
}

void Grid::Search(const Point3 &point, std::vector<pair<unsigned,double> >& nbs, bool selfcounting){
    nbs.clear();

    unsigned cell_index[3];
    Locate(point, cell_index);

    vector<unsigned> nb_cells;
    GetNeighborCells(cell_index, nb_cells);
    for(unsigned i=0; i<nb_cells.size(); i++){
        vector<unsigned>& nb_cell = p_indices[ nb_cells[i] ];
        for(unsigned j=0; j<nb_cell.size(); j++){
            unsigned nb_point_id = nb_cell[j];
            const Point3& nb_point = p_points[ nb_point_id ];
            double sqaured_dist = (point-nb_point).squared_length();
            if(!selfcounting && sqaured_dist<1e-9f) continue;

            if(sqaured_dist < squared_cell_size){
                nbs.push_back(make_pair(nb_point_id, sqrtf(sqaured_dist)));
            }
        }
    }
}

void Grid::Search(unsigned index, std::vector<pair<unsigned,double> >& nbs, bool selfcounting){
    nbs.clear();

    unsigned cell_index[3];
    const Point3& point =  p_points[index];
    Locate(point, cell_index);
    vector<unsigned> nb_cells;
    GetNeighborCells(cell_index, nb_cells);
    for(unsigned i=0; i<nb_cells.size(); i++){
        vector<unsigned>& nb_cell = p_indices[ nb_cells[i] ];
        for(unsigned j=0; j<nb_cell.size(); j++){
            unsigned nb_point_id = nb_cell[j];
            if(!selfcounting && nb_point_id==index) continue;
            double sqaured_dist = (point-p_points[ nb_point_id ]).squared_length();

            if(sqaured_dist < squared_cell_size){
                nbs.push_back(make_pair(nb_point_id, sqrtf(sqaured_dist)));
            }
        }
    }
}


void Grid::Locate(const Point3 &pos, unsigned cell_index[]){
    cell_index[0] = unsigned(floorf((pos[0]-bbox[0])/cell_size));
    cell_index[0] = min(cell_index[0], xdim-1);
    cell_index[1] = unsigned(floorf((pos[1]-bbox[2])/cell_size));
    cell_index[1] = min(cell_index[1], ydim-1);
    cell_index[2] = unsigned(floorf((pos[2]-bbox[4])/cell_size));
    cell_index[2] = min(cell_index[2], zdim-1);
}

unsigned Grid::IndexToId(unsigned x, unsigned y, unsigned z){
    assert(x<xdim && y<ydim && z<zdim);
    return x*yzdim+y*zdim+z;
}

void Grid::GetNeighborCells(unsigned cell_index[3], vector<unsigned>& cells){
    cells.clear();
    int x=cell_index[0], y=cell_index[1], z=cell_index[2];
    assert(x<xdim && y<ydim && z<zdim);
    for(int i=-1; i<2; i++){
        int xi=x+i;
        for(int j=-1; j<2; j++){
            int yj=y+j;
            for(int k=-1; k<2; k++){
                int zk=z+k;
                if(xi>=0 && xi<xdim && yj>=0 && yj<ydim && zk>=0 && zk<zdim)
                    cells.push_back(IndexToId(xi, yj, zk));
            }
        }
    }
}
