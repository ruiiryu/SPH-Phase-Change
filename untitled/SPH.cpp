#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <time.h>
#include <complex>

#include <GL/gl.h>
#include <GL/glut.h>

#include "SPH.h"
#include "vec3.h"
#include "gui.h"

//water properties from http://www.engineeringtoolbox.com/water-thermal-properties-d_162.html
//http://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
//and http://www1.lsbu.ac.uk/water/water_properties.html
//https://wiki.blender.org/index.php/Dev:Ref/Release_Notes/2.66/Physics
bool saveframe = false;
unsigned particle_num;
char compare_particle[256];
unsigned num_comp = 3;
double f_ps=0;
double f_vis=0;
double acc_aGravity = 0;
SPHParticleSystem::SPHParticleSystem(vector<Particle> &particles , double xmin, double xmax,
                                     double ymin, double ymax,
                                     double zmin, double zmax, Point3 center,double ixmin, double ixmax,
                                     double iymin, double iymax,
                                     double izmin, double izmax){
    cout<<"Init SPH system"<<std::endl;
    time_current = 0.0f;
    time_end = 1.0f;
    time_delta = 0.01f;// (seconds per frame)
    frame = unsigned(floorf(time_current/time_delta+0.5f));
    list_particles = particles;
    bbox[0] = xmin;     bbox[1] = xmax;
    bbox[2] = ymin;     bbox[3] = ymax;
    bbox[4] = zmin;     bbox[5] = zmax;

    particle_num = list_particles.size();
    num_ice = list_particles.size();
    num_w = 0;

    /*       y
            |
            |
            |
    z-------\
             \
              \x


                 p7         p3


                     p8             p4
                 p5         p1


                     p6             p2

                    Natual      unit            Modify      unit
   rest density :   1000        kg/m**3
   viscosity    :   1.52        kg/ms
   k            :   3000        m/s
*/
    p_grid = NULL;
    InitSPH();
}

void SPHParticleSystem::InitSPH(){
    h = 0.3f;//unit: m
    //h = 0.1f;
    hh = h*h;
    gamma = 1.0f;
    sigma = 0.0f;
    zeta = 1.0f;
    reden=1000.0f;//rest density unit: kg/m^3
    //reden=125f;
    vis = 1.52f;//viscosity(0.89 centiPoise) kg/ms
    //vis_ice = 0.89;
    k = 1484.0f;//(stiffness) speed of sound unit: m/s //Individual Gas Constant - R 461.5 (J/kg K)

    max_acceleration = 300.0f;
    max_acceleration_sqrt = sqrtf(max_acceleration);
    max_velocity = max_acceleration*time_delta;
    max_velocity_sqrt = sqrtf(max_velocity);
    damp = 0.6f;
    ambient_gravity = Vector3(0.0f, -10.0f, 0.0f);

    T_room = 298.15f;

    //for Heat Tranfer
    d=0.1f;//voxel width
    c_ice=0.00118f;//ice thermal diffusion (mm^2/s) 0.001181976 m^2/s
    c_w=0.00014f;//water thermal diffusion (mm^2/s) 0.00013916 m^2/s
    hc=2.22f;//thermal conductivity of ice (W/mK)
    sb_cons=0.0000000567f;//Stefan Boltzmann constant
    k_w=0.5f;//interfacial tension coefficient for w-w particle (N/m)
    k_ice=0.069f;//interfacial tension coefficient for w-ice particle
    ch = 2.108f;//specific heat capacity of the ice (kJ/kgK)
    LH=33.4;//Latent Heat of Melting (kJ/kg)

    kn=0.05;//
    uk=0.02;//0.02-0.09 ice-ice

    kernel_poly6 = 315.0f/(64.0f*PI*powf(h, 9.0f));
    kernel_poly6_diff = -945.0f/(32.0f*PI*powf(h, 9.0f));
    kernel_poly6_2diff = 945.0f/(32.0f*PI*powf(h, 6.0f));
    kernel_spiky_diff = -45.0f/(PI*powf(h, 6.0f));
    kernel_laplacian = 45.0f/(PI*powf(h, 6.0f));

    for(unsigned i = 0; i<particle_num ; i++){
        Particle& p = list_particles[i];
        p.id = i;
        p.area = 6*4*p.radius*p.radius;
        //p.mass = p.density/particle_num*(p.radius*p.radius*p.radius*8);
    }

    sprintf(compare_particle,
            "TperFrame/SPH_PTh%.2f-reden%.2f-vis%.2f-k%.2f-damp%.2f-Cice%.5f-CW%.5f-hc%.2f-kw%.2f-kice%.3f-ch%.3f-LH%.2f.txt",
            h,reden,vis,k,damp,c_ice,c_w,hc,k_w,k_ice,ch,LH);

}

bool SPHParticleSystem::savefile(){
    std::ofstream out;
    try{
        if(frame == 0){
            out.open(compare_particle, std::ofstream::out);
        }else{
            out.open(compare_particle, std::ofstream::out | std::ofstream::app);
        }

    }catch(int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
        return false;
    }
    //for(unsigned i=0; i<particle_num; i++){
        Particle& curr = list_particles[10];
        out << "Frame: " << frame << std::endl;
        out << "particle number: " << 10 << std::endl;
        out << "pos: " << curr.pos[0] << " " << curr.pos[1] << " " << curr.pos[2] << std::endl;
        out << "nbs-size: " << curr.nbs.size() << std::endl;
        out << "mass: " << curr.mass << std::endl;
        out << "area: " << curr.area << std::endl;
        out << "density: " << curr.density << std::endl;
        out << "pressure: " << curr.pressure << std::endl;
        out << "acc: " << curr.acc[0] << " " << curr.acc[1] << " " << curr.acc[2] << std::endl;
        out << "vel: " << curr.velocity[0] << " " << curr.velocity[1] << " " << curr.velocity[2] << " " <<std::endl;
        out << "Temperature: " << curr.T[0] << std::endl;
        out << "Heat: " << curr.heatvalue << std::endl;
        out << "------------------------------------------------" << std::endl;

    //}
    out.close();
    return true;
}

double SPHParticleSystem::SPHStep(){
    clock_t start = clock();
    if(p_grid != NULL) delete p_grid;
    p_grid = new Grid(bbox, h, list_particles);
    SearchParticleNeighbors();
    //UpdateParticleEnergy();
    UpdateParticlePressure();

    UpdateParticleAcceleration();
    UpdateParticlePosition();

    //save file
    char buffer[256];
    sprintf(buffer, "video/ice_%05d.tga", frame);


    if(savefile()){
    }else cout<<"unsucessful save compare"<<std::endl;

    if(saveframe){
    // we will store the image data here
      unsigned char *pixels;
      // the thingy we use to write files
      FILE * shot;
      // we get the width/height of the screen into this array
      int screenStats[4];

      // get the width/height of the window
      glGetIntegerv(GL_VIEWPORT, screenStats);

      // generate an array large enough to hold the pixel data
      // (width*height*bytesPerPixel)
      pixels = new unsigned char[screenStats[2]*screenStats[3]*3];
      // read in the pixel data, TGA's pixels are BGR aligned
      glReadPixels(0, 0, screenStats[2], screenStats[3], GL_BGR,
                                       GL_UNSIGNED_BYTE, pixels);


      // open the file for writing. If unsucessful, return 1
      if((shot=fopen(buffer, "wb"))==NULL){
          cout<<"unsucessful image"<<std::endl;
      }

      // this is the tga header it must be in the beginning of
      // every (uncompressed) .tga
      unsigned char TGAheader[12]={0,0,2,0,0,0,0,0,0,0,0,0};
      // the header that is used to get the dimensions of the .tga
      // header[1]*256+header[0] - width
      // header[3]*256+header[2] - height
      // header[4] - bits per pixel
      // header[5] - ?
      unsigned char header[6]={((int)(screenStats[2]%256)),
                       ((int)(screenStats[2]/256)),
                       ((int)(screenStats[3]%256)),
               ((int)(screenStats[3]/256)),24,0};

      // write out the TGA header
      fwrite(TGAheader, sizeof(unsigned char), 12, shot);
      // write out the header
      fwrite(header, sizeof(unsigned char), 6, shot);
      // write the pixels
      fwrite(pixels, sizeof(unsigned char),
                     screenStats[2]*screenStats[3]*3, shot);

      // close the file
      fclose(shot);
      // free the memory
      delete [] pixels;
    }

    frame++;

    f_ps = 0;
    f_vis = 0;
    acc_aGravity = 0;
    clock_t end = clock();
    return double(end - start)/CLOCKS_PER_SEC;
}

void SPHParticleSystem::SearchParticleNeighbors(){
    cout<<"SearchParticleNeighbors."<<std::endl;
    for(unsigned i=0; i<particle_num; i++){
        p_grid->Search(i, list_particles[i].nbs);
    }
}

void SPHParticleSystem::UpdateParticleTemperature(unsigned index, unsigned c){
        Particle& particle = list_particles[index];
        if(c == 1){
            double cd = 0;
            if(particle.type == Particle::ICE){
                cd = c_ice;
            }else if(particle.type == Particle::WATER){
                cd = c_w;
            }else{
                cd = 0;
            }
            for(unsigned j=0; j<particle.nbs.size(); j++){
                Particle& nb_particle = list_particles[particle.nbs[j].first];
                double distance = particle.nbs[j].second;
                //particle-particle
                //T-delta = cd*massJ/densityJ*(Tj-Ti)*kernel_laplacian
                particle.T[1] += cd*nb_particle.mass/nb_particle.density*
                        (nb_particle.T[0]-particle.T[0])*kernel_laplacian*(h-distance);
            }
        }

        //particle-air
        //Qi = h(T air − T i )δA
        //find not touch air
        double Q;

        for(unsigned j=0; j<particle.nbs.size(); j++){
            Particle& nb_particle = list_particles[particle.nbs[j].first];
            if(particle.pos[0]>nb_particle.pos[0]){//-x
                particle.nbs_notAir[0]++;
            }else if(particle.pos[0]<nb_particle.pos[0]){//x
                particle.nbs_notAir[1]++;
            }
            if(particle.pos[1]>nb_particle.pos[1]){//-y
                particle.nbs_notAir[2]++;
            }else if(particle.pos[1]<nb_particle.pos[1]){//y
                particle.nbs_notAir[3]++;
            }
            if(particle.pos[2]>nb_particle.pos[2]){//-z
                particle.nbs_notAir[4]++;
            }else if(particle.pos[2]<nb_particle.pos[2]){//z
                particle.nbs_notAir[5]++;
            }
        }

        double count = 0;
        for(unsigned k=0; k<6; k++){
            if(particle.nbs_notAir[k]<45){
                count++;
            }
        }
        cout <<" particle: "<< index <<" air count: "<< count << std::endl;

        Q = hc*(T_room-(particle.T[0]+particle.T[1]))*(count/6.0f)*particle.area;
        if(c == 1){
            particle.T[1] += Q/(ch*particle.mass);
            particle.T[0] += particle.T[1];
        }
        particle.heatvalue += Q;

}

void SPHParticleSystem::UpdateParticleEnergy(){
    for(unsigned i=0; i<particle_num; i++){
        Particle& particle = list_particles[i];
        if(particle.T[0]>273.15f && particle.type == Particle::ICE && particle.heatvalue>LH){
            particle.type = Particle::WATER;
            listWaterparticle.push_back(i);
            particle.color[0] = 0.0f;
            particle.color[1] = 1.0f;
            particle.color[2] = 0.0f;
            UpdateParticleTemperature(i,1);
        }else if(particle.T[0]>=273.15f && particle.type == Particle::ICE){
            particle.color[0] = 0.0f;
            particle.color[1] = 1.0f;
            particle.color[2] = 1.0f;
            UpdateParticleTemperature(i,0);
        }else{
            UpdateParticleTemperature(i,1);
        }
    }

}

void SPHParticleSystem::UpdateParticlePressure(){
    cout<<"UpdateParticlePressure."<<std::endl;
    for(unsigned i=0; i<particle_num; i++){
        Particle& particle = list_particles[i];
        particle.density = 0.0f;
        for(unsigned j=0; j<particle.nbs.size(); j++){
            unsigned nb_id = particle.nbs[j].first;
            double distance = particle.nbs[j].second;
            double c = hh-distance*distance;
            particle.density += list_particles[nb_id].mass*kernel_poly6*(c*c*c);
        }
        particle.density = max(particle.density, reden);
        particle.pressure = k*(particle.density-reden);
        if(particle.pressure < 0)
            particle.pressure*=zeta;
    }
}

void SPHParticleSystem::UpdateParticleAcceleration(){
    cout<<"UpdateParticleAcceleration."<<std::endl;
    for(unsigned i=0; i<particle_num; i++){
        Particle& particle = list_particles[i];
        particle.acc = Vector3();
        for(unsigned j=0; j<particle.nbs.size(); j++){
            Particle& nb_particle = list_particles[particle.nbs[j].first];
            double distance = particle.nbs[j].second;
            //acceleration from pressure
            double c1 = -0.5f*nb_particle.mass*(particle.pressure+nb_particle.pressure);
            double c2 = kernel_spiky_diff*(h-distance)*(h-distance);
            double c3 = 1.0f/nb_particle.density;
            particle.acc += (particle.pos-nb_particle.pos).normalize()*(c1*c2*c3/particle.mass);

            //acceleration from viscosity
            double c4 = vis*kernel_laplacian*(h-distance)*nb_particle.mass/(nb_particle.density/**particle.density*/);
            particle.acc += (nb_particle.velocity-particle.velocity)*(c4/particle.mass);

        }

        particle.acc += ambient_gravity;

    }
}

void SPHParticleSystem::UpdateParticlePosition(){
    cout<<"UpdateParticlePosition."<<std::endl;
    for(unsigned i=0; i<particle_num; i++){
        Particle& particle = list_particles[i];

        Vector3 prev_velocity = particle.velocity_delta;
        particle.pre_pos = particle.pos;
        particle.velocity_delta += particle.acc*time_delta;
        if(particle.velocity_delta.length() > max_velocity)
            particle.velocity_delta = particle.velocity_delta.normalize() * max_velocity;

        particle.pos += particle.velocity_delta*time_delta;

        EnforceVelocityConstraints(i);

        particle.velocity = (prev_velocity + particle.velocity_delta)*0.5f;

    }


}

void SPHParticleSystem::EnforceVelocityConstraints(unsigned index){
    Particle& particle = list_particles[index];
    while(particle.pos[0]<bbox[0] || particle.pos[1]<bbox[2] || particle.pos[2]<bbox[4]
        || particle.pos[0]>bbox[1] || particle.pos[1]>bbox[3] || particle.pos[2]>bbox[5]){
        for(unsigned i=0; i<3; i++){
            unsigned ii=i+i;
            if(particle.pos[i]<bbox[ii]){
                double pd = particle.pos[i];
                double dist_diff = bbox[ii]-pd;
                double time_diff = dist_diff/fabsf(particle.velocity_delta[i]);
                particle.velocity_delta[i]=-particle.velocity_delta[i];
                particle.velocity_delta = particle.velocity_delta * damp;
                particle.pos[i] = bbox[ii]+fabsf(particle.velocity_delta[i])*time_diff;
            }
        }

        for(unsigned i=0; i<3; i++){
            unsigned ii=i+i+1;
            if(particle.pos[i]>bbox[ii]){
                double pd = particle.pos[i];
                double dist_diff = pd-bbox[ii];
                double time_diff = dist_diff/fabsf(particle.velocity_delta[i]);
                particle.velocity_delta[i]=-particle.velocity_delta[i];
                particle.velocity_delta = particle.velocity_delta * damp;
                particle.pos[i] = bbox[ii]-fabsf(particle.velocity_delta[i])*time_diff;
            }
        }
    }
}
