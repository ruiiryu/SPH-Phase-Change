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
#include <GL/glext.h>
#include "SPH.h"
#include "vec3.h"
#include "gui.h"

//water properties from http://www.engineeringtoolbox.com/water-thermal-properties-d_162.html
//http://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
//and http://www1.lsbu.ac.uk/water/water_properties.html
//https://en.wikipedia.org/wiki/Speed_of_sound
bool saveframe = false;
unsigned particle_num;
char compare_particle[256];
vector<unsigned> less_nds;
SPHParticleSystem::SPHParticleSystem(vector<Particle> &particles , double xmin, double xmax,
                                     double ymin, double ymax,
                                     double zmin, double zmax, Point3 center,double ixmin, double ixmax,
                                     double iymin, double iymax,
                                     double izmin, double izmax){
    cout<<"Init SPH system"<<std::endl;
    time_current = 0.0f;
    time_end = 1.0f;
    time_delta = 0.01f;// (0.0005 seconds per frame)
    frame = unsigned(floorf(time_current/time_delta+0.5f));
    list_particles = particles;
    bbox[0] = xmin;     bbox[1] = xmax;
    bbox[2] = ymin;     bbox[3] = ymax;
    bbox[4] = zmin;     bbox[5] = zmax;

    particle_num = list_particles.size();
    num_ice = list_particles.size();
    num_w = 0;
	cout << ixmin << " " << ixmax << " : " << iymin << " " << iymax << " : " << izmin << " " << izmax << std::endl;
	volume = (ixmax-ixmin)*(iymax-iymin)*(izmax-izmin);

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

*/
    p_grid = NULL;
    less_nds.clear();
    InitSPH();
}

void SPHParticleSystem::InitSPH(){
    //
    if(infname == 0){
        h = 0.3f;
    }else{
        h = 0.2f;
    }
    hh = h*h;
	h2 = 2*h;
	double h4 = 4.0f*h;
    gamma = 1.0f;
    sigma = 0.0f;
    zeta = 1.0f;

    //rest density unit: kg/m^3 20C
    //reden=0.9983f;//g/cm^3
	reden=0.9983*h2*h2*h2;//g/cm^3
	//reden = 1000.0f;
    //vis = 1.0f;//viscosity kg/ms 10.0g/cms 20C
    vis = 0.005f;//g/cm(0.0005)s
	//vis = 10.0f;
    //k = 343.2f;//(stiffness) speed of sound unit: m/s 34320 cm/s 20C
	k = 17.16f;//cm/(0.0005)s
	//k = 34320.0f;

    damp = 0.6f;
    //ambient_gravity = Vector3(0.0f, -10.0f, 0.0f);//
	//ambient_gravity = Vector3(0.0f, -1000.0f, 0.0f);//
	ambient_gravity = Vector3(0.0f, -5.0f, 0.0f);


    kernel_poly6 = 315.0f/(64.0f*PI*powf(h, 9.0f));
    kernel_spiky_diff = -45.0f/(PI*powf(h, 6.0f));
    kernel_laplacian = 45.0f/(PI*powf(h, 6.0f));

	double mass = volume*reden/particle_num;
	cout << volume << " : " << mass << std::endl;
    
    for(unsigned i = 0; i<particle_num ; i++){
        Particle& p = list_particles[i];
        p.id = i;
        p.area = 6*4*p.radius*p.radius;
        p.mass = mass;
        //p.velocity = ambient_gravity*time_delta;
        //p.velocity_delta = ambient_gravity*time_delta;
    }

    sprintf(compare_particle,
            "TperFrame/SPH_PTScnce%dh%.1f-reden%.4f-vis%.3f-k%.2f-damp%.1f-mass%.3f.txt",
            infname,h,reden,vis,k,damp,mass);
	//sprintf(compare_particle,
    //        "TperFrame/SPH_PT-ParticleNum0-h%.1f-reden%.4f-vis%.3f-k%.2f-damp%.1f.txt",
    //        h,reden,vis,k,damp,c_ice,c_w,hc,k_w,k_ice,ch,LH);

}

bool SPHParticleSystem::savefile(){
    std::ofstream out;
    try{
        if(frame == 0){
            out.open(compare_particle, std::ofstream::out);
            out << "kernel_poly6: " << kernel_poly6 << std::endl;
            out << "kernel_spiky_diff: " << kernel_spiky_diff << std::endl;
            out << "kernel_laplacian: " << kernel_laplacian << std::endl;
        }else{
            out.open(compare_particle, std::ofstream::out | std::ofstream::app);
        }

    }catch(int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
        return false;
    }
    out << "Frame: " << frame << std::endl;
    for(unsigned i=0; i<particle_num; i++){
        Particle& curr = list_particles[i];//less_nds[i]
        out << "particle number: " << i << std::endl;
        out << "nbs-size: " << curr.nbs.size() << std::endl;
        out << "mass: " << curr.mass << std::endl;
        //out << "area: " << curr.area << std::endl;
        out << "density: " << curr.density << std::endl;
        out << "pressure: " << curr.pressure << std::endl;
        out << "fp: " << curr.force_pressure[0] << " " << curr.force_pressure[1] << " " << curr.force_pressure[2] << std::endl;
        out << "fv: " << curr.force_vis[0] << " " << curr.force_vis[1] << " " << curr.force_vis[2] << " " <<std::endl;
        out << "acc: " << curr.acc[0] << " " << curr.acc[1] << " " << curr.acc[2] << " " <<std::endl;
        out << "vel_delta: " << curr.velocity_delta[0] << " " << curr.velocity_delta[1] << " " << curr.velocity_delta[2] << " " <<std::endl;
        out << "pos: " << curr.pos[0] << " " << curr.pos[1] << " " << curr.pos[2] << std::endl;
        //out << "Temperature: " << curr.T[0] << std::endl;
        //out << "Heat: " << curr.heatvalue << std::endl;

    }
    out << "------------------------------------------------" << std::endl;
    out.close();
    return true;
}

double SPHParticleSystem::SPHStep(){
    clock_t start = clock();
    if(p_grid != NULL) delete p_grid;
    p_grid = new Grid(bbox, h, list_particles);
    SearchParticleNeighbors();
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
    less_nds.clear();
    clock_t end = clock();
    return double(end - start)/CLOCKS_PER_SEC;
}

//Search index of neighbor particle from paper "Analysis of the clustering Properties of Hilbert space-Filling Curve"
void SPHParticleSystem::SearchParticleNeighbors(){
    cout<<"SearchParticleNeighbors."<<std::endl;
    for(unsigned i=0; i<particle_num; i++){
        p_grid->Search(i, list_particles[i].nbs);
    }

}

//Use Equation from paper "Particle-based Fluid Simulation for Interaction Applications" 
void SPHParticleSystem::UpdateParticlePressure(){
    cout<<"UpdateParticlePressure."<<std::endl;
	//Equation 3. and kernel_poly6 is Smoothing kernel Poly6 in Equation 20.
    for(unsigned i=0; i<particle_num; i++){
        Particle& particle = list_particles[i];
        particle.density = 0.0f;
        for(unsigned j=0; j<particle.nbs.size(); j++){
            unsigned nb_id = particle.nbs[j].first;
            double distance = particle.nbs[j].second;
            double c = hh-distance*distance;
            particle.density += list_particles[nb_id].mass*kernel_poly6*(c*c*c);
        }
        if(particle.density <= 0.0f){
           particle.density = particle.mass;
        }
		//Equation 12.
        particle.pressure = k*(particle.density-reden);

    }
}
//Use Equation from paper "Particle-based Fluid Simulation for Interaction Applications"
void SPHParticleSystem::UpdateParticleAcceleration(){
    cout<<"UpdateParticleAcceleration."<<std::endl;
    for(unsigned i=0; i<particle_num; i++){
        Particle& particle = list_particles[i];
        particle.acc = Vector3();
        particle.force_vis = Vector3();
        particle.force_pressure = Vector3();
        for(unsigned j=0; j<particle.nbs.size(); j++){
            Particle& nb_particle = list_particles[particle.nbs[j].first];
            double distance = particle.nbs[j].second;
            //acceleration from pressure calculated by equation 10.
            double c1 = -0.5f*nb_particle.mass*(particle.pressure+nb_particle.pressure);
			//kernel_spiky_diff is gradient Smoothing kernel Spiky in Equation 21.
            double c2 = kernel_spiky_diff*(h-distance)*(h-distance);
            double c3 = 1.0f/nb_particle.density;
            particle.force_pressure += (particle.pos-nb_particle.pos).normalize()*c1*c2*c3;
			
            //acceleration from viscosity calculated by equation 14.
			//kernel_laplacian is laplacian Smoothing kernel viscosity in Equation 22.
            double c4 = vis*kernel_laplacian*(h-distance)*nb_particle.mass/(nb_particle.density);
            particle.force_vis += (nb_particle.velocity-particle.velocity)*c4;

        }
		//acceleration total calculated by equation 8.
        particle.acc += particle.force_pressure/particle.density;
        particle.acc += particle.force_vis/particle.density;
        particle.acc += ambient_gravity;

    }
}
//Use Equation from paper "Particle-based Fluid Simulation for Interaction Applications"
//In Section 3.6 Simulation
void SPHParticleSystem::UpdateParticlePosition(){
    cout<<"UpdateParticlePosition."<<std::endl;
	//Integration equation 8. use Leap-Frog scheme.
	//http://www.ch.embnet.org/MD_tutorial/pages/MD.Part1.html
    for(unsigned i=0; i<particle_num; i++){
        Particle& particle = list_particles[i];
		particle.pre_velocity = particle.velocity_delta;
		particle.pre_pos = particle.pos;
		
		particle.velocity_delta += particle.acc*time_delta;
		particle.pos += particle.velocity_delta*time_delta;
		particle.velocity = (particle.pre_velocity + particle.velocity_delta)*0.5f;
        
        wall_collide(i);

        

    }


}


void SPHParticleSystem::wall_collide(unsigned index){
    Particle& particle = list_particles[index];
    while(particle.pos[0]<bbox[0]
          || particle.pos[1]<bbox[2]
          || particle.pos[2]<bbox[4]
          || particle.pos[0]>bbox[1]
          || particle.pos[1]>bbox[3]
          || particle.pos[2]>bbox[5]){
        for(unsigned i=0; i<3; i++){
            unsigned ii=i+i;
            if(particle.pos[i]<bbox[ii]){
                double pd = particle.pos[i];
                double dist_diff = bbox[ii]-pd;
                double time_diff = dist_diff/fabsf(particle.velocity_delta[i]);
                particle.velocity_delta[i] = fabsf(particle.velocity_delta[i]);
                particle.velocity_delta = particle.velocity_delta * damp;
				//particle.pre_velocity[i] = -particle.pre_velocity[i];
                particle.pos[i] = bbox[ii]+fabsf(particle.velocity_delta[i])*time_diff;
            }
        }
        for(unsigned i=0; i<3; i++){
            unsigned ii=i+i+1;
            if(particle.pos[i]>bbox[ii]){
                double pd = particle.pos[i];
                double dist_diff = pd-bbox[ii];
                double time_diff = dist_diff/fabsf(particle.velocity_delta[i]);
                particle.velocity_delta[i] = -fabsf(particle.velocity_delta[i]);
                particle.velocity_delta = particle.velocity_delta * damp;
				//particle.pre_velocity[i] = -particle.pre_velocity[i];
                particle.pos[i] = bbox[ii]-fabsf(particle.velocity_delta[i])*time_diff;
            }
        }



    }


}


