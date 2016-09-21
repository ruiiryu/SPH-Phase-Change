#include <Windows.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;
#include <GL/gl.h>
#include <GL/glut.h>
#include "gui.h"
#include "particle.h"
#include "SPH.h"

const char* test_scene;
int infname;
float bbox[6];
vector<Particle> particles;
bool simulation_pause = true;
Point3 center;
double frame_time = 0;
SPHParticleSystem* fluid;
int win_id;
double xmin = 1e10;     double ymin = 1e10;     double zmin = 1e10;
double xmax = -1e10;    double ymax = -1e10;    double zmax = -1e10;
void timer(int value){
    if(simulation_pause == false){
        frame_time = ((SPHParticleSystem*)fluid)->SPHStep();
        if(frame_time > 10000){
            simulation_pause = true;
        }
        glutPostRedisplay();
        glutTimerFunc(50, timer, value+1);

    }
    else{
        glutPostRedisplay();
        glutTimerFunc(50, timer, value);
    }
}

void InitGL(int argc, char** argv){

    extern Camera camera;
    camera.m_zoom=(bbox[4]+bbox[5])/8;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH);
    glutInitWindowSize(800.0f, 600.0f);

    win_id=glutCreateWindow("SPH Phase Transition Demo");
    glutDisplayFunc(display);
    glutMouseFunc(mousefunc);
    glutMotionFunc(mousemotion);
    glutPassiveMotionFunc(mousemotion2);
    glutKeyboardFunc(keyfunc);
    glutReshapeFunc(changesize);
    glutTimerFunc(50, timer, 0);

    set_lighting();
    glutMainLoop();
}



bool createParticle(vector<Particle> &particles, double r, const char* scene, int box){
    std::ifstream sc(scene, std::ios::in | std::ios::binary);
    int count=0;
    while (!sc.eof()) {
        sc>>bbox[count];
        cout<<"bbox" << bbox[count]<<std::endl;
        count++;

    }
    sc.close();
    count = 0;
	xmin = 10000.0;		ymin = 10000.0;		zmin = 10000.0;
	xmax = -10000.0;	ymax = -10000.0;	zmax = -10000.0;
    double r2 = r*2;
	double r4 = r*4;
	if(box == 0){//water drop
        for(double i= bbox[0]; i<bbox[1]; i=i+r2){//-4 4 = 8
            for(double j= bbox[2]; j<bbox[2]+0.5f; j=j+r2){//-4 -3.5 = 0.5
                for(double k= bbox[4]; k<bbox[5]; k=k+r2){//-4 4 = 8
                    if(count%2 == 0){
                        particles.push_back(Particle(i, j, k,r));
                    }else{
                        particles.push_back(Particle(i+r, j+r, k+r,r));
                    }
                    count++;
                }
            }
        }
		xmin = bbox[0];		ymin = bbox[2];			zmin = bbox[4];
		xmax = bbox[1];		ymax = bbox[2]+0.5f;	zmax = bbox[5];
		//create water droplet
		for(double i= -1.0f; i<1.0f; i=i+r2){//-1 1 = 2
            for(double j= 0.0f; j<bbox[3]/4.0f; j=j+r2){//0 1 = 1
                for(double k= -1.0f; k<1.0f; k=k+r2){//-1 1 = 2
                    if(count%2 == 0){
                        particles.push_back(Particle(i, j, k,r));
                    }else{
                        particles.push_back(Particle(i+r, j+r, k+r,r));
                    }
                    count++;
                }
            }
        }
    }else if(box == 1){//cube water at top
        for(double i= -1.0; i<1; i=i+r+r){
            for(double j= 3.0; j<bbox[3]; j=j+r+r){
                for(double k= -1.0; k<1; k=k+r+r){
                    if(count%2 == 0){
                        particles.push_back(Particle(i, j, k,r));
                    }else{
                        particles.push_back(Particle(i+r, j, k+r,r));
                    }
					if(k < zmin){
						zmin = k;
					}else if(k > zmax){
						zmax = k;
					}
                    count++;
                }
				if(j < ymin){
						ymin = j;
				}else if(j > ymax){
						ymax = j;
				}
            }
			if(i < xmin){
				xmin = i;
			}else if(i > xmax){
				xmax = i;
			}
        }
		xmax += r+r;	ymax += r+r;	zmax += r+r;
    }
	
    return true;

}

int main(int argc, char *argv[])
{
    if(argc > 1){
        test_scene = argv[1];
        infname = atoi(argv[2]);
    }

    double radius = 0.05f;

    if(createParticle(particles, radius, test_scene, infname)){
        cout<<"create particle success." << particles.size()<<std::endl;
    }else{
        cout<<"create particle fail."<<std::endl;
    }

    cout<<"stat SPH"<<std::endl;
    fluid = new SPHParticleSystem(particles, bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5], center
            ,xmin,xmax,ymin,ymax,zmin,zmax);
    int count = 0;
    InitGL(argc, argv);
    count++;

    return 0;
}
