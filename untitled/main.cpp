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
        count++;

    }
    sc.close();

    if(box == 0){//rectangle water at back
        for(int i= bbox[0]*10; i<bbox[1]*10; i++){
            for(int j= bbox[2]*10; j<bbox[2]*5; j++){
                for(int k= bbox[4]*10; k<bbox[4]*8.5; k++){
                    particles.push_back(Particle(i*0.1, j*0.1, k*0.1,r));
                }
            }
        }
    }else if(box == 1){//cube water at center
        for(int i= -10; i<10; i++){
            for(int j= -10; j<10; j++){
                for(int k= -10; k<10; k++){
                    particles.push_back(Particle(i*0.1, (j*0.1)+bbox[2]+1, k*0.1,r));
                }
            }
        }
    }
    return true;

}

int main(int argc, char *argv[])
{
    if(argc > 1){
        test_scene = argv[1];
        infname = atoi(argv[2]);
    }

    double radius = 0.04f;

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
