/*
 * gui_util.cpp
 *
 *  Created on: Aug 4, 2010
 *      Author: leiwang
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
using namespace std;


#include "ArcBall.h"
#include "SPH.h"
#include "gui.h"

Matrix4fT   Transform   = {  1.0f,  0.0f,  0.0f,  0.0f,				// NEW: Final Transform
0.0f,  1.0f,  0.0f,  0.0f,
0.0f,  0.0f,  1.0f,  0.0f,
0.0f,  0.0f,  0.0f,  1.0f };

Matrix3fT   LastRot     = {  1.0f,  0.0f,  0.0f,					// NEW: Last Rotation
0.0f,  1.0f,  0.0f,
0.0f,  0.0f,  1.0f };

Matrix3fT   ThisRot     = {  1.0f,  0.0f,  0.0f,					// NEW: This Rotation
0.0f,  1.0f,  0.0f,
0.0f,  0.0f,  1.0f };
ArcBall_t ArcBall(800.0f, 600.0f);				                // NEW: ArcBall Instance
Point2fT MousePt;												// NEW: Current Mouse Point
bool        isClicked  = false;										// NEW: Clicking The Mouse?
bool        isRClicked = false;										// NEW: Clicking The Right Mouse Button?
bool        isMClicked=false;
//bool        isDragging = false;					                    // NEW: Dragging The Mouse?
int show_axis=1;
extern float bbox[6];
extern float pcbox[6];//edit
///////////////////////////////////////////

Camera::Camera(){
    m_xpos = 0.0f;
    m_ypos = 0.0f;
    m_zoom = 1.0f;
    m_lastMouseX=m_lastMouseY=-1;
}

Camera camera;

void DrawString(float x, float y, float z, const char* string)
{
    glRasterPos3f(x, y, z);
    const char *s;
    for(s = string; *s; s++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *s);
}

void DrawTitleAtBottomCenter(const char* title, const float* title_color)
{
    if(title_color!=NULL)
        glColor3fv(title_color);
    DrawString(15, 2.5, 0, title);
}

void DrawTitleAtTopLeft(const char* title, const float* title_color)
{
    if(title_color!=NULL)
        glColor3fv(title_color);
    DrawString(2.5, 35, 0, title);
}

void DrawScaleAtRightCenter(float min_scale, float max_scale, const float* title_color, int color_scale_num, int str_scale_num)
{
    int i=0;
    const float scale_left=34;
    const float scale_right=35;
    const float scale_top=32;
    const float scale_bottom=8;
    float scale_y=scale_top;
    int scale_num=color_scale_num;
    float scale_dec=(scale_top-scale_bottom)/scale_num;
    float s;
    float color[4];

    glDisable(GL_CULL_FACE);
    glBegin(GL_QUADS);
    for(i=0;i<scale_num;i++){
        val_to_color(scale_y-scale_dec/2, scale_bottom, scale_top, color);
        glColor3fv(color);
        glVertex2f(scale_left, scale_y);
        glVertex2f(scale_right, scale_y);
        scale_y-=scale_dec;
        glVertex2f(scale_right, scale_y);
        glVertex2f(scale_left, scale_y);
    }
    glEnd();
    glEnable(GL_CULL_FACE);

    scale_y=scale_top;
    scale_num=str_scale_num;
    scale_dec=(scale_top-scale_bottom)/(scale_num+1);
    char buffer[128];
    if(title_color!=NULL)
        glColor3fv(title_color);
    scale_y=scale_top-scale_dec;
    for(i=0;i<scale_num;i++){
        s=min_scale+(max_scale-min_scale)*(scale_y-scale_bottom)/(scale_top-scale_bottom);
        sprintf(buffer, "%4.3e", s);
        DrawString(scale_right+0.5, scale_y, 0, buffer);
        scale_y-=scale_dec;
    }
}

void draw_bbox(){
    glBegin(GL_LINES);
    glColor3f(1.0f, 1.0f, 1.0f);
    //1
    glVertex3f(bbox[0], bbox[2], bbox[4]);
    glVertex3f(bbox[1], bbox[2], bbox[4]);
    //2
    glVertex3f(bbox[0], bbox[2], bbox[4]);
    glVertex3f(bbox[0], bbox[3], bbox[4]);
    //3
    glVertex3f(bbox[0], bbox[2], bbox[4]);
    glVertex3f(bbox[0], bbox[2], bbox[5]);
    //4
    glVertex3f(bbox[0], bbox[3], bbox[4]);
    glVertex3f(bbox[0], bbox[3], bbox[5]);
    //5
    glVertex3f(bbox[0], bbox[3], bbox[4]);
    glVertex3f(bbox[1], bbox[3], bbox[4]);
    //6
    glVertex3f(bbox[1], bbox[3], bbox[4]);
    glVertex3f(bbox[1], bbox[3], bbox[5]);
    //7
    glVertex3f(bbox[1], bbox[3], bbox[5]);
    glVertex3f(bbox[1], bbox[2], bbox[5]);
    //8
    glVertex3f(bbox[1], bbox[3], bbox[5]);
    glVertex3f(bbox[0], bbox[3], bbox[5]);
    //9
    glVertex3f(bbox[0], bbox[3], bbox[5]);
    glVertex3f(bbox[0], bbox[2], bbox[5]);
    //10
    glVertex3f(bbox[0], bbox[2], bbox[5]);
    glVertex3f(bbox[1], bbox[2], bbox[5]);
    //11
    glVertex3f(bbox[1], bbox[2], bbox[4]);
    glVertex3f(bbox[1], bbox[3], bbox[4]);
    //12
    glVertex3f(bbox[1], bbox[2], bbox[4]);
    glVertex3f(bbox[1], bbox[2], bbox[5]);

    glEnd();
}
/*
void draw_pcbox(){
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    //1
    glVertex3f(pcbox[0], pcbox[2], pcbox[4]);
    glVertex3f(pcbox[1], pcbox[2], pcbox[4]);
    //2
    glVertex3f(pcbox[0], pcbox[2], pcbox[4]);
    glVertex3f(pcbox[0], pcbox[3], pcbox[4]);
    //3
    glVertex3f(pcbox[0], pcbox[2], pcbox[4]);
    glVertex3f(pcbox[0], pcbox[2], pcbox[5]);
    //4
    glVertex3f(pcbox[0], pcbox[3], pcbox[4]);
    glVertex3f(pcbox[0], pcbox[3], pcbox[5]);
    //5
    glVertex3f(pcbox[0], pcbox[3], pcbox[4]);
    glVertex3f(pcbox[1], pcbox[3], pcbox[4]);
    //6
    glVertex3f(pcbox[1], pcbox[3], pcbox[4]);
    glVertex3f(pcbox[1], pcbox[3], pcbox[5]);
    //7
    glVertex3f(pcbox[1], pcbox[3], pcbox[5]);
    glVertex3f(pcbox[1], pcbox[2], pcbox[5]);
    //8
    glVertex3f(pcbox[1], pcbox[3], pcbox[5]);
    glVertex3f(pcbox[0], pcbox[3], pcbox[5]);
    //9
    glVertex3f(pcbox[0], pcbox[3], pcbox[5]);
    glVertex3f(pcbox[0], pcbox[2], pcbox[5]);
    //10
    glVertex3f(pcbox[0], pcbox[2], pcbox[5]);
    glVertex3f(pcbox[1], pcbox[2], pcbox[5]);
    //11
    glVertex3f(pcbox[1], pcbox[2], pcbox[4]);
    glVertex3f(pcbox[1], pcbox[3], pcbox[4]);
    //12
    glVertex3f(pcbox[1], pcbox[2], pcbox[4]);
    glVertex3f(pcbox[1], pcbox[2], pcbox[5]);

    glEnd();
}
*/
void DrawTitleAndScale(const char* title//, const char* titlepos
                       ,const float* title_color, bool drawscale, float max_scale,
                                             float min_scale, int color_scale_num, int str_scale_num)
{
    glMatrixMode(GL_PROJECTION); //change to Ortho view
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, 40, 0, 40);

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();

    if(title_color!=NULL)
        glColor3fv(title_color);

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);

    DrawTitleAtBottomCenter(title, title_color);
    //DrawTitleAtTopLeft(titlepos, title_color);

    if(drawscale==true){
       DrawScaleAtRightCenter(min_scale, max_scale, title_color, color_scale_num, str_scale_num);
    }

    glPopAttrib();

    glEnable(GL_LIGHTING);
    glPopMatrix();

    //pop GL_PROJECTION
    glMatrixMode(GL_PROJECTION); //change to Pers view
    glPopMatrix();

    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();
}

void draw_obj(){
    extern SPHParticleSystem* fluid;

    const vector<Particle>& particles = fluid->GetParticles();
    unsigned particle_num = particles.size();
    bbox[0] = fluid->bbox[0];   bbox[1] = fluid->bbox[1];
    bbox[2] = fluid->bbox[2];   bbox[3] = fluid->bbox[3];
    bbox[4] = fluid->bbox[4];   bbox[5] = fluid->bbox[5];

    glDisable(GL_LIGHTING);
    draw_bbox();
    glEnable(GL_LIGHTING);


    //pcbox[0] = ; pcbox[1] = ; pcbox[2] = ;
    //pcbox[3] = ; pcbox[4] = ; pcbox[5] = ;

    //draw particle as a sphere
    if(particle_num < 10000){
        GLUquadric *pQuad=gluNewQuadric();

        for(unsigned i=0; i<particle_num; i++){
            glPushMatrix();
            glColor3f(particles[i].color[0], particles[i].color[1], particles[i].color[2]);
            glTranslatef(particles[i].pos[0], particles[i].pos[1], particles[i].pos[2]);
            gluSphere(pQuad, particles[i].radius, 8, 8);
            glPopMatrix();
        }
    }else{
        //draw particle as dot
        glDisable(GL_LIGHTING);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        for(unsigned i=0; i<particle_num; i++){
            glColor3f(particles[i].color[0], particles[i].color[1], particles[i].color[2]);
            glVertex3fv(particles[i].pos.v);
        }
        glEnd();
        glEnable(GL_LIGHTING);
    }


    vector<Point3> isosurface;

    glColor3f(0.5f, 0.5f, 0.9f);
    glDisable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for(unsigned i=0; i<isosurface.size(); i+=3){
        glVertex3f(isosurface[i].x(), isosurface[i].y(), isosurface[i].z());
        glVertex3f(isosurface[i+1].x(), isosurface[i+1].y(), isosurface[i+1].z());
        glVertex3f(isosurface[i+2].x(), isosurface[i+2].y(), isosurface[i+2].z());
    }
    glEnd();
    glEnable(GL_LIGHTING);


    float title_color[3]={0,0,0};
    char title[512];
    extern double frame_time;
    sprintf(title, "SPH  - %s. Frame : %d at rate %.1f fps", "Phase Change", fluid->GetCurrentFrame(), 1.0f/frame_time);

    DrawTitleAndScale(title//,titlepos
                      , title_color, false, Particle::max_temperature, Particle::min_temperature, 64, 4);
}

void draw_axes()
{
    glDisable(GL_LIGHTING);

    glMatrixMode(GL_PROJECTION); //change to Ortho view
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0,20,0,20);

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();

    glTranslatef(1.1f,1,0);
    glMultMatrixf(Transform.M);

    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex3f(0,0,0);
    glVertex3f(1,0,0);
    glColor3f(0,1,0);
    glVertex3f(0,0,0);
    glVertex3f(0,1,0);
    glColor3f(0,0,1);
    glVertex3f(0,0,0);
    glVertex3f(0,0,1);
    glEnd();

    //SetFont(GLUT_BITMAP_HELVETICA_12);
    glColor3f(0.0, 0.0, 1.0);
    DrawString(0, 0, 1, "Z");
    glColor3f(0.0, 1.0, 0.0);
    DrawString(0, 1, 0, "Y");
    glColor3f(1.0, 0.0, 0.0);
    DrawString(1, 0, 0, "X");

    glPopMatrix();

    //pop GL_PROJECTION
    glMatrixMode(GL_PROJECTION); //change to Pers view
    glPopMatrix();

    glMatrixMode( GL_MODELVIEW );
    glEnable(GL_LIGHTING);
}

void set_lighting(){

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glEnable(GL_NORMALIZE);

    glEnable(GL_COLOR_MATERIAL);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    //glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    GLfloat light_ambient[] = {.3,.3,.3,1};
    GLfloat light_diffuse[] = {.3,.3,.3,1};
    GLfloat light_specular[] = {.8,.8,.8,1};
    GLfloat light_pos[] = {0,15,0,1};

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);

    glEnable(GL_LIGHT1);
    GLfloat light_diffuse2[] = {.2,.2,.2,1};
    GLfloat light_specular2[] = {.5,.5,.5,1};
    GLfloat light_pos2[] = {0,15,-2,0};
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse2);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular2);
    glLightfv(GL_LIGHT1, GL_POSITION, light_pos2);

    glEnable(GL_LIGHT2);
    GLfloat light_diffuse3[] = {.2,.2,.2,1};
    GLfloat light_specular3[] = {.5,.5,.5,1};
    GLfloat light_pos3[] = {15,2,0,0};
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse2);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular2);
    glLightfv(GL_LIGHT2, GL_POSITION, light_pos3);

    GLfloat mat_specular[] = {.3, .3, .3, 1};
    GLfloat mat_shine[] = {50};

    glMaterialfv(GL_FRONT/*_AND_BACK*/, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT/*_AND_BACK*/, GL_SHININESS, mat_shine);
}

void display()
{
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.8,0.8,0.8,1);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    glTranslatef( camera.m_xpos, -camera.m_ypos, camera.m_zoom );
    glMultMatrixf(Transform.M);
    glTranslatef(-(bbox[0]+bbox[1])/2, -(bbox[2]+bbox[3])/2, (bbox[4]+bbox[5])/2);

    if(show_axis){
        glDisable(GL_LIGHTING);
        draw_axes();
        glEnable(GL_LIGHTING);
    }

    glShadeModel(GL_SMOOTH);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(2,2);
    draw_obj();

    glutSwapBuffers();
}


void mousefunc(int button, int state, int x, int y)
{
    if(state == GLUT_UP){
        if (button == GLUT_LEFT_BUTTON){
            isClicked   = false;
        }
        else if (button == GLUT_RIGHT_BUTTON){
            isRClicked  = false;
        }
        else if (button == GLUT_MIDDLE_BUTTON){
            isMClicked=false;
        }
    }
    else if(state == GLUT_DOWN){
        if (button == GLUT_LEFT_BUTTON){
            isClicked   = true;
            LastRot = ThisRot;										// Set Last Static Rotation To Last Dynamic One
            ArcBall.click(&MousePt);
        }
        else if (button == GLUT_RIGHT_BUTTON){
            isRClicked  = true;
            camera.m_lastMouseY = y;
        }
        else if (button == GLUT_MIDDLE_BUTTON){
            isMClicked=true;
            camera.m_lastMouseX = x;
            camera.m_lastMouseY = y;
        }
    }
}

/////////////////////////////////////////////////////////

void mousemotion(int x, int y)
{
    MousePt.s.X = (GLfloat)x;
    MousePt.s.Y = (GLfloat)y;

    if(isClicked){
        Quat4fT ThisQuat;

        ArcBall.drag(&MousePt, &ThisQuat);								// Update End Vector And Get Rotation As Quaternion
        Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);		// Convert Quaternion Into Matrix3fT
        Matrix3fMulMatrix3f(&ThisRot, &LastRot);				// Accumulate Last Rotation Into This One
        Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);	// Set Our Final Transform's Rotation From This One

        // Left mouse button is being
        // pressed. Rotate the camera.
        if ( camera.m_lastMouseX != -1 )
        {

            //camera.m_yrot += y - camera.m_lastMouseY;
            //camera.m_xrot += x - camera.m_lastMouseX;
            // Redraw the viewport.
            glutPostRedisplay();
        }
        camera.m_lastMouseX = x;
        camera.m_lastMouseY = y;
    }
    else if(isMClicked){
        // Middle mouse button is being
        // pressed. Pan the camera.
        if ( camera.m_lastMouseX != -1 )
        {
            camera.m_xpos += (x - camera.m_lastMouseX) * (bbox[5] - bbox[4])/500;
            camera.m_ypos += (y - camera.m_lastMouseY) * (bbox[5] - bbox[4])/500;
            // Redraw the viewport.
            glutPostRedisplay();
        }
        camera.m_lastMouseX = x;
        camera.m_lastMouseY = y;
    }
    else if(isRClicked){
        // Right mouse button is being
        // pressed. Zoom the camera.
        if ( camera.m_lastMouseY != -1 )
        {
            camera.m_zoom += (bbox[5] - bbox[4])/1000*(y - camera.m_lastMouseY);
            // Redraw the viewport.
            glutPostRedisplay();
        }
        camera.m_lastMouseY = y;
    }
}

void mousemotion2(int x, int y)
{
    MousePt.s.X = (GLfloat)x;
    MousePt.s.Y = (GLfloat)y;

    camera.m_lastMouseX = -1;
    camera.m_lastMouseY = -1;
}

void Quit(){
    glutDestroyWindow ( win_id );
    exit (0);
}

///////////////////////////////////////////////////////////////


void keyfunc(unsigned char key, int x, int y)
{
    extern bool simulation_pause;
    extern bool saveframe;
    switch(key) {
        case 'n': show_axis=1-show_axis;
            break;
        case 'p': simulation_pause=!simulation_pause;
            break;
        case 's': saveframe=!saveframe;
            break;
        case 'q': Quit();
            break;
    }
    glutPostRedisplay();
}
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
//	When changing the window size, some work has to be done.
//	Say the viewport, and trackball parameters
void changesize(int width, int height)
{
    glViewport (0, 0, (GLsizei)(width), (GLsizei)(height));				// Reset The Current Viewport

    glMatrixMode (GL_PROJECTION);										// Select The Projection Matrix
    glLoadIdentity ();													// Reset The Projection Matrix
    gluPerspective (45.0f, (GLfloat)(width)/(GLfloat)(height),			// Calculate The Aspect Ratio Of The Window
                    0.1f, 5000.0f);
                                                    // Reset The Modelview Matrix
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    ArcBall.setBounds((GLfloat)width, (GLfloat)height);                 //*NEW* Update mouse bounds for arcball
    glutPostRedisplay();
}

void val_to_color(float val, float min, float max, float color[3]){
    if(max<=min){
        color[0]=color[2]=0.0f;
        color[1]=1.0f;
    }
    float s=(val-min)/(max-min);
    if(s<0 || s>1) color[0]=color[1]=color[2]=0.0f;
    if(s<=0.25f){
        color[0]=0.0f;
        color[1]=s*4;
        color[2]=1.0f;
    }
    else if(s<=0.5f){
        color[0]=0.0f;
        color[1]=1.0f;
        color[2]=1.0f-(s-0.25f)*4;
    }
    else if(s<=0.75f){
        color[0]=(s-0.5f)*4;
        color[1]=1.0f;
        color[2]=0.0f;
    }
    else{
        color[0]=1.0f;
        color[1]=1.0f-(s-0.75f)*4;
        color[2]=0.0f;
    }
}
