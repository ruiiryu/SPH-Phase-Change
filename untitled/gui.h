#ifndef GUI_H_
#define GUI_H_

#include <GL/gl.h>
#include <GL/glut.h>



struct Camera{
public:
    Camera();
    float m_zoom,m_xpos,m_ypos;
    int	  m_lastMouseX,m_lastMouseY;
};

extern float bbox[6];
extern int win_id;
void display();
void display_using_VBO();
void mousefunc(int button, int state, int x, int y);
void mousemotion(int x, int y);
void mousemotion2(int x, int y);
void keyfunc(unsigned char key, int x, int y);
void changesize(int width, int height);
void set_lighting();
void DrawString(float x, float y, float z, const char* string);
void DrawTitleAndScale(const char* title, const float* title_color, bool drawscale, float max_scale,
                       float min_scale, int color_scale_num, int str_scale_num);
void globalSetDrawScaleFunc(void (*pfuncDrawScale) (float max_scale, float min_scale, const float* title_color, int color_scale_num, int str_scale_num));
void DrawScaleAtRightCenter(float min_scale, float max_scale, const float* title_color, int color_scale_num, int str_scale_num);
void DrawScaleAtLeftCenter(float min_scale, float max_scale, const float* title_color, int color_scale_num, int str_scale_num);
void DrawScaleAtBottomCenter(float min_scale, float max_scale, const float* title_color, int color_scale_num, int str_scale_num);
void DrawScaleAtTopCenter(float min_scale, float max_scale, const float* title_color, int color_scale_num, int str_scale_num);
void DrawScaleAtVoid(float min_scale, float max_scale, const float* title_color, int color_scale_num, int str_scale_num);

void globalSetDrawTitleFunc(void (*pfuncDrawTitle) (const char* title, const float* title_color));
void DrawTitleAtBottomCenter(const char* title, const float* title_color);
void DrawTitleAtTopCenter(const char* title, const float* title_color);
void DrawTitleAtVoid(const char* title, const float* title_color);

void val_to_color(float val, float min, float max, float color[3]);



#endif
