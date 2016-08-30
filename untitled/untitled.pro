#-------------------------------------------------
#
# Project created by QtCreator 2016-04-28T20:35:10
#
#-------------------------------------------------

QT       += core gui opengl

QT       -= gui

TARGET = untitled
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    particle.cpp \
    SPH.cpp \
    grid.cpp \
    vec3.cpp \
    gui.cpp \
    ArcBall.cpp

HEADERS += \
    particle.h \
    SPH.h \
    grid.h \
    vec3.h \
    gui.h \
    ArcBall.h

# this is the important part
unix|win32: LIBS += -lGLU
unix|win32: LIBS += -lGL
unix|win32: LIBS += -lglut
