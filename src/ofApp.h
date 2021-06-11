#pragma once

#include "ofMain.h"

class ofApp : public ofBaseApp {

    public:

    void setup();
    void draw();
    void update();  

    void keyPressed(int key);

    int w;
    int h;

    int imgw;
    int imgh;

    float seed;

    ofEasyCam cam;

    ofFbo fbo;
    ofFbo fbo1;    

    ofPixels px;
    ofImage img;
    string screen_size; 

    ofShader render;
    ofShader buffer;

    ofFileDialogResult res;
    filesystem::path path; 

    string src_frag;
    string src_vert;


};
