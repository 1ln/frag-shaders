#pragma once

#include "ofMain.h"

class ofApp : public ofBaseApp {

    public:

    void setup();
    void draw();
    void update();  

    void keyPressed(int key);

    void mouseMoved(int x,int y);   

    int w;
    int h;
    ofFbo b;

    int imgw;
    int imgh;
    ofPixels px;
    ofImage img;
    string screen_size; 

    ofPlanePrimitive plane;

    ofShader shader;      

    ofVec2f mouse;

    ofFileDialogResult res;
    filesystem::path path; 

    string frag;
    string src_frag;
    string src_vert;


};
