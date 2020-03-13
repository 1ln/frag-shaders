#pragma once

#include "ofMain.h"

class ofApp : public ofBaseApp {

    public:

    void setup();
    void draw();
    void update();  

    void keyPressed(int key);
    void mouseMoved(int x,int y);   
    void mousePressed(int x, int y,int button);
    void mouseReleased(int x,int y,int button);
 
    int w;
    int h;
       
    ofEasyCam cam;

    ofFbo buf;

    ofTexture ntex;
    ofTexture rtex;

    ofShader shader;
    ofImage img;      

    ofNode light;

    float trace_dist;
    float eps;
       
    int steps;
    bool cam_orbit;

    ofVec2f mouse;

};

