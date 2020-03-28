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

    ofxPanel gui;

    ofxFloatSlider trace_dist;
    ofxFloatSlider eps;  
    ofxIntSlider steps;
   
    ofxIntSlider aa;
  
    ofxIntSlider shad_steps;
    ofxFloatSlider shad_eps;
    ofxFloatSlider shad_min;
    ofxFloatSlider shad_max;
    ofxFloatSlider shad_k;
    
    ofxColor dif;
    ofxColor amb;
    ofxColor spe;

    ofxVec3Slider light_pos;
    ofxFloatSlider light_intensity;
    ofxToggle light_display;

    ofxVec3Slider plane_offset;
    ofxVec3Slider plane_orient;
    ofxToggle plane_display;


    bool cam_orbit;

    ofVec2f mouse;

};

