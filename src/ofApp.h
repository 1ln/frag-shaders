#pragma once

#include "ofMain.h"
#include "ofxGui.h"

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
   
    ofxVec3Slider cam_pos;
    ofxToggle cam_orbit;

    ofxIntSlider aa;
  
    ofxIntSlider shad_steps;
    ofxFloatSlider shad_eps;
    ofxFloatSlider shad_min;
    ofxFloatSlider shad_max;
    ofxFloatSlider shad_k;
    
    ofxColorSlider dif;
    ofxColorSlider amb;
    ofxColorSlider spe;

    ofxVec3Slider light_pos;
    ofxFloatSlider light_intensity;
    ofxToggle light_display;
    ofxToggle light_cam;

    ofxVec3Slider plane_offset;
    ofxVec3Slider plane_orient;
    ofxToggle plane_display;

    ofVec2f mouse;

};

