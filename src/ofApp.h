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
    void windowResized(int w,int h); 

    int w;
    int h;
       
    ofEasyCam cam;

    ofFbo s;

    ofTexture ntex;
    ofTexture rtex;

    ofShader shader;
    ofImage img;      

    ofxPanel gui;

    ofParameter<float> dist;
    ofParameter<float> eps;  
    ofParameter<int> steps;
   
    ofParameter<glm::vec3> cam_pos;
    ofParameter<bool> cam_orbit;

    ofParameter<int> aa;
  
    ofParameter<int> shad_steps;
    ofParameter<float> shad_eps;
    ofParameter<float> shad_min;
    ofParameter<float> shad_max;
    ofParameter<float> shad_k;
    
    ofParameter<ofColor> dif;
    ofParameter<ofColor> amb;
    ofParameter<ofColor> spe;

    ofParameter<glm::vec3> light_pos;
    ofParameter<float> light_intensity;
    ofParameter<bool> light_display;
    ofParameter<bool> light_cam;

    ofParameter<glm::vec3> plane_offset;
    ofParameter<glm::vec3> plane_orient;
    ofParameter<bool> plane_display;

    ofParameter<string> screen_size;
    ofParameter<string> fps_counter;
   
    
    ofVec2f mouse;

};

