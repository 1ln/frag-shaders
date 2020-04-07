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
 
    ofPlanePrimitive plane;      
    ofEasyCam cam;

    ofFbo s;

    ofTexture ntex;
    ofTexture rtex;

    ofShader shader;
    ofImage img;      

    ofxPanel gui;

    ofParameterGroup render_grp;
    ofParameterGroup light_grp;
    ofParameterGroup reflect_grp;
    ofParameterGroup shad_grp;

    ofParameter<float> dist;
    ofParameter<float> eps;  
    ofParameter<int> steps;
    ofParameter<float> hash;   
    ofParameter<int> aa;
  
    ofParameter<int> reflect_steps;
    ofParameter<float> reflect_dist;

    ofParameter<int> shad_steps;
    ofParameter<float> shad_dist;
    ofParameter<float> shad_k;

    ofParameter<glm::vec3> col;
    ofParameter<glm::vec3> dif;
    ofParameter<glm::vec3> amb;
    ofParameter<glm::vec3> spe;
    ofParameter<glm::vec3> light_pos;
    ofParameter<float> light_intensity;

    ofParameter<string> screen_size;
    ofParameter<string> fps_counter;
   
    bool gui_hide;   
    ofVec2f mouse;

};
