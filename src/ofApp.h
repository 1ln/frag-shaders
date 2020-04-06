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
    ofParameterGroup info_grp;
    ofParameterGroup sdf_grp;
    ofParameterGroup light_grp;
    ofParameterGroup reflect_grp;
    ofParameterGroup sine_grp;
    ofParameterGroup fnoise_grp;
    ofParameterGroup shad_grp;
    ofParameterGroup plane_grp;

    ofParameter<float> dist;
    ofParameter<float> eps;  
    ofParameter<int> steps;
    ofParameter<float> hash;   
    ofParameter<int> aa;
  
    ofParameter<int> reflect_steps;
    ofParameter<float> reflect_eps;
    ofParameter<float> reflect_min;
    ofParameter<float> reflect_max;

    ofParameter<int> shad_steps;
    ofParameter<float> shad_eps;
    ofParameter<float> shad_min;
    ofParameter<float> shad_max;
    ofParameter<float> shad_k;

    ofParameter<glm::vec3> col;
    ofParameter<glm::vec3> dif;
    ofParameter<glm::vec3> amb;
    ofParameter<glm::vec3> spe;
    ofParameter<glm::vec3> light_pos;
    ofParameter<float> light_intensity;
    ofParameter<bool> light_display;
    ofParameter<bool> light_cam;

    ofParameter<float> plane_offset;
    ofParameter<bool> plane_display;

    ofParameter<bool> sine_display;
    ofParameter<float> sine_height;   
    ofParameter<float> sine_offset;

    ofParameter<bool> fnoise_display;
    ofParameter<int> fnoise_octaves;
    ofParameter<float> fnoise_frequency;
    ofParameter<float> fnoise_offset;

    ofParameter<bool> sphere;
    ofParameter<bool> box;
    ofParameter<bool> torus;
    ofParameter<bool> octahedron;
    ofParameter<bool> prism;
 
    ofParameter<string> screen_size;
    ofParameter<string> fps_counter;
   
    bool gui_hide;   
    ofVec2f mouse;

};
