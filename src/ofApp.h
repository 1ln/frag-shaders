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
    void windowResized(int w,int h); 

    void printInfo();

    int w;
    int h;

    int imgw;
    int imgh;
    ofFbo imgbuff;
    ofPixels px;
    ofImage img;

    string screen_size; 

    ofPlanePrimitive plane;

    ofBoxPrimitive box;
    bool unit_cube;
    
    ofEasyCam cam;

    ofLight light;

    ofFbo db;
    ofFbo::Settings db_settings;

    ofShader shader;      

    bool info;
   
    ofVec2f mouse;
    bool mouse_released_left;
    bool mouse_pressed_left;
    
    ofFileDialogResult res;
    filesystem::path path; 

    string frag;
    string src_frag;
    string src_vert;


};
