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
    void mouseScrolled(int x,int y,float scrollX,float scrollY);

    void windowResized(int w,int h); 

    void printInfo();

    int w;
    int h;
    ofFbo b;

    int imgw;
    int imgh;
    ofFbo ib;
    ofPixels px;
    ofImage img;
    string screen_size; 

    ofFbo d;

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
    bool mouse_le;
    bool mouse_ri;
    
    ofVec2f scroll;

    ofFileDialogResult res;
    filesystem::path path; 

    string frag;
    string src_frag;
    string src_vert;


};
