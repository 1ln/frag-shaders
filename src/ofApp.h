#pragma once

#include "ofMain.h"
#include "ofxGui.h"

/*
#include "hash.h"
#include "noise.h"

#include "gridBoxes.h"
#include "concentricCircles.h"
#include "stackedCircles.h"
*/

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

    string screen_size; 

    ofPlanePrimitive plane;
    ofBoxPrimitive box;
    
    ofEasyCam cam;

    ofLight light;

    ofFbo db;
    ofFbo::Settings db_settings;

    ofShader shader;
    ofImage img;      

    bool info;
    bool unit_cube;
   
    ofVec2f mouse;

    ofFileDialogResult res;
    filesystem::path path; 

    string frag;
    string src_frag;
    string src_vert;


};
