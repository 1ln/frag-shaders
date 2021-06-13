#pragma once

#include "ofMain.h"
#include "ofxGui.h"

class ofApp : public ofBaseApp {

    public:

    void setup();
    void draw();
    void update();  

    void keyPressed(int key);

    int w;
    int h;

    int imgw;
    int imgh;

    ofxPanel gui;
 
    ofParameter<int> seed;
    ofParameter<bool> info;

    bool hide_panel;

    ofEasyCam cam;

    ofFbo fbo;
    ofFbo fbo1;    

    ofPixels px;
    ofImage img;

    ofParameter<string> screen_size; 
    ofParameter<string> fps;
    ofParameter<string> frame;

    ofShader render;
    ofShader buffer;

    ofFileDialogResult res;
    filesystem::path path; 

    string src_frag;
    string src_vert;


};
