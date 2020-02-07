#pragma once

#include "ofMain.h"

class ofApp : public ofBaseApp {

    public:

    void setup();
    void draw();
    void update();

    ofShader shader;
    
    int w;
    int h;

};

