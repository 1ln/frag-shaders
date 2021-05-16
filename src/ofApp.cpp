#include "ofApp.h"

void ofApp::setup() {

ofEnableDepthTest();

frag = "frag.frag"; 

path = filesystem::path("../../src");
shader.load(path/"render.vert",path/frag);

w = 895; 
h = 645;

ofSetWindowShape(w,h);

b.allocate(w,h);
b.getTexture().getTextureData().bFlipTexture = true;

plane.set(w,h);
plane.setPosition(w/2,h/2,0);
plane.setResolution(2,2);

mouse.x = 0;
mouse.y = 0;
 
}

void ofApp::draw() {

shader.begin();

shader.setUniform1f("time",ofGetElapsedTimeMillis()/1000.);
shader.setUniform1f("frame",ofGetLastFrameTime());
shader.setUniform2f("resolution",w,h);
shader.setUniform2f("mouse",mouse.x,mouse.y);

plane.set(w,h);
plane.setPosition(w/2,h/2,0);

plane.draw();

shader.end();


}

void ofApp::update() {

w = ofGetWidth();
h = ofGetHeight();

}

void ofApp::keyPressed(int key) {

    if(key == 'g') {
        img.grabScreen(0,0,w,h);
        img.save("../../frag.png");

    }

    if(key == 'f') {
      
        res = ofSystemLoadDialog("Load fragment shader",false,"../../src");
      
        if(res.bSuccess) {
            src_frag = res.getPath();
        }
        
        shader.load(path/"render.vert",src_frag);
    } 
}

void ofApp::mouseMoved(int x,int y) {

    mouse.x = x;
    mouse.y = y;

}
