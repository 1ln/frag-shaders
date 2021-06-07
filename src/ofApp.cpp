#include "ofApp.h"

void ofApp::setup() {

ofDisableDepthTest();

frag = "frag.frag"; 

path = filesystem::path("../../src");
shader.load(path/"render.vert",path/frag);
buffer.load(path/"render.vert",path/buffer.glsl);

w = ofGetWidth(); 
h = ofGetHeight();

cam.setPosition(glm::vec3(1.0));
cam.setLookAt(glm::vec3(0.0));



seed = ofRandom(0,1);
 
}

void ofApp::draw() {

shader.begin();


shader.setUniform1f("time",ofGetElapsedTimeMillis()/1000.);
shader.setUniform1f("frame",ofGetLastFrameTime());
shader.setUniform2f("resolution",w,h);
shader.setUniform1f("seed",seed);
shader.setUniform3f("camPos",glm::vec3(cam.getPosition()));

ofDrawRectangle(0,0,w,h);

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
