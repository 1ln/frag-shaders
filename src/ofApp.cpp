#include "ofApp.h"

void ofApp::setup() {

shader.load("render.vert","nanfloor.frag");

w = ofGetWidth();
h = ofGetHeight();

cam_orbit = true;

steps = 64;
eps = 0.0001;
trace_dist = 500.0;

cam.setPosition(glm::vec3(0.0,0.0,5.0));
cam.lookAt(glm::vec3(0.0));

mouse.x = 0;
mouse.y = 0;

}

void ofApp::draw() {

if(cam_orbit == true) {
cam.begin();
cam.end();
}

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());
shader.setUniform1i("u_steps",steps);
shader.setUniform1f("u_trace_dist",trace_dist);
shader.setUniform1f("u_eps",eps);
shader.setUniform2f("u_res",w,h);
shader.setUniform3f("u_cam_pos",glm::vec3(cam.getPosition()));
shader.setUniform3f("u_cam_target",glm::vec3(cam.getTarget().getPosition()));
shader.setUniform2f("u_mouse_pos",mouse.x,mouse.y);

ofDrawRectangle(0,0,w,h);

shader.end();


}

void ofApp::update() {

w = ofGetWidth();
h = ofGetHeight();

}

void ofApp::keyPressed(int key) {
    if(key == 's') {
        img.grabScreen(0,0,w,h);
        img.save("source.png");
    }
}

void ofApp::mousePressed(int x,int y,int button) {
}

void ofApp::mouseReleased(int x,int y,int button) {
}

void ofApp::mouseMoved(int x,int y) {
    mouse.x = x;
    mouse.y = y;
}
