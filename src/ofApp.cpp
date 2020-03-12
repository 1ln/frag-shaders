#include "ofApp.h"

void ofApp::setup() {

shader.load("render.vert","render.frag");

cam_orbit = true;

steps = 64;
eps = 0.0001;
trace_dist = 500.0;

cam.setTarget(glm::vec3(0.0));

mouse.x = 0;
mouse.y = 0;

//light.setPosition(ofVec3f(10.0));

}

void ofApp::draw() {

cam.setPosition(0.0,0.0,-5.0);

if(cam_orbit == true) {
cam.begin();
cam.end();
}

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());
shader.setUniform1i("u_steps",steps);
shader.setUniform1f("u_trace_dist",trace_dist);
shader.setUniform1f("u_eps",eps);
shader.setUniform2f("u_res",ofGetWidth(),ofGetHeight());
shader.setUniform3f("u_cam_pos",cam.getPosition().x,cam.getPosition().y,cam.getPosition().z);
shader.setUniform2f("u_mouse_pos",mouse.x,mouse.y);
//shader.setUniform3fv("u_light_pos",ofVec3f(light.getPosition()));

ofDrawRectangle(0,0,ofGetWidth(),ofGetHeight());

shader.end();


}

void ofApp::update() {
}

void ofApp::keyPressed(int key) {
    if(key == 's') {
        img.grabScreen(0,0,ofGetWidth(),ofGetHeight());
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






