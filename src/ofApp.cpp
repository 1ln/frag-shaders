#include "ofApp.h"

void ofApp::setup() {

shader.load("render.vert","render.frag");

w = ofGetWidth();
h = ofGetHeight();

}

void ofApp::draw() {

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());
shader.setUniform2f("u_res",w,h);

ofDrawRectangle(0,0,w,h);

shader.end();


}

void ofApp::update() {
}
