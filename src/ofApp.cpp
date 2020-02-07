#include "ofApp.h"

void ofApp::setup() {

shader.load("render.vert","render.frag");

w = 512;
h = 512;

}

void ofApp::draw() {

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimef());
shader.setUniform2f("u_resolution",w,h);

ofDrawRectangle(0,0,w,h);

shader.end();


}

void ofApp::update() {
}
