#include "ofApp.h"

void ofApp::setup() {

ofDisableDepthTest();

path = filesystem::path("../../src");
render.load(path/"render.vert",path/render.glsl);
buffer.load(path/"render.vert",path/buffer.glsl);

w = ofGetWidth(); 
h = ofGetHeight();

ofFbo::Settings settings;
settings.width = w;
settings.height = h;

fbo.allocate(settings);
fbo1.allocate(w,h);
fbo1.getTexture().getTextureData().bFlipTexture = true;

cam.setPosition(glm::vec3(1.0));
cam.setLookAt(glm::vec3(0.0));

seed = ofRandom(0,1);
 
}

void ofApp::draw() {

fbo.begin();
ofClear();
buffer.begin();
buffer.setUniform2f("resolution",w,h);
buffer.setUniform1f("seed",seed);
buffer.end();
fbo.end();

fbo1.begin();
render.begin();

cam.begin();
render.setUniform1f("time",ofGetElapsedTimeMillis()/1000.);
render.setUniform1f("frame",ofGetLastFrameTime());
render.setUniform2f("resolution",w,h);
render.setUniform1f("seed",seed);
render.setUniform3f("camPos",glm::vec3(cam.getPosition()));
render.setUniformTexture("tex",fbo.getTexture(0),1);

ofDrawRectangle(0,0,w,h);
cam.end();


render.end();
fbo1.end();
fbo1.draw(0,0,w,h);

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
        
        render.load(path/"render.vert",src_frag);
    } 
}
