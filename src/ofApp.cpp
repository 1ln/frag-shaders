#include "ofApp.h"

void ofApp::setup() {

ofDisableDepthTest();

path = filesystem::path("../../src");
render.load(path/"render.vert",path/"render.glsl");
buffer.load(path/"render.vert",path/"buffer.glsl");

w = ofGetWidth(); 
h = ofGetHeight();

ofFbo::Settings settings;
settings.width = w;
settings.height = h;

fbo.allocate(settings);
fbo1.allocate(w,h);
fbo1.getTexture().getTextureData().bFlipTexture = true;

cam.setPosition(glm::vec3(1.0));
cam.lookAt(glm::vec3(0.0));

hide_panel = true;
gui.setup(); 

gui.add(seed.set("seed",(int)(ofRandom(0,1)*1000000)));
gui.add(info.set("info",false));

key_up = false;
key_down = false;
key_right = false;
key_left = false;
key_x = false;
key_z = false;

}

void ofApp::draw() {

fbo.begin();
ofClear(0);
buffer.begin();
buffer.setUniform2f("resolution",w,h);
buffer.setUniform1i("seed",seed);
buffer.end();
fbo.end();

fbo1.begin();
render.begin();

cam.begin();
cam.end();

render.setUniform1f("time",ofGetElapsedTimeMillis()/1000.);
render.setUniform1f("frame",ofGetLastFrameTime());
render.setUniform2f("resolution",w,h);
render.setUniform1i("seed",seed);
render.setUniform3f("camPos",glm::vec3(cam.getPosition()));
render.setUniform1i("up",key_up);
render.setUniform1i("down",key_down);
render.setUniform1i("right",key_right);
render.setUniform1i("left",key_left);
render.setUniform1i("key_x",key_x);
render.setUniform1i("key_z",key_z);
render.setUniformTexture("tex",fbo.getTexture(0),1);



ofDrawRectangle(0,0,w,h);

render.end();
fbo1.end();
fbo1.draw(0,0,w,h);

if(!hide_panel) {
gui.draw();
} 

}

void ofApp::update() {

w = ofGetWidth();
h = ofGetHeight();

}

void ofApp::keyPressed(int key) {

    if(key == OF_KEY_UP) {
        key_up = true;
    }
 
    if(key == OF_KEY_DOWN) {
        key_down = true;
    }

    if(key == OF_KEY_RIGHT) {
        key_right = true;
    }

    if(key == OF_KEY_LEFT) {
        key_left = true;
    }

    if(key == 'x') {
        key_x = true;
    }

    if(key == 'z') {
        key_z = true;
    }

    if(key == 'g') {
        img.grabScreen(0,0,w,h);
        img.save("../../frag.png");

    }

    if(key == 'h') {
    hide_panel = !hide_panel;
    }    

    if(key == 'f') {
      
        res = ofSystemLoadDialog("Load fragment shader",false,"../../src");
      
        if(res.bSuccess) {
            src_frag = res.getPath();
        }
        
        render.load(path/"render.vert",src_frag);
    } 
}
