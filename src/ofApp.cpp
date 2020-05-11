#include "ofApp.h"

void ofApp::setup() {

ofEnableDepthTest();

frag = "frag.frag"; 

path = filesystem::path("../../src");
shader.load(path/"render.vert",path/frag);

w = 512; 
h = 512;

ofSetWindowShape(w,h);

db_settings.width = w;
db_settings.height = h;
db_settings.internalformat = GL_RGBA; 
db_settings.useDepth = true;
db_settings.depthStencilAsTexture = true; 

db.allocate(db_settings);

db.begin();
ofClear(0);
db.end();

plane.set(w,h);
plane.setPosition(w/2,h/2,0);
plane.setResolution(2,2);

cam.setPosition(glm::vec3(1.0,2.0,5.0));
cam.lookAt(glm::vec3(0.0));
cam.setNearClip(.01);

light.setup();
light.setPosition(glm::vec3(0,5,10));

info = false;
unit_cube = false;

mouse.x = 0;
mouse.y = 0;

mouse_released_left = false;
mouse_pressed_left = false;
 
}

void ofApp::draw() {

shader.begin();
cam.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());

shader.setUniform2f("u_res",w,h);

shader.setUniform3f("u_cam_pos",glm::vec3(cam.getPosition()));
shader.setUniform3f("u_cam_tar",glm::vec3(cam.getTarget().getPosition())); 

shader.setUniform2f("u_mouse_pos",mouse.x,mouse.y);

shader.setUniform1f("u_mouse_released_left",mouse_released_left); 
shader.setUniform1f("u_mouse_pressed_left",mouse_pressed_left);

cam.end();

plane.set(w,h);
plane.setPosition(w/2,h/2,0);

plane.draw();

shader.end();

if(unit_cube) {
ofNoFill();
ofDrawBox(glm::vec3(0),1.);
}

}

void ofApp::update() {

w = ofGetWidth();
h = ofGetHeight();

if(info) {
    printInfo();
}

}

void ofApp::printInfo() {

   cout << "Frame Rate: " + ofToString(ofGetFrameRate()) << endl;
   cout << "Camera Position : " + ofToString(cam.getPosition()) << endl;

}

void ofApp::windowResized(int w,int h) {

screen_size = ofToString(w) + "," + ofToString(h);

}

void ofApp::keyPressed(int key) {

    if(key == 's') {
        img.grabScreen(0,0,w,h);
        img.save("../../ " + frag + ".png");
    }

    if(key == 'f') {
      
        res = ofSystemLoadDialog("Load fragment shader",false,"../../src");
      
        if(res.bSuccess) {
            src_frag = res.getPath();
        }
        
        shader.load(path/"render.vert",src_frag);
    } 
}

void ofApp::mousePressed(int x,int y,int button) {

    if(button == 0) { 
    mouse_pressed_left = true;  
    }
}

void ofApp::mouseReleased(int x,int y,int button) {

    if(button == 0) {
    mouse_released_left = true;
    }

} 

void ofApp::mouseMoved(int x,int y) {

    mouse.x = x;
    mouse.y = y;

}
