#include "ofApp.h"

void ofApp::setup() {

ofEnableDepthTest();

frag = "frag.frag"; 

path = filesystem::path("../../src");
shader.load(path/"render.vert",path/frag);

w = 512; 
h = 512;

imgw = 2500;
imgh = 2500;

ofSetWindowShape(w,h);

b.allocate(w,h);
b.getTexture().getTextureData().bFlipTexture = true;

d.allocate(16,16);
d.getTexture().disableMipmap();
d.getTexture().setTextureMinMagFilter(GL_NEAREST,GL_NEAREST);

ib.allocate(imgw,imgh);

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

scroll.x = 0;
scroll.y = 0;

mouse_ri = false;
mouse_le = false;
 
}

void ofApp::draw() {

b.begin();
ofClear(0);

shader.begin();
cam.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());
shader.setUniform1f("u_dtime",ofGetLastFrameTime());

shader.setUniform1i("u_hour",ofGetHours());
shader.setUniform1i("u_minutes",ofGetMinutes());
shader.setUniform1i("u_seconds",ofGetSeconds());

shader.setUniform2f("u_res",w,h);

shader.setUniform3f("u_cam_pos",glm::vec3(cam.getPosition()));
shader.setUniform3f("u_cam_tar",glm::vec3(cam.getTarget().getPosition())); 

shader.setUniform2f("u_mouse_pos",mouse.x,mouse.y);
shader.setUniform1i("u_mouse_le",mouse_le); 
shader.setUniform1i("u_mouse_ri",mouse_ri);
shader.setUniform1f("u_sclx",scroll.x);
shader.setUniform1f("u_scly",scroll.y);

shader.setUniformTexture("dtex",d.getTexture(0),0);

cam.end();

plane.set(w,h);
plane.setPosition(w/2,h/2,0);

plane.draw();

shader.end();

b.end();
b.draw(0,0,w,h);

d.begin();
d.end();

if(unit_cube) {
ofNoFill();
ofDrawBox(glm::vec3(0),1.);
}

}

void ofApp::update() {

//imgbuff.readToPixels(px);
//ofSaveImage(px,"../.. " + frag + ".png"); 

w = ofGetWidth();
h = ofGetHeight();

if(info) {
    printInfo();
}

}

void ofApp::printInfo() {

   cout << "Frame Rate: " + ofToString(ofGetFrameRate()) << endl;
   cout << "Camera Position : " + ofToString(cam.getPosition()) << endl;
   cout << "Screen size : " + screen_size << endl;

}

void ofApp::windowResized(int w,int h) {

screen_size = ofToString(w) + "," + ofToString(h);

}

void ofApp::keyPressed(int key) {

    if(key == 'i') { info = true; }

    if(key == 'g') {
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

    if(button == 0) { mouse_le = true; }
    if(button == 1) { mouse_ri = true; }

}

void ofApp::mouseReleased(int x,int y,int button) {

    if(button == 0) { mouse_le = false; }
    if(button == 1) { mouse_ri = false; }

} 

void ofApp::mouseScrolled(int x,int y,float scrollX,float scrollY) {

    scroll.x = scrollX;
    scroll.y = scrollY;

}


void ofApp::mouseMoved(int x,int y) {

    mouse.x = x;
    mouse.y = y;

}
