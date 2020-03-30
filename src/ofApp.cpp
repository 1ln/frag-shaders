#include "ofApp.h"

void ofApp::setup() {

shader.load("render.vert","frag.frag");

w = ofGetWidth();
h = ofGetHeight();

plane.set(w,h);
plane.setPosition(w/2,h/2,0);
plane.setResolution(2,2);

cam_orbit = true;

gui_hide = false;

cam.setPosition(glm::vec3(0.0,0.0,5.0));
cam.lookAt(glm::vec3(0.0));

mouse.x = 0;
mouse.y = 0;

gui.setup();

gui.add(steps.set("Steps",100,0,1000));
gui.add(eps.set("Epsilon",0.0001,0.0001,0.01));
gui.add(dist.set("Distance",500,0,10000));

gui.add(aa.set("Anti Alias",1,0,3));

gui.add(shad_steps.set("Shadow Steps",25,0,1000));
gui.add(shad_eps.set("Shadow Epsilon",0.001,0.0,1.));
gui.add(shad_min.set("Shadow Minimum",0.005,0.,100.));
gui.add(shad_max.set("Shadow Maximum",0.005,0.,100.));
gui.add(shad_k.set("Shadow K",10.,0.,100.));

dif_listener = addListener(this,&ofApp::updateDif);
amb_listener = addListener(this,&ofApp::updateAmb);
spe_listener = addListener(this,&ofApp::updateSpe);

gui.add(dif.set("Diffuse",ofColor(255.,0.,0.),ofColor(0.),ofColor(255.)));
gui.add(amb.set("Ambient",ofColor(.05),ofColor(0.),ofColor(255.)));
gui.add(spe.set("Specular",ofColor(1.),ofColor(0.),ofColor(255.)));

gui.add(light_pos.set("Light Position",glm::vec3(1.0),glm::vec3(-1.),glm::vec3(1.)));
gui.add(light_intensity.set("Light Intensity",10.0,0.0,100.0));
gui.add(light_cam.set("Light Camera Follow",false));
 
gui.add(plane_offset.set("Plane offset",glm::vec3(0.0),glm::vec3(-10.),glm::vec3(10.0)));
gui.add(plane_orient.set("Plane orient",glm::vec3(0.0,1.0,0.0),glm::vec3(0.),glm::vec3(1.)));
gui.add(plane_display.set("Plane Display",false));

gui.add(screen_size.set("Screen Size",""));


}

void ofApp::draw() {

if(cam_orbit == true) {
cam.begin();
cam.end();
}

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());

shader.setUniform1i("u_steps",steps);
shader.setUniform1f("u_dist",dist);
shader.setUniform1f("u_eps",eps);

shader.setUniform1i("u_aa",aa);

shader.setUniform1i("u_shad_steps",shad_steps);
shader.setUniform1f("u_shad_eps",shad_eps);
shader.setUniform1f("u_shad_min",shad_min);
shader.setUniform1f("u_shad_max",shad_max);
shader.setUniform1f("u_shad_k",shad_k); 

shader.setUniform3f("u_light_pos",light_pos);
shader.setUniform1f("u_light_intensity",light_intensity);
shader.setUniform1i("u_light_cam",light_cam);

shader.setUniform3f("u_dif",dif);
shader.setUniform3f("u_amb",amb);
shader.setUniform3f("u_spe",spe);

shader.setUniform1i("u_plane_display",plane_display);
shader.setUniform3f("u_plane_orient",plane_orient);
shader.setUniform3f("u_plane_offset",plane_offset);

shader.setUniform2f("u_res",w,h);

shader.setUniform3f("u_cam_pos",glm::vec3(cam.getPosition()));
shader.setUniform3f("u_cam_tar",glm::vec3(cam.getTarget().getPosition()));

shader.setUniform2f("u_mouse_pos",mouse.x,mouse.y);

plane.draw();
//ofDrawRectangle(0,0,w,h);

shader.end();

if(!gui_hide) {
gui.draw();
}

}

void ofApp::update() {

w = ofGetWidth();
h = ofGetHeight();

}

void ofApp::windowResized(int w,int h) {

screen_size = ofToString(w) + "," + ofToString(h);
}

void ofApp::keyPressed(int key) {

    if(key == 's') {
        img.grabScreen(0,0,w,h);
        img.save("source.png");
    }

    if(key == 'h') {
        gui_hide = !gui_hide;
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

void ofApp::updateDif(ofColor & col) {
    dif = col;
}

void ofApp::updateAmb(ofColor & col) {
    amb = col;
}

void ofApp::updateSpe(ofColor & col) {
    spe = col;
}






