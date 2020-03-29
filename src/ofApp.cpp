#include "ofApp.h"

void ofApp::setup() {

shader.load("render.vert","frag.frag");

w = ofGetWidth();
h = ofGetHeight();

cam_orbit = true;

cam.setPosition(glm::vec3(0.0,0.0,5.0));
cam.lookAt(glm::vec3(0.0));

mouse.x = 0;
mouse.y = 0;

gui.setup();

gui.add(steps.setup("Steps",100,0,1000));
gui.add(eps.setup("Epsilon",0.001,0.0000001,1.));
gui.add(trace_dist.setup("Trace Distance",500,0,10000));

gui.add(aa.setup("Anti Alias",1,0,3));

gui.add(shad_steps.setup("Shadow Steps",25,0,1000));
gui.add(shad_eps.setup("Shadow Epsilon",0.001,0.0000001,1.));
gui.add(shad_min.setup("Shadow Minimum",0.005,0.,100.));
gui.add(shad_max.setup("Shadow Maximum",0.005,0.,100.));
gui.add(shad_k.setup("Shadow K",10.,0.,100.));

gui.add(dif.setup("Diffuse",ofColor(255.,0.,0.),ofColor(0.),ofColor(255.)));
gui.add(amb.setup("Ambient",ofColor(.05),ofColor(0.),ofColor(255.)));
gui.add(spe.setup("Specular",ofColor(1.),ofColor(0.),ofColor(255.)));

gui.add(light_pos.setup("Light Position",glm::vec3(0.0,0.0,5.0));
gui.add(light_intensity.setup("Light Intensity",10.0,0.0,100.0));
gui.add(light_cam.setup("Light Camera Follow",false));
 
gui.add(plane_offset.setup("Plane offset",glm::vec3(0.0,0.0,5.0));
gui.add(plane_orient.setup("Plane orient",glm::vec3(0.0,1.0,0.0));
gui.add(plane_display.setup("Plane Display",false));

}

void ofApp::draw() {

if(cam_orbit == true) {
cam.begin();
cam.end();
}

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());

shader.setUniform1i("u_steps",steps);
shader.setUniform1f("u_dist",trace_dist);
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

ofDrawRectangle(0,0,w,h);

shader.end();

gui.draw();


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
