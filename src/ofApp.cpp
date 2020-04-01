#include "ofApp.h"

void ofApp::setup() {

shader.load("render.vert","frag.frag");

w = ofGetWidth();
h = ofGetHeight();

s.allocate(w,h,GL_RGBA);

s.begin();
ofClear(0);
s.end();

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

render_grp.setName("Renderer");
render_grp.add(steps.set("Steps",250,0,2500));
render_grp.add(eps.set("Epsilon",0.001,0.00001,.01));
render_grp.add(dist.set("Distance",500,0,2500));
render_grp.add(aa.set("Anti Alias",1,0,3));

info_grp.setName("Information");
//info_grp.add(hash.set("Hash",hash);
//info_grp.add(cam_pos.set("Camera position",glm::vec3(cam.getPosition()),""));
//info_grp.add(screen_size.set("Screen Size",""));
//info_grp.add(fps_counter.set("Frame per second",""));

//cam.add(cam_fov.set("Camera fov",45,0,2500));
//cam.add(cam_orbit.set("Camera orbit",true));

shad_grp.setName("Shadow");
shad_grp.add(shad_steps.set("Shadow Steps",25,0,500));
shad_grp.add(shad_eps.set("Shadow Epsilon",0.001,0.00001,0.01));
shad_grp.add(shad_min.set("Shadow Minimum",0.005,0.0,10.));
shad_grp.add(shad_max.set("Shadow Maximum",10.,0.,100.));
shad_grp.add(shad_k.set("Shadow K",10.,0.,100.));

light_grp.setName("Lighting");
light_grp.add(dif.set("Diffuse",glm::vec3(1.,0.,0.),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(amb.set("Ambient",glm::vec3(.05),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(spe.set("Specular",glm::vec3(1.),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(light_pos.set("Light Position",glm::vec3(1.0),glm::vec3(-1.),glm::vec3(1.)));
light_grp.add(light_intensity.set("Light Intensity",10.0,0.0,100.0));
light_grp.add(light_cam.set("Light Camera Follow",false));
 
plane_grp.setName("Plane Reflector");
plane_grp.add(plane_offset.set("Plane offset",glm::vec3(0.0),glm::vec3(-10.),glm::vec3(10.0)));
plane_grp.add(plane_orient.set("Plane orient",glm::vec3(0.0,1.0,0.0),glm::vec3(0.),glm::vec3(1.)));
plane_grp.add(plane_display.set("Plane Display",false));

gui.add(render_grp);
gui.add(info_grp);
gui.add(shad_grp);
gui.add(light_grp);
gui.add(plane_grp);

//gui.getGroup("Camera").minimize();
gui.getGroup("Shadow").minimize();
gui.getGroup("Lighting").minimize();
gui.getGroup("Plane Reflector").minimize();



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

plane.set(w,h);
plane.setPosition(w/2,h/2,0);
plane.draw();
//ofDrawRectangle(0,0,w,h);

shader.end();

if(!gui_hide) {
gui.draw();
}

}

void ofApp::update() {

fps_counter = ofToString(ofGetFrameRate());

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






