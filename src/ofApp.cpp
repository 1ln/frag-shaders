#include "ofApp.h"

void ofApp::setup() {

path = filesystem::path("../../src");
shader.load(path/"render.vert",path/"frag.frag");

w = ofGetWidth();
h = ofGetHeight();

s.allocate(w,h,GL_RGBA);

s.begin();
ofClear(0);
s.end();

plane.set(w,h);
plane.setPosition(w/2,h/2,0);
plane.setResolution(2,2);

cam.setPosition(glm::vec3(1.0,2.0,5.0));
cam.lookAt(glm::vec3(0.0));
cam.setNearClip(.01);

light.setup();
light.setPosition(glm::vec3(0,5,10));

gui_display = false;
info = false;
unit_cube = false;

mouse.x = 0;
mouse.y = 0;

eps = 0.0001;
hash = 0.155126;

gui.setup();

render_grp.setName("Renderer");
render_grp.add(steps.set("Steps",250,0,2500));
render_grp.add(dist.set("Distance",500,0,2500));
render_grp.add(aa.set("Anti Alias",1,0,3));

reflect_grp.setName("Reflection");
reflect_grp.add(reflect_steps.set("Reflect Steps",25,0,500));
reflect_grp.add(reflect_dist.set("Reflect Distance",10.,0.,100.));

shad_grp.setName("Shadow");
shad_grp.add(shad_steps.set("Shadow Steps",100.,0,500));
shad_grp.add(shad_dist.set("Shadow Distance",33.,1.,100.));
shad_grp.add(shad_k.set("Shadow K",10.,0.,100.));

light_grp.setName("Lighting");
light_grp.add(col.set("Color",glm::vec3(0.,1.,0.),glm::vec3(0.),glm::vec3(1.))); 
light_grp.add(dif.set("Diffuse",glm::vec3(.5),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(amb.set("Ambient",glm::vec3(.05),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(spe.set("Specular",glm::vec3(1.),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(light_pos.set("Light Position",light.getPosition(),glm::vec3(-100.),glm::vec3(100.)));
light_grp.add(light_intensity.set("Light Intensity",10.0,0.0,100.0));

gui.add(render_grp);
gui.add(reflect_grp);
gui.add(shad_grp);
gui.add(light_grp);

gui.getGroup("Reflection").minimize();
gui.getGroup("Shadow").minimize();
gui.getGroup("Lighting").minimize();

}

void ofApp::draw() {

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());

shader.setUniform1i("u_steps",steps);
shader.setUniform1f("u_dist",dist);
shader.setUniform1f("u_eps",eps);

shader.setUniform1i("u_aa",aa);

shader.setUniform1i("u_refl_steps",reflect_steps);
shader.setUniform1f("u_refl_dist",reflect_dist);

shader.setUniform1i("u_shad_steps",shad_steps);
shader.setUniform1f("u_shad_dist",shad_dist);
shader.setUniform1f("u_shad_k",shad_k); 

shader.setUniform3f("u_light_pos",light_pos);
shader.setUniform1f("u_light_intensity",light_intensity);

shader.setUniform3f("u_col",col);
shader.setUniform3f("u_dif",dif);
shader.setUniform3f("u_amb",amb);
shader.setUniform3f("u_spe",spe);

shader.setUniform2f("u_res",w,h);

shader.setUniform3f("u_cam_pos",glm::vec3(cam.getPosition()));
shader.setUniform3f("u_cam_tar",glm::vec3(cam.getTarget().getPosition()));

shader.setUniform2f("u_mouse_pos",mouse.x,mouse.y);

plane.set(w,h);
plane.setPosition(w/2,h/2,0);
plane.draw();
shader.end();

cam.begin();

if(unit_cube) {
ofNoFill();
ofDrawBox(glm::vec3(0),1.);
}

cam.end(); 

if(gui_display) {
gui.draw();
}

}

void ofApp::update() {

w = ofGetWidth();
h = ofGetHeight();

light.setPosition(light_pos);

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
        img.save("../../source.png");
    }

    if(key == 'f') {
      
        res = ofSystemLoadDialog("Load fragment shader",false,"../../src");
      
        if(res.bSuccess) {
            src_frag = res.getPath();
        }
        
        shader.load(path/"render.vert",src_frag);
    } 

    if(key == 'h') {
        gui_display = !gui_display;
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






