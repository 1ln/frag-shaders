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

cam.setPosition(glm::vec3(0.0,0.0,5.0));
cam.lookAt(glm::vec3(0.0));

gui_hide = false;

mouse.x = 0;
mouse.y = 0;

gui.setup();

render_grp.setName("Renderer");
render_grp.add(steps.set("Steps",250,0,2500));
render_grp.add(eps.set("Epsilon",0.001,0.00001,.01));
render_grp.add(dist.set("Distance",500,0,2500));
render_grp.add(aa.set("Anti Alias",1,0,3));
render_grp.add(hash.set("Hash",.55121));

sdf_grp.setName("Signed Distance");
sdf_grp.add(sphere.set("Sphere",true));
sdf_grp.add(box.set("Box",false));
sdf_grp.add(torus.set("Torus",false));
sdf_grp.add(octahedron.set("Octahedron",false));
sdf_grp.add(prism.set("Prism",false));

info_grp.setName("Information");
info_grp.add(screen_size.set("Screen Size",""));
info_grp.add(fps_counter.set("Frame per second",""));

reflect_grp.setName("Reflection");
reflect_grp.add(reflect_steps.set("Reflect Steps",25,0,500));
reflect_grp.add(reflect_eps.set("Reflect Espilon",0.001,0.00001,0.01));
reflect_grp.add(reflect_min.set("Reflect Minimum",2.,0.,100.));
reflect_grp.add(reflect_max.set("Reflect Maximum",10.,0.,100.));

shad_grp.setName("Shadow");
shad_grp.add(shad_steps.set("Shadow Steps",25,0,500));
shad_grp.add(shad_eps.set("Shadow Epsilon",0.001,0.00001,0.01));
shad_grp.add(shad_min.set("Shadow Minimum",2.,0.0,100.));
shad_grp.add(shad_max.set("Shadow Maximum",10.,0.,100.));
shad_grp.add(shad_k.set("Shadow K",10.,0.,100.));

light_grp.setName("Lighting");
light_grp.add(col.set("Color",glm::vec3(0.,1.,0.),glm::vec3(0.),glm::vec3(1.))); 
light_grp.add(dif.set("Diffuse",glm::vec3(.5),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(amb.set("Ambient",glm::vec3(.05),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(spe.set("Specular",glm::vec3(1.),glm::vec3(0.),glm::vec3(1.)));
light_grp.add(light_pos.set("Light Position",glm::vec3(10.0),glm::vec3(-100.),glm::vec3(100.)));
light_grp.add(light_intensity.set("Light Intensity",10.0,0.0,100.0));
light_grp.add(light_cam.set("Light Camera Follow",false));
 
plane_grp.setName("Plane Reflector");
plane_grp.add(plane_offset.set("Plane Offset",5,0,100));
plane_grp.add(plane_display.set("Plane Display",false));

fnoise_grp.setName("Fractional Noise"); 
fnoise_grp.add(fnoise_display.set("Display",false));
fnoise_grp.add(fnoise_octaves.set("Octaves",2,1,8));
fnoise_grp.add(fnoise_frequency.set("Frequency",.5,0.,.5));
fnoise_grp.add(fnoise_offset.set("Offset",.25,0.,1.));

sine_grp.setName("Sine Distortion");
sine_grp.add(sine_display.set("Sine Display",false));
sine_grp.add(sine_offset.set("Sine Offset",0.025,0,1));
sine_grp.add(sine_height.set("Sine Height",10,0,100));

gui.add(render_grp);
gui.add(sdf_grp);
gui.add(info_grp);
gui.add(reflect_grp);
gui.add(shad_grp);
gui.add(light_grp);
gui.add(plane_grp);
gui.add(sine_grp);
gui.add(fnoise_grp);

gui.getGroup("Signed Distance").minimize();
gui.getGroup("Information").minimize();
gui.getGroup("Reflection").minimize();
gui.getGroup("Shadow").minimize();
gui.getGroup("Lighting").minimize();
gui.getGroup("Plane Reflector").minimize();
gui.getGroup("Sine Distortion").minimize();
gui.getGroup("Fractional Noise").minimize();

}

void ofApp::draw() {

cam.begin();
cam.end();

shader.begin();

shader.setUniform1f("u_time",ofGetElapsedTimeMillis());

shader.setUniform1i("u_steps",steps);
shader.setUniform1f("u_dist",dist);
shader.setUniform1f("u_eps",eps);

shader.setUniform1i("u_aa",aa);

shader.setUniform1i("u_refl_steps",reflect_steps);
shader.setUniform1f("u_refl_eps",reflect_eps);
shader.setUniform1f("u_refl_min",reflect_min);
shader.setUniform1f("u_refl_max",reflect_max);

shader.setUniform1i("u_shad_steps",shad_steps);
shader.setUniform1f("u_shad_eps",shad_eps);
shader.setUniform1f("u_shad_min",shad_min);
shader.setUniform1f("u_shad_max",shad_max);
shader.setUniform1f("u_shad_k",shad_k); 

shader.setUniform3f("u_light_pos",light_pos);
shader.setUniform1f("u_light_intensity",light_intensity);
shader.setUniform1i("u_light_cam",light_cam);

shader.setUniform3f("u_col",col);
shader.setUniform3f("u_dif",dif);
shader.setUniform3f("u_amb",amb);
shader.setUniform3f("u_spe",spe);

shader.setUniform1i("u_sphere",sphere);
shader.setUniform1i("u_box",box);
shader.setUniform1i("u_torus",torus);
shader.setUniform1i("u_octahedron",octahedron);
shader.setUniform1i("u_prism",prism);

shader.setUniform1i("u_plane_display",plane_display);
shader.setUniform1f("u_plane_offset",plane_offset);

shader.setUniform1i("u_sine_display",sine_display);
shader.setUniform1f("u_sine_offset",sine_offset);
shader.setUniform1f("u_sine_height",sine_height);

shader.setUniform1i("u_fnoise_display",fnoise_display);
shader.setUniform1i("u_fnoise_octaves",fnoise_octaves);
shader.setUniform1f("u_fnoise_frequency",fnoise_frequency);
shader.setUniform1f("u_fnoise_offset",fnoise_offset);

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






