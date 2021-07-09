#version 330 core     

//dolson
//2020

out vec4 FragColor;

uniform vec2 resolution;
uniform float time;
uniform int seed; 

#define F x*x*(x*5.-4.)


#ifdef PCG
uint h(uint p) {
    uint state = 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}
#else
float h(float p) {
    return fract(sin(p)*43578.5453+seed);
}
#end

float n(vec2 x) {
    float i = floor(x);
    float f = fract(x);
    
    vec2 a,b,c,d;
    a = vec2(h(i));
    b = vec2(h(i+vec2(0.,1.));
    c = vec2(h(i+vec2(1.,0.));
    d = vec2(h(i+vec2(1.,1.));
     
    vec2 u = vec2(F);

    return mix(a,b,u.x) + (c-a)*u.y*(1.-u.x) + (d-b) * u.x*u.y;
}

float f6(vec2 p) {

    float f = 1.;

    f += 0.5      * n(p); p = m*p*2.01;
    f += 0.25     * n(p); p = m*p*2.02;
    f += 0.125    * n(p); p = m*p*2.04;
    f += 0.0625   * n(p); p = m*p*2.03;
    f += 0.0325   * n(p); p = m*p*2.06;
    f += 0.015625 * n(p);
    return f/0.92;    

}


vec3 color(float e) {
return vec3(.5) + vec3(.45,.33,.12)*cos(radians(180.) *
                  (vec3(1.)*e+vec3(0.,.1,.25)));
}
 
void main() {
 
    vec2 uv = (gl_FragCoord.xy-.5*resolution)/resolution.y;  
     
    float scl = 10.;
    uv *= scl;
    
    vec2 loc = floor(uv);
    c = vec3(h(loc.x+h(loc.y)));
    
    uv.y += time*.005;

    c = mix(c,vec3(1.,0.,0.),.5);

    FragColor = vec4(c,1.);

}
