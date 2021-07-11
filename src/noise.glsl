#version 330 core     

//dolson
//2020

out vec4 FragColor;

uniform vec2 resolution;
uniform float time;
uniform int seed; 

#define F f*f*(f*5.-4.)

uint pcg(uint p) {
    uint state = p* 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

float h(vec2 p) {
    return fract(sin(p.x+fract(sin(p.y+float(seed))))*float(seed));
}

float n(vec2 x) {
    vec2 i = floor(x);
    vec2 f = fract(x);

    float a,b,c,d;
    a = h(i);
    b = h(i+vec2(1.,0.));
    c = h(i+vec2(0.,1.));
    d = h(i+vec2(1.,1.));

    vec2 u = vec2(F);

    return mix(a,b,u.x) + (c-a)*u.y*(1.-u.x) + (d-b) * u.x*u.y;
}

mat2 m = mat2(0.8,0.6,-0.6,0.8);
float f(vec2 p) {

    float f = 1.;

    f += 0.5      * n(p); p = p*2.01*m;
    f += 0.25     * n(p); p = p*2.03*m;
    f += 0.125    * n(p); p = p*2.05*m;
    f += 0.0625   * n(p); p = p*2.07*m;
    f += 0.03125  * n(p); p = p*2.09*m;
    f += 0.015625 * n(p); p = p*2.11*m;
    f += 0.007825 * n(p); p = p*2.13*m;
    f += 0.003925 * n(p);
    return f/0.92;    

}

vec3 color(float e) {
return vec3(.5) + vec3(.45,.33,.12)*cos(radians(180.) *
                  (vec3(1.)*e+vec3(0.,.1,.25)));
}

void main() {
 
    vec2 uv = (gl_FragCoord.xy-.5*resolution)/resolution.y;  

    uv.y += time * .05;    
     
    float scl = 10.;
    vec2 loc = floor(uv);

    vec3 c;
    c = vec3(color(f(uv+f(uv+f(uv)))));
    FragColor = vec4(c,1.);

}
