//Signed Distance with raymarching

// dolson,2019

precision mediump float;

uniform float u_hash;

uniform vec2 u_mouse;
uniform int u_mouse_pressed;

uniform vec2 u_resolution;
uniform float u_time;

const float E   =  2.7182818;
const float PI  =  radians(180.0); 
const float PHI =  (1.0 + sqrt(5.0)) / 2.0;
const float nhash = .5125094;

float hash(float p) {
    return fract(sin(p) * nhash * 43758.5453); 
}

float hash21(vec2 p) {
    return fract(sin(dot(p.xy,vec2(12.9898,78.233))) * 43758.5453123);
}
  

float cell(vec2 x,float iterations,int type) {
 
    x *= iterations;

    vec2 p = floor(x);
    vec2 f = fract(x);
 
    float min_dist = 1.0;
    
    for(int i = -1; i <= 1; i++) {
        for(int j = -1; j <= 1; j++) {
            
                vec2 b = vec2(float(i),float(j));
                vec2 r = hash21( p + b );
                
                vec2 diff = (b + r - f);

                float d = length(diff);

                    if(type == 0) { 
                        min_dist = min(min_dist,d);
                    }
 
                    if(type == 1) {
                        min_dist = min(min_dist,abs(diff.x)+abs(diff.y)+abs(diff.z));
                    }

                    if(type == 2) {
                        min_dist = min(min_dist,max(abs(diff.x),max(abs(diff.y),abs(diff.z))));
                    }

            }
        }
    }
 
    return min_dist;  

}

float noise(vec2 x) {

vec2 p = floor(x);
vec2 f = fract(x);

f = f*f*(3.-2.*f);
float n = p.x + p.y * 157.;


return mix(mix(hash(n + 0.), 
               hash(n + 1.),f.x), 
       mix(hash(n + 157.), 
               hash(n + 158.),f.x),f.y);
}


float fractional(vec2 x,float h) {

    float t = 0.0;

    float g = exp2(-h); 

    float a = 0.5;
    float f = 1.0;

    for(int i = 0; i < 5; i++) {
 
    t += a * noise(f * x); 
    f *= 2.0; 
    a *=  g;  
    
    }    

    return t;
}

float sin2(vec2 p,float h) {
    
    return sin(p.x*h) * sin(p.y*h);
}

float envImpulse(float x,float k) {

    float h = k * x;
    return h * exp(1.0 - h);
}

float envStep(float x,float k,float n) {

    return exp(-k * pow(x,n));
}

float cubicImpulse(float x,float c,float w) {

    x = abs(x - c);
    if( x > w) { return 0.0; }
    x /= w;
    return 1.0 - x * x  * (3.0 - 2.0 * x);

}

float sincPhase(float x,float k) {

    float a = PI * (k * x - 1.0);
    return sin(a)/a;
}

void main() {

vec2 uv = gl_FragCoord.xy / u_resolution.xy;
uv.x *= u_resolution.x / u_resolution.y;

vec3 col = vec3(0.0);
col += fractional(uv * 10.,.5);

gl_FragColor = vec4(col,1.0);

}
