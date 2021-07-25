#version 330 core     

//Dan Olson

out vec4 FragColor;

uniform vec2 resolution;
uniform float time;

uniform int seed;

const int steps = 250;
float eps = 0.00001;
float dmin = 0.001;
float dmax = 24;
const int aa = 2;
 
const int shsteps = 45; 
float shblur = 100.;
float shmax = 3.; 

vec3 ta = vec3(0.,0.5,3.);
vec3 ta2 = vec3(.5,-2.,10.);
vec3 ta3 = vec3(1.,.5,-5.);

const float PI2 = radians(180.)*2.;

float h11(float p) {
    uvec2 n = uint(int(p)) * uvec2(uint(int(seed)),2531151992.0);
    uint h = (n.x ^ n.y) * uint(int(seed));
    return float(h) * (1./float(0xffffffffU));
}

mat2 rot(float a) {
    return mat2(cos(a),-sin(a),sin(a),cos(a));
}

mat3 rotAxis(vec3 axis,float theta) {

axis = normalize(axis);

    float c = cos(theta);
    float s = sin(theta);

    float oc = 1.0 - c;

    return mat3(
 
        oc * axis.x * axis.x + c, 
        oc * axis.x * axis.y - axis.z * s,
        oc * axis.z * axis.x + axis.y * s, 
    
        oc * axis.x * axis.y + axis.z * s,
        oc * axis.y * axis.y + c, 
        oc * axis.y * axis.z - axis.x * s,

        oc * axis.z * axis.x - axis.y * s,
        oc * axis.y * axis.z + axis.x * s, 
        oc * axis.z * axis.z + c);

}

mat3 camEuler(float yaw,float pitch,float roll) {

     vec3 f = -normalize(vec3(sin(yaw),sin(pitch),cos(yaw)));
     vec3 r = normalize(cross(f,vec3(0.0,1.0,0.0)));
     vec3 u = normalize(cross(r,f));

     return rotAxis(f,roll) * mat3(r,u,f);
}

vec2 opu(vec2 d1,vec2 d2) {
    return (d1.x < d2.x) ? d1 : d2;
} 

float smou(float d1,float d2,float k) {

    float h = clamp(0.5 + 0.5 * (d2-d1)/k,0.0,1.0);
    return mix(d2,d1,h) - k * h * (1.0 - h);
}

float smod(float d1,float d2,float k) {

    float h = clamp(0.5 - 0.5 * (d2+d1)/k,0.0,1.0);
    return mix(d2,-d1,h) + k * h * (1.0 - h);
}

float smoi(float d1,float d2,float k) {

    float h = clamp(0.5 + 0.5 * (d2-d1)/k,0.0,1.0);
    return mix(d2,d1,h) + k * h * (1.0 - h);

}

vec2 scene(vec3 p) {

    vec2 res = vec2(1.,0.);

    float d = 0.; 

    vec3 q = p; 
    q.xz *= rot(time*.01);   
    d = -(length(q)-21.5);

    res = opu(res,vec2(d,25.)); 
    res = opu(res,vec2(length(p-ta)-.25,1.));
    res = opu(res,vec2(length(p-ta2)-1.,2.));
    res = opu(res,vec2(length(p-ta3)-.5,3.));

    return res;

}

vec2 rayScene(vec3 ro,vec3 rd) {
    
    float d = -1.0;
    float s = dmin;
    float e = dmax;  

    for(int i = 0; i < steps; i++) {

        vec3 p = ro + s * rd;
        vec2 dist = scene(p);
   
        if(abs(dist.x) < eps || e <  dist.x ) { break; }
        s += dist.x;
        d = dist.y;

        }
 
        if(e < s) { d = -1.0; }
        return vec2(s,d);

}

float expStep(float x,float k) {
    return exp((x*k)-k);
}

float shadow(vec3 ro,vec3 rd) {

    float res = 1.0;
    float t = 0.005;
    float ph = 1e10;
    
    for(int i = 0; i < shsteps; i++ ) {
        
        float h = scene(ro + rd * t  ).x;

        float y = h * h / (2. * ph);
        float d = sqrt(h*h-y*y);         
        res = min(res,shblur * d/max(0.,t-y));
        ph = h;
        t += h;
    
        if(res < eps || t > shmax) { break; }

        }

        return clamp(res,0.0,1.0);

}

vec3 calcNormal(vec3 p) {

    vec2 e = vec2(1.0,-1.0) * eps;

    return normalize(vec3(
    vec3(e.x,e.y,e.y) * scene(p + vec3(e.x,e.y,e.y)).x +
    vec3(e.y,e.x,e.y) * scene(p + vec3(e.y,e.x,e.y)).x +
    vec3(e.y,e.y,e.x) * scene(p + vec3(e.y,e.y,e.x)).x + 
    vec3(e.x,e.x,e.x) * scene(p + vec3(e.x,e.x,e.x)).x

    ));
    
}

vec3 renderScene(vec3 ro,vec3 rd) {

vec2 d = rayScene(ro, rd);

vec3 col = vec3(1.) * max(0.,rd.y);

if(d.y >= 0.) {

vec3 p = ro + rd * d.x;
vec3 n = calcNormal(p);

if(d.y == 1.) {
    col = vec3(.5);
}

if(d.y == 2.) {
    col = vec3(1.,.5,.25);
}

if(d.y == 3.) {
    col = vec3(1.);
}

if(d.y == 25.) {
    col = vec3(.5);
}

vec3 l = normalize(ta-vec3(.25));
float rad = dot(rd,l);
col += col * vec3(.05,.25,.16)*expStep(rad,75.);
col += col * vec3(.5,.25,.1)*expStep(rad,100.);
col += col * vec3(.25)*expStep(rad,25.);
col += col * vec3(1.,.5,.5)*expStep(rad,100.);

vec3 la = normalize(ta2-vec3(2.,10.,-5.));
float rad1 = dot(rd,la);
col += col * vec3(.25,.89,.01)*expStep(rad1,25.);
col += col * vec3(1.)*expStep(rad1,112.);
col += col * vec3(.25)*expStep(rad1,131.);
col += col * vec3(.5)*expStep(rad1,16.);

vec3 lb = normalize(ta3-vec3(-5.));
float rad2 = dot(rd,lb);
col += col * vec3(1.,.12,.5)*expStep(rad2,100.);
col += col * vec3(.5,.25,.95)*expStep(rad2,50.);
col += col * vec3(.05,.35,.0025)*expStep(rad2,45.);
col += col * vec3(1.)*expStep(rad2,25.);

vec3 h = normalize(l - rd);
vec3 r = reflect(rd,n);}

float amb = clamp(0.5 + 0.5 * n.y,0.,1.);

float dif = clamp(dot(n,l),0.0,1.0);

float spe = pow(clamp(dot(n,h),0.0,1.0),16.)
* dif * (.04 + 0.9 * pow(clamp(1. + dot(h,rd),0.,1.),5.));

float fre = pow(clamp(1. + dot(n,rd),0.0,1.0),2.0);
float ref = smoothstep(-.2,.2,r.y);

vec3 linear = vec3(0.);

dif *= shadow(p,l);
ref *= shadow(p,r);

linear += .05+ dif * vec3(0.5);
linear += .02+ amb * vec3(.005);

col = col * linear;
col += spe * vec3(.025); 
col = mix(col,vec3(1.),1.-exp(-0.00001 * d.x*d.x*d.x)); 

}

return col;
}

void main() {
 
vec3 color = vec3(0.);
vec3 ro = vec3(0.);

for(int k = 0; k < aa; ++k) {
    for(int l = 0; l < aa; ++l) {

    vec2 o = vec2(float(l),float(k)) / float(aa) - .5;

    vec2 uv = (2. * (gl_FragCoord.xy + o) -
    resolution.xy) / resolution.y; 

    mat3 cm = camEuler(time*PI2*0.01,0.45,0.);
    vec3 rd = cm * normalize(vec3(uv.xy,2.));

    vec3 render = renderScene(ro,rd);   

    render = pow(render,vec3(.4545));
    color += render;
    }
}

color /= float(aa*aa);
FragColor = vec4(color,1.0);

}
