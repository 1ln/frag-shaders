#version 430 core

// dolson,2019

out vec4 FragColor;

uniform vec2 resolution;
uniform float time;

#define SEED 1523595.
#define EPS 0.0001
#define EXPONENTIAL
#define BLEND_K 2.

float hash(float p) {
    return fract(sin(p) * 4358.5453);
}

float hash(vec2 p) {
   return fract(sin(dot(p.xy,vec2(12.9898,78.233)))*4358.5353);
}

vec2 mod289(vec2 p) { return p - floor(p * (1. / 289.)) * 289.; }
vec3 mod289(vec3 p) { return p - floor(p * (1. / 289.)) * 289.; }
vec3 permute(vec3 p) { return mod289(((p * 34.) + 1.) * p); } 

float ns2(vec2 p) {

    const float k1 = (3. - sqrt(3.))/6.;
    const float k2 = .5 * (sqrt(3.) -1.);
    const float k3 = -.5773;
    const float k4 = 1./41.;

    const vec4 c = vec4(k1,k2,k3,k4);
    
    vec2 i = floor(p + dot(p,c.yy));
    vec2 x0 = p - i + dot(i,c.xx);
  
    vec2 i1;
    i1 = (x0.x > x0.y) ? vec2(1.,0.) : vec2(0.,1.);
    vec4 x12 = x0.xyxy + c.xxzz;
    x12.xy -= i1;

    i = mod289(i);
    
    vec3 p1 = permute(permute(i.y + vec3(0.,i1.y,1.))
        + i.x + vec3(0.,i1.x,1.));
  
    p1 = permute(mod289(p1 + vec3(float(SEED))));

    vec3 m = max(.5 - 
    vec3(dot(x0,x0),dot(x12.xy,x12.xy),dot(x12.zw,x12.zw)),0.);
    m = m * m; 
    m = m * m;

    vec3 x = fract(p1 * c.www) - 1.;
    vec3 h = abs(x) - .5;
    vec3 ox = floor(x + .5);
    vec3 a0 = x - ox; 
    m *= 1.792842 - 0.853734 * (a0 * a0 + h * h);
     
    vec3 g;
    g.x = a0.x * x0.x + h.x * x0.y;
    g.yz = a0.yz * x12.xz + h.yz * x12.yw;
    return 130. * dot(m,g);
}

float f(vec2 x,int octaves) {

    float f = 0.;

    for(int i = 1; i < octaves; i++) {
 
    float e = pow(2.,float(i));
    float s = (1./e);
    f += ns2(x*e)*s;   
    
    }    

    return f * .5 + .5;
}

float envImp(float x,float k) {

    float h = k * x;
    return h * exp(1.0 - h);
}

vec3 fmCol(float t,vec3 a,vec3 b,vec3 c,vec3 d) {
    
    return a + b * cos( (radians(180)*2.0) * (c * t + d));
}

float easeOut4(float t) {

    return -1.0 * t * (t - 2.0);

}

float easeInOut3(float t) {

    if((t *= 2.0) < 1.0) {
        return 0.5 * t * t * t;
    } else { 
        return 0.5 * ((t -= 2.0) * t * t + 2.0);

    }
}

mat2 rot2(float a) {

    float c = cos(a);
    float s = sin(a);
    
    return mat2(c,-s,s,c);
}

vec2 opu(vec2 d1,vec2 d2) {

    return (d1.x < d2.x) ? d1 : d2;
} 

#ifdef POLYNOMIAL
float smin(float d1,float d2,float k) {

    float h = clamp(0.5 + 0.5 * (d2-d1)/k,0.0,1.0);
    return mix(d2,d1,h) - k * h * (1.0 - h);
}
#endif

#ifdef EXPONENTIAL
float smin(float d1,float d2,float k) {
    float res = exp2(-k * d1) + exp2(-k * d2);
    return -log2(res)/k;
}
#endif

#ifdef POWER
float smin(float d1,float d2,float k) {
     d1 = pow(d1,k);
     d2 = pow(d2,k);
     return pow((d1*d2) / (d1+d2),1./k);
}  
#endif

vec2 blend(vec2 d1,vec2 d2) {

    float d = smin(d1.x,d2.x,BLEND_K);
    float m = mix(d1.y,d2.y,clamp(d1.x-d,0.,1.));
    return vec2(d,m);
}

float plane(vec3 p,vec4 n) {

    return dot(p,n.xyz) + n.w;
}

vec2 scene(vec3 p) {

    vec2 res = vec2(1.,0.);
    vec3 pl = p;
    vec3 q = p;    

    p.y += ns2(p.xz * .005 + f(p.xz * .025,6) * .125) * 10.;

    vec2 h = vec2(plane(p,vec4(0.,1.,0.,1.)),25.);
    vec2 l = vec2(pl.y + 1.,45.);

    res = blend(l,h);
    return res;

}

vec2 rayScene(vec3 ro,vec3 rd) {
    
    float d = -1.0;
    float s = 0.;
    float e = 100.;  

    for(int i = 0; i < 10; i++) {

        vec3 p = ro + s * rd;
        vec2 dist = scene(p);
   
        if(abs(dist.x) < EPS || e <  dist.x ) { break; }
        s += dist.x;
        d = dist.y;

        }
 
        if(e < s) { d = -1.0; }
        return vec2(s,d);

}

vec3 fog(vec3 col,vec3 fcol,float fdist,float fdens,float y) {
    float fdep = 1. - exp(-fdist * pow(fdens,y));
    return mix(col,fcol,fdep);
}

vec3 scatter(vec3 col,float distance,float density,vec3 rd,vec3 ld) {
    float fog  = 1. - exp(-distance * density);
    float light = max(dot(rd,ld),0.);
    vec3 fog_col = mix(vec3(.5,.6,.7),vec3(.6,.5,.1),pow(light,8.));
    return mix(col,fog_col,fog);
}

float shadow(vec3 ro,vec3 rd ) {

    float res = 1.0;
    float t = 0.005;
    float ph = 1e10;
    
    for(int i = 0; i < 105; i++ ) {
        
        float h = scene(ro + rd * t  ).x;

        float y = h * h / (2. * ph);
        float d = sqrt(h*h-y*y);         
        res = min(res,64. * d/max(0.,t-y));
        ph = h;
        t += h;
    
        if(res < EPS || t > 10.) { break; }

        }

        return clamp(res,0.0,1.0);

}

vec3 calcNormal(vec3 p) {

    vec2 e = vec2(1.0,-1.0) * EPS;

    return normalize(vec3(
    vec3(e.x,e.y,e.y) * scene(p + vec3(e.x,e.y,e.y)).x +
    vec3(e.y,e.x,e.y) * scene(p + vec3(e.y,e.x,e.y)).x +
    vec3(e.y,e.y,e.x) * scene(p + vec3(e.y,e.y,e.x)).x + 
    vec3(e.x,e.x,e.x) * scene(p + vec3(e.x,e.x,e.x)).x

    ));
    
}

vec3 rayCamDir(vec2 uv,vec3 camPosition,vec3 camTarget,float fPersp) {

     vec3 camForward = normalize(camTarget - camPosition);
     vec3 camRight = normalize(cross(vec3(0.0,1.0,0.0),camForward));
     vec3 camUp = normalize(cross(camForward,camRight));


     vec3 vDir = normalize(uv.x * camRight + 
                 uv.y * camUp + camForward * fPersp);  

     return vDir;
}

vec3 render(vec3 ro,vec3 rd) {
 
vec2 d = rayScene(ro, rd);

vec3 col = vec3(1.) - max(rd.y,0.);

if(d.y >= 0.) {

vec3 p = ro + rd * d.x;
vec3 n = calcNormal(p);
vec3 l = normalize(vec3(10.,15.,15.));
vec3 h = normalize(l - rd);
vec3 r = reflect(rd,n);

col = .2+.2*sin(2.*d.y*vec3(.45));

float amb = sqrt(clamp(0.5 + 0.5 * n.y,0.0,1.0));
float dif = clamp(dot(n,l),0.0,1.0);

float spe = pow(clamp(dot(n,h),0.0,1.0),16.) * dif * (.04 + 0.9 * pow(clamp(1. + dot(h,rd),0.,1.),5.));

float ind = clamp(dot(n,normalize(h*vec3(-1.,0.,1.))),0.,1.);

float fre = pow(clamp(1. + dot(n,rd),0.0,1.0),2.0);
float ref = smoothstep(-.2,.2,r.y);
vec3 linear = vec3(0.);

dif *= shadow(p,l);
ref *= shadow(p,r);

linear += dif * vec3(1.5,1.,1.25);
linear += amb * vec3(0.02,0.4,0.1);
linear += ref * vec3(0.05,0.01,0.04);
linear += fre * vec3(0.04,0.01,0.05);
linear += ind * vec3(0.04,0.5,0.33);

col = col * linear;
col += 5. * spe * vec3(0.06,0.005,0.004 );

col = scatter(col,.00025,d.x*d.x,rd,l);

}

return col;
}

void main() {
 
vec3 color = vec3(0.);

vec3 cam_tar = vec3(0.);
vec3 cam_pos = vec3(25.,15.,35.);

cam_tar.z += time;
cam_pos += cam_tar;

vec2 uv = (2. * gl_FragCoord.xy - resolution.xy) / resolution.y; 

vec3 dir = rayCamDir(uv,cam_pos,cam_tar,1.); 
color = render(cam_pos,dir);  
color = pow(color,vec3(.4545));      
FragColor = vec4(color,1.0);

}
