#version 330 core

//dan olson
//2020

uniform vec2 resolution;
uniform float time;

float hash(float p) {
    uvec2 n = uint(int(p)) * uvec2(1391674541U,2531151992.0 * 21885);
    uint h = (n.x ^ n.y) * 1391674541U;
    return float(h) * (1.0/float(0xffffffffU));
}
 
float noise(vec3 x) {

    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f * f * (3.0 - 2.0 * f);
    float n = p.x + p.y * 157.0 + 113.0 * p.z;

    return mix(mix(mix(hash(  n +   0.0) , hash(   n +   1.0)  ,f.x),
                   mix(hash(  n + 157.0) , hash(   n + 158.0)   ,f.x),f.y),
               mix(mix(hash(  n + 113.0) , hash(   n + 114.0)   ,f.x),
                   mix(hash(  n + 270.0) , hash(   n + 271.0)   ,f.x),f.y),f.z);
}


float f3(vec3 x,int octaves,float h) {

    float t = 0.0;

    float g = exp2(-h); 

    float a = 0.5;
    float f = 1.0;

    for(int i = 0; i < octaves; i++) {
 
    t += a * noise(f * x); 
    f *= 2.0; 
    a *=  g;  
    
    }    

    return t;
}

mat2 rot2(float a) {

    float c = cos(a);
    float s = sin(a);
    
    return mat2(c,-s,s,c);
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

vec2 res = vec2(1.0,0.0);

vec3 q = vec3(p);

vec4 v = vec4(0.,1.,-1.,0.5 + sqrt(1.25));
v /= length(v.zw);

float d;
d = abs(dot(p,v.xyw))-1.;
d = max(d,abs(dot(p,v.ywx))-1.);
d = max(d,abs(dot(p,v.wxy))-1.);
d = max(d,abs(dot(p,v.xzw))-1.);
d = max(d,abs(dot(p,v.zwx))-1.);
d = max(d,abs(dot(p,v.wxz))-1.);

res = opu(res,vec2(d,2.));

float pl = dot(q,vec3(0.,1.,0.)) + 1.;
res = opu(res,vec2(pl,1.));  


return res;

}

vec2 rayScene(vec3 ro,vec3 rd) {
    
    float depth = 0.0;
    float d = -1.0;

    for(int i = 0; i < 512; i++) {

        vec3 p = ro + depth * rd;
        vec2 dist = scene(p);
   
        if(abs( dist.x) < 0.0001 || 2500. <  dist.x ) { break; }
        depth += dist.x;
        d = dist.y;

        }
 
        if(2500. < depth) { d = -1.0; }
        return vec2(depth,d);

}

vec3 fog(vec3 c,vec3 fc,float b,float distance) {
    float depth = 1. - exp(-distance *b);
    return mix(c,fc,depth);
}

float reflection(vec3 ro,vec3 rd,float dmin,float dmax ) {

    float depth = dmin;
    float h = 0.;

    for(int i = 0; i < 100; i++ ) {
        h = scene(ro + rd * depth).x;

        if(abs( h) < 0.0001 || depth > dmax ) { break; }
        
        depth += h;
    }
    
    return depth;
}

float shadow(vec3 ro,vec3 rd,float dmin,float dmax) {

    float res = 1.0;
    float t = dmin;
    float ph = 1e10;
    
    for(int i = 0; i < 150; i++ ) {
        
        float h = scene(ro + rd * t  ).x;

        float y = h*h / (2. * ph);
        float d = sqrt(h*h-y*y);         
        res = min(res,10.*d/max(0.,t-y)); 
        ph = h;
        t += h;
        
        if(res < 0.0 || t > dmax ) { break; }

        }

        return clamp(res,0.0,1.0);

}

vec3 calcNormal(vec3 p) {

    vec2 e = vec2(1.0,-1.0) * 0.0001;

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

     vec3 vDir = normalize(uv.x * 
     camRight + uv.y * camUp + camForward * fPersp);  

     return vDir;
}

vec3 light(vec3 ro,vec3 rd,vec3 n,vec3 l,vec2 d) { 

vec3 lig_dir = l - rd;
float lig_dist = max(length(lig_dir),0.001 );
lig_dir /= lig_dist;

float at = 1. / (1. + lig_dist * .2 + lig_dist *lig_dist * .1 );  
float dif = max(dot(n,lig_dir),0.0);

float spe = pow(max(dot(reflect(-lig_dir,n),-rd),0.),8.);  

vec3 col = vec3(1.);

if(d.y == 1.) {
col = vec3(.5,.4,.5);
}

if(d.y == 2.) {
col = vec3(1.,.5,1.);
}

col = (col * (dif + 100.  ) + vec3(1.  ) * spe * 2.) * at;
 
return col;

}

void main() {
 
vec3 color = vec3(0.);

vec3 cam_target = vec3(0.0);
vec3 cam_pos = vec3(2.,1.,2.);

vec3 lig_pos = vec3(5.,10.,10.);

vec2 uv = -1. + 2. * gl_FragCoord.xy / resolution.xy ;
uv.x *= resolution.x / resolution.y; 

vec3 direction = rayCamDir(uv,cam_pos,cam_target,1.); 
vec2 d = rayScene(cam_pos,direction);

cam_pos += direction * d.x;
vec3 n = calcNormal(cam_pos);
       
color = light(cam_pos,direction,n,lig_pos,d);
       
float sh = shadow(cam_pos,normalize(lig_pos),0.005,250.);
direction = reflect(direction,n);
d.x += reflection(cam_pos + direction  * 0.01,direction,0.005,25.);
       
cam_pos += direction  * d.x;
n = calcNormal(cam_pos);
         
color += light(cam_pos,direction,n,lig_pos,d) * .25  ;
color *= sh + d.x * .005;

color = pow(color,vec3(.4545)); 
gl_FragColor = vec4(color,1.0);  

}
