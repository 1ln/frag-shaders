#version 330 core     

//Dan Olson

out vec4 FragColor;

uniform vec2 resolution;
uniform float time;
uniform int seed;

float eps = 0.0001;

float h11(float p) {
    uvec2 n = uint(int(p)) * uvec2(uint(int(seed)),2531151992.0);
    uint h = (n.x ^ n.y) * uint(int(seed));
    return float(h) * (1./float(0xffffffffU));
}

float n2(vec2 x) { 

    vec2 p = floor(x);
    vec2 f = fract(x);

    f = f * f * (3.0 - 2.0 * f);  
    float n = p.x + p.y * 57.;  

    return mix(mix(h11(n+0.),h11(n+1.),f.x),
               mix(h11(n+57.),h11(n+58.),f.x),f.y);  
}

float f2(vec2 x) {

    float s = 0.;
    float h = exp2(-.5);     
    float f = 1.;
    float a = 0.5;

    for(int i = 1; i < 6; i++) {
 
        s += a * n2(f * x);
        f *= 2.;
        a *= h;
    }    

    return s;
}

float dd(vec2 p) {

   vec2 q = vec2(f2(p+vec2(3.,0.5)),
                 f2(p+vec2(1.,2.5)));

   vec2 r = vec2(f2(p + 4. * q + vec2(7.5,4.35)),
                 f2(p + 4. * q + vec2(5.6,2.2))); 

   return f2(p + 4. * r);
}

mat2 rot2(float a) {

    float c = cos(a);
    float s = sin(a);
    
    return mat2(c,-s,s,c);
}

mat3 camOrthographic(vec3 ro,vec3 ta,float r) {
     
     vec3 w = normalize(ta - ro); 
     vec3 p = vec3(sin(r),cos(r),0.);           
     vec3 u = normalize(cross(w,p)); 
     vec3 v = normalize(cross(u,w));

     return mat3(u,v,w); 
} 

vec3 repLimit(vec3 p,float c,vec3 l) {
  
    vec3 q = p - c * clamp( floor((p/c)+0.5) ,-l,l);
    return q; 
}

vec2 opu(vec2 d1,vec2 d2) {
    return (d1.x < d2.x) ? d1 : d2;
} 

float box(vec3 p,vec3 b) {

    vec3 d = abs(p) - b;
    return length(max(d,0.0)) + min(max(d.x,max(d.y,d.z)),0.0);
}

vec2 scene(vec3 p) {

    vec2 res = vec2(1.,0.);

    float d = 0.;     
    float t = time*.05;  
    
    vec3 q = p;
    vec3 pl = p; 

    p.xz *= rot2(t);
    p.zy *= rot2(t);

    float scale = 0.05;
    q = repLimit(p/scale,4.,vec3(5.))*scale;  
    d = box(q/scale,vec3(1.))*scale;
    res = opu(res,vec2(d,0.25)); 

    res = opu(res,
    vec2(pl.y+1.,0.45));
    return res;

}

vec2 rayScene(vec3 ro,vec3 rd) {
    
    float d = -1.0;
    float s = 0.;
    float e = 100.;  

    for(int i = 0; i < 125; i++) {

        vec3 p = ro + s * rd;
        vec2 dist = scene(p);
   
        if(abs(dist.x) < eps || e <  dist.x ) { break; }
        s += dist.x;
        d = dist.y;

        }
 
        if(e < s) { d = -1.0; }
        return vec2(s,d);

}

float shadow(vec3 ro,vec3 rd) {

    float res = 1.0;
    float t = 0.005;
    float ph = 1e10;
    
    for(int i = 0; i < 125; i++ ) {
        
        float h = scene(ro + rd * t  ).x;

        float y = h * h / (2. * ph);
        float d = sqrt(h*h-y*y);         
        res = min(res,255. * d/max(0.,t-y));
        ph = h;
        t += h;
    
        if(res < eps || t > 12.) { break; }

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

vec3 col = vec3(1.);

if(d.y >= 0.) {

vec3 p = ro + rd * d.x;
vec3 n = calcNormal(p);
vec3 l = normalize(vec3(10.));
vec3 h = normalize(l - rd);
vec3 r = reflect(rd,n);

if(d.y == 0.25) {
    col = vec3(1.,0.,0.);
    col += dd(p.xz);
} else {
    col = vec3(.5);
}

float amb = clamp(0.5 + 0.5 * n.y,0.,1.);
float dif = clamp(dot(n,l),0.0,1.0);
float spe = pow(clamp(dot(n,h),0.0,1.0),16.)
* dif * (.04 + 0.9 * pow(clamp(1. + dot(h,rd),0.,1.),5.));

float fre = pow(clamp(1. + dot(n,rd),0.0,1.0),2.0);
float ref = smoothstep(-.2,.2,r.y);

vec3 linear = vec3(0.);

dif *= shadow(p,l);
ref *= shadow(p,r);

linear += dif * vec3(1.);
linear += amb * vec3(0.002,0.01,0.05);
linear += fre * vec3(.05);

col = col * linear;
col += spe * vec3(1.); 
}

return col;
}

void main() {
 
vec3 color = vec3(0.);
vec3 ro = vec3(0.,3.,.5)*2.;
vec3 ta = vec3(0.0);

const int aa = 2;
for(int k = 0; k < aa; ++k) {
    for(int l = 0; l < aa; ++l) {

    vec2 o = vec2(float(l),float(k)) / float(aa) - .5;

    vec2 uv = (2. * (gl_FragCoord.xy + o) -
    resolution.xy) / resolution.y; 

    mat3 cm = camOrthographic(ro,ta,0.);
    vec3 rd = cm * normalize(vec3(uv.xy,2.));
    vec3 render = renderScene(ro,rd);    

    color += render;
    }

color /= float(aa*aa);
color = pow(color,vec3(.4545));

FragColor = vec4(color,1.0);
}

}
