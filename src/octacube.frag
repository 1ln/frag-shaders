#version 150

// dolson,2019

out vec4 out_FragColor; 

uniform vec2 u_res;
uniform float u_time;

uniform vec2 u_mouse_pos; 
uniform int u_mouse_pressed_left;
uniform int u_mouse_released_left;

const float E   =  2.7182818;
const float PI  =  radians(180.0);
const float PI2 =  PI * 2.; 
const float PHI =  (1.0 + sqrt(5.0)) / 2.0;

//15551*89491 = 1391674541
float hash(float p) {
    uvec2 n = uint(int(p)) * uvec2(1391674541U,2531151992.0);
    uint h = (n.x ^ n.y) * 1391674541U;
    return float(h) * (1.0/float(0xffffffffU));
}

mat4 rotAxis(vec3 axis,float theta) {

axis = normalize(axis);

    float c = cos(theta);
    float s = sin(theta);

    float oc = 1.0 - c;

    return mat4( 
        oc * axis.x * axis.x + c, 
        oc * axis.x * axis.y - axis.z * s, 
        oc * axis.z * axis.x + axis.y * s, 0.0,

        oc * axis.x * axis.y + axis.z * s, 
        oc * axis.y * axis.y + c, 
        oc * axis.y * axis.z - axis.x * s, 0.0,

        oc * axis.z * axis.x - axis.y * s, 
        oc * axis.y * axis.z + axis.x * s, 
        oc * axis.z * axis.z + c, 0.0,
        0.0,0.0,0.0,1.0);

}

vec2 opu(vec2 d1,vec2 d2) {

    return (d1.x < d2.x) ? d1 : d2;
} 

float smou(float d1,float d2,float k) {

    float h = clamp(0.5 + 0.5 * (d2-d1)/k,0.0,1.0);
    return mix(d2,d1,h) - k * h * (1.0 - h);
}

float box(vec3 p,vec3 b) {

    vec3 d = abs(p) - b;
    return length(max(d,0.0)) + min(max(d.x,max(d.y,d.z)),0.0);
}

float octahedron(vec3 p,float s) {

    p = abs(p);

    float m = p.x + p.y + p.z - s;
    vec3 q;

    if(3.0 * p.x < m) {
       q = vec3(p.x,p.y,p.z);  
    } else if(3.0 * p.y < m) {
       q = vec3(p.y,p.z,p.x); 
    } else if(3.0 * p.z < m) { 
       q = vec3(p.z,p.x,p.y);
    } else { 
       return m * 0.57735027;
    }

    float k = clamp(0.5 *(q.z-q.y+s),0.0,s);
    return length(vec3(q.x,q.y-s+k,q.z - k)); 
} 

vec2 scene(vec3 p) { 

vec2 res = vec2(1.0,0.0);

vec2 m = u_mouse_pos / u_res.xy;

mat4 mx = rotAxis(vec3(1.,0.,0.),m.x * PI2); 
mat4 my = rotAxis(vec3(0.,1.,0.),m.y * PI2); 

//if(u_mouse_pressed_left == 1 && u_mouse_released_left != 0) {
p = (vec4(p,1.) * mx * my).xyz; 
//}

res = opu(res,vec2(smou(
    octahedron(p,1.),
    box(p,vec3(.5 )),
    .01),2.));

return res;

}

vec2 rayScene(vec3 ro,vec3 rd) {
    
    float depth = 0.0;
    float d = -1.0;

    for(int i = 0; i < 100; i++) {

        vec3 p = ro + depth * rd;
        vec2 dist = scene(p);
   
        if(abs( dist.x) < 0.001 || 500. <  dist.x ) { break; }
        depth += dist.x;
        d = dist.y;

        }
 
        if(500. < depth) { d = -1.0; }
        return vec2(depth,d);

}

float shadow(vec3 ro,vec3 rd ) {

    float res = 1.0;
    float t = 0.005;
    float ph = 1e10;
    
    for(int i = 0; i < 35; i++ ) {
        
        float h = scene(ro + rd * t  ).x;

        float y = h * h / (2. * ph);
        float d = sqrt(h*h-y*y);         
        res = min(res,10. * d/max(0.,t-y));
        ph = h;
        t += h;
    
        if(res < 0.001 || t > 5. ) { break; }

        }

        return clamp(res,0.0,1.0);

}

vec3 calcNormal(vec3 p) {

    vec2 e = vec2(1.0,-1.0) * 0.001;

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


     vec3 vDir = normalize(uv.x * camRight + uv.y * camUp + camForward * fPersp);  

     return vDir;
}

vec3 render(vec3 ro,vec3 rd) {

float t = u_time;
 
vec2 d = rayScene(ro, rd);

vec3 cf = vec3(1.);                         
vec3 col = cf - max(rd.y,0.);

if(d.y >= 0.) {

vec3 p = ro + rd * d.x;
vec3 n = calcNormal(p);
vec3 l = normalize( vec3(0.,10.,0.));
vec3 h = normalize(l - rd);
vec3 r = reflect(rd,n);
float amb = sqrt(clamp(0.5 + 0.5 * n.y,0.0,1.0));
float dif = clamp(dot(n,l),0.0,1.0);
float spe = pow(clamp(dot(n,h),0.0,1.0),16.) * dif * (.04 + 0.9 * pow(clamp(1. + dot(h,rd),0.,1.),5.));
float fre = pow(clamp(1. + dot(n,rd),0.0,1.0),2.0);
float ref = smoothstep(-.2,.2,r.y);
vec3 linear = vec3(0.);

dif *= shadow(p,l);
ref *= shadow(p,r);

linear += dif * vec3(0.,1.,0.);
linear += amb * vec3(0.02);
linear += ref * vec3(1.);
linear += fre * vec3(1.);

if(d.y == 1.) {
col += vec3(.5);
}

col = col * linear;
col += 5. * spe;
col = mix(col,cf,1. - exp(-0.0001 * d.x * d.x * d.x));

}

return col;
}

void main() {
 
vec3 out_color = vec3(0.);
int aa = 2;

vec3 cam_target = vec3(0.);
vec3 cam_pos = vec3(.0,.0,2.);

for(int k = 0; k < aa; k++ ) {
   for(int l = 0; l < aa; l++) {
   
       vec2 o = vec2(float(k),float(l)) / float(aa) * .5;

       vec2 uv = (gl_FragCoord.xy+o) / u_res.xy;
       uv = uv * 2. - vec2(1.); 
       uv.x *= u_res.x/u_res.y; 

       vec3 direction = rayCamDir(uv,cam_pos,cam_target,1.); 
       vec3 color = render(cam_pos,direction);
      
       out_color += color;  

   }
    
   out_color /= float(aa*aa);
   out_color = pow(out_color,vec3(.4545));      
   out_FragColor = vec4(out_color,1.0);
 
}

}
