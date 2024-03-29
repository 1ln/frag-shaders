#version 150

// dolson,2019

//precision mediump float;

out vec4 out_FragColor; 

uniform vec2 u_res;
uniform float u_time;
uniform vec3 u_cam_pos;

//uniform sampler2D u_noise_tex;

const float E   =  2.7182818;
const float PI  =  radians(180.0); 
const float PHI =  (1.0 + sqrt(5.0)) / 2.0;

//15551*89491 = 1391674541
float hash(float p) {
    uvec2 n = uint(int(p)) * uvec2(1391674541U,2531151992.0);
    uint h = (n.x ^ n.y) * 1391674541U;
    return float(h) * (1.0/float(0xffffffffU));
}

vec3 hash3(vec3 p) {
   uvec3 h = uvec3(ivec3(  p)) *  uvec3(1391674541U,2531151992.0,2860486313U);
   h = (h.x ^ h.y ^ h.z) * uvec3(1391674541U,2531151992U,2860486313U);
   return vec3(h) * (1.0/float(0xffffffffU));

}
 
float cell(vec3 x,float iterations,int type) {
 
    x *= iterations;

    vec3 p = floor(x);
    vec3 f = fract(x);
 
    float min_dist = 1.0;
    
    for(int i = -1; i <= 1; i++) {
        for(int j = -1; j <= 1; j++) {
            for(int k = -1; k <= 1; k++) { 

                vec3 b = vec3(float(k),float(j),float(i));
                vec3 r = hash3( p + b );
                
                vec3 diff = (b + r - f);

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


float fractal(vec3 x,int octaves,float h) {

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

float sin3(vec3 p,float h) {
    
    return sin(p.x*h) * sin(p.y*h) * sin(p.z*h);
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

vec3 fmCol(float t,vec3 a,vec3 b,vec3 c,vec3 d) {
    
    return a + b * cos( (PI*2.0) * (c * t + d));
}

float easeIn4(float t) {
    return t * t;
}

float easeOut4(float t) {
    return -1.0 * t * (t - 2.0);
}

float easeInOut4(float t) {
    if((t *= 2.0) < 1.0) {
        return 0.5 * t * t;
    } else {
        return -0.5 * ((t - 1.0) * (t - 3.0) - 1.0);
    }
}

float easeIn3(float t) {
    return t * t * t;
}

float easeOut3(float t) {
    return (t = t - 1.0) * t * t + 1.0;
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

mat4 rotAxis(vec3 axis,float theta) {

axis = normalize(axis);

    float c = cos(theta);
    float s = sin(theta);

    float oc = 1.0 - c;

    return mat4( 
        oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, 0.0,
        oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, 0.0,
        oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c, 0.0,
        0.0,0.0,0.0,1.0);

}

mat4 translate(vec3 p) {
 
    return mat4(
        vec4(1,0,0,p.x),
        vec4(0,1,0,p.y),
        vec4(0,0,1,p.z),
        vec4(0,0,0,1)  
);
}

vec3 repeatLimit(vec3 p,float c,vec3 l) {
  
    vec3 q = p - c * clamp( floor((p/c)+0.5) ,-l,l);
    return q; 
}

vec3 repeat(vec3 p,vec3 s) {
   
    vec3 q = mod(p,s) - 0.5 * s;
    return q;
}

vec3 pmod(inout vec3 p,vec3 s) {
    vec3 c = floor(( p + s * 0.5 ) / s );
    p = mod(p + s * 0.5,s ) - s * 0.5;
    return c;
}

float refx(vec3 p) {
    
   return p.x = abs(p.x);
} 

vec2 opu(vec2 d1,vec2 d2) {

    return (d1.x < d2.x) ? d1 : d2;
} 

float opu(float d1,float d2) {
    
    return min(d1,d2);
}

float opi(float d1,float d2) {

    return max(d1,d2);
}

float opd(float d1,float d2) {

    return max(-d1,d2);
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

float layer(float d,float h) {

    return abs(d) - h;
}

float sphere(vec3 p,float r) { 
     
    return length(p) - r;
}

float nsphere(vec3 p,float r) {

    return abs(length(p)-r);
}

float ellipsoid(vec3 p,vec3 r) {

    float k0 = length(p/r); 
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float cone(vec3 p,vec2 c) {

    float q = length(p.xy);
    return dot(c,vec2(q,p.z));
}

float roundCone(vec3 p,float r1,float r2,float h) {

    vec2 q = vec2(length(vec2(p.x,p.z)),p.y);
    float b = (r1-r2)/h;
    float a = sqrt(1.0 - b*b);
    float k = dot(q,vec2(-b,a));

    if( k < 0.0) return length(q) - r1;
    if( k > a*h) return length(q - vec2(0.0,h)) - r2;

    return dot(q,vec2(a,b)) - r1;
}

float solidAngle(vec3 p,vec2 c,float ra) {
    
    vec2 q = vec2(length(vec2(p.x,p.z)),p.y);
    float l = length(q) - ra;
    float m = length(q - c * clamp(dot(q,c),0.0,ra));
    return max(l,m * sign(c.y * q.x - c.x * q.y));
}

float link(vec3 p,float le,float r1,float r2) {

    vec3 q = vec3(p.x,max(abs(p.y) -le,0.0),p.z);
    return length(vec2(length(q.xy)-r1,q.z)) - r2;
}

float plane(vec3 p,vec4 n) {

    return dot(p,n.xyz) + n.w;
}

float capsule(vec3 p,vec3 a,vec3 b,float r) {

    vec3 pa = p - a;
    vec3 ba = b - a;
    float h = clamp(dot(pa,ba)/dot(ba,ba),0.0,1.0);
    return length(pa - ba * h) - r;
} 

float prism(vec3 p,vec2 h) {

    vec3 q = abs(p);
    return max(q.z - h.y,max(q.x * 0.866025 + p.y * 0.5,-p.y) - h.x * 0.5); 
}

float box(vec3 p,vec3 b) {

    vec3 d = abs(p) - b;
    return length(max(d,0.0)) + min(max(d.x,max(d.y,d.z)),0.0);
}

float roundbox(vec3 p,vec3 b,float r) {

    vec3 q = abs(p) - b;
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}

float torus(vec3 p,vec2 t) {

    vec2 q = vec2(length(vec2(p.x,p.z)) - t.x,p.y);
    return length(q) - t.y; 
}

float cylinder(vec3 p,float h,float r) {
    
    float d = length(vec2(p.x,p.z)) - r;
    d = max(d, -p.y - h);
    d = max(d, p.y - h);
    return d; 
}

float hexPrism(vec3 p,vec2 h) {
 
    const vec3 k = vec3(-0.8660254,0.5,0.57735);
    p = abs(p); 
    p.xy -= 2.0 * min(dot(k.xy,p.xy),0.0) * k.xy;
 
    vec2 d = vec2(length(p.xy - vec2(clamp(p.x,-k.z * h.x,k.z * h.x),h.x)) * sign(p.y-h.x),p.z-h.y);
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
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

float crossbox(vec3 p,float l,float d) {

    float b0 = box(p.xyz,vec3(l,d,d));
    float b1 = box(p.yzx,vec3(d,l,d));
    float b2 = box(p.zxy,vec3(d,d,l));
    
    return min(b0,min(b1,b2));
} 

vec2 scene(vec3 p) { 

float s  = 0.0001;
float t  = u_time;

vec2 res = vec2(1.0,0.0);

vec3 q = vec3(p);

res = opu(res,vec2(plane(p +vec3(0.,0.,1.),vec4(0.,0.,1.,0.)) ,1.)); 
res = opu(res,vec2(smod(cylinder(p,1.,.5),box(p,vec3(1.)),.5),2.)   );

res = opu(res,vec2(1.,0.));

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


     vec3 vDir = normalize(uv.x * camRight + uv.y * camUp + camForward * fPersp);  

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
col = vec3(1.);
}

if(d.y == 2.) {
col = vec3(1.,0.,0.);
}

col = (col * (dif + 100.  ) + vec3(1.  ) * spe * 2.) * at;
 
return col;

}

void main() {
 
vec3 out_color = vec3(0.);
int aa = 1;

vec3 cam_target = vec3(0.0);
vec3 cam_pos = vec3(.0);
cam_pos = vec3(0.,-5.,.50);
//cam_pos = u_cam_pos;

vec3 lig_pos = vec3(5.,10.,10.);

for(int k = 0; k < aa; k++ ) {

   for(int l = 0; l < aa; l++) {
   
       vec2 o = vec2(float(k),float(l)) / float(aa) * .5;

       vec2 uv = (gl_FragCoord.xy+o) / u_res.xy;
       uv = uv * 2. - vec2(1.); 
       uv.x *= u_res.x/u_res.y; 

       vec3 direction = rayCamDir(uv,cam_pos,cam_target,1.); 

       vec2 d = rayScene(cam_pos,direction);

       cam_pos += direction * d.x;
       vec3 n = calcNormal(cam_pos);
       
       vec3 color = light(cam_pos,direction,n,lig_pos,d);
       
       float sh = shadow(cam_pos,normalize(lig_pos),0.005,250.);
       direction = reflect(direction,n);
       d.x += reflection(cam_pos + direction  * 0.01,direction,0.005,25.);
       
       cam_pos += direction  * d.x;
       n = calcNormal(cam_pos);
         
       color += light(cam_pos,direction,n,lig_pos,d) * .25  ;
       color *= sh + d.x * .005;

       color = pow(color,vec3(.4545)); 
       out_color += color;

   }

   out_color /= float(aa*aa);      
   out_FragColor = vec4(out_color,1.0);
   
   }

}
