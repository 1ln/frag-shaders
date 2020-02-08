uniform mat4 modelViewProjectionMatrix;
varying vec4 position;

void main() {

gl_Position = modelViewProjectionMatrix * position;

}
