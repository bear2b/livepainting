var canvas=document.createElement("canvas");
document.body.appendChild(canvas)
canvas.setAttribute("id","hello")
var _scene = new THREE.Scene();
var camera = new THREE.PerspectiveCamera(75,window.innerWidth/window.innerHeight,0.1,1000);
var renderer = new THREE.WebGLRenderer({canvas:canvas});
renderer.setSize(window.innerWidth,window.innerHeight);
document.body.appendChild(renderer.domElement);


window.addEventListener('resize',function(){
  var width = window.innerWidth;
  var height = window.innerHeight;
  renderer.setSize(width,height);
  camera.aspect = width/height;
  camera.updateProjectionMatrix();
})

//creates cube
// Create the Geometry passing the size
var geometry = new THREE.BoxGeometry( 1, 1, 1 );
// Create the Material passing the color
var material = new THREE.MeshBasicMaterial( { color: "#433F81" } );
// Create the Mesh
var cube = new THREE.Mesh( geometry, material );
_scene.add( cube );
_scene.add(new THREE.AmbientLight(0xcccccc));
const light = new THREE.DirectionalLight( 0xffffff, 5.0 );
light.position.set(10,10,10);
light.target.position.set(cube.position);
_scene.add(light);

camera.position.set(0,5,0);
camera.lookAt(cube.position);

const controls = new THREE.OrbitControls(camera, canvas);
controls.target.set(0, 0, 0);
controls.update();

var update = function()
{
  cube.rotation.x += 0.01;
  cube.rotation.y += 0.01;
  cube.rotation.z += 0.01;
//  controls.update();
};


//drawscene
var render = function()
{
  renderer.render(_scene,camera);
};

//run GameLoop
var GameLoop = function()
{
  requestAnimationFrame(GameLoop);
  update();
  render();
};



GameLoop();
