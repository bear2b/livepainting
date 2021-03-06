

var scene = new THREE.Scene();
scene.background = new THREE.Color('blue');
var camera = new THREE.PerspectiveCamera(75,window.innerWidth/window.innerHeight,0.1,1000 );
var renderer = new THREE.WebGLRenderer({antialias: true});

renderer.setSize(window.innerWidth,window.innerHeight);
document.body.appendChild(renderer.domElement);
var geometry = new THREE.BoxGeometry();
var material = new THREE.MeshStandardMaterial( { color: 0x00ff00 } );
var cube = new THREE.Mesh( geometry, material );
//scene.add(new THREE.AmbientLight(0xcccccc));
const light = new THREE.DirectionalLight( 0xffffff, 3.0 );
light.position.set(1,-1,0);
scene.add(light);
//scene.add(cube);
camera.position.set( 1, 0, 0 );
var loader = new THREE.GLTFLoader();
//loader.load("model/singe.glb",function(objet){
  //scene.add(object);
//})
loader.load('./model/singe.glb',

    // called when the resource is loaded
    function (gltf) {
      console.log(gltf);
      var scale =0.1;
      const singe = gltf.scene;

      //singe.rotation.set(0,1.5708,0)
      singe.position.set(0,0,0);
      singe.scale.set(scale,scale,scale);

      singe.traverse( ( child ) => {
    if ( child.name== "Cube004" ) {
      console.log("Singe trouvé!");
      let texture = new THREE.TextureLoader().load("./model/singe.jpg");
      child.rotation.set(-1.57,0,0)
      child.position.set(0,0,0);
      child.scale.set(scale,scale,scale);
      texture.encoding = THREE.sRGBEncoding;
      //texture.flipY = false;
      //texture.anisotropy = 16;
      child.material=  new THREE.MeshPhongMaterial({color: 0xff0000,emissive: 0x2a0000,shininess: 10,specular: 0xffffff});
      child.material = new THREE.MeshBasicMaterial( { map: texture } );
      child.material.needsUpdate = true;


     //console.log(child.material)
     scene.add(child);
   }
 });



});
var controls = new THREE.OrbitControls(camera, renderer.domElement);



camera.position.z = 5;

var update = function()
{
  //cube.rotation.x += 0.01;
//  cube.rotation.y += 0.01;
//  cube.rotation.z += 0.01;
  controls.update();
};
//drawscene
var render = function()
{
renderer.render(scene,camera);
};
//run GameLoop
var GameLoop = function()
{
  requestAnimationFrame(GameLoop);
  update();
  render();
};
function onWindowResize() {

  // set the aspect ratio to match the new browser window aspect ratio
  camera.aspect = window.innerWidth / window.innerHeight;


  // update the camera's frustum
  camera.updateProjectionMatrix();

  // update the size of the renderer AND the canvas
  renderer.setSize( window.innerWidth , window.innerHeight );

}

window.addEventListener( 'resize', onWindowResize );

GameLoop();
