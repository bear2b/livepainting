

var scene = new THREE.Scene();
scene.background = new THREE.Color('blue');
var camera = new THREE.PerspectiveCamera(75,window.innerWidth/window.innerHeight,0.1,1000 );
var renderer = new THREE.WebGLRenderer({antialias: true});

renderer.setSize(window.innerWidth,window.innerHeight);
document.body.appendChild(renderer.domElement);
var geometry = new THREE.BoxGeometry();
var material = new THREE.MeshStandardMaterial( { color: 0x00ff00 } );
var cube = new THREE.Mesh( geometry, material );
var alight=new THREE.AmbientLight(0x404040,0.1);
window.light=alight;
scene.add(alight);
const light = new THREE.DirectionalLight( 0xffffff, 1.0 );
light.position.set(1,-1,0);
scene.add(light); 

const light2 = new THREE.DirectionalLight( 0xffffff, 1.0 );
light2.position.set(1,1,0);
scene.add(light2); 
//scene.add(cube);
camera.position.set( 1, 0, 0 );
var turtle, mixer;
var loader = new THREE.GLTFLoader();
var clock = new THREE.Clock();
var newSkin 

loader.load('./model/turtle.glb',

    // called when the resource is loaded
    function (gltf) {
      console.log(gltf);
      var scale =3;
      const turtleScene = gltf.scene;
      mixer= new THREE.AnimationMixer(turtleScene);
      gltf.animations.forEach((clip) => {mixer.clipAction(clip).play(); });

      
      turtleScene.position.set(0,0,0);
      turtleScene.scale.set(scale,scale,scale);
      console.log("loaded turtle!")
      turtleScene.traverse( ( child ) => {
    if ( child.name== "wugui" ) {
      console.log("wugui trouvé!");
      turtle=child.children[0];
      scene.add(turtle);
      turtle.scale.set(0.1,0.1,0.1)
      turtle.rotation.set(0,3.14,0)
      window.scene=scene;
   }
 });



},function ( xhr ) {
  console.log( (xhr.loaded / xhr.total * 100) + '% loaded' );
},

// onError callback
function ( err ) {
  console.error( 'An error happened' +err);
}

);

function changeSkin(){
  var textureLoader = new THREE.TextureLoader();
	textureLoader.load( "./model/turtle.png", function ( map ) {
    console.log("new texture")
    map.flipY=false;
    turtle.traverse(child=>{
      
      if(child.animations!=undefined){
        console.log(child.name)
      }
      if(child.material!=undefined && child.material.map != undefined){
        child.material.map = map;
        child.material.needsUpdate = true;
      }
    })
		
    
	});
}


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
  var delta = clock.getDelta();
  if(mixer!=undefined){
    mixer.update( delta )
  }

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
