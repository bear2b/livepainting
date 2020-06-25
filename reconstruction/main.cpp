#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/imgcodecs.hpp"
#include <opencv2/opencv.hpp>
#include "maillage.hpp"
#include "LaplacianMesh.hpp"
#include "eigen3/Eigen/Eigen"
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include "Visualization.hpp"
#include "DefinedMacros.hpp"
using namespace arma;


std::string datafolder = "./data/";



int main() {

//chargement camera

mat camintr(3,3);
std::ifstream fcamintr("./data/cam.intr"); //"./im_corners.txt"};""
if(fcamintr.is_open()){
    float val;
       std::cout<<"fichier camera intr ouvert"<<std::endl;
    for( int i= 0; i<3; i++ ){
        for (int j=0; j<3;j++){
            fcamintr>>val;
            camintr(i,j)= val;
        }
    }
}
else{
    std::cout<<"erreur"<<std::endl;
}


std::cout<<"Matrice Camera intrasèque"<< std::endl;
std::cout<<camintr<<std::endl;

mat camext(3,4);
std::ifstream fcamext("./data/cam.ext"); //"./im_corners.txt"};""
if(fcamintr.is_open()){

    float val;
       std::cout<<"fichier camera extr ouvert"<<std::endl;
    for( int i= 0; i<3; i++ ){
        for (int j=0; j<4;j++){
            fcamext>>val;
            camext(i,j)= val;
        }
    }
}
else{
    std::cout<<"erreur"<<std::endl;
}
std::cout<<"Matrice Camera extr"<< std::endl;
std::cout<<camext<<std::endl;
 Camera camera = Camera(camintr,camext);

//initialisation mesh de reference
//initialisation maillage de référence
arma::mat ptsmesh;
arma::umat ptsconnect;

//chargement coordonnées point
std::ifstream fpts("./data/mesh.pts");
if(fpts.is_open()){
    double val;
    //std::cout<<"fichier coordonées ouvert"<<std::endl;
    //lecture du nombre de points
    int sizepts=0;
    char ligne[50];
    while(fpts.getline(ligne,50)){
        sizepts +=1;
    }
    //remontée au début du fichier

    fpts.clear();
    fpts.seekg(1);
    ptsmesh = mat(sizepts,3);
    for(int i=0;i<sizepts;i++){
        for( int j= 0; j<3; j++ ){
            fpts>>val;
            ptsmesh(i,j)=val;
            }
    }
    }

std::cout<<ptsmesh<<std::endl;

//chargement connectivité
std::ifstream ftri("./data/mesh.tri");
if(ftri.is_open()){
    int val;
    //std::cout<<"fichier connectivité ouvert"<<std::endl;
    //lecture du nombre de points
    int sizetri=0;
    char ligne[50];
    while(ftri.getline(ligne,50)){
        sizetri +=1;
    }
    std::cout<<sizetri<<std::endl;
    ftri.clear();
    ftri.seekg(0);
    ptsconnect = umat(sizetri,3);
    for(int i=0;i<sizetri;i++){
        for( int j= 0; j<3; j++ ){
           ftri >>val;
           ptsconnect(i,j)=val;

           }

    }
    }
std::cout<<ptsconnect<<std::endl;



//chargement des points de controle
arma::urowvec ptsctrl;

std::ifstream fctrl("./data/ControlPointIDs.txt");
int val;
std::cout<<"fichier points de controle ouvert"<<std::endl;
//lecture du nombre de points
int sizectrl=0;
while(fctrl>>val){
    sizectrl +=1;
fctrl.clear();
fctrl.seekg(0);
ptsctrl = urowvec(sizectrl);
for(int i=0;i<sizectrl;i++){
       fctrl >>val;
       ptsctrl(i)=val;

       }
}

//initialisation mesh laplacien
LaplacianMesh * mesh =  new LaplacianMesh();
mesh->Load(ptsmesh,ptsconnect);
mesh->SetCtrlPointIDs(ptsctrl);
//std::cout<<mesh->GetFacets()<<std::endl;
//passage des coordonnées du mesh dans le ref de la caméra
//mesh->TransformToCameraCoord(camera);
//calcul des matrices A et P
mesh->ComputeAPMatrices();
mesh->computeFacetNormalsNCentroids();

//initialisation du mesh résultat
LaplacianMesh resmesh= *mesh;

// FIN FONCTION INIT CODE DE BASE

// Definition de l'image de base
cv::Mat refImg= cv::imread("./data/template_print.jpg");
//Chargement coin de l'image
mat camcoin(4,3);
std::ifstream fcoin("./data/im_corners.txt");
if(fcoin.is_open()){
    float val;
       std::cout<<"fichier coins ouvert"<<std::endl;
    for( int i= 0; i<4; i++ ){
        for (int j=0; j<3;j++){
            fcoin>>val;
            camcoin(i,j)= val;
        }
    }
}
int pad=0;
cv::Point topLeft    (camcoin(0,0) - pad, camcoin(0,1) - pad);
cv::Point topRight   (camcoin(1,0) + pad, camcoin(1,1) - pad);
cv::Point bottomRight(camcoin(2,0) + pad, camcoin(2,1) + pad);
cv::Point bottomLeft (camcoin(3,0) - pad, camcoin(3,1) + pad);
//
// MANQUANT initialisation pointmatcher
//

//fin init
//std::cout<<resmesh.GetEdges()<<std::endl;
//visualisation du maillage sur l'image de référence
 Visualization::DrawAQuadrangle(refImg, topLeft, topRight, bottomRight, bottomLeft, TPL_COLOR);
Visualization::DrawProjectedMesh(refImg, resmesh, camera, CV_RGB(255,0,0));
cv::imshow("Reference image", refImg);
cv::waitKey(0);




}





