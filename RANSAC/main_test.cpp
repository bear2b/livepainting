
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/imgproc.hpp"
// #include <opencv2/core/eigen.hpp>
#include "eigen3/Eigen/Geometry"
#include "./src/pattern.h"
#define PI 3.141592
#include "ceres/ceres.h"
#include "ceres/rotation.h"
//#include "./src/patterndetector.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace Eigen;
#include "ceres/ceres.h"
#include "glog/logging.h"
typedef Eigen::NumTraits<double> EigenDouble;
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Matrix<double, 3, 3> Mat3;
typedef Eigen::Matrix<double, 2, 1> Vec2;
typedef Eigen::Matrix<double, Eigen::Dynamic,  8> MatX8;
typedef Eigen::Vector3d Vec3;

class SolverOptions {
    // Default settings for homography estimation which should be suitable
    // for a wide range of use cases.
public:
    SolverOptions()
        :  max_num_iterations(500000),
          expected_average_symmetric_distance(1.f) {}
    // Maximal number of iterations for the refinement step.

    int max_num_iterations;
    // Expected average of symmetric geometric distance between
    // actual destination points and original ones transformed by
    // estimated homography matrix.
    //
    // Refinement will finish as soon as average of symmetric
    // geometric distance is less or equal to this value.
    //
    // This distance is measured in the same units as input points are.
    double expected_average_symmetric_distance;
};
//données physique de la page
// transforme la position

//deformation de la feuille
//parametres:
//L: longueur de base de la feuille
//l: abscisse d'où la déformation commence
//def: deformation de la page (def represente le décalage en abscisse entre la deformée et la feuille de référence, qui cause la déformation)
Matrix<double,4,1> transformation2(Matrix<double,3,1> position,double l,double def){
    double x = position[0];
    double z= position[2];
    double newz;
    double newx;
    if(x<l){
        double h= 54/7.*def/l;
        newz = h/(l*l)*x*x*x-2*h/l*x*x+h*x;
        newx = x*(l-def)/l;
    }
    else{
        newz=position[2];
        newx = x-def;
    }
    return Matrix<double,4,1>(newx,position[1],newz,1);
}
float distance(cv::Point2f refpts,cv::Point2f newpts){
    return sqrt(pow((refpts.x-newpts.x),2)+pow((refpts.y-newpts.y),2));

}

//void transformation_inv(Matrix3d position){
//    return position/2;
//}
//calcul de la matrice extr F à partir des rotation Rx, Ry, Rz, et des translation tx,ty,tz
Matrix<double,4,4> calcul_F(double Rx,double Ry, double Rz,double tx,double ty,double tz){
    //init F
    Matrix<double,4,4> F;
    F.row(3).setZero();
    Matrix<double,3,1>  translation;
    translation<<tx,ty,tz;
    double cosX = std::cos(Rx);
    double sinX = std::sin(Rx);
    double cosY= std::cos(Ry);
    double sinY = std::sin(Ry);
    double cosZ = std::cos(Rz);
    double sinZ = std::sin(Rz);
    Matrix<double,3,3>  rotx;
    rotx<<1,0,0,
            0,cosX,-sinX,
            0,sinX,cosX;
    Matrix<double,3,3>  roty;
    roty<<cosY,0,-sinY,
            0,1,0,
            sinY,0,cosY;
    Matrix<double,3,3>  rotz;
    rotz<<cosZ,-sinZ,0,
            sinZ,cosZ,0,
            0,0,1;
    Matrix<double,3,3>  r;
    r= rotx*roty*rotz;
    //ajout rotation totale dans F
    for(int i=0;i<3;i++){
        for(int j=0; j<3;j++){
            F(i,j)=r(i,j);
        }
    }
    //ajout translation dans F
    for(int k=0;k<3;k++){
        F(k,3)=translation(k);
    }
    F(3,3)=1;
    return F;
}
// calcul de la matrice intrinseque
Matrix<double,3,3> calcul_K(double fu,double fv,double cu,double cv){
    Matrix<double,3,3> K;
    K<<fu,0,cu,
            0,fu,cv,
            0,0,1;
    return K;
}
Matrix<double,3,4> from4x4to3x4 (Matrix<double,4,4> F){
    Matrix<double,3,4> mat3x4;
    mat3x4.setZero();
    mat3x4(0,0)=1;
    mat3x4(1,1)=1;
    mat3x4(2,2)=1;
    Matrix<double,3,4> F3x4(3,4);
    F3x4=mat3x4*F;
    return F3x4;
}

Matrix<double,3,4> calcul_P(double Rx,double Ry, double Rz,double tx,double ty,double tz,double fu,double fv,double cu,double cv){
    Matrix<double,4,4> F = calcul_F( Rx, Ry, Rz, tx, ty, tz);
    Matrix<double,3,3> K =calcul_K(fu,fv,cu,cv);
    Matrix<double,3,4> F3x4 =from4x4to3x4(F);
    return K*F3x4;
}
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//                          FONCTIONS POUR LE SOLVER
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
template<typename T>
Matrix<T,3,1>  RotateAndTranslatePoint(const T* rotation, const T* translation, const Matrix<T,4,1> point){
    Matrix<T,3,1>  result;
    //definition K
    Matrix<T,3,3> K;
    T fu = T(424);
    T cu = T(0);
    T cv = T(0);
    K<<fu,T(0),cu,
            T(0),fu,cv,
            T(0),T(0),T(1);
    T cosX = cos(rotation[0]);
    T sinX = sin(rotation[0]);
    T cosY = cos(rotation[1]);
    T sinY = sin(rotation[1]);
    T cosZ = cos(rotation[2]);
    T sinZ = sin(rotation[2]);

    Matrix<T,3,3>  rotx;
    rotx<<T(1),T(0),T(0),
            T(0),T(cosX),-sinX,
            T(0),sinX,cosX;
    Matrix<T,3,3>  roty;
    roty<<cosY,T(0),-sinY,
            T(0),T(1),T(0),
            sinY,T(0),cosY;
    Matrix<T,3,3>  rotz;
    rotz<<cosZ,-sinZ,T(0),
            sinZ,cosZ,T(0),
            T(0),T(0),T(1);
    Matrix<T,3,3>  r;
    r= rotx*roty*rotz;

    //definition F3x4
    Matrix<T,3,4> F;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            F(i,j)=r(i,j);
        }
        F(i,3)=translation[i];
    }
    result=K*F*point;

    // std::cout<<"point trans "<<result<<std::endl;
    return result;
    //    std::cout<<K*F<<std::endl;

    //    std::cout<<"point trans "<<result<<std::endl;
}






Eigen::Matrix<double, 2, 1> transform(Eigen::Matrix<double, 2, 1> trainpt,double rotation[3],double translation[3],double L,double def){
    //    std::cout<<"translation"<<std::endl;
    //    std::cout<<translation[0]<<"   "<<translation[1]<<"   "<<translation[2]<<std::endl;
    //    std::cout<<"rotation"<<std::endl;
    //    std::cout<<rotation[0]<<"   "<<rotation[1]<<"   "<<rotation[2]<<std::endl;
    Matrix<double,3,1> position(trainpt[0],trainpt[1],1);
    //    std::cout<<"l: "<<L<<" "<<"def: "<<def<<std::endl;
    Matrix<double,4,1> pt= transformation2(position,L/2,def);
    Matrix<double,3,4> P = calcul_P(rotation[0],rotation[1],rotation[2],translation[0],translation[1],translation[2],424,424,0,0);
    //std::cout<<"P"<<std::endl;
    //std::cout<<P<<std::endl;
    Eigen::Matrix<double, 3, 1> res3d=P*pt;
    return Eigen::Matrix<double,2,1>(res3d[0]/res3d[2],res3d[1]/res3d[2]);


}
Eigen::Matrix<double, 2, 1> transPoint(const double * param,
                                       const double * x1){//,
    // T backward_error[2]) {
    typedef Eigen::Matrix<double, 3, 1> Vec3;
    double rotation[3];
    double translation[3];
    double L;
    double def;
    for (int i=0;i>3;i++){
        rotation[i]=param[i];
        translation[i]=param[i+3];
    }
    L = param[6];
    def =param[7];
    Eigen::Matrix<double,2,1> x(x1[0],x1[1]);
    Eigen::Matrix<double,2,1> x1_trans =transform(x,rotation,translation,L,def);
    return x1_trans;
}
struct transformationFunctor {
    bool operator()(const double* parametres,const double* x, double * value) const {
        Eigen::Matrix<double, 2, 1> x1_trans;
        //param<<parametres[0],parametres[1],parametres[2],parametres[3],
        //         parametres[4],parametres[5],parametres[6],parametres[7];
        //  x1<<x[0]<<x[1];
        x1_trans= transPoint(parametres,x) ;
        value[0]=x1_trans[0];
        value[1]=x1_trans[1];
        return true;
    }
};
class TestFunctor {
public:
    TestFunctor(const Vec &x,
                const Vec &y)
        : x_(x), y_(y) {  compute_pt_trans.reset(new ceres::CostFunctionToFunctor<2,8,2>(
                                                     new ceres::NumericDiffCostFunction<transformationFunctor,
                                                     ceres::CENTRAL,
                                                     2,
                                                     8,
                                                     2>(
                                                         new transformationFunctor)));}
    template<typename T>
    bool operator()(const T* parametres, T* residuals) const {
        std::cout<<&parametres<<std::endl;
        T * value;
        T * x1;
        x1[0]=T(x_[0]);
        x1[1]=T(x_[1]);

        (*compute_pt_trans)(parametres,x1,value);
        Eigen::Matrix<T, 2,1>xt(value[0],value[1]);
        Eigen::Matrix<T, 2,1> x2(y_(0), y_(1));
        residuals[0]= xt[0]-x2[0];
        residuals[1]= xt[1]-x2[1];
        return true;
    }
    const Vec2 x_;
    const Vec2 y_;
    std::unique_ptr<ceres::CostFunctionToFunctor<2,8,2> > compute_pt_trans;
};


class distanceFunctor {
public:
    distanceFunctor(const Vec &x,
                    const Vec &y,
                    const double &Long)
        : x_(x), y_(y),L_(Long) {}
    template<typename T>
    bool operator()(const T* parametres, T* residuals) const {
        T rotation[3] = {parametres[0],parametres[1],parametres[2]};
        T  translation[3] = {parametres[3],parametres[4],parametres[5]};
        T def =parametres[6];
        T L = T(L_) ;

        //std::cout<< rotation <<std::endl;

        //coordonnées du point de reference
        T x = T(x_[0]);
        T y = T(x_[1]);
        T newz;
        T newx;
        T l=L/T(2);
        if(x<l){
            //std::cout<<"point tranformé"<<std::endl;
            T h= T(54/7.)*def/l;
            newz = h/(l*l)*x*x*x-T(2)*h/l*x*x+h*x;
            newx = x*(l-def)/l;
        }
        else{
            newz=T(0);
            newx = x-def;
        }
        Matrix<T,4,1> point(newx,y,newz,T(1));

        Matrix<T,3,1> transpoint;
        //std::cout<<"point :"<<point<<std::endl;
        //        std::cout<< "parametres initiaux"<<std::endl;
        //        std::cout<< "rotations"<<std::endl;
        //        std::cout<< "Rx: "<<parametres[0]<< " Ry: "<<parametres[1]<< " Rz: "<<parametres[2]<<std::endl;
        //        std::cout<< "Translations"<<std::endl;
        //        std::cout<< "Tx: "<<parametres[3]<< " Ty: "<<parametres[4]<< " Tz: "<<parametres[5]<<std::endl;
        //        std::cout<< "Parametres déformation"<<std::endl;

        //        std::cout<< "Longueur: "<<parametres[6]<<" deformation: "<<parametres[7] <<std::endl;
        transpoint = RotateAndTranslatePoint<T>(rotation,translation,point);
        //std::cout<<"transpoint: "<< transpoint<<std::endl;
        T new_x1[2];
        if(transpoint[2]==T(0)){
            new_x1[0] =transpoint[0];
            new_x1[1]=transpoint[1];

        }
        else{

            new_x1[0] =transpoint[0]/transpoint[2];
            new_x1[1] =transpoint[1]/transpoint[2];

        }
        //std::cout<<"x      :"<< x<<" "<< y<<std::endl;
        //std::cout<<"point  :"<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
        //std::cout<<"pttrans:"<< newx<<" "<<y <<" "<< newz<<std::endl;
        //std::cout<<" "<<std::endl;
        //std::cout<<"pt rot : "<< transpoint[0]<<" "<< transpoint[1]<<std::endl;
        //std::cout<<"new_x  : "<< new_x1[0]<<" "<< new_x1[1]<<std::endl;
        //        T new_x= T(y_[0])-new_x1[0];
        //        T new_y= T(y_[1])-new_x1[1];
        //       residuals[0]=new_x*new_x+new_y*new_y;
        residuals[0]=T(y_[0])-new_x1[0];
        residuals[1]=T(y_[1])-new_x1[1];
        // std::cout<<"residuals: "<<residuals[0]<<std::endl;//" "<<residuals[1]
        // conditions sur les variables
//        if((parametres[6]<T(0))){//||(parametres[6]>L/T(2)))){
//            return false;
//        }
        return true;
    }

    const Vec2 x_;
    const Vec2 y_;
    const double L_;
};





Vec2 testdistance( double L,
                   const Eigen::Matrix<double,1,7>   param,
                   const Eigen::Matrix<double, 2, 1> &x1, //train point
                   const Eigen::Matrix<double, 2, 1> &x2){ //query point){//,
    // T backward_error[2]) {
    // std::cout<<"param"<<std::endl;
    //    std::cout<<param<<std::endl;

    typedef Eigen::Matrix<double, 3, 1> Vec3;
    double rotation[3];
    double translation[3];
    double def;
    //  std::cout<<"rotation:           translation:"<<std::endl;
    for (int i=0;i<3;i++){
        rotation[i]=param[i];
        translation[i]=param[i+3];
        //  std::cout<<rotation[i]<<"                 "<<translation[i]<<std::endl;
    }

    def = param[6];
    // def =param[7];
    Eigen::Matrix<double,2,1> x1_trans =transform(x1,rotation,translation,L,def);
   // std::cout<<"point transformé  : "<<x1_trans<<std::endl;
    // distance[0] =x2(0)-x1_trans(0);
    //distance[1]= x2(1)-x1_trans(1);
    //std::cout<<x2-x1_trans<<std::endl;
    return x2-x1_trans;
}
double Distancetransfo( double L,const Eigen::Matrix<double, 1, 7> &param,const  Eigen::Matrix<double, 2, 1> &x1,
                        const  Eigen::Matrix<double, 2, 1> &x2) {
    Vec2 distance;
    distance =testdistance(L,param,x1,x2);
    return distance.norm();

}
double Distancetransfo2(  double L, Eigen::Matrix<double, 7, 1> &param,const  Eigen::Matrix<double, 2, 1> &x1,
                          const  Eigen::Matrix<double, 2, 1> &x2) {
    Vec2 distance;
    distance =testdistance(L,param,x1,x2);
    return distance.squaredNorm();
}
class testconditionfinCallback : public ceres::IterationCallback {
public:
    testconditionfinCallback( double L, Eigen::Matrix<double, 7, 1> *param,const Mat &x1, const Mat &x2,const SolverOptions &options): options_(options), x1_(x1), x2_(x2),param_(param),L_(L) {}



    virtual ceres::CallbackReturnType operator()(
            const ceres::IterationSummary& summary) {
        // If the step wasn't successful, there's nothing to do.
        if (!summary.step_is_successful) {
            return ceres::SOLVER_CONTINUE;
        }
        // Calculate average of symmetric geometric distance.
        double average_distance = 0.0;
        for (int i = 0; i < x1_.cols(); i++) {

//            std::cout<<"point d'origine   : "<<x1_.col(i)<<std::endl;
//            std::cout<<"point but         : "<<x2_.col(i)<<std::endl;
            double dist =Distancetransfo2(L_,*param_,
                                          x1_.col(i),
                                          x2_.col(i));

//            std::cout<<"   distance       : "<<dist<<std::endl;



            average_distance += dist;
        }
        average_distance /= x1_.cols();
        //std::cout<<"distance moyenne : "<<average_distance<<std::endl
       // std::cout<<"nombre pts: "<< x1_.cols();
       // std::cout<<"distance moyenne: "<<average_distance<<" "<<"précision limite: "<<1 /*options_.expected_average_symmetric_distance*/ <<std::endl;
        if (average_distance <= 1.f){//options_.expected_average_symmetric_distance) {
            //std::cout<<"distance ok"<<std::endl;
            return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
        }

        return ceres::SOLVER_CONTINUE;
    }
private:
    const SolverOptions &options_;
    const Mat &x1_;
    const Mat &x2_;
    Eigen::Matrix<double, 7,1 >  *param_;
    double L_;
};

bool testmodelfrom2DFromCorrespondences( Mat x1tot, Mat x2tot,const double wT,const Mat &x1,const Mat &x2, double parametresguess[7],const SolverOptions &options,Matrix<double,7,1> * result) {
    assert(2 == x1.rows());
    assert(5 <= x1.cols());
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());
    double L = wT;
    // Step 1: Algebraic homography estimation.
    // Assume algebraic estimation always succeeds.
    //j'etablis le guess

    //definition des parametres de base
    //    std::cout<< "parametres initiaux"<<std::endl;
    //    std::cout<< "rotations"<<std::endl;
    //    std::cout<< "Rx: "<<parametresguess[0]<< " Ry: "<<parametresguess[1]<< " Rz: "<<parametresguess[2]<<std::endl;
    //    std::cout<< "Translations"<<std::endl;
    //    std::cout<< "Tx: "<<parametresguess[3]<< " Ty: "<<parametresguess[4]<< " Tz: "<<parametresguess[5]<<std::endl;
    //    std::cout<< "Parametres déformation"<<std::endl;

    //    std::cout<< "Longueur: "<<parametresguess[6]<<" deformation: "<<parametresguess[7] <<std::endl;

    //    for (int i =0;i<6;i++){
    //        parametresguess[i]+= (std::rand()/(double)RAND_MAX )*parametresguess[i];
    //    }
    *result<<parametresguess[0],parametresguess[1],parametresguess[2],parametresguess[3],parametresguess[4],parametresguess[5],parametresguess[6];
    std::cout<<" "<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<< "parametres déformés"<<std::endl;
    std::cout<< "rotations"<<std::endl;
    std::cout<< "Rx: "<<parametresguess[0]<< " Ry: "<<parametresguess[1]<< " Rz: "<<parametresguess[2]<<std::endl;
    std::cout<< "Translations"<<std::endl;
    std::cout<< "Tx: "<<parametresguess[3]<< " Ty: "<<parametresguess[4]<< " Tz: "<<parametresguess[5]<<std::endl;
    std::cout<< "Parametres déformation"<<std::endl;

    std::cout<< "Longueur: "<<parametresguess[6]<<" deformation: "<<parametresguess[7] <<std::endl;
    std::cout<<*result<<std::endl;
    // Step 2: Refine matrix using Ceres minimizer.
    ceres::Problem problem;
    //for (int i = 0; i <8 ; i++) {
    for (int i = 0; i < x1.cols(); i++) {
        distanceFunctor
                *parametres_cost_function =
                new distanceFunctor(x1.col(i),x2.col(i),wT);
        problem.AddResidualBlock(
                    new ceres::AutoDiffCostFunction<
                    distanceFunctor,
                    2,  // num_residuals
                    7>(parametres_cost_function),
                    NULL,
                    result->data());
    }
    // Configure the solve.
    ceres::Solver::Options solver_options;
    solver_options.linear_solver_type = ceres::DENSE_SCHUR;
    solver_options.minimizer_progress_to_stdout = true;
    solver_options.max_num_iterations = options.max_num_iterations;
    solver_options.update_state_every_iteration = true;
    // Terminate if the average symmetric distance is good enough.
    testconditionfinCallback callback(L, result,x1tot, x2tot, options);
    solver_options.callbacks.push_back(&callback);
    // Run the solve.
    ceres::Solver::Summary summary;

    ceres::Solve( solver_options, &problem, &summary);
    //std::cout<<summary.FullReport()<<std::endl;
    return summary.IsSolutionUsable();

}

void affiche_avec_param( double L,int wT,int hT,Matrix<double,7,1> param){

    cv::Mat Queryimg = cv::imread("./data/nat.png",cv::COLOR_BGR2GRAY);
    double ratioW= 1.*Queryimg.size().width/wT;
    double ratioH= 1.*Queryimg.size().height/hT;
    double cu = 0;//imrefcolor.size().width/2;
    double cv = 0;//imrefcolor.size().height/2;
    double fu = 424;
    double Rx= param[0];
    double Ry= param[1];
    double Rz= param[2];
    double tx= param[3];
    double ty =param[4];
    double tz= param[5];
    double def =param[6];
    Matrix<double,4,4> F = calcul_F(Rx,Ry,Rz,tx,ty,tz);
    Matrix<double,3,3> K = calcul_K(fu,fu,cu,cv);

    //definition des points de départs
    std::vector<Matrix<double,3,1>> originpts;
    for(int i= 0; i<wT; i=i+100){
        for (int j= 0; j<hT;j=j+100){
            Matrix<double,3,1> pts(i,j,0);
            originpts.push_back(pts);
        }
    }

    //passage en 3x4de F
    Matrix<double,3,4> F3x4= from4x4to3x4(F);
    //cv::Mat img(wT,hT/2, CV_8UC3, cv::Scalar(1000,1000, 1000));
    const int nbpts= originpts.size();
    Matrix<double,4,1> pos;
    Matrix<double,4,1> refpos;

    for(int i = 0; i<nbpts; i++){

        //std::cout<<"pt de base : "<<originpts[i]<<std::endl;
        pos=transformation2(originpts[i],L/2,def);
        //  std::cout<<"pos : "<<pos<<std::endl;
        refpos(0)= originpts[i][0];
        refpos(1)=originpts[i][1];
        refpos(2)=originpts[i][2];
        refpos(3)=1;
        Matrix<double,4,1> ref;
        Matrix<double,3,1> res;
        ref = refpos;
        res = K*F3x4*pos;
        cv::Point2d pt(res(0)/res(2),res(1)/res(2));
        //  std::cout<<"Point : "<<pt<<std::endl;
        cv::Point2d refpt(ref(0),ref(1));//(ref(0)/ref(2),ref(1)/ref(2));*ratioW
        // std::cout<<refpt<<std::endl;
        cv::drawMarker(Queryimg,pt,cv::Scalar(1000,0,10));

        //cv::drawMarker(Queryimg,refpt,cv::Scalar(0,0,1000));
    }
    //cv::imwrite("testres.jpg",img);
    //cv::imshow("test",trainimg);
    cv::namedWindow("test",cv::WINDOW_NORMAL);
    //cv::resizeWindow("test",wT,hT);
    //cv::destroyAllWindows();
    cv::imshow("test",Queryimg);
    cv::imwrite("./finalmarion1.jpg",Queryimg);
    // cv::waitKey(0);
}

void affiche_avec_param2(double L,int wT,int hT,std::vector<int> inliers,Matrix<double,7,1> &param){

    cv::Mat Queryimg = cv::imread("./data/nat.png",cv::COLOR_BGR2GRAY);
    double ratioW= 1.*Queryimg.size().width/wT;
    double ratioH= 1.*Queryimg.size().height/hT;
    double cu = 0;//imrefcolor.size().width/2;
    double cv = 0;//imrefcolor.size().height/2;
    double fu = 424;
    double Rx= param[0];
    double Ry= param[1];
    double Rz= param[2];
    double tx= param[3];
    double ty =param[4];
    double tz= param[5];
    double def =param[6];
    Matrix<double,4,4> F = calcul_F(Rx,Ry,Rz,tx,ty,tz);
    Matrix<double,3,3> K = calcul_K(fu,fu,cu,cv);

    //definition des points de départs
    std::vector<Matrix<double,3,1>> originpts;
    for(int i= 0; i<wT; i=i+100){
        for (int j= 0; j<hT;j=j+100){
            Matrix<double,3,1> pts(i,j,0);
            originpts.push_back(pts);
        }
    }

    //passage en 3x4de F
    Matrix<double,3,4> F3x4= from4x4to3x4(F);
    //cv::Mat img(wT,hT/2, CV_8UC3, cv::Scalar(1000,1000, 1000));
    const int nbpts= originpts.size();
    Matrix<double,4,1> pos;
    Matrix<double,4,1> refpos;

    for(int i = 0; i<nbpts; i++){

        //std::cout<<"pt de base : "<<originpts[i]<<std::endl;
        pos=transformation2(originpts[i],L/2,def);
        //  std::cout<<"pos : "<<pos<<std::endl;
        refpos(0)= originpts[i][0];
        refpos(1)=originpts[i][1];
        refpos(2)=originpts[i][2];
        refpos(3)=1;
        Matrix<double,4,1> ref;
        Matrix<double,3,1> res;
        ref = refpos;
        res = K*F3x4*pos;
        cv::Point2d pt(res(0)/res(2),res(1)/res(2));
        cv::Point2d refpt(ref(0),ref(1));//(ref(0)/ref(2),ref(1)/ref(2));*ratioW

        cv::drawMarker(Queryimg,pt,cv::Scalar(1000,0,10));

        //cv::drawMarker(Queryimg,refpt,cv::Scalar(0,0,1000));
    }
    //cv::imwrite("testres.jpg",img);
    //cv::imshow("test",trainimg);
    cv::namedWindow("test",cv::WINDOW_NORMAL);
    //cv::resizeWindow("test",wT,hT);
    //cv::destroyAllWindows();
    cv::imshow("test",Queryimg);
    cv::imwrite("./finalmarion1.jpg",Queryimg);
    cv::waitKey(0);
}
void affiche_match_et_set_avec_param( double L, std::vector<int> inliers,std::vector<cv::Point2f>  ptref, int wT,int hT,Matrix<double,7,1> param){
    //std::cout<<"in"<<std::endl;
    cv::Mat Queryimg = cv::imread("./data/nat.png",cv::COLOR_BGR2GRAY);
    cv::Mat trainimg = cv::imread("./data/test2.jpg",cv::COLOR_BGR2GRAY);
    double cu = 0;//imrefcolor.size().width/2;
    double cv = 0;//imrefcolor.size().height/2;
    double fu = 424;
    double Rx= param[0];
    double Ry= param[1];
    double Rz= param[2];
    double tx= param[3];
    double ty =param[4];
    double tz= param[5];
    double def =param[6];
    Matrix<double,4,4> F = calcul_F(Rx,Ry,Rz,tx,ty,tz);
    Matrix<double,3,3> K = calcul_K(fu,fu,cu,cv);
    Matrix<double,3,4> F3x4= from4x4to3x4(F);
    Matrix<double,4,1> pos;
    Matrix<double,4,1> refpos;
    //definition des points de l'image de ref correspondant aux inliers
    std::vector<Matrix<double,3,1>> originpts;
    for(int i= 0; i<inliers.size(); i++){
        Matrix<double,3,1> pts(ptref[inliers[i]].x,ptref[inliers[i]].y,0);
        originpts.push_back(pts);
    }
    std::vector<Matrix<double,3,1>> setoriginpts;
    for(int i= 0; i<=wT; i=i+100){
        for (int j= 0; j<=hT;j=j+100){
            Matrix<double,3,1> pts(i,j,0);
            setoriginpts.push_back(pts);
        }
    }
    //affichage inliers
    for(int i = 0; i<originpts.size(); i++){
        //std::cout<<"pt de base : "<<originpts[i]<<std::endl;
        pos=transformation2(originpts[i],L/2,def);
        //  std::cout<<"pos : "<<pos<<std::endl;
        refpos(0)= originpts[i][0];
        refpos(1)=originpts[i][1];
        refpos(2)=originpts[i][2];
        refpos(3)=1;
        Matrix<double,4,1> ref;
        Matrix<double,3,1> res;
        ref = refpos;
        res = K*F3x4*pos;
        cv::Point2d pt(res(0)/res(2),res(1)/res(2));//
        cv::Point2d refpt(ref(0),ref(1));//(ref(0)/ref(2),ref(1)/ref(2));

        // std::cout<<refpt<<std::endl;

        cv::drawMarker(Queryimg,pt,cv::Scalar(0,0,1000));//pt2[inliers[i]]

        cv::drawMarker(trainimg,refpt,cv::Scalar(1000,0,0));



        //        std::cout<<"point no: "<<inliers[i]<<std::endl;
        //        std::cout<<"pt x: "<<pt2[inliers[i]].x<<" y:"<<pt2[inliers[i]].y<<std::endl;

        //        std::cout<<"point reprojeté: "<<std::endl;
        //        std::cout<<"pt x: "<<pt.x<<" y:"<<pt.y<<std::endl;
        //        std::cout<<"distance        : "<<std::sqrt((pt2[inliers[i]].x-pt.x)*(pt2[inliers[i]].x-pt.x)+(pt2[inliers[i]].y-pt.y)*(pt2[inliers[i]].y-pt.y))<<std::endl;
        //        std::cout<<"distance solveur: "<<inliers_dist[i]<<std::endl;
        //        std::cout<<" "<<std::endl;

    }
    for(int i = 0; i<setoriginpts.size(); i++){

        //std::cout<<"pt de base : "<<originpts[i]<<std::endl;
        pos=transformation2(setoriginpts[i],L/2,def);
        //  std::cout<<"pos : "<<pos<<std::endl;
        refpos(0)= setoriginpts[i][0];
        refpos(1)=setoriginpts[i][1];
        refpos(2)=setoriginpts[i][2];
        refpos(3)=1;
        Matrix<double,4,1> ref;
        Matrix<double,3,1> res;
        ref = refpos;
        res = K*F3x4*pos;
        cv::Point2d pt(res(0)/res(2),res(1)/res(2));
        //       std::cout<<"Point : "<<pt<<std::endl;
        cv::Point2d refpt(ref(0),ref(1));//(ref(0)/ref(2),ref(1)/ref(2));*ratioW
        // std::cout<<refpt<<std::endl;
        cv::drawMarker(Queryimg,pt,cv::Scalar(1000,0,0));

        cv::drawMarker(trainimg,refpt,cv::Scalar(0,0,1000));
    }
    cv::imshow("test1",trainimg);
    cv::namedWindow("test2",cv::WINDOW_NORMAL);
    cv::resizeWindow("test2",wT,hT);
    cv::imshow("test2",Queryimg);

    cv::waitKey(0);
    cv::destroyAllWindows();
    cv::imwrite("./data/affichage_set_et_inliers_TRAIN.png",trainimg);
    cv::imwrite("./data/affichage_set_et_inliers_QUERY.png",Queryimg);

}
void affiche_matches_avec_param(  double L,std::vector<int> inliers,std::vector<double> inliers_dist,std::vector<cv::Point2f> pt2,std::vector<cv::Point2f>  ptref, int wT,int hT,Matrix<double,7,1> param){
    //std::cout<<"in"<<std::endl;
    cv::Mat Queryimg = cv::imread("./data/nat.png",cv::COLOR_BGR2GRAY);
    cv::Mat trainimg = cv::imread("./data/test2.jpg",cv::COLOR_BGR2GRAY);
    double ratioW= 1.*Queryimg.size().width/wT;
    double ratioH= 1.*Queryimg.size().height/hT;
    double cu = 0;//imrefcolor.size().width/2;
    double cv = 0;//imrefcolor.size().height/2;
    double fu = 424;
    double Rx= param[0];
    double Ry= param[1];
    double Rz= param[2];
    double tx= param[3];
    double ty =param[4];
    double tz= param[5];

    double def =param[6];
    Matrix<double,4,4> F = calcul_F(Rx,Ry,Rz,tx,ty,tz);
    Matrix<double,3,3> K = calcul_K(fu,fu,cu,cv);
    //definition des points de départs
    std::vector<Matrix<double,3,1>> originpts;

    for(int i= 0; i<inliers.size(); i++){
        Matrix<double,3,1> pts(ptref[inliers[i]].x,ptref[inliers[i]].y,0);
        originpts.push_back(pts);
    }

    //passage en 3x4de F
    Matrix<double,3,4> F3x4= from4x4to3x4(F);
    //cv::Mat img(wT,hT/2, CV_8UC3, cv::Scalar(1000,1000, 1000));
    const int nbpts= originpts.size();
    Matrix<double,4,1> pos;
    Matrix<double,4,1> refpos;


    for(int i = 0; i<nbpts; i++){
        //std::cout<<"pt de base : "<<originpts[i]<<std::endl;
        pos=transformation2(originpts[i],L/2,def);
        //  std::cout<<"pos : "<<pos<<std::endl;
        refpos(0)= originpts[i][0];
        refpos(1)=originpts[i][1];
        refpos(2)=originpts[i][2];
        refpos(3)=1;
        Matrix<double,4,1> ref;
        Matrix<double,3,1> res;
        ref = refpos;
        res = K*F3x4*pos;
        cv::Point2d pt(res(0)/res(2),res(1)/res(2));//
        cv::Point2d refpt(ref(0),ref(1));//(ref(0)/ref(2),ref(1)/ref(2));

        // std::cout<<refpt<<std::endl;

        cv::drawMarker(Queryimg,pt,cv::Scalar(0,0,1000));//pt2[inliers[i]]

        cv::drawMarker(trainimg,refpt,cv::Scalar(0,0,1000));

        cv::imshow("test1",trainimg);
        cv::namedWindow("test2",cv::WINDOW_NORMAL);
        cv::resizeWindow("test2",wT,hT);
        cv::imshow("test2",Queryimg);

        // cv::waitKey(0);
        cv::destroyAllWindows();

        //        std::cout<<"point no: "<<inliers[i]<<std::endl;
        //        std::cout<<"pt x: "<<pt2[inliers[i]].x<<" y:"<<pt2[inliers[i]].y<<std::endl;

        //        std::cout<<"point reprojeté: "<<std::endl;
        //        std::cout<<"pt x: "<<pt.x<<" y:"<<pt.y<<std::endl;
        //        std::cout<<"distance        : "<<std::sqrt((pt2[inliers[i]].x-pt.x)*(pt2[inliers[i]].x-pt.x)+(pt2[inliers[i]].y-pt.y)*(pt2[inliers[i]].y-pt.y))<<std::endl;
        //        std::cout<<"distance solveur: "<<inliers_dist[i]<<std::endl;
        //        std::cout<<" "<<std::endl;

    }
    cv::imwrite("./data/affichageinliers_train.png",trainimg);
    cv::imwrite("./data/affichageinliers_query.png",Queryimg);

}

void affiche_matches_avec_param2( double L, std::vector<int> inliers, std::vector<double> inliers_dist,std::vector<cv::Point2f> pt2,std::vector<cv::Point2f>  ptref, int wT,int hT,Matrix<double,7,1> &param){

    cv::Mat Queryimg = cv::imread("./data/nat.png",cv::COLOR_BGR2GRAY);
    cv::Mat trainimg = cv::imread("./data/test2.jpg",cv::COLOR_BGR2GRAY);
    double ratioW= 1.*Queryimg.size().width/wT;
    double ratioH= 1.*Queryimg.size().height/hT;

    //definition des points de départs
    std::vector<Matrix<double,2,1>> originpts;
    std::vector<Matrix<double,2,1>> originquerypts;
    // std::cout<<inliers.data()<<std::endl;
    for(int i= 0; i<inliers.size(); i++){
        Matrix<double,2,1> pts(ptref[inliers[i]].x,ptref[inliers[i]].y);
        Matrix<double,2,1> ptsquery(pt2[inliers[i]].x,pt2[inliers[i]].y);

        originpts.push_back(pts);
        originquerypts.push_back(ptsquery);
    }

    const int nbpts= originpts.size();
    for(int i = 0; i<nbpts; i++){
        double parametre[7]= {param[0],param[1],param[2],param[3],param[4],param[5],param[6]};
        double point[2]= {originpts[i][0],originpts[i][1]};
        //con mais obligatoire pour utiliser transPoint
        Eigen::Matrix<double, 2,1> res= transPoint(parametre,point);
        Eigen::Matrix<double, 2,1> ref = originpts[i];
        cv::Point2d pt(res(0)*ratioW,res(1)*ratioH);//
        cv::Point2d refpt(ref(0),ref(1));//(ref(0)/ref(2),ref(1)/ref(2));
        Eigen::Matrix<double, 2,1> dist = originquerypts[i]-res;
        // std::cout<<refpt<<std::endl;

        cv::drawMarker(Queryimg,pt,cv::Scalar(0,0,1000));//pt2[inliers[i]]

        cv::drawMarker(trainimg,refpt,cv::Scalar(0,0,1000));

        cv::imshow("test1",trainimg);
        cv::namedWindow("test2",cv::WINDOW_NORMAL);
        cv::resizeWindow("test2",wT,hT);
        cv::imshow("test2",Queryimg);

        cv::waitKey(0);
        cv::destroyAllWindows();

        //        std::cout<<"point no: "<<inliers[i]<<std::endl;
        //        std::cout<<"pt x: "<<pt2[i].x<<" y:"<<pt2[i].y<<std::endl;

        //        std::cout<<"point reprojeté: "<<std::endl;
        //        std::cout<<"pt x: "<<pt.x<<" y:"<<pt.y<<std::endl;
        //        std::cout<<"distance: "<<dist.norm()<<std::endl;
        //        std::cout<<" "<<std::endl;

    }

}
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//                      MAIN
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int nombre_parametre = 7;
int main(int argc, char **argv)
{


    //initialisation variable
    std::string imRef("./data/test2.jpg");     //nom de l'image de reference
    std::string imquery("./data/nat.png");



    //definition des images
    cv::Mat trainimg = cv::imread(imRef,cv::COLOR_BGR2GRAY);
    cv::Mat Queryimg = cv::imread(imquery,cv::COLOR_BGR2GRAY);
    //definition width and height
    int wT=trainimg.size().width;
    int hT=trainimg.size().height;
    float wQ=Queryimg.size().width;
    float hQ=Queryimg.size().height;
    //definition des variables auxiliaires
    int nbit =0; //nombre iteration
    int maxNbInliers = 0;


    //copy calculate pos

    std::cout<<"initialisation"<<std::endl;
    //initialisation F
    double cu = 0;//imrefcolor.size().width/2;
    double cv = 0;//imrefcolor.size().height/2;
    double fu = 424;
    double Rx= 180./180*PI;
    double Ry= -90./180*PI;
    double Rz= 0./180*PI;
    double tx= 300;
    double ty = 300;
    double tz=1000;
    double def =0;

    Matrix<double,4,4> F = calcul_F(Rx,Ry,Rz,tx,ty,tz);

    Matrix<double,3,3> K = calcul_K(fu,fu,cu,cv);

    //definition des points de départs
    std::vector<Matrix<double,3,1>> originpts;
    for(int i= 0; i<trainimg.size().width; i=i+100){
        for (int j= 0; j<trainimg.size().height;j=j+100){
            Matrix<double,3,1> pts(i,j,0);
            originpts.push_back(pts);
        }
    }

    //passage en 3x4de F
    Matrix<double,3,4> F3x4= from4x4to3x4(F);
    std::cout<<K*F3x4<<std::endl;
    std::vector<cv::Point2f> respts;
    cv::Mat img(wT,hT/2, CV_8UC3, cv::Scalar(1000,1000, 1000));
    std::cout<<trainimg.size().width<<" "<< trainimg.size().height<<std::endl;
    const int nbpts= originpts.size();
    std::cout<<nbpts<<std::endl;
    Matrix<double,4,1> pos;
    Matrix<double,4,1> refpos;
    float L =trainimg.size().width;
    //variables solveur
    Matrix<double,2,70> X1,X2;
    Matrix<double,7,1>  result;
    double parametresguess[7]={Rx,Ry,Rz,tx,ty,tz,def};
    SolverOptions options;
    options.expected_average_symmetric_distance = 0.02;


    for(int i = 0; i<nbpts; i++){

        //std::cout<<"pt de base : "<<originpts[i]<<std::endl;
        pos=transformation2(originpts[i],L/2,def);
        //  std::cout<<"pos : "<<pos<<std::endl;
        refpos(0)= originpts[i][0];
        refpos(1)=originpts[i][1];
        refpos(2)=originpts[i][2];
        refpos(3)=1;
        Matrix<double,4,1> ref;
        Matrix<double,3,1> res;
        ref = refpos;
        res = K*F3x4*pos;
        cv::Point2d pt(res(0)/res(2),res(1)/res(2));
        cv::Point2d refpt(ref(0),ref(1));//(ref(0)/ref(2),ref(1)/ref(2));
        X1(0,i)=refpt.x;
        X1(1,i)=refpt.y;
        X2(0,i)=pt.x;
        X2(1,i)=pt.y;
        // std::cout<<refpt<<std::endl;
        cv::drawMarker(img,pt,cv::Scalar(1000,0,0));
        cv::drawMarker(img,refpt,cv::Scalar(0,0,1000));
    }
    //cv::imwrite("testres.jpg",img);
    //cv::imshow("test",trainimg);
    cv::namedWindow("trans_théorique",cv::WINDOW_NORMAL);
    cv::resizeWindow("trans_théorique",wT,hT);
    cv::imshow("trans_théorique",img);
    // cv::waitKey(0);
    //
    //detection et calcul de l'homographie
    cv::Ptr<cv::Feature2D> detector = cv::ORB::create(500);
    std::vector<cv::KeyPoint> keypoints2,kpref;
    cv::Mat descriptors2, descriptors,queryimggray,trainimggray;
    std::vector<cv::DMatch> matches;
    //passage en gris
    cv::cvtColor(trainimg,trainimggray,cv::COLOR_BGR2GRAY);
    //    capture>>imcam;
    cv::cvtColor(Queryimg,queryimggray,cv::COLOR_BGR2GRAY);

    std::cout<<"detection"<< std::endl;
    //detection des points de reference
    detector->detectAndCompute(trainimggray, cv::Mat(), kpref, descriptors);
    detector->detectAndCompute(queryimggray, cv::Mat(), keypoints2, descriptors2);
    std::cout<<"comparaison"<< std::endl;
    // comparaison des points
    cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("BruteForce-Hamming");
    matcher->match(descriptors2,descriptors,matches,cv::Mat());
    //classement des matches les plus serieux
    std::sort(matches.begin(),matches.end());


    //suppression des plus faux
    matches.erase(matches.begin()+matches.size()*0.20f, matches.end());
    //cv::waitKey(0);
    //affichage des matches
    cv::Mat imMatches;
    cv::drawMatches( Queryimg,keypoints2,trainimg,kpref , matches, imMatches);
    // cv::drawMatches(imcam, keypoints2, patref->data,patref->keypoints, matches, imMatches);
    cv::imwrite("./matchesnew.jpg", imMatches);

    std::vector<cv::Point2f> pt2, ptref;
    for(int i = 0; i<matches.size();i++){
        pt2.push_back( keypoints2[ matches[i].queryIdx ].pt );
        ptref.push_back( kpref[ matches[i].trainIdx ].pt );
    }
    cv::Mat imred;
    cv::Mat homo;
    //calcul de l'homographie via openCV
    homo = cv::findHomography( pt2, ptref, cv::RANSAC );
    std::cout<<"calcul homographie ok"<<std::endl;
    int nb_inliers=0;
    Matrix<double,7,1> tempresult;

    //int ptrefsize= 70;
    int nb_iteration=100;
    srand(time(NULL));
    std::vector<int> inliers;
    std::vector<int> inliers_temp;
    std::vector<double> inliers_dist;
    std::vector<double> inliers_dist_temp;
    int nb_ptsref=7;
   const int nb_ptmodel=X1.cols();//ptref.size();//X1.cols();
    Mat x1(2, nb_ptsref);
    Mat x2(2, nb_ptsref);



    //      int it_high_match[10];
    //      //selection du matches le plus haut
    //      double coord_max[10] = {hT,hT,hT,hT,hT,hT,hT,hT,hT,hT};
    //for(int i  =0; i<ptref.size();i++ ){
    //    for (int j =0; j<10;j++){

    //    if (ptref[i].y<coord_max[j]){
    //        coord_max[j]=ptref[i].y;
    //        it_high_match[j]=i;
    //        break;
    //    }
    //}

    //}
    //std::cout<<"match le plus haut"<<std::endl;
    //for (int j =0; j<10;j++){

    //std::cout<<it_high_match[j]<<std::endl;
    //}
    int left,right,high,low;
    double maxleft,maxright,maxhigh,maxlow;
    maxleft=wT;
    maxright=0;
    maxhigh=hT;
    maxlow=0;
    for (int it =0; it <ptref.size();it++){
        if (ptref[it].y<maxhigh){
            maxhigh =ptref[it].y;
            high=it;
        }
        if ((ptref[it].y>maxlow)&&(ptref[it].x>wT/2)){
            maxlow =ptref[it].y;
            low=it;
        }
        if (ptref[it].x>maxright){
            maxright =ptref[it].x;
            right=it;
        }
        if ((ptref[it].x<maxleft)){
            maxleft =ptref[it].x;
            left=it;
        }
    }
    //       std::cout<<"high  : "<<high<< " "<< ptref[high]<<std::endl;
    //       std::cout<<"low   : "<<low<<" "<< ptref[low]<<std::endl;
    //       std::cout<<"right : "<<right<<" "<< ptref[right]<<std::endl;
    //       std::cout<<"left  : "<<left<<" "<< ptref[left]<<std::endl;
    cv::drawMarker(trainimg,ptref[high],cv::Scalar(0,1000,10));
    cv::drawMarker(trainimg,ptref[low],cv::Scalar(0,1000,10));
    cv::drawMarker(trainimg,ptref[left],cv::Scalar(0,1000,10));
    cv::drawMarker(trainimg,ptref[right],cv::Scalar(0,1000,10));
    cv::imshow("max",trainimg);
    //       cv::waitKey(0);


std::cout<<"max ok ok"<<std::endl;
std::cout<<"nb matches: "<< nb_ptmodel<<std::endl;

    //conversion pt ref et pt2 en mat
    //const int si= nb_ptmodel;
    Matrix<double,2,100> mat_ptref,mat_pt2;
    for (int i=0;i<ptref.size();i++){
        mat_ptref(0,i)=ptref[i].x;
        mat_ptref(1,i)=ptref[i].y;
        mat_pt2(0,i)=pt2[i].x;
        mat_pt2(1,i)=pt2[i].y;

    }
std::cout<<"conversion ok"<<std::endl;




    for (int ite=0;ite<nb_iteration;ite++){
        double param_guess[7]={Rx,Ry,Rz,tx,ty,tz,def};//{1.57,0,0,300,300,300,61.4};
        //   tempresult<<0,0,0,300,100,250,wT/2,0;
        std::cout<<"iteration no "<<ite<<std::endl;

        bool test = true;

        int points[nb_ptsref];
        while(test){
            for (int ke=0;ke<nb_ptsref;ke++){
                points[ke]= rand() % nb_ptmodel;
                //std::cout<<points[ke]<<std::endl;
            }
            //            int k = rand() %10;
            //            points[nb_ptsref-1]=high;
            //            points[nb_ptsref-2]=low;
            //            points[nb_ptsref-3]=left;
            //            points[nb_ptsref-4]=right;
            test=false;
            for( int a =0;a<nb_ptsref;a++){
                for (int b = 0;b<nb_ptsref;b++ ){
                    if(a!=b){
                        if (points[a]==points[b]){
                            test=true;
                        }
                    }

                }
                if(test){
                    break;
                }
            }
            // bool cond[nb_ptsref] =false;
            test=false;
            for( int a =0;a<nb_ptsref;a++){
                for (int b = 0;b<nb_ptsref;b++ ){
                    if(a!=b){
                        if (distance(ptref[points[a]],ptref[points[b]])<0){
                            test=true;
                        }
                    }
                }
                if(test){
                    break;
                }
            }
        }
        //assert
        for (int itee=0;itee<nb_ptsref;itee++){
            //   std::cout<<"point no: "<<points[itee]<<" : "<< ptref[points[itee]]<<std::endl;
        }
        int nb_inliers_temp=0;
        for (int itpt=0;itpt<nb_ptsref;itpt++){

//            x1(0, itpt) = ptref[points[itpt]].x;
//            x1(1, itpt) = ptref[points[itpt]].y;
//            x2(0, itpt) = pt2[points[itpt]].x;
//            x2(1, itpt) = pt2[points[itpt]].y;
                                        x1(0, itpt) = X1(0,points[itpt]);
                                        x1(1, itpt) = X1(1,points[itpt]);
                                        x2(0, itpt) = X2(0,points[itpt]);
                                        x2(1, itpt) = X2(1,points[itpt]);
        }

        testmodelfrom2DFromCorrespondences(X2,X1,wT,x2,x1,param_guess,options,&tempresult);//(mat_ptref,mat_pt2
        std::cout<<" solve ok"<<std::endl;

        inliers_temp.clear();
        for (int it=0;it<nb_ptmodel;it++){
            double threshold= 5;
           // Matrix<double,1,2> x1_(ptref[it].x,ptref[it].y);  //
            //Matrix<double,1,2> x2_(pt2[it].x,pt2[it].y); // //

                       Matrix<double,1,2> x1_(X1(0,it),X2(1,it));
                        Matrix<double,1,2> x2_(X2(0,it),X2(1,it));
            //verif inliers
            double dist = Distancetransfo2(L,tempresult,x1_,x2_);

            //std::cout<<it<<std::endl;
            if (dist<=threshold){
                bool cond = true;
                //                for( auto test : inliers_temp){
                //                    if(distance(pt2[test],pt2[it])<5){
                //                        cond =false;
                //                        break;
                //                    }
                //                }
                if(cond){
                    nb_inliers_temp++;
                    inliers_temp.push_back(it);
                    inliers_dist_temp.push_back(dist);
                }
                //std::cout<<"point no : "<<it<<std::endl;
                // std::cout<<"distance : "<<dist<<std::endl;

            }

        }
        std::cout<<"nombre inliers: "<<nb_inliers_temp<<std::endl;
        if(nb_inliers_temp>nb_inliers){
            inliers.clear();
            inliers_dist.clear();
            inliers_dist=inliers_dist_temp;
            nb_inliers=nb_inliers_temp;
            inliers=inliers_temp;
            result=tempresult;
            // std::cout<<"inliers "<<inliers.<<std::endl;
            affiche_avec_param(L,wT, hT,result);

        }
        //        std::cout<<"tempresult"<<std::endl;
        //        std::cout<<tempresult<<std::endl;
        //        std::cout<<"result"<<std::endl;
        //        std::cout<<result<<std::endl;
        //        std::cout<<" "<<std::endl;
        std::cout<<std::flush;
        std::cout<<"nombre inliers max: "<<nb_inliers<<std::endl;
        //        std::cout<< "parametres"<<std::endl;
        //        std::cout<< "rotations"<<std::endl;
        //        std::cout<< "Rx: "<<result[0]<< " Ry: "<<result[1]<< " Rz: "<<result[2]<<std::endl;
        //        std::cout<< "Translations"<<std::endl;
        //        std::cout<< "Tx: "<<result[3]<< " Ty: "<<result[4]<< " Tz: "<<result[5]<<std::endl;
        //        std::cout<< "Parametres déformation"<<std::endl;
        //        std::cout<< "Longueur: "<<L<<" deformation: "<<result[6] <<std::endl;
    }





    //    for(int i =0; i<8;i++){
    //        std::cout<<parametresguess[i]<<std::endl;
    //    }
    //testmodelfrom2DFromCorrespondences(X1,X2,parametresguess,options,&result);
    //cv::warpPerspective(Queryimg, imred, homo, trainimg.size());
    std::cout<<"fin fonction"<<std::endl;
    std::cout<<result<<std::endl;
    for (int i= 0;i<70;i++){
        cv::Point2d pt(X2(0,i),X2(1,i));
        cv::drawMarker(trainimg,pt,cv::Scalar(0,1000,1000));
    }
        cv::imshow("test_reconstruction",trainimg);
         cv::waitKey(0);
    //std::string outFilename("./finalmario2.jpg");
    //cv::imwrite(outFilename,Queryimg );
    //affiche_matches_avec_param(L,inliers,inliers_dist,pt2,ptref,wT,hT,result);
    affiche_match_et_set_avec_param(L,inliers,ptref,wT,hT,result);
    
}
