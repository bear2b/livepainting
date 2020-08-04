
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
namespace {
// This structure contains options that controls how the homography
// estimation operates.
//
// Defaults should be suitable for a wide range of use cases, but
// better performance and accuracy might require tweaking.
struct EstimateHomographyOptions {
    // Default settings for homography estimation which should be suitable
    // for a wide range of use cases.
    EstimateHomographyOptions()
        :  max_num_iterations(50),
          expected_average_symmetric_distance(1e-16) {}
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
// Calculate symmetric geometric cost terms:
//
// forward_error = D(H * x1, x2)
// backward_error = D(H^-1 * x2, x1)
//
// Templated to be used with autodifferentiation.
template <typename T>
void SymmetricGeometricDistanceTerms(const Eigen::Matrix<T, 3, 3> &H,
                                     const Eigen::Matrix<T, 2, 1> &x1,
                                     const Eigen::Matrix<T, 2, 1> &x2,
                                     T forward_error[2],
                                     T backward_error[2]) {
    typedef Eigen::Matrix<T, 3, 1> Vec3;
    Vec3 x(x1(0), x1(1), T(1.0));
    Vec3 y(x2(0), x2(1), T(1.0));
    Vec3 H_x = H * x;
    Vec3 Hinv_y = H.inverse() * y;
    H_x /= H_x(2);
    Hinv_y /= Hinv_y(2);
    forward_error[0] = H_x(0) - y(0);
    forward_error[1] = H_x(1) - y(1);
    backward_error[0] = Hinv_y(0) - x(0);
    backward_error[1] = Hinv_y(1) - x(1);
}
// Calculate symmetric geometric cost:
//
//   D(H * x1, x2)^2 + D(H^-1 * x2, x1)^2
//
double SymmetricGeometricDistance(const Mat3 &H,
                                  const Vec2 &x1,
                                  const Vec2 &x2) {
    Vec2 forward_error, backward_error;
    SymmetricGeometricDistanceTerms<double>(H,
                                            x1,
                                            x2,
                                            forward_error.data(),
                                            backward_error.data());
    return forward_error.squaredNorm() +
            backward_error.squaredNorm();
}
// A parameterization of the 2D homography matrix that uses 8 parameters so
// that the matrix is normalized (H(2,2) == 1).
// The homography matrix H is built from a list of 8 parameters (a, b,...g, h)
// as follows
//
//         |a b c|
//     H = |d e f|
//         |g h 1|
//
template<typename T = double>
class Homography2DNormalizedParameterization {
public:
    typedef Eigen::Matrix<T, 8, 1> Parameters;     // a, b, ... g, h
    typedef Eigen::Matrix<T, 3, 3> Parameterized;  // H
    // Convert from the 8 parameters to a H matrix.
    static void To(const Parameters &p, Parameterized *h) {
        *h << p(0), p(1), p(2),
                p(3), p(4), p(5),
                p(6), p(7), 1.0;
    }
    // Convert from a H matrix to the 8 parameters.
    static void From(const Parameterized &h, Parameters *p) {
        *p << h(0, 0), h(0, 1), h(0, 2),
                h(1, 0), h(1, 1), h(1, 2),
                h(2, 0), h(2, 1);
    }
};
// 2D Homography transformation estimation in the case that points are in
// euclidean coordinates.
//
//   x = H y
//
// x and y vector must have the same direction, we could write
//
//   crossproduct(|x|, * H * |y| ) = |0|
//
//   | 0 -1  x2|   |a b c|   |y1|    |0|
//   | 1  0 -x1| * |d e f| * |y2| =  |0|
//   |-x2  x1 0|   |g h 1|   |1 |    |0|
//
// That gives:
//
//   (-d+x2*g)*y1    + (-e+x2*h)*y2 + -f+x2          |0|
//   (a-x1*g)*y1     + (b-x1*h)*y2  + c-x1         = |0|
//   (-x2*a+x1*d)*y1 + (-x2*b+x1*e)*y2 + -x2*c+x1*f  |0|
//
bool Homography2DFromCorrespondencesLinearEuc(
        const Mat &x1,
        const Mat &x2,
        Mat3 *H,
        double expected_precision) {
    assert(2 == x1.rows());
    assert(4 <= x1.cols());
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());
    int n = x1.cols();
    MatX8 L = Mat::Zero(n * 3, 8);
    Mat b = Mat::Zero(n * 3, 1);
    for (int i = 0; i < n; ++i) {
        int j = 3 * i;
        L(j, 0) =  x1(0, i);             // a
        L(j, 1) =  x1(1, i);             // b
        L(j, 2) =  1.0;                  // c
        L(j, 6) = -x2(0, i) * x1(0, i);  // g
        L(j, 7) = -x2(0, i) * x1(1, i);  // h
        b(j, 0) =  x2(0, i);             // i
        ++j;
        L(j, 3) =  x1(0, i);             // d
        L(j, 4) =  x1(1, i);             // e
        L(j, 5) =  1.0;                  // f
        L(j, 6) = -x2(1, i) * x1(0, i);  // g
        L(j, 7) = -x2(1, i) * x1(1, i);  // h
        b(j, 0) =  x2(1, i);             // i
        // This ensures better stability
        // TODO(julien) make a lite version without this 3rd set
        ++j;
        L(j, 0) =  x2(1, i) * x1(0, i);  // a
        L(j, 1) =  x2(1, i) * x1(1, i);  // b
        L(j, 2) =  x2(1, i);             // c
        L(j, 3) = -x2(0, i) * x1(0, i);  // d
        L(j, 4) = -x2(0, i) * x1(1, i);  // e
        L(j, 5) = -x2(0, i);             // f
//        std::cout<<"x1 : x: "<<x1(0,i)<<" y: "<<x1(1,i)<<std::endl;
//        std::cout<<"x2 : x: "<<x2(0,i)<<" y: "<<x2(1,i)<<std::endl;

    }
//    std::cout<<"matrice L "<<L<<std::endl;
//    std::cout<<"vecteur b "<<b<<std::endl;
    // Solve Lx=B
    const Vec h = L.fullPivLu().solve(b);
    //std::cout<<"vecteur h: "<<h<<std::endl;
    Homography2DNormalizedParameterization<double>::To(h, H);
    return (L * h).isApprox(b, expected_precision);
}
// Cost functor which computes symmetric geometric distance
// used for homography matrix refinement.
class HomographySymmetricGeometricCostFunctor {
public:
    HomographySymmetricGeometricCostFunctor(const Vec2 &x,
                                            const Vec2 &y)
        : x_(x), y_(y) { }
    template<typename T>
    bool operator()(const T* homography_parameters, T* residuals) const {
        typedef Eigen::Matrix<T, 3, 3> Mat3;
        typedef Eigen::Matrix<T, 2, 1> Vec2;
        Mat3 H(homography_parameters);
        Vec2 x(T(x_(0)), T(x_(1)));
        Vec2 y(T(y_(0)), T(y_(1)));
        SymmetricGeometricDistanceTerms<T>(H,
                                           x,
                                           y,
                                           &residuals[0],
                                           &residuals[2]);
        return true;
    }
    const Vec2 x_;
    const Vec2 y_;
};
// Termination checking callback. This is needed to finish the
// optimization when an absolute error threshold is met, as opposed
// to Ceres's function_tolerance, which provides for finishing when
// successful steps reduce the cost function by a fractional amount.
// In this case, the callback checks for the absolute average reprojection
// error and terminates when it's below a threshold (for example all
// points < 0.5px error).
class TerminationCheckingCallback : public ceres::IterationCallback {
public:
    TerminationCheckingCallback(const Mat &x1, const Mat &x2,
                                const EstimateHomographyOptions &options,
                                Mat3 *H)
        : options_(options), x1_(x1), x2_(x2), H_(H) {}
    virtual ceres::CallbackReturnType operator()(
            const ceres::IterationSummary& summary) {
        // If the step wasn't successful, there's nothing to do.
        if (!summary.step_is_successful) {
            return ceres::SOLVER_CONTINUE;
        }
        // Calculate average of symmetric geometric distance.
        double average_distance = 0.0;
        for (int i = 0; i < x1_.cols(); i++) {
            average_distance += SymmetricGeometricDistance(*H_,
                                                           x1_.col(i),
                                                           x2_.col(i));
        }
        average_distance /= x1_.cols();
        if (average_distance <= options_.expected_average_symmetric_distance) {
            return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
        }
        return ceres::SOLVER_CONTINUE;
    }
private:
    const EstimateHomographyOptions &options_;
    const Mat &x1_;
    const Mat &x2_;
    Mat3 *H_;
};
bool EstimateHomography2DFromCorrespondences(
        const Mat &x1,
        const Mat &x2,
        const EstimateHomographyOptions &options,
        Mat3 *H) {
    assert(2 == x1.rows());
    assert(4 <= x1.cols());
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());
    // Step 1: Algebraic homography estimation.
    // Assume algebraic estimation always succeeds.
    Homography2DFromCorrespondencesLinearEuc(x1,
                                             x2,
                                             H,
                                             EigenDouble::dummy_precision());
    std::cout<< "Estimated matrix after algebraic estimation:\n" << *H<<std::endl;
    // Step 2: Refine matrix using Ceres minimizer.
    ceres::Problem problem;
    for (int i = 0; i < x1.cols(); i++) {
        HomographySymmetricGeometricCostFunctor
                *homography_symmetric_geometric_cost_function =
                new HomographySymmetricGeometricCostFunctor(x1.col(i),
                                                            x2.col(i));
        problem.AddResidualBlock(
                    new ceres::AutoDiffCostFunction<
                    HomographySymmetricGeometricCostFunctor,
                    4,  // num_residuals
                    9>(homography_symmetric_geometric_cost_function),
                    NULL,
                    H->data());
    }
    // Configure the solve.
    ceres::Solver::Options solver_options;
    solver_options.linear_solver_type = ceres::DENSE_QR;
    solver_options.max_num_iterations = options.max_num_iterations;
    solver_options.update_state_every_iteration = true;
    // Terminate if the average symmetric distance is good enough.
    TerminationCheckingCallback callback(x1, x2, options, H);
    solver_options.callbacks.push_back(&callback);
    // Run the solve.
    ceres::Solver::Summary summary;
    ceres::Solve(solver_options, &problem, &summary);
    std::cout<<summary.FullReport()<<std::endl;
    return summary.IsSolutionUsable();

}
}  // namespace libmv

//données physique de la page
// transforme la position

//deformation de la feuille
//parametres:
//L: longueur de base de la feuille
//l: abscisse d'où la déformation commence
//def: deformation de la page (def represente le décalage en abscisse entre la deformée et la feuille de référence, qui cause la déformation)
Vector4f transformation2(cv::Point3f position,float l,float def){
    double x = position.x;
    double z= position.z;
    double newz;
    double newx;
    if(x<l){
        float h= 54/7*def/l;
        newz = h/(l*l)*x*x*x-2*h/l*x*x+h*x;
        newx = x*(l-def)/l;
    }
    else{
        newz=position.z;
        newx = x-def;
    }
    return Vector4f(newx,position.y,newz,1);
}
float distance(cv::Point2f refpts,cv::Point2f newpts){
    return pow((refpts.x-newpts.x),2)+pow((refpts.y-newpts.y),2);

}

//void transformation_inv(Matrix3d position){
//    return position/2;
//}
//calcul de la matrice extr F à partir des rotation Rx, Ry, Rz, et des translation tx,ty,tz
MatrixXf calcul_F(float Rx,float Ry, float Rz,float tx,float ty,float tz){
    //init F
    MatrixXf F(4,4);
    F.row(3).setZero();
    VectorXf translation(3);
    translation<<tx,ty,tz;
    float cosX = std::cos(Rx);
    float sinX = std::sin(Rx);
    float cosY= std::cos(Ry);
    float sinY = std::sin(Ry);
    float cosZ = std::cos(Rz);
    float sinZ = std::sin(Rz);
    MatrixXf rotx(3,3);
    rotx<<1,0,0,
            0,cosX,-sinX,
            0,sinX,cosX;
    MatrixXf roty(3,3);
    roty<<cosY,0,-sinY,
            0,1,0,
            sinY,0,cosY;
    MatrixXf rotz(3,3);
    rotz<<cosZ,-sinZ,0,
            sinZ,cosZ,0,
            0,0,1;
    MatrixXf r(3,3);
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
MatrixXf calcul_K(float fu,float fv,float cu,float cv){
    MatrixXf K(3,3);
    K<<fu,0,cu,
            0,fu,cv,
            0,0,1;
    return K;
}
MatrixXf from4x4to3x4(MatrixXf F){
    MatrixXf mat3x4(3,4);
    mat3x4.setZero();
    mat3x4(0,0)=1;
    mat3x4(1,1)=1;
    mat3x4(2,2)=1;
    MatrixXf F3x4(3,4);
    F3x4=mat3x4*F;
    return F3x4;
}



int main(int argc, char **argv)
{

    //initialisation variable
    std::string imRef("./data/test2.jpg");     //nom de l'image de reference
    std::string imquery("./data/nat.png");



    //definition des images
    cv::Mat trainimg = cv::imread(imRef,cv::COLOR_BGR2GRAY);
    cv::Mat Queryimg = cv::imread(imquery,cv::COLOR_BGR2GRAY);
    //definition width and height
    float wT=trainimg.size().width;
    float hT=trainimg.size().height;
    float wQ=Queryimg.size().width;
    float hQ=Queryimg.size().height;
    //definition des variables auxiliaires
    int nbit =0; //nombre iteration
    int maxNbInliers = 0;


    //copy calculate pos
    std::cout<<"initialisation"<<std::endl;
    //initialisation F
    float cu = 0;//imrefcolor.size().width/2;
    float cv = 0;//imrefcolor.size().height/2;
    float fu = 424;
    for(int k = 0; k<10 ;k++){
        float Rx= k/180*PI;
        float Ry= 0/180*PI;
        float Rz= 0/180*PI;
        float tx= 500;
        float ty = 500;
        float tz=1000;
        MatrixXf F = calcul_F(Rx,Ry,Rz,tx,ty,tz);
        MatrixXf K = calcul_K(fu,fu,cu,cv);
        //definition des points de départs
        std::vector<cv::Point3d> originpts;
        for(int i= 0; i<trainimg.size().width; i=i+100){
            for (int j= 0; j<trainimg.size().height;j=j+100){
                cv::Point3d pts(i,j,0);
                originpts.push_back(pts);
            }
        }

        //passage en 3x4de F
        MatrixXf F3x4= from4x4to3x4(F);
        std::vector<cv::Point2f> respts;
        cv::Mat img(trainimg.size().width,trainimg.size().height, CV_8UC3, cv::Scalar(1000,1000, 1000));
        int nbpts= originpts.size();
        Vector4f pos(4);
        Vector4f refpos(4);
        for(int i = 0; i<nbpts; i++){
            float L =trainimg.size().width;
            std::cout<<"pt de base : "<<originpts[i]<<std::endl;
            pos=transformation2(originpts[i],L/2,k*L/40);
            std::cout<<"pos : "<<pos<<std::endl;
            refpos(0)= originpts[i].x;
            refpos(1)=originpts[i].y;
            refpos(2)=originpts[i].z;
            refpos(3)=1;
            VectorXf ref(3);
            VectorXf res(3);
            ref = K*F3x4*refpos;
            res = K*F3x4*pos;
            cv::Point2f pt(res(0)/res(2),res(1)/res(2));
            cv::Point2f refpt(ref(0)/ref(2),ref(1)/ref(2));
            cv::drawMarker(img,pt,cv::Scalar(100,0,0));
        }
        //cv::imwrite("testres.jpg",img);
        //cv::imshow("test",trainimg);
        cv::imshow("test",img);
        //cv::waitKey(0);
    }


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
    std::cout<<"ok"<<std::endl;

    //suppression des plus faux
    matches.erase(matches.begin()+matches.size()*0.15f, matches.end());

    //affichage des matches
    cv::Mat imMatches;
    cv::drawMatches( Queryimg,keypoints2,trainimg,kpref , matches, imMatches);
    //cv::drawMatches(imcam, keypoints2, patref->data,patref->keypoints, matches, imMatches);
    cv::imwrite("./matches.jpg", imMatches);

    std::vector<cv::Point2f> pt2, ptref;
    for(int i = 0; i<matches.size();i++){
        pt2.push_back( keypoints2[ matches[i].queryIdx ].pt );
        ptref.push_back( kpref[ matches[i].trainIdx ].pt );
        }

        cv::Mat imred;
        cv::Mat homo;
        //calcul de l'homographie via openCV
        homo = cv::findHomography( pt2, ptref, cv::RANSAC );

        //calcul homographie via Solveur
        Mat x1(2, ptref.size());
        for (int i = 0; i <ptref.size(); ++i) {
            x1(0, i) = ptref[i].x;
            x1(1, i) = ptref[i].y;
        }
        Mat3 homography_matrix;

        // This matrix has been dumped from a Blender test file of plane tracking.
//        homo.at(0,0),homo.at(0,1),homo.at(0,2),
//                homo.at(1,0),homo.at(1,1),homo.at(1,2)
//                homo.at(2,0),homo.at(2,1),homo.at(2,2);
        homography_matrix <<20.3715, -0.461057, -111.964454,
                0.0,       0.617589, -192.379252,
                0.0,      -0.000983,    1.0;
        //cv::cv2eigen(homo,homography_matrix);
        Mat x2 = x1;
        for (int i = 0; i < x2.cols(); ++i) {

            //Vec3 homogenous_x1 = Vec3(x1(0, i), x1(1, i), 1.0);
            //Vec3 homogenous_x1 = Vec3(pt2[i].x, pt2[i].y, 1.0);
//            Vec3 homogenous_x2 = homography_matrix * homogenous_x1;
//            x2(0, i) = homogenous_x2[0];
//            x2(1, i) = homogenous_x2[1];
            x2(0, i) =pt2[i].x;
            x2(1, i) =  pt2[i].y;
        }
        Mat3 estimated_matrix;
        EstimateHomographyOptions options;
        options.expected_average_symmetric_distance = 0.02;
        EstimateHomography2DFromCorrespondences(x1, x2, options, &estimated_matrix);
        // Normalize the matrix for easier comparison.
        estimated_matrix /= estimated_matrix(2 ,2);
        std::cout << "Original matrix:\n" << homography_matrix << "\n";
        std::cout << "Estimated matrix:\n" << estimated_matrix << "\n";
//        cv::Mat homography(3,3,6);
//        for(int i=0;i<3;i++){
//           for (int j=0;j<3;j++){
//               homography.at<float>(i,j)=estimated_matrix(i,j);
//        }
//        }
//          std::cout<<"homo: "<<homo<<std::endl;
//          std::cout<<"homography: "<<std::endl;
//        std::cout<<homography<<std::endl;
          //cv::imshow("image nat",Queryimg);
          cv::warpPerspective(Queryimg, imred, homo, trainimg.size());

            //cv::waitKey(0);
             std::string outFilename("./final.jpg");
             cv::imwrite(outFilename, imred);

}


//NEW CODE

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
using namespace Eigen;
#include "glog/logging.h"
typedef Eigen::NumTraits<double> EigenDouble;
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Matrix<double, 3, 3> Mat3;
typedef Eigen::Matrix<double, 2, 1> Vec2;
typedef Eigen::Matrix<double, Eigen::Dynamic,  8> MatX8;
typedef Eigen::Vector3d Vec3;

class EstimateHomographyOptions {
    // Default settings for homography estimation which should be suitable
    // for a wide range of use cases.
public:
    EstimateHomographyOptions()
        :  max_num_iterations(50),
          expected_average_symmetric_distance(1e-16) {}
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
namespace {
//fonction pour le solveur
Eigen::Matrix<double, 2, 1> transform(Eigen::Matrix<double, 2, 1> trainpt,double rotation[3],double translation[3],double L,double def){

    Matrix<double,3,1> position(trainpt[0],trainpt[1],1);
    Matrix<double,4,1> pt= transformation2(position,L,def);
    Matrix<double,3,4> P = calcul_P(rotation[0],rotation[1],rotation[2],translation[0],translation[1],translation[2],444,444,0,0);
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
void testdistance(const Eigen::Matrix<double,1,8>   &param,
                  const Eigen::Matrix<double, 2, 1> &x1, //train point
                  const Eigen::Matrix<double, 2, 1> &x2, //query poi
                  double distance[2]){//,
    // T backward_error[2]) {
    typedef Eigen::Matrix<double, 3, 1> Vec3;
    double rotation[3];
    double translation[3];
    double L;
    double def;
    for (int i=0;i>3;i++){
        rotation[i]=param(i);
        translation[i]=param(i+3);
    }
    L = param[6];
    def =param[7];
    Eigen::Matrix<double,2,1> x1_trans =transform(x1,rotation,translation,L,def);
    distance[0] =x2(0)-x1_trans(0);
    distance[1]= x2(1)-x1_trans(1);
}
double Distancetransfo(const Eigen::Matrix<double, 1, 8> &param,const  Eigen::Matrix<double, 2, 1> &x1,
                       const  Eigen::Matrix<double, 2, 1> &x2) {
    Vec2 distance;
    testdistance(param,x1,x2,distance.data());
    return distance.squaredNorm();
    //                         x1,
    //                         x2,
    //                         forward_error.data());
    //    return forward_error.squaredNorm();
}

struct transformationFunctor {
    bool operator()(const double* parametres,const double* x, double * value) const {
        const Eigen::Matrix<double, 1,8> param;
        const Eigen::Matrix<double, 2, 1> x1;
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
        : x_(x), y_(y) {  compute_distortion.reset(new ceres::CostFunctionToFunctor<1, 1>(
                                                       new ceres::NumericDiffCostFunction<transformationFunctor,
                                                       ceres::CENTRAL,
                                                       ,
                                                       1>(
                                                           new transformationFunctor)));}
    template<typename T>
    bool operator()(const T* parametres, T* residuals) const {
        T * value;
        T * x1;
        x1[0]=T(x_[0]);
        x1[1]=T(x_[1]);

        (*compute_distortion)(parametres,x1,value);
        Eigen::Matrix<T, 2,1>xt(value[0],value[1]);
        Eigen::Matrix<T, 2,1> x2(y_(0), y_(1));
        residuals[0]= xt[0]-x2[0];
        residuals[1]= xt[1]-x2[1];
        return true;
    }
    const Vec2 x_;
    const Vec2 y_;
    std::unique_ptr<ceres::CostFunctionToFunctor<1, 1> > compute_distortion;
};

class testconditionfinCallback : public ceres::IterationCallback {
public:
    testconditionfinCallback(const Eigen::Matrix<double, 1, 8> &param,const Mat &x1, const Mat &x2,
                             const EstimateHomographyOptions &options): options_(options), x1_(x1), x2_(x2),param_(param) {}
    virtual ceres::CallbackReturnType operator()(
            const ceres::IterationSummary& summary) {
        // If the step wasn't successful, there's nothing to do.
        if (!summary.step_is_successful) {
            return ceres::SOLVER_CONTINUE;
        }
        // Calculate average of symmetric geometric distance.
        double average_distance = 0.0;
        for (int i = 0; i < x1_.cols(); i++) {
            average_distance += Distancetransfo(param_,
                                                x1_.col(i),
                                                x2_.col(i));
        }
        average_distance /= x1_.cols();
        if (average_distance <= options_.expected_average_symmetric_distance) {
            return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
        }
        return ceres::SOLVER_CONTINUE;
    }
private:
    const EstimateHomographyOptions &options_;
    const Mat &x1_;
    const Mat &x2_;
    Eigen::Matrix<double, 1, 8> param_;
};

bool testmodelfrom2DFromCorrespondences(const Mat &x1,const Mat &x2,const double parametresguess[8],const EstimateHomographyOptions &options,Matrix<double,1,8> *result) {
    assert(2 == x1.rows());
    assert(5 <= x1.cols());
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());
    // Step 1: Algebraic homography estimation.
    // Assume algebraic estimation always succeeds.
    //j'etablis le guess

    //definition des parametres de base
    std::cout<< "parametres initiaux"<<std::endl;
    std::cout<< "rotations"<<std::endl;
    std::cout<< "Rx: "<<parametresguess[0]<< " Ry: "<<parametresguess[1]<< " Rz: "<<parametresguess[2]<<std::endl;
    std::cout<< "Translations"<<std::endl;
    std::cout<< "Tx: "<<parametresguess[3]<< " Ty: "<<parametresguess[4]<< " Tz: "<<parametresguess[5]<<std::endl;
    std::cout<< "Parametres déformation"<<std::endl;

    std::cout<< "Longueur: "<<parametresguess[6]<<" deformation: "<<parametresguess[7] <<std::endl;
    // Step 2: Refine matrix using Ceres minimizer.
    ceres::Problem problem;
    for (int i = 0; i < x1.cols(); i++) {
        TestFunctor
                *parametres_cost_function =
                new TestFunctor(x1.col(i),x2.col(i));
        problem.AddResidualBlock(
                    new ceres::AutoDiffCostFunction<
                    TestFunctor,
                    2,  // num_residuals
                    8>(parametres_cost_function),
                    NULL,
                    result->data());
    }
    // Configure the solve.
    ceres::Solver::Options solver_options;
    solver_options.linear_solver_type = ceres::DENSE_QR;
    solver_options.max_num_iterations = options.max_num_iterations;
    solver_options.update_state_every_iteration = true;
    // Terminate if the average symmetric distance is good enough.
    testconditionfinCallback callback( *result,x1, x2, options);
    solver_options.callbacks.push_back(&callback);
    // Run the solve.
    ceres::Solver::Summary summary;
    ceres::Solve(solver_options, &problem, &summary);
    std::cout<<summary.FullReport()<<std::endl;
    return summary.IsSolutionUsable();

}


// Calculate symmetric geometric cost:
//
//   D(H * x1, x2)^2 + D(H^-1 * x2, x1)^2
//





// Termination checking callback. This is needed to finish the
// optimization when an absolute error threshold is met, as opposed
// to Ceres's function_tolerance, which provides for finishing when
// successful steps reduce the cost function by a fractional amount.
// In this case, the callback checks for the absolute average reprojection
// error and terminates when it's below a threshold (for example all
// points < 0.5px error).

//données physique de la page
// transforme la position

//deformation de la feuille
//parametres:
//L: longueur de base de la feuille
//l: abscisse d'où la déformation commence
//def: deformation de la page (def represente le décalage en abscisse entre la deformée et la feuille de référence, qui cause la déformation)

float distance(cv::Point2f refpts,cv::Point2f newpts){
    return pow((refpts.x-newpts.x),2)+pow((refpts.y-newpts.y),2);

}


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
    float cu = 0;//imrefcolor.size().width/2;
    float cv = 0;//imrefcolor.size().height/2;
    float fu = 424;
    float Rx= 90/180*PI;
    float Ry= 0/180*PI;
    float Rz= 0/180*PI;
    float tx= 500;
    float ty = 500;
    float tz=1000;
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
    std::vector<cv::Point2f> respts;
    cv::Mat img(trainimg.size().width,trainimg.size().height, CV_8UC3, cv::Scalar(1000,1000, 1000));
    int nbpts= originpts.size();
    Matrix<double,4,1> pos;
    Matrix<double,4,1> refpos;
    for(int i = 0; i<nbpts; i++){
        float L =trainimg.size().width;
        std::cout<<"pt de base : "<<originpts[i]<<std::endl;
        pos=transformation2(originpts[i],L/2,5*L/40);
        std::cout<<"pos : "<<pos<<std::endl;
        refpos(0)= originpts[i][0];
        refpos(1)=originpts[i][1];
        refpos(2)=originpts[i][2];
        refpos(3)=1;
        Matrix<double,4,1> ref;
       Matrix<double,3,1> res;
        ref = refpos;
        res = K*F3x4*pos;
        cv::Point2d pt(res(0)/res(2),res(1)/res(2));
        cv::Point2d refpt(ref(0)/ref(2),ref(1)/ref(2));
        cv::drawMarker(img,pt,cv::Scalar(100,0,0));
        cv::drawMarker(img,refpt,cv::Scalar(0,100,0));
    }
    //cv::imwrite("testres.jpg",img);
    //cv::imshow("test",trainimg);
    cv::imshow("test",img);
     cv::waitKey(0);



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
    std::cout<<"ok"<<std::endl;

    //suppression des plus faux
    matches.erase(matches.begin()+matches.size()*0.15f, matches.end());

    //affichage des matches
    cv::Mat imMatches;
    cv::drawMatches( Queryimg,keypoints2,trainimg,kpref , matches, imMatches);
    //cv::drawMatches(imcam, keypoints2, patref->data,patref->keypoints, matches, imMatches);
    cv::imwrite("./matches.jpg", imMatches);

    std::vector<cv::Point2f> pt2, ptref;
    for(int i = 0; i<matches.size();i++){
        pt2.push_back( keypoints2[ matches[i].queryIdx ].pt );
        ptref.push_back( kpref[ matches[i].trainIdx ].pt );
    }

    cv::Mat imred;
    cv::Mat homo;
    //calcul de l'homographie via openCV
    homo = cv::findHomography( pt2, ptref, cv::RANSAC );

    //calcul homographie via Solveur
    Mat x1(2, ptref.size());
    for (int i = 0; i <ptref.size(); ++i) {
        x1(0, i) = ptref[i].x;
        x1(1, i) = ptref[i].y;
    }
    Mat3 homography_matrix;
    Mat x2 = x1;
    for (int i = 0; i < x2.cols(); ++i) {
        x2(0, i) =pt2[i].x;
        x2(1, i) =  pt2[i].y;
    }
    Mat3 estimated_matrix;
    EstimateHomographyOptions options;
    options.expected_average_symmetric_distance = 0.02;


    //definition des points
    //definition de x1, cad les points de références
    std::vector<Vector3d> pts;
    std::vector<Vector3d> ptstrans;
    //Matrix<double,2,hT*wT/100> x1,x2;
    for(int i= 0; i<wT; i=i+wT/10){
        for (int j= 0; j<hT;j=j+hT/10){
            Vector3d pts(i,j,0);
           // pts.push_back(pts);
        }
    }






    //testmodelfrom2DFromCorrespondences(x1, x2, options, &estimated_matrix);
    // Normalize the matrix for easier comparison.
    estimated_matrix /= estimated_matrix(2 ,2);
    std::cout << "Original matrix:\n" << homography_matrix << "\n";
    std::cout << "Estimated matrix:\n" << estimated_matrix << "\n";
    //        cv::Mat homography(3,3,6);
    //        for(int i=0;i<3;i++){
    //           for (int j=0;j<3;j++){
    //               homography.at<float>(i,j)=estimated_matrix(i,j);
    //        }
    //        }
    //          std::cout<<"homo: "<<homo<<std::endl;
    //          std::cout<<"homography: "<<std::endl;
    //        std::cout<<homography<<std::endl;
    //cv::imshow("image nat",Queryimg);
    cv::warpPerspective(Queryimg, imred, homo, trainimg.size());

    //cv::waitKey(0);
    std::string outFilename("./final.jpg");
    cv::imwrite(outFilename, imred);

}
}



