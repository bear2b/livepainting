#ifndef PATTERN_H
#define PATTERN_H

#import "opencv2/opencv.hpp"
class Pattern
{
public:
    Pattern(cv::Mat image);
    cv::Size size;
    cv::Mat data;
    std::vector<cv::KeyPoint> keypoints;
    cv::Mat descriptors;
    std::vector<cv::Point2f>  points2d;
    std::vector<cv::Point3f>  points3d;
    void buildPatternFromImage(const cv::Mat& image, Pattern& pattern) const;

};

#endif // PATTERN_H
