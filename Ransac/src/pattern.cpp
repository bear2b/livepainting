#include "pattern.h"

Pattern::Pattern(cv::Mat image)
{
this->data=image.clone();
this->size=cv::Size(image.size().width,image.size().height);
this->points2d=
}
