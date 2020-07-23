#ifndef MVDESIGN_H
#define MVDESIGN_H

#include<iostream>
#include<opencv2/opencv.hpp>
#include<opencv2/core/core.hpp>

#define pi 3.141593
const double ANGLE_ARANGE = 60;
const double START_ANGLE = -30;
const double DANGLE = 5;
const double ANGLE_NUMBLE = ANGLE_ARANGE / DANGLE;

using namespace cv;
using namespace std;

void rotate_arbitrarily_angle(Mat& src, Mat& dst, float angle);
void Image_Mean_Variance(Mat src);
float NCC(int r, int c, Mat TemplateImage, Mat TargetImage, float TemplateMean, float TemplateVariance);
float tempData[2] = { 0 };
float sqrt3(float x);		//开平方根快速算法


#endif // !MVDESIGN_H

