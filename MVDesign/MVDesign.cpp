#include "MVDesign.h"

int main(void)
{
	Mat tempImage[4], tarImage[4];

	//读图
	Mat temp = imread("../images/pattern.bmp");
	Mat tar = imread("../images/IMAGEB12.bmp");
	//转为灰度图
	cvtColor(temp, tempImage[0], COLOR_BGR2GRAY);
	cvtColor(tar, tarImage[0], COLOR_BGR2GRAY);

	//降采样
	pyrDown(tarImage[0], tarImage[1], Size(tarImage[0].cols / 2, tarImage[0].rows / 2));
	pyrDown(tarImage[1], tarImage[2], Size(tarImage[1].cols / 2, tarImage[1].rows / 2));
	pyrDown(tarImage[2], tarImage[3], Size(tarImage[2].cols / 2, tarImage[2].rows / 2));

	pyrDown(tempImage[0], tempImage[1], Size(tempImage[0].cols / 2, tempImage[0].rows / 2));
	pyrDown(tempImage[1], tempImage[2], Size(tempImage[1].cols / 2, tempImage[1].rows / 2));
	pyrDown(tempImage[2], tempImage[3], Size(tempImage[2].cols / 2, tempImage[2].rows / 2));

	int maxRow = 0, maxCol = 0;
	float maxNCC = 0, ncc = 0, maxAngle = 0;
	Mat tempRotate;		//旋转之后的图像

	double start = static_cast<double>(getTickCount());

	/* 先算金字塔顶端 */
	for (int i = 0; i < ANGLE_NUMBLE; i++) {
		//旋转模板
		rotate_arbitrarily_angle(tempImage[3], tempRotate, i * DANGLE);
		//求模板均值和方差
		Image_Mean_Variance(tempRotate);
		//求(r,c)NCC值
		for (int r = 0; r < tarImage[3].rows - tempRotate.rows; r++) {
			for (int c = 0; c < tarImage[3].cols - tempRotate.cols; c++) {
				ncc = NCC(r, c, tempRotate, tarImage[3], tempData[0], tempData[1]);
				if (ncc > maxNCC) {
					maxNCC = ncc;
					maxRow = r;
					maxCol = c;
					maxAngle = i * DANGLE;
				}
			}
		}
	}
	/* 在依次往上走 */
	//后面加temp都是暂存量
	float maxAngletemp, step = DANGLE;
	int maxRowtemp, maxColtemp;
	for (int t = 2; t >= 0; t--)
	{
		maxAngletemp = maxAngle;
		maxRow *= 2;
		maxCol *= 2;
		maxNCC = 0;
		maxRowtemp = maxRow;
		maxColtemp = maxCol;
		step /= 5;

		if (1)
		{
			for (int i = -4; i < 4; i++) {
				rotate_arbitrarily_angle(tempImage[t], tempRotate, maxAngle + i * step);
				//cout << i << "    ";
				Image_Mean_Variance(tempRotate);
				//cout << *tempData << "   " << *(tempData + 2) << endl;
				for (int r = -1; r <= 1; r++) {
					//cout << r << endl;
					for (int c = -1; c <= 1; c++) {
						//cout << maxRow * 2 << endl;
						ncc = NCC(maxRowtemp + r, maxColtemp + c, tempRotate, tarImage[t], tempData[0], tempData[1]);
						if (ncc > maxNCC) {
							maxNCC = ncc;
							maxRow = maxRowtemp + r;
							maxCol = maxColtemp + c;
							maxAngle = maxAngletemp + i * step;
						}
					}
				}
			}
		}

	}
	cout << "坐标(" << maxCol << "," << maxRow << ")   角度" << maxAngle << endl;
	double time = ((double)getTickCount() - start) / getTickFrequency();
	cout << "所用时间:" << time << "s" << endl;

	/* 画矩形 */
	int centerCol = maxCol + tempRotate.cols / 2, centerRow = maxRow + tempRotate.rows / 2;
	float len = sqrt(tempImage[0].cols * tempImage[0].cols + tempImage[0].rows * tempImage[0].rows) / 2;
	float theta0 = atan2f(tempImage[0].rows, tempImage[0].cols);
	float phi1 = theta0 + maxAngle * pi / 180, phi2 = theta0 - maxAngle * pi / 180;
	float deltaCol = len * cosf(phi1), deltaRow = len * sinf(phi2);

	line(tarImage[0], Point(centerCol - len * cosf(phi1), centerRow + len * sinf(phi1)), Point(centerCol + len * cosf(phi2), centerRow + len * sinf(phi2)), Scalar(255, 0, 0), 2);
	line(tarImage[0], Point(centerCol + len * cosf(phi2), centerRow + len * sinf(phi2)), Point(centerCol + len * cosf(phi1), centerRow - len * sinf(phi1)), Scalar(255, 0, 0), 2);
	line(tarImage[0], Point(centerCol + len * cosf(phi1), centerRow - len * sinf(phi1)), Point(centerCol - len * cosf(phi2), centerRow - len * sinf(phi2)), Scalar(255, 0, 0), 2);
	line(tarImage[0], Point(centerCol - len * cosf(phi2), centerRow - len * sinf(phi2)), Point(centerCol - len * cosf(phi1), centerRow + len * sinf(phi1)), Scalar(255, 0, 0), 2);

	circle(tarImage[0], Point(centerCol, centerRow), 3, Scalar(255, 0, 0), -1);	//矩形中心点

	imshow("result", tarImage[0]);
	imshow("template", tempImage[0]);

	waitKey(0);
	return 0;
}

//求图像的均值和方差
void Image_Mean_Variance(Mat src)
{
	int n = src.cols * src.rows;   //像素总数
	//均值
	float mt = 0;
	for (int r = 0; r < src.rows; r++) {
		for (int c = 0; c < src.cols; c++) {
			mt += src.at<uchar>(r, c);
		}
	}
	mt /= n;
	//方差
	float st2 = 0;
	for (int r = 0; r < src.rows; r++) {
		for (int c = 0; c < src.cols; c++) {
			st2 += (src.at<uchar>(r, c) - mt) * (src.at<uchar>(r, c) - mt);
		}
	}
	st2 /= n;

	tempData[0] = mt;
	tempData[1] = st2;
}

//(r,c)的NCC值
float NCC(int r, int c, Mat TemplateImage, Mat TargetImage, float TemplateMean, float TemplateVariance)
{
	float NCC = 0, mf = 0, sf2 = 0;
	int Num = TemplateImage.rows * TemplateImage.cols;

	if (r > TargetImage.rows - TemplateImage.rows || c > TargetImage.cols - TemplateImage.cols) {
		return 0;
	}
	//计算(r,c)的mf
	for (int u = 0; u < TemplateImage.rows; u++) {
		uchar* pTI = TargetImage.ptr<uchar>(r + u);
		for (int v = 0; v < TemplateImage.cols; v++) {
			mf += pTI[c + v];
		}
	}
	mf /= Num;

	//计算(r,c)的sf2
	for (int u = 0; u < TemplateImage.rows; u++) {
		uchar* pTI = TargetImage.ptr<uchar>(r + u);
		for (int v = 0; v < TemplateImage.cols; v++) {
			sf2 += (pTI[c + v] - mf) * (pTI[c + v] - mf);
		}
	}
	sf2 /= Num;
	//cout << "sf2    " << sf2 << endl;
	float sf = sqrt3(sf2);

	float st = sqrt3(TemplateVariance);

	//计算(r,c)的NCC
	for (int u = 0; u < TemplateImage.rows; u++) {
		uchar* pTempI = TemplateImage.ptr<uchar>(u);
		uchar* pTargetI = TargetImage.ptr<uchar>(r + u);
		for (int v = 0; v < TemplateImage.cols; v++) {
			NCC += ((pTempI[v] - TemplateMean) / st) * ((pTargetI[v + c] - mf) / sf);
		}
	}
	NCC /= Num;
	return NCC;
}

//旋转图像
void rotate_arbitrarily_angle(Mat& src, Mat& dst, float angle)
{
	float radian = (float)(angle / 180.0 * CV_PI);

	//填充图像
	int maxBorder = (int)(max(src.cols, src.rows) * 1.414); //即为sqrt(2)*max
	int dx = (maxBorder - src.cols) / 2;
	int dy = (maxBorder - src.rows) / 2;
	copyMakeBorder(src, dst, dy, dy, dx, dx, BORDER_CONSTANT);

	//旋转
	Point2f center((float)(dst.cols / 2), (float)(dst.rows / 2));
	Mat affine_matrix = getRotationMatrix2D(center, angle, 1.0);//求得旋转矩阵
	warpAffine(dst, dst, affine_matrix, dst.size());

	//计算图像旋转之后包含图像的最大的矩形
	float sinVal = abs(sin(radian));
	float cosVal = abs(cos(radian));
	Size targetSize((int)(src.cols * cosVal + src.rows * sinVal),
		(int)(src.cols * sinVal + src.rows * cosVal));

	//剪掉多余边框
	int x = (dst.cols - targetSize.width) / 2;
	int y = (dst.rows - targetSize.height) / 2;
	Rect rect(x, y, targetSize.width, targetSize.height);
	dst = Mat(dst, rect);
}

float Q_rsqrt(float number)
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y = number;
	i = *(long*)&y;  // evil floating point bit level hacking
	i = 0x5f3759df - (i >> 1); // what the fuck?
	y = *(float*)&i;
	y = y * (threehalfs - (x2 * y * y)); // 1st iteration
	// y  = y * ( threehalfs - ( x2 * y * y ) ); // 2nd iteration, this can be removed

#ifndef Q3_VM
#ifdef __linux__
	assert(!isnan(y));
#endif
#endif
	return y;
}

/* 开平方根快速算法 */
float sqrt3(float x)
{
	if (x == 0) return 0;
	float result = x;
	float xhalf = 0.5f * result;
	int i = *(int*)&result;

	i = 0x5f3759df - (i >> 1);
	result = *(float*)&i;
	result = result * (1.5f - xhalf * result * result);
	result = result * (1.5f - xhalf * result * result);
	return 1.0f / result;
}