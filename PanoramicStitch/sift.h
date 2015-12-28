#ifndef _SIFT_H__
#define _SIFT_H__

#include <opencv2/opencv.hpp>

using namespace cv;


#ifndef _KEY_POINT_CLASS
#define _KEY_POINT_CLASS

class Key_point
{
public:
	int octave;        //keypoint所在的组
	int layer;         //keypoint所在的组中的层
	float scale;       //keypoint的尺度
	float angle;       //keypoint的主方向所在直方图
	Point location;    //keypoint在图片中的具体位置
	float* hists;      //keypoint的灰度直方图
	float* magnitudes;  //用来描述keypoint特征的4*4*8的模
	Mat* img;           //keypoint所在图片的指针
public:
	Key_point();
	~Key_point();
	bool Compare(Key_point& point_from_another_img, float& result);  // 使用欧式距离比较两者之间的差异
};

#endif

#ifndef _DOG_PYRAMID_CLASS_
#define _DOG_PYRAMID_CLASS_

class DoGPyramid
{
// Attributes
public:
	Mat*		origin_image;         // 原始图片
	int			n_octave;             // DOG金字塔的组数
	double*		sigma_layer;          // 每一层的sigma作为一个数组
	Vector<Mat> Gauss_layer;          // 高斯金字塔的每一组中每一层
	Vector<Mat> DoG_layer;            // DoG金字塔的每一组中每一层
	Vector<Key_point> keypoints;
// Functions
public:
	DoGPyramid();
	~DoGPyramid();
	void Initialize(Mat& dest);
	void BuildGaussPyramid();  // 高斯金字塔
	void Build();              // DoG金字塔
	int GetKeypoint();         // 得到关键点
	bool AccurateLocalizeAndEmitCusrvature(int octave, int layer, int row, 
		                                    int col, Key_point& key_point);   // 关键点精确定位和消除边缘效应
	float CalcOrientationForPoint(Key_point& p);   // 计算梯度直方图，确定梯度方向
	void CalcKeyPointDescriptor();
};

#endif
#endif