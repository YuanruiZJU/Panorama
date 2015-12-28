#ifndef _SIFT_H__
#define _SIFT_H__

#include <opencv2/opencv.hpp>

using namespace cv;


#ifndef _KEY_POINT_CLASS
#define _KEY_POINT_CLASS

class Key_point
{
public:
	int octave;        //keypoint���ڵ���
	int layer;         //keypoint���ڵ����еĲ�
	float scale;       //keypoint�ĳ߶�
	float angle;       //keypoint������������ֱ��ͼ
	Point location;    //keypoint��ͼƬ�еľ���λ��
	float* hists;      //keypoint�ĻҶ�ֱ��ͼ
	float* magnitudes;  //��������keypoint������4*4*8��ģ
	Mat* img;           //keypoint����ͼƬ��ָ��
public:
	Key_point();
	~Key_point();
	bool Compare(Key_point& point_from_another_img, float& result);  // ʹ��ŷʽ����Ƚ�����֮��Ĳ���
};

#endif

#ifndef _DOG_PYRAMID_CLASS_
#define _DOG_PYRAMID_CLASS_

class DoGPyramid
{
// Attributes
public:
	Mat*		origin_image;         // ԭʼͼƬ
	int			n_octave;             // DOG������������
	double*		sigma_layer;          // ÿһ���sigma��Ϊһ������
	Vector<Mat> Gauss_layer;          // ��˹��������ÿһ����ÿһ��
	Vector<Mat> DoG_layer;            // DoG��������ÿһ����ÿһ��
	Vector<Key_point> keypoints;
// Functions
public:
	DoGPyramid();
	~DoGPyramid();
	void Initialize(Mat& dest);
	void BuildGaussPyramid();  // ��˹������
	void Build();              // DoG������
	int GetKeypoint();         // �õ��ؼ���
	bool AccurateLocalizeAndEmitCusrvature(int octave, int layer, int row, 
		                                    int col, Key_point& key_point);   // �ؼ��㾫ȷ��λ��������ԵЧӦ
	float CalcOrientationForPoint(Key_point& p);   // �����ݶ�ֱ��ͼ��ȷ���ݶȷ���
	void CalcKeyPointDescriptor();
};

#endif
#endif