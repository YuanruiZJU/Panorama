#include "sift.h"
#include <iostream>
using namespace cv;
using namespace std;

int main(int argc, char* argv[]) {
	Mat image = imread("D://img2.jpg"); 
	DoGPyramid test;
	test.origin_image = &image;
	//test.InitializeImg();
	Mat dest;
	test.Build();
	test.GetKeypoint();
	
	dest = test.DoG_layer[6];
	namedWindow("hello"); // ����һ������Ϊ "hello" �Ĵ���
	imshow("hello", image); // �ڴ��� "hello" ����ʾͼƬ
	waitKey(0); // �ȴ��û����¼���
	destroyWindow("hello"); // ���ٴ��� "hello"
	return 0;
}
