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
	namedWindow("hello"); // 创建一个标题为 "hello" 的窗口
	imshow("hello", image); // 在窗口 "hello" 中显示图片
	waitKey(0); // 等待用户按下键盘
	destroyWindow("hello"); // 销毁窗口 "hello"
	return 0;
}
