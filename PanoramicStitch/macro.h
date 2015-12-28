#ifndef _MACRO_H__
#define _MACRO_H__

// The image's initial sigma, we assume the default value 0.5 
#define IMG_INITIAL_SIGMA      0.5f

// The initial sigma that we use to process the img
#define INITIAL_SIGMA_PROCESS  1.6f

// The kernel width of the guass kernel
#define GAUSS_KERNEL_WIDTH     6

// The width of subdomain of the keypoint descriptor
#define DESC_SUBDOMAIN_WIDTH   4

// The number of orientations for a seed point for a descriptor 
#define DESC_NUM_ORIENT        8

// If the distance of two keypoint descriptor, assume that keypoints are the same point
#define IMAGE_MATCH_EBSILON    1.0f

// The number of layers in a octave
#define NUMBER_LAYER_OCTAVE    3

// To initial the image, double size the image for -1 octave
#define DOUBLE_IMAGE_SIZE      true

// The Border reserved to conviniencely compute the descriptor
#define IMG_BORDER             5

// To accurately localize the keypoint, we use an iterative algorithm
#define ITER_STEPS            5

// 关键点精确定位所求偏移量的阀值
#define OFFSET_UPPER_LIMIT    0.5f

// The Contrast threshold 
#define CONTRAST_THRESHOLD     0.04f

// The curvature threshold    
#define CUSRVATURE_THRESHOLD  10.f

// 计算梯度直方图时的副方向与主方向比重的阀值
#define WEIGHT_LOWER_LIMIT  0.8f

// 对于关键点的梯度直方图的数量
#define KEYPOINT_ORIENT_HIST_BINS 36

// sigma的倍数
#define KEYPOINT_ORIENT_HIST_SIGMA 1.5f


#endif