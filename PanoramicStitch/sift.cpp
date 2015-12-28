#include "macro.h"
#include "sift.h"

static bool isExtrema(const Mat& cur, const Mat& prev, const Mat& next, int row, int col, float threshold)
{
	float cur_point = cur.at<float>(row, col);
	if (abs(cur_point) < threshold)
		return false;
	bool More = false;
	bool Less = false;
	for (int i = row - 1; i < row + 2; i++)
	{
		for (int j = col - 1; j < col + 2; j++)
		{
			if (cur_point > cur.at<float>(i, j) || cur_point > prev.at<float>(i, j) || cur_point > next.at<float>(i, j))
				More = true;
			if (cur_point < cur.at<float>(i, j) || cur_point < prev.at<float>(i, j) || cur_point < next.at<float>(i, j))
				Less = true;
			if (More && Less)
				return false;
		}
	}
	if (cur_point >0 && More || cur_point <0 && Less)
	    return true;
	return false;
}


/*******************************************************************/
/*************         Key_point class       ***********************/
/*******************************************************************/
Key_point::Key_point()
{
	octave = 0;
	layer = 0;
	scale = 0;
	angle = 0;
	location = Point(0, 0);

	hists = NULL;
	magnitudes = NULL;
}

Key_point::~Key_point()
{
	if (hists != NULL)
		delete(hists);
	if (magnitudes != NULL)
		delete(magnitudes);
	hists = NULL;
	magnitudes = NULL;
}

bool Key_point::Compare(Key_point& point_from_another_img, float& result)
{
	Key_point& p = point_from_another_img;

	// if one of two points is not keypoint, returen false
	if (p.img == NULL || p.hists == NULL || p.magnitudes == NULL)
	{
		std::cout << "The point you compared is NOT keypoint!" << std::endl;
		return false;
	}

	if (img == NULL || hists == NULL || magnitudes == NULL)
	{
		std::cout << "The point is NOT keypoint!" << std::endl;
		return false;
	}

	// if p is not from another image, return false
	if (img == p.img)
	{
		std::cout << "The point you compared is from the same image." << std::endl;
		return false;
	}

	// compute the Euclidean distance
	int size_mag = 4 * DESC_SUBDOMAIN_WIDTH * DESC_NUM_ORIENT;
	float distance = 0;
	for (int i = 0; i < size_mag; i++)
		distance += (p.magnitudes[i] - magnitudes[i]) * (p.magnitudes[i] - magnitudes[i]);
	result = sqrtf(distance);
	return true;
}

/**********************************************************************/
/*********************  DoGPyramid class  *****************************/
/**********************************************************************/
DoGPyramid::DoGPyramid()
{
	origin_image = NULL;
	n_octave = 0;
	sigma_layer = NULL;
}

DoGPyramid::~DoGPyramid()
{
	origin_image = NULL;

	if (sigma_layer != NULL)
		delete(sigma_layer);
	sigma_layer = NULL;
}

void DoGPyramid::Initialize(Mat& dest)
{
	if (origin_image == NULL)
	{
		std::cout << "The image is NULL." << std::endl;
		return;
	}
	
	// compute the number of octave
	n_octave = cvRound(std::log((double)std::min(origin_image->cols, origin_image->rows)) /
		std::log(2.0) - 2) + 1;

	if (!DOUBLE_IMAGE_SIZE)
		n_octave = n_octave - 1;
	
	if (n_octave <= 0)
	{
		std::cout << "The image is not valid." << std::endl;
		return;
	}

	int gauss_layer_total = NUMBER_LAYER_OCTAVE + 3;
	
	// prepare the sigma for every layer in the gauss pyramid
	sigma_layer = new double[gauss_layer_total];
	sigma_layer[0] = INITIAL_SIGMA_PROCESS;
	float tempK = powf(2.0f, 1.0f / NUMBER_LAYER_OCTAVE);
	for (int i = 1; i < gauss_layer_total; i++)
	{
		float temp2 = sqrtf(powf(tempK, (float)i) - powf(tempK, (float)(i-1)));
		sigma_layer[i] = temp2 * sigma_layer[0];
	}
	
	// convert the image to gray image
	dest = Mat(origin_image->size(), CV_32FC1);
	
	if (origin_image->channels() == 3 || origin_image->channels() == 4)
	    cvtColor(*origin_image, dest, COLOR_BGR2GRAY);
	dest.convertTo(dest, DataType<float>::type, 1.0f, 0.0f);
//	normalize(dest, dest, 1.0, 0.0, CV_MINMAX);

	// if need double size image, double size it
	if (DOUBLE_IMAGE_SIZE)
	{
		resize(dest, dest, Size(origin_image->cols * 2, origin_image->rows * 2),
			0, 0, INTER_LINEAR);
	}
}

void DoGPyramid::BuildGaussPyramid()
{
	// Build the first layer of the first octave of the pyramid
	Mat layer0;
	Initialize(layer0);
    
	double first_sigma = DOUBLE_IMAGE_SIZE ? sigma_layer[0] * sigma_layer[0] - IMG_INITIAL_SIGMA * IMG_INITIAL_SIGMA * 4 :
		sigma_layer[0] * sigma_layer[0] - IMG_INITIAL_SIGMA * INITIAL_SIGMA_PROCESS;
	GaussianBlur(layer0, layer0, Size(), first_sigma, first_sigma);
	Gauss_layer.push_back(layer0);

	// Build another layer and another octave
	int gauss_layer_number = NUMBER_LAYER_OCTAVE + 3;

	for (int layer = 1; layer < gauss_layer_number; layer++)
	{
		Mat currentLayer;
		GaussianBlur(Gauss_layer[layer - 1], currentLayer, Size(), 
			sigma_layer[layer], sigma_layer[layer]);
		Gauss_layer.push_back(currentLayer);
	}

	for (int octave = 1; octave < n_octave; octave++)
	{
		for (int layer = 0; layer < gauss_layer_number; layer++)
		{
			Mat currentLayer;
			if (layer == 0)
			{
				Gauss_layer[(octave-1) * gauss_layer_number + NUMBER_LAYER_OCTAVE].copyTo(currentLayer);
				resize(currentLayer, currentLayer, Size(currentLayer.cols/2, currentLayer.rows/2), INTER_LINEAR);
				Gauss_layer.push_back(currentLayer);
			}
			else
			{
				GaussianBlur(Gauss_layer[octave*gauss_layer_number + layer - 1], currentLayer, Size(),
					sigma_layer[layer], sigma_layer[layer]);
				Gauss_layer.push_back(currentLayer);
			}
		}
	}
}

void DoGPyramid::Build()
{
	BuildGaussPyramid();

	int dog_layer_number = NUMBER_LAYER_OCTAVE + 2;
	int gauss_layer_number = NUMBER_LAYER_OCTAVE + 3;
	for (int octave = 0; octave < n_octave; octave++)
	{
		for (int layer = 0; layer < dog_layer_number; layer++)
		{
			const Mat& src1 = Gauss_layer[octave * gauss_layer_number + layer];
			const Mat& src2 = Gauss_layer[octave * gauss_layer_number + layer+ 1];
			Mat dest;
			subtract(src2, src1, dest, noArray(), DataType<float>::type);
			DoG_layer.push_back(dest);
		}
	}
}



int DoGPyramid::GetKeypoint()
{
	int number = 0;
	float threshold = cvFloor(0.5f * CONTRAST_THRESHOLD/NUMBER_LAYER_OCTAVE);
	int layer_number = NUMBER_LAYER_OCTAVE + 1;
	int dog_layer_number = NUMBER_LAYER_OCTAVE + 2;

	for (int octave = 0; octave < n_octave; octave ++)
	{
		for (int layer = 1; layer < layer_number; layer++)
		{
			Mat& current_img = DoG_layer[octave * dog_layer_number + layer];
			Mat& prev_img = DoG_layer[octave * dog_layer_number +layer - 1];
			Mat& next_img = DoG_layer[octave * dog_layer_number + layer + 1];

			for (int i = IMG_BORDER; i < current_img.rows - IMG_BORDER; i++)
			{
				for (int j = IMG_BORDER; j < current_img.cols - IMG_BORDER; j++)
				{
					Key_point p;
					if (isExtrema(current_img, prev_img, next_img, i, j, threshold) &&
						AccurateLocalizeAndEmitCusrvature(octave, layer, i,j, p))
					{
						float max_hist = CalcOrientationForPoint(p);

						float limit = WEIGHT_LOWER_LIMIT * max_hist;
						for (int t = 0; t < KEYPOINT_ORIENT_HIST_BINS; t++)
						{
							int left_t = t == 0 ? KEYPOINT_ORIENT_HIST_BINS - 1 : t - 1;
							int right_t = t == KEYPOINT_ORIENT_HIST_BINS - 1 ? 0 : t + 1;
							if (p.hists[t] > p.hists[left_t] && p.hists[t] > p.hists[right_t] && p.hists[t] >= limit)
							{
								Key_point temp_key;
								temp_key.octave = p.octave;
								float bin = t + 0.5 * (p.hists[left_t] - p.hists[right_t]) / 
									(p.hists[left_t] + p.hists[right_t] - 2 * p.hists[t]);
								if (bin < 0)
									bin += KEYPOINT_ORIENT_HIST_BINS;
								if (bin >= KEYPOINT_ORIENT_HIST_BINS)
									bin -= KEYPOINT_ORIENT_HIST_BINS;
								temp_key.angle = 360.0f - 360.0 / KEYPOINT_ORIENT_HIST_BINS *bin;
								if (std::abs(temp_key.angle - 360) < FLT_EPSILON)
									temp_key.angle = 0.f;
								temp_key.layer = p.layer;
								temp_key.img = p.img;
								temp_key.location = p.location;
								temp_key.scale = p.scale;
								keypoints.push_back(temp_key);
								std::cout << "octave:" << p.octave << " layer:" << p.layer << " scale:" << p.scale << " angle:" << temp_key.angle << std::endl;
								number++;
							}
						}
   					}
				}
			}
		}
	}
	std::cout << number << std::endl;
	return number;
}


bool DoGPyramid::AccurateLocalizeAndEmitCusrvature(int octave, int layer, int row, int col, Key_point& p)
{
	float drr, dcc, dss, drc, drs, dsc;
	float dr, dc, ds;
	int dog_layer_number = NUMBER_LAYER_OCTAVE + 2;
	int i;
	Vec3f offset;
	float limit_cusrvature = (CUSRVATURE_THRESHOLD + 1) * (CUSRVATURE_THRESHOLD + 1) / CUSRVATURE_THRESHOLD;
	for (i = 0; i < ITER_STEPS; i++)
	{
		const Mat& cur = DoG_layer[dog_layer_number * octave + layer];
		const Mat& prev = DoG_layer[dog_layer_number * octave + layer - 1];
		const Mat& next = DoG_layer[dog_layer_number * octave + layer + 1];
		float cur_val_double = cur.at<float>(row, col) * 2;
		dr = (cur.at<float>(row + 1, col) - cur.at<float>(row - 1, col)) * 0.5f;
		dc = (cur.at<float>(row, col + 1) - cur.at<float>(row, col - 1)) * 0.5f;
		ds = (next.at<float>(row, col) - prev.at<float>(row, col)) * 0.5f;

		Vec3f dD = { dr, dc, ds };

		drr = cur.at<float>(row + 1, col) + cur.at<float>(row - 1, col) - cur_val_double;
		dcc = cur.at<float>(row, col + 1) + cur.at<float>(row, col - 1) - cur_val_double;
		dss = next.at<float>(row, col) + prev.at<float>(row, col) - cur_val_double;
		drc = (cur.at<float>(row - 1, col - 1) + cur.at<float>(row + 1, col + 1) -
			cur.at<float>(row - 1, col + 1) - cur.at<float>(row + 1, col - 1))*0.25f;
		drs = (prev.at<float>(row - 1, col) + next.at<float>(row + 1, col) -
			prev.at<float>(row + 1, col) - next.at<float>(row - 1, col))*0.25f;
		dsc = (prev.at<float>(row, col - 1) + next.at<float>(row, col + 1) -
			prev.at<float>(row, col + 1) - next.at<float>(row, col - 1))*0.25f;

		Matx33f D_mat    = { drr, drc, drs,
			                 drc, dcc, dsc,
			                 dsc, drs, dss };
		offset = D_mat.solve(dD, DECOMP_LU);

		offset = -offset;
		if (offset[0] < 0.5f && offset[1] < 0.5f && offset[2] < 0.5f)
			break;

		if (std::abs(offset[0]) > (float)(INT_MAX / 3) ||
			std::abs(offset[1]) > (float)(INT_MAX / 3) ||
			std::abs(offset[2]) > (float)(INT_MAX / 3))
			return false;
		

		row += cvRound(offset[0]);
		col += cvRound(offset[1]);
		layer += cvRound(offset[2]);
		
		if (layer < 1 || layer > dog_layer_number - 2 ||
			row < IMG_BORDER || row > cur.rows - IMG_BORDER ||
			col < IMG_BORDER || col > cur.cols - IMG_BORDER)
			return false;
	}
	
	if (i >= ITER_STEPS)
		return false;

	Matx31f dD = { dr, dc, ds };
	Matx31f dX = { offset[0], offset[1], offset[2] };
	const Mat& img = DoG_layer[dog_layer_number * octave + layer];
	float Dx = img.at<float>(row, col) + dD.dot(dX) * 0.5f;

	if (std::abs(Dx) < CONTRAST_THRESHOLD / NUMBER_LAYER_OCTAVE * 255)
		return false;

	float tr = drr + dcc;
	float det = drr * dcc - drc * drc;
	if (det <= 0 || tr/det >= limit_cusrvature )
		return false;

	p.img = origin_image;
	p.layer = layer;
	p.location = Point(col, row);
	p.octave = octave;
	p.scale = INITIAL_SIGMA_PROCESS * powf(2.0f, (float)layer / NUMBER_LAYER_OCTAVE);
	return true;
}

float DoGPyramid::CalcOrientationForPoint(Key_point& p)
{
	p.hists = new float[KEYPOINT_ORIENT_HIST_BINS];
	int layer = p.layer;
	int octave = p.octave;
	int gauss_layer_number = NUMBER_LAYER_OCTAVE + 3;
	const Mat& img = Gauss_layer[octave * gauss_layer_number + layer];
	double scale = p.scale;
	double sigma = KEYPOINT_ORIENT_HIST_SIGMA * scale;
	int radius = cvRound(3.0f * sigma);
	int col = p.location.x;
	int row = p.location.y;
	double exp_scale = -1.0 / (2.0 * sigma*sigma);
	int number_bin = KEYPOINT_ORIENT_HIST_BINS;
		
	int len = (2 * radius + 1)*(2 * radius + 1);
	float* exps = new float[len];
	float* mag = new float[len];
	float* theta = new float[len];
	float* dr_array = new float[len];
	float* dc_array = new float[len];
	float* temp_hist = new float[number_bin + 4];
	
	int k = 0;
	for (int i = -radius; i <= radius; i++)
	{
		int tempRow = row + i;
		if (tempRow <= 0 || tempRow >= img.rows - 1)
			continue;
		for (int j = -radius; j <= radius; j++)
		{			
			int tempCol = col + j;
			if (tempCol <= 0 || tempCol >= img.cols - 1)
				continue;
			float dr = img.at<float>(tempRow + 1, tempCol) - img.at<float>(tempRow-1, tempCol);
			float dc = img.at<float>(tempRow, tempCol + 1) - img.at<float>(tempRow, tempCol - 1);
			exps[k] = ((double)(i*i + j*j)) * exp_scale;
			mag[k] = sqrtf(dr * dr + dc * dc);
			dr_array[k] = dr;
			dc_array[k] = dc;
			k++;
		}
	}

	len = k;
	exp(exps, exps, len);
	fastAtan2(dr_array, dc_array, theta, len, true);



	for (int i = 0; i < number_bin+4; i++)
		temp_hist[i] = 0.0f;
	for (int i = 0; i < len; i++)
	{
		int ori_hist_bin = cvRound(number_bin / 360.0f * theta[i]);
		if (ori_hist_bin < 0)
			ori_hist_bin += number_bin;
		if (ori_hist_bin >= number_bin)
			ori_hist_bin -= number_bin;
		temp_hist[ori_hist_bin + 2] += exps[i] * mag[i];
	}
	
	
	//Use linear interpolation to smooth the histogram
	temp_hist[number_bin + 2] = temp_hist[2];
	temp_hist[number_bin + 3] = temp_hist[3];
	temp_hist[0] = temp_hist[number_bin];
	temp_hist[1] = temp_hist[number_bin + 1];
	for (int i = 2; i < number_bin + 2; i++)
	{
		p.hists[i - 2] = (temp_hist[i - 2] + temp_hist[i + 2])*1.0f / 16.0f +
			(temp_hist[i - 1] + temp_hist[i + 1])*4.0f / 16.0f +
			temp_hist[i] * 6.0f / 16.0f;
	}
	
	float maxval = p.hists[0];
	for (int i = 1; i < number_bin; i++)
		if (p.hists[i] > maxval)
			maxval = p.hists[i];
	

	delete[] exps;
	delete[] mag;
	delete[] dr_array;
	delete[] dc_array;
	delete[] theta;
	delete[] temp_hist;

	return maxval;
}

void DoGPyramid::CalcKeyPointDescriptor()
{
	int keypsize = this->keypoints.size();
	int d = DESC_SUBDOMAIN_WIDTH;
	double temp = 1.5f * sqrtf(2.0) * (d + 1);
	for (int i = 0; i < keypsize; i++)
	{
		Key_point& p = this->keypoints[i];
		p.magnitudes = new float[128];
		for (int j = 0; j < 128; i++)
			p.magnitudes[j] = 0.0f;
		int octave = p.octave;
		int layer = p.layer;
		int row = p.location.y;
		int col = p.location.x;
		int gauss_layer_number = NUMBER_LAYER_OCTAVE + 3;
		Mat& img = Gauss_layer[octave * gauss_layer_number + layer];
		float angle = p.angle;
		float scale = p.scale;
		float cost = cosf(angle * (float)(CV_PI / 180.0f)) / 3.0f / scale;
		float sint = sinf(angle * (float)(CV_PI / 180.0f)) / 3.0f / scale;
		float bins_per_rad = DESC_NUM_ORIENT / 360.0f;
		float exp_scale = -1.f / (d * d * 0.5f);

		int radius = cvRound(temp * scale);
		radius = std::min(radius, (int)sqrt((double)img.cols*img.cols + img.rows * img.rows));

		int len = (2 * radius + 1) * (2 * radius + 1);
		float* dx_array = new float[len];
		float* dy_array = new float[len];
		float* x_bin_array = new float[len];
		float* y_bin_array = new float[len];
		float* exps = new float[len];
		float* theta = new float[len];
		float* mag = new float[len];
		int k = 0;
		for (int x = -radius; x <= radius; x++)
		{
			for (int y = -radius; y <= radius; y++)
			{
				// 旋转之后的坐标，并且放到d*d的区域中
				float x_rot = cost * (float)x - sint * (float)y;
				float y_rot = sint * (float)x + cost * (float)y;

				// 将上述坐标平移到第一象限
				float x_bin = x_rot + d / 2.0f - 0.5f;
				float y_bin = y_rot + d / 2.0f - 0.5f;
				int r = row + y;
				int c = col + x;

				// 
				if (x_bin > 0 && x_bin < d && y_bin > 0 && y_bin < d &&
					r > 0 && r < row - 1 && c > 0 && c < col - 1)
				{
					float dy = img.at<float>(r + 1, c) - img.at<float>(r - 1, c);
					float dx = img.at<float>(r, c + 1) - img.at<float>(r, c - 1);
					dx_array[k] = dx;
					dy_array[k] = dy;
					x_bin_array[k] = x_bin;
					y_bin_array[k] = y_bin;
					exps[k] = (x*x + y*y) * exp_scale;
					mag[k] = sqrtf(dx * dx + dy * dy);
					k++;
				}
			}
		}

		// Compute the angles
		len = k;
		fastAtan2(dy_array, dx_array, theta, len, true);
		exp(exps, exps, len);

		//Compute the hists
		int n = DESC_NUM_ORIENT;
		int hist_len = (d + 2)*(d + 2)*(n + 2);
		float* temp_hist = new float[hist_len];

		for (int j = 0; j < hist_len; j++)
		{
			temp_hist[j] = 0;
		}
		for (int j = 0; j < len; j++)
		{
			float rbin = y_bin_array[j];
			float cbin = x_bin_array[j];
			float ori_bin = (theta[j] - angle) * bins_per_rad;
			float magni = exps[k] * mag[k];

			int row0 = cvFloor(rbin);
			int col0 = cvFloor(cbin);
			int ori0 = cvFloor(ori_bin);

			rbin = rbin - row0;
			cbin = cbin - col0;
			ori_bin = ori_bin - ori0;

			if (ori0 >= n)
				ori0 -= n;
			if (ori0 < 0)
				ori0 += n;

			int index = ((row0 + 1) * (d + 2) + col0 + 1) * (n + 2) + ori0;

			float mag0 = magni * rbin;
			float mag1 = magni - mag0;
			float mag00 = mag0 * cbin;
			float mag01 = mag0 - mag00;
			float mag10 = mag1 * cbin;
			float mag11 = mag1 - mag10;
			float mag000 = mag0 * ori_bin;
			float mag001 = mag000 - mag001;
			float mag010 = mag01 * ori_bin;
			float mag011 = mag01 - mag010;
			float mag100 = mag10 * ori_bin;
			float mag101 = mag10 - mag100;
			float mag110 = mag11 * ori_bin;
			float mag111 = mag11 - mag110;

			temp_hist[index] = 
		}

		delete[] dx_array;
		delete[] dy_array;
		delete[] exps;
		delete[] theta;
		delete[] mag;
		delete[] temp_hist;
	}
}