#include "Application.h"
#include "qt_opengl_framework.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <assert.h>
using namespace std;

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene(void)
{

	ui_instance = Qt_Opengl_Framework::getInstance();

}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage(QString filePath)
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath)
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB(void)
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (!img_data)
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0; j < img_width; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j * 4), rgb + (out_offset + j * 3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0; i < 3; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();

	unsigned char Rdiscrete = 256 / 8;
	unsigned char Gdiscrete = 256 / 8;
	unsigned char Bdiscrete = 256 / 4;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			img_data[offset_rgba + rr] = img_data[offset_rgba + rr] / Rdiscrete * Rdiscrete;
			img_data[offset_rgba + gg] = img_data[offset_rgba + gg] / Gdiscrete * Gdiscrete;
			img_data[offset_rgba + bb] = img_data[offset_rgba + bb] / Bdiscrete * Bdiscrete;

			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();

	// 陣列順序依序為 r g b count
	std::vector<std::vector<unsigned short>> allColor(32768, std::vector<unsigned short>(4, 0));
	// 執行 uniform quantilization 每個顏色取5bit
	unsigned char Rdiscrete = 256 - 8;
	unsigned char Gdiscrete = 256 - 8;
	unsigned char Bdiscrete = 256 - 8;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			rgb[offset_rgb + rr] = rgb[offset_rgb + rr] & Rdiscrete;
			rgb[offset_rgb + gg] = rgb[offset_rgb + gg] & Gdiscrete;
			rgb[offset_rgb + bb] = rgb[offset_rgb + bb] & Bdiscrete;

			std::vector<unsigned short> & ref = allColor[rgb[offset_rgb + rr] * 128 + rgb[offset_rgb + gg] * 4 + rgb[offset_rgb + bb] / 8];

			ref[rr] = rgb[offset_rgb + rr];
			ref[gg] = rgb[offset_rgb + gg];
			ref[bb] = rgb[offset_rgb + bb];
			++ref[3];
		}
	}
	// 依顏色使用順序進行排序，並砍掉多出256以外的顏色
	sort(allColor.begin(), allColor.end(), [](std::vector<unsigned short> &item1, std::vector<unsigned short> &item2) {
		return item1[3] > item2[3];
	});
	if (allColor.size() > 256u) {
		allColor.erase(allColor.begin() + 256, allColor.end());
	}

	// 依每個pixel 來算最接近的距離並且替換成最接近的顏色
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int MAX_DISTANCE = INT_MAX, index = -1;
			for (unsigned int k = 0; k < allColor.size(); ++k )
			{
				auto & color = allColor[k];
				int distance = int(pow(color[rr] - img_data[offset_rgba + rr], 2.f) + pow(color[gg] - img_data[offset_rgba + gg], 2.f) + pow(color[bb] - img_data[offset_rgba + bb], 2.f));
				if (MAX_DISTANCE > distance) {
					MAX_DISTANCE = distance;
					index = k;
				}
			}
			
			img_data[offset_rgba + rr] = allColor[index][rr];
			img_data[offset_rgba + gg] = allColor[index][gg];
			img_data[offset_rgba + bb] = allColor[index][bb];
		}
	}


	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			if (gray / 256 > 0.5)
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = 255;
			}
			else
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = 0;
			}

			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{

	unsigned char *rgb = this->To_RGB();
	srand(time(NULL));
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			float intesity = (gray / 256) + (((float)rand() / (float)RAND_MAX*0.4) - 0.2);

			if (intesity > 0.5)
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = 255;
			}
			else
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = 0;
			}
		}
	}


	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();


}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	Gray();
	Dither_Color();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{


	unsigned char *rgb = this->To_RGB();

	float sum = 0;
	float ct = 0;

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			sum += gray / 256;
			ct++;
		}
	}

	float avg = sum / ct;

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			if ((gray / 256) > avg)
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = 255;
			}
			else
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = 0;
			}

			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float mask[4][4] = { 0.7059, 0.3529, 0.5882, 0.2353,
				0.0588, 0.9412, 0.8235, 0.4118,
				0.4706, 0.7647,0.8824, 0.1176,
				0.1765, 0.5294, 0.2941, 0.6471 };

			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			if (gray / 256 >= mask[j % 4][i % 4]) {
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = 255;
			}
			else {
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = gray;
			}

			img_data[offset_rgba + aa] = WHITE;
		}
	}


	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int offset;
			float val;
			float r, g, b;

			float oldR = img_data[offset_rgba + rr];
			float oldG = img_data[offset_rgba + gg];
			float oldB = img_data[offset_rgba + bb];

			float Rdiscrete = 255.0 / 7;
			float Gdiscrete = 255.0 / 7;
			float Bdiscrete = 255.0 / 7;

			float newR = round(oldR / Rdiscrete) * Rdiscrete;
			float newG = round(oldG / Gdiscrete) * Gdiscrete;
			float newB = round(oldB / Bdiscrete) * Bdiscrete;

			img_data[offset_rgba + rr] = newR;
			img_data[offset_rgba + gg] = newG;
			img_data[offset_rgba + bb] = newB;

			float qtErrR = oldR - newR;
			float qtErrG = oldG - newG;
			float qtErrB = oldB - newB;

			if ((j + 1) < img_width)
			{
				offset = i * img_width * 4 + (j + 1) * 4;
				val = 7.0 / 16.0;
				r = img_data[offset + rr] + qtErrR * val;
				g = img_data[offset + gg] + qtErrG * val;
				b = img_data[offset + bb] + qtErrB * val;
				if (r > 255)r = 255; if (r < 0)r = 0;
				if (g > 255)g = 255; if (g < 0)g = 0;
				if (b > 255)b = 255; if (b < 0)b = 0;
				img_data[offset + rr] = r;
				img_data[offset + gg] = g;
				img_data[offset + bb] = b;
			}

			if ((i + 1) < img_height && (j - 1) >= 0)
			{
				offset = (i + 1) * img_width * 4 + (j - 1) * 4;
				val = 3.0 / 16.0;
				r = img_data[offset + rr] + qtErrR * val;
				g = img_data[offset + gg] + qtErrG * val;
				b = img_data[offset + bb] + qtErrB * val;
				if (r > 255)r = 255; if (r < 0)r = 0;
				if (g > 255)g = 255; if (g < 0)g = 0;
				if (b > 255)b = 255; if (b < 0)b = 0;
				img_data[offset + rr] = r;
				img_data[offset + gg] = g;
				img_data[offset + bb] = b;
			}

			if ((i + 1) < img_height)
			{
				offset = (i + 1) * img_width * 4 + j * 4;
				val = 5.0 / 16.0;
				r = img_data[offset + rr] + qtErrR * val;
				g = img_data[offset + gg] + qtErrG * val;
				b = img_data[offset + bb] + qtErrB * val;
				if (r > 255)r = 255; if (r < 0)r = 0;
				if (g > 255)g = 255; if (g < 0)g = 0;
				if (b > 255)b = 255; if (b < 0)b = 0;
				img_data[offset + rr] = r;
				img_data[offset + gg] = g;
				img_data[offset + bb] = b;
			}

			if ((i + 1) < img_height && (j + 1) < img_width)
			{
				offset = (i + 1) * img_width * 4 + (j + 1) * 4;
				val = 1.0 / 16.0;
				r = img_data[offset + rr] + qtErrR * val;
				g = img_data[offset + gg] + qtErrG * val;
				b = img_data[offset + bb] + qtErrB * val;
				if (r > 255)r = 255; if (r < 0)r = 0;
				if (g > 255)g = 255; if (g < 0)g = 0;
				if (b > 255)b = 255; if (b < 0)b = 0;
				img_data[offset + rr] = r;
				img_data[offset + gg] = g;
				img_data[offset + bb] = b;
			}

		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering(double filter[][5])
{
	unsigned char *rgb = this->To_RGB();



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::filtering(double **filter, int n)
{
	bool edgeFlag = false;
	if (n == -5)
	{
		edgeFlag = true;
		n = 5;
	}
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int startPixel = offset_rgb - 2 * img_width * 3 - 2 * 3;

			int sum = 0;

			if (!edgeFlag)
				for (int x = 0; x < n; x++)
				{
					for (int y = 0; y < n; y++)
					{
						sum += filter[x][y];
					}
				}
			else
				sum = 256;


			for (int k = 0; k < 3; k++)
			{
				int sttemp = startPixel + k;
				double colorSum = 0;
				for (int x = 0; x < n; x++)
				{
					for (int y = 0; y < n; y++)
					{
						int pixelAt = sttemp + x * img_width * 3 + y * 3;
						if (pixelAt < 0 || pixelAt>img_height*img_width * 3)
							colorSum += 0;
						else
							colorSum += rgb[pixelAt] * filter[x][y];
					}
				}
				if (colorSum < 0) {
					colorSum = 0;
				}
				img_data[offset_rgba + k] = colorSum / sum;
			}

			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	double ** filter;
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			filter[i][j] = 1;

	filtering(filter, 5);

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	double ** filter;
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));
	for (int i = 0; i < 5; i++)
	{
		int n;
		if (i == 1 || i == 3)
			n = 2;
		else if (i == 2)
			n = 3;
		else
			n = 1;
		for (int j = 0; j < 5; j++)
		{
			int m;
			if (j == 1 || j == 3)
				m = 2;
			else if (j == 2)
				m = 3;
			else
				m = 1;

			filter[i][j] = n * m;
		}
	}

	filtering(filter, 5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	double ** filter;
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));

	for (int i = 0; i < 5; i++) {

		int tmpS = 1, tmpE = 1;
		for (int j = 0; j < i; j++)
			tmpS *= (5 - 1 - j);
		for (int j = 1; j <= i; j++)
			tmpE *= j;

		filter[0][i] = tmpS / tmpE;
		filter[i][0] = tmpS / tmpE;
	}
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			if (j != 0 && i != 0)
				filter[i][j] = filter[0][j] * filter[i][0];
		}
	}

	filtering(filter, 5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////



void Application::Filter_Gaussian_N(unsigned int N)
{
	// testing for Error
	unsigned char * currentRGB = this->To_RGB();
	ofstream out;
	out.open("Gaussian_error.csv");


	double ** filter;
	filter = (double **)malloc(N * sizeof(double *));
	for (int i = 0; i < N; i++)
		filter[i] = (double *)malloc(N * sizeof(double));

	for (int i = 0; i < N; i++) {

		int tmpS = 1, tmpE = 1;
		for (int j = 0; j < i; j++)
			tmpS *= (N - 1 - j);
		for (int j = 1; j <= i; j++)
			tmpE *= j;

		filter[0][i] = tmpS / tmpE;
		filter[i][0] = tmpS / tmpE;
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j != 0 && i != 0)
				filter[i][j] = filter[0][j] * filter[i][0];
		}
	}

	filtering(filter, N);

	// testing for Error
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			float error = sqrt(
				pow(img_data[offset_rgba + rr] - currentRGB[offset_rgb + rr], 2.f) +
				pow(img_data[offset_rgba + gg] - currentRGB[offset_rgb + gg], 2.f) +
				pow(img_data[offset_rgba + bb] - currentRGB[offset_rgb + bb], 2.f)
			);
			out << error << ',';
		}
		out << endl;
	}
	out.close();
}




///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	double ** filter;
	filter = (double **)malloc(5 * sizeof(double *));
	for (int i = 0; i < 5; i++)
		filter[i] = (double *)malloc(5 * sizeof(double));


	for (int i = 0; i < 5; i++) {

		int tmpS = 1, tmpE = 1;
		for (int j = 0; j < i; j++)
			tmpS *= (5 - 1 - j);
		for (int j = 1; j <= i; j++)
			tmpE *= j;

		filter[0][i] = tmpS / tmpE;
		filter[i][0] = tmpS / tmpE;
	}

	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			if (j != 0 && i != 0)
				filter[i][j] = filter[0][j] * filter[i][0];
		}
	}

	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			filter[i][j] = -filter[i][j];
		}
	}


	filter[2][2] = 220;

	filtering(filter, -5);

	/*	unsigned char *rgb = this->To_RGB();

		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				int offset_rgb = i * img_width * 3 + j * 3;
				int offset_rgba = i * img_width * 4 + j * 4;

				for (int k = 0; k < 3; k++) {
					if(img_data[offset_rgba + k] >10)
					img_data[offset_rgba + k] +=10;
				}
			}
		}

		delete[] rgb;
		mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
		renew();*/
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	unsigned char *rgb = this->To_RGB();

	this->Filter_Edge();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			for (int k = 0; k < 3; k++) {
				double rgbTemp = rgb[offset_rgb + k];
				double imgDataTemp = img_data[offset_rgba + k];

				if (imgDataTemp + rgbTemp > 255) {
					img_data[offset_rgba + k] = 255;
				}
				else {
					img_data[offset_rgba + k] = imgDataTemp + rgbTemp;
				}

			}
		}
	}


	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////

Stroke applyFilter(unsigned char *img_data_rgba, int img_width, int img_height, int x, int y, std::vector<std::vector<double>> &filter)
{
	int filterSizeY = (int)filter.size();
	int filterSizeX = (int)filter.front().size();
	Stroke result(0,0,0,0,0,0,0);
	int sum = 0;
	for (std::vector<double> &a : filter)
	{
		for (double &b : a)
		{
			sum += b;
		}
	}
	if (x >= 300 && y >= 200) {
		int a = 000;
	}

	const int rr = 2, gg = 1, bb = 0, aa = 3;
	int r = 0, g = 0, b = 0, a = 0;
	for (int i = -filterSizeY / 2; i < filterSizeY - filterSizeY / 2; i++)
	{
		for (int j = -filterSizeX / 2; j < filterSizeX - filterSizeX / 2; j++)
		{
			int rgb_offset = (i + y) * img_width * 4 + (j + x) * 4;
			if (i + y >= img_height || i + y < 0 || j + x >= img_width || j + x < 0)
			{
				continue;
			} 
			else
			{
				r += img_data_rgba[rgb_offset + rr] * filter[i + filterSizeY / 2][j + filterSizeX / 2];
				g += img_data_rgba[rgb_offset + gg] * filter[i + filterSizeY / 2][j + filterSizeX / 2];
				b += img_data_rgba[rgb_offset + bb] * filter[i + filterSizeY / 2][j + filterSizeX / 2];
				a += img_data_rgba[rgb_offset + aa] * filter[i + filterSizeY / 2][j + filterSizeX / 2];
			}
		}
	}

	result.r = r/sum;
	result.g = g/sum;
	result.b = b/sum;
	result.a = a/sum;
	return result;
}
Stroke applyFilterRGB(unsigned char *img_data_rgb, int img_width, int img_height, int x, int y, std::vector<std::vector<int>> &filter)
{
	int filterSizeY = (int)filter.size();
	int filterSizeX = (int)filter.front().size();
	Stroke result(0, 0, 0, 0, 0, 0, 0);
	int sum = 0;
	for (std::vector<int> &a : filter)
	{
		for (int &b : a)
		{
			sum += b;
		}
	}

	const int rr = 2, gg = 1, bb = 0;
	int r = 0, g = 0, b = 0;
	for (int i = -filterSizeY / 2; i < filterSizeY - filterSizeY / 2; i++)
	{
		for (int j = -filterSizeX / 2; j < filterSizeX - filterSizeX / 2; j++)
		{
			int rgb_offset = (i + y) * img_width * 3 + (j + x) * 3;
			if (i + y >= img_height || i + y < 0 || j + x >= img_width || j + x < 0)
			{
				continue;
			}
			else
			{
				r += img_data_rgb[rgb_offset + rr] * filter[i + filterSizeY / 2][j + filterSizeX / 2];
				g += img_data_rgb[rgb_offset + gg] * filter[i + filterSizeY / 2][j + filterSizeX / 2];
				b += img_data_rgb[rgb_offset + bb] * filter[i + filterSizeY / 2][j + filterSizeX / 2];
			}
		}
	}

	result.r = r / sum;
	result.g = g / sum;
	result.b = b / sum;
	return result;
}
Stroke applyFilterRGB(unsigned char *img_data_rgb, int img_width, int img_height, int x, int y, std::vector<std::vector<double>> &filter)
{
	int filterSizeY = (int)filter.size();
	int filterSizeX = (int)filter.front().size();
	Stroke result(0, 0, 0, 0, 0, 0, 0);
	double sum = 0;
	for (std::vector<double> &a : filter)
	{
		for (double &b : a)
		{
			sum += b;
		}
	}

	const int rr = 2, gg = 1, bb = 0;
	double r = 0, g = 0, b = 0;
	for (int i = -filterSizeY / 2; i < filterSizeY - filterSizeY / 2; i++)
	{
		for (int j = -filterSizeX / 2; j < filterSizeX - filterSizeX / 2; j++)
		{
			int rgb_offset = (i + y) * img_width * 3 + (j + x) * 3;
			if (i + y >= img_height || i + y < 0 || j + x >= img_width || j + x < 0)
			{
				continue;
			}
			else
			{
				r += (double)img_data_rgb[rgb_offset + rr] * (double)filter[i + filterSizeY / 2][j + filterSizeX / 2];
				g += (double)img_data_rgb[rgb_offset + gg] * (double)filter[i + filterSizeY / 2][j + filterSizeX / 2];
				b += (double)img_data_rgb[rgb_offset + bb] * (double)filter[i + filterSizeY / 2][j + filterSizeX / 2];
			}
		}
	}
	assert(255.f >= r / sum);
	assert(255.f >= g / sum);
	assert(255.f >= b / sum);
	result.r = unsigned char (r / sum);
	result.g = unsigned char (g / sum);
	result.b = unsigned char (b / sum);
	return result;
}

void Application::Half_Size()
{
	unsigned char *rgb = this->To_RGB();

	int new_img_width = img_width  / 2;
	int new_img_height = img_height / 2;
	unsigned char *new_img_data = new unsigned char[new_img_width * new_img_height * 4];

	std::vector<std::vector<double>> filter = { {1,2,1},{2,4,2},{1,2,1} };
	for (int i = 0; i < new_img_height; i++)
	{
		for (int j = 0; j < new_img_width; j++)
		{
			int new_offset_rgba = i * new_img_width * 4 + j * 4;
			Stroke myStroke = applyFilter(img_data,img_width,img_height, 2 *j, 2 * i, filter);
			new_img_data[new_offset_rgba + rr] = myStroke.r;
			new_img_data[new_offset_rgba + gg] = myStroke.g;
			new_img_data[new_offset_rgba + bb] = myStroke.b;
			new_img_data[new_offset_rgba + aa] = myStroke.a;
		}
	}

	delete[] rgb;
	// delete[] img_data;
	img_data = new_img_data;

	mImageDst = QImage(img_data, new_img_width, new_img_height, QImage::Format_ARGB32);
	img_width = new_img_width;
	img_height = new_img_height;
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	unsigned char *rgb = this->To_RGB();

	int new_img_width = img_width * 2;
	int new_img_height = img_height * 2;
	unsigned char *new_img_data = new unsigned char[new_img_width * new_img_height * 4];

	std::vector<std::vector<double>> filter1 = { {1,2,1},{2,4,2}, {1,2,1} };
	std::vector<std::vector<double>> filter2 = { {1,3,3,1},{3,9,9,3},{3,9,9,3},{1,3,3,1} };
	std::vector<std::vector<double>> filter3 = { {1,2,1},{3,6,3},{3,6,3},{1,2,1} };
	for (int i = 0; i < new_img_height; i++)
	{
		for (int j = 0; j < new_img_width; j++)
		{
			int new_offset_rgba = i * new_img_width * 4 + j * 4;
			Stroke myStroke;
			if (i % 2 == 0 && j % 2 == 0)
			{
				myStroke = applyFilter(img_data, img_width, img_height, j / 2, i / 2, filter1);
			}
			else if (i % 2 == 1 && j % 2 == 1)
			{
				Stroke tmp1 = applyFilter(img_data, img_width, img_height, j / 2 - 1, i / 2 - 1, filter2);
				Stroke tmp2 = applyFilter(img_data, img_width, img_height, j / 2, i / 2, filter2);
				Stroke tmp3 = applyFilter(img_data, img_width, img_height, j / 2 + 1, i / 2 + 1, filter2);
				Stroke tmp4 = applyFilter(img_data, img_width, img_height, j / 2 + 2, i / 2 + 2, filter2);
				myStroke.r = tmp1.r / 4 + tmp2.r / 4 + tmp3.r / 4 + tmp4.r / 4;
				myStroke.g = tmp1.g / 4 + tmp2.g / 4 + tmp3.g / 4 + tmp4.g / 4;
				myStroke.b = tmp1.b / 4 + tmp2.b / 4 + tmp3.b / 4 + tmp4.b / 4;
				myStroke.a = tmp1.a / 4 + tmp2.a / 4 + tmp3.a / 4 + tmp4.a / 4;

			}
			else
			{
				Stroke tmp1 = applyFilter(img_data, img_width, img_height, j / 2 - 1, i / 2 - 1, filter3);
				Stroke tmp2 = applyFilter(img_data, img_width, img_height, j / 2, i / 2, filter3);
				Stroke tmp3 = applyFilter(img_data, img_width, img_height, j / 2 + 1, i / 2 + 1, filter3);
				Stroke tmp4 = applyFilter(img_data, img_width, img_height, j / 2 + 2, i / 2 + 2, filter3);
				myStroke.r = tmp1.r / 4 + tmp2.r / 4 + tmp3.r / 4 + tmp4.r / 4;
				myStroke.g = tmp1.g / 4 + tmp2.g / 4 + tmp3.g / 4 + tmp4.g / 4;
				myStroke.b = tmp1.b / 4 + tmp2.b / 4 + tmp3.b / 4 + tmp4.b / 4;
				myStroke.a = tmp1.a / 4 + tmp2.a / 4 + tmp3.a / 4 + tmp4.a / 4;
			}
			new_img_data[new_offset_rgba + rr] = myStroke.r;
			new_img_data[new_offset_rgba + gg] = myStroke.g;
			new_img_data[new_offset_rgba + bb] = myStroke.b;
			new_img_data[new_offset_rgba + aa] = myStroke.a;
		}
	}

	delete[] rgb;
	img_data = new_img_data;

	mImageDst = QImage(img_data, new_img_width, new_img_height, QImage::Format_ARGB32);
	this->img_width = new_img_width;
	this->img_height = new_img_height;
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize(float scale)
{
	unsigned char *rgb = this->To_RGB();

	int new_img_width = img_width * scale;
	int new_img_height = img_height * scale;
	unsigned char *new_img_data = new unsigned char[new_img_width * new_img_height * 4];

	std::vector<std::vector<double>> filter = { {1,2,1},{2,4,2},{1,2,1} };

	for (int i = 0; i < new_img_height; i++)
	{
		for (int j = 0; j < new_img_width; j++)
		{
			int new_offset_rgba = i * new_img_width * 4 + j * 4;

			int src_j = img_width * j / new_img_width;
			int src_i = img_height * i / new_img_height;

			Stroke mystroke = applyFilter(img_data,img_width,img_height,src_j,src_i,filter);
			new_img_data[new_offset_rgba + rr] = mystroke.r;
			new_img_data[new_offset_rgba + gg] = mystroke.g;
			new_img_data[new_offset_rgba + bb] = mystroke.b;
			new_img_data[new_offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	img_data = new_img_data;

	mImageDst = QImage(img_data, new_img_width, new_img_height, QImage::Format_ARGB32);
	img_width = new_img_width;
	img_height = new_img_height;
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////

void Application::Rotate(float angleDegrees)
{
	unsigned char *rgb = this->To_RGB();
	double radian = angleDegrees * 3.1415926 / 180.f;

	int center_x = img_width / 2;
	int center_y = img_height / 2;
	img_data = new unsigned char[img_height * img_width * 4];

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{

			float delta_j = (j - center_x) * cos(radian) - (center_y - i) * sin(radian);
			float delta_i = (j - center_x) * sin(radian) + (center_y - i) * cos(radian);
			int src_j = center_x + delta_j;
			int src_i = center_y - delta_i;

			int dst_offset_rgba = i * img_width * 4 + j * 4;
			int src_offset_rgb = src_i * img_width * 3 + src_j * 3;


			if (src_j >= img_width || src_j < 0 || src_i >= img_height || src_i < 0)
			{
				img_data[dst_offset_rgba + rr] = 0;
				img_data[dst_offset_rgba + gg] = 0;
				img_data[dst_offset_rgba + bb] = 0;
				img_data[dst_offset_rgba + aa] = WHITE;
			}
			else 
			{
				img_data[dst_offset_rgba + rr] = rgb[src_offset_rgb + rr];
				img_data[dst_offset_rgba + gg] = rgb[src_offset_rgb + gg];
				img_data[dst_offset_rgba + bb] = rgb[src_offset_rgb + bb];
				img_data[dst_offset_rgba + aa] = WHITE;

			}
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge(QString filePath)
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image(int tMethod)
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	// img_data over img_data2
	if (img_height == img_height2 && img_width == img_width2)
	{
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				int offset_rgba = i * img_width * 4 + j * 4;
				// img_data[offset_rgba + k]  + img_data2[offset_rgba + k] * (img_data[offset_rgba + aa] - img_data[offset_rgba + aa]) / 255;
				// img_data[offset_rgba + k]  + img_data2[offset_rgba + k] * (255 - img_data[offset_rgba + aa]) / 255;
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = img_data[offset_rgba + k] * img_data[offset_rgba + aa] / 255 + img_data2[offset_rgba + k] * (255 - img_data[offset_rgba + aa]) / 255;
				img_data[offset_rgba + 3] = WHITE;
			}
		}
		mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
		renew();
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	// img_data in img_data2
	if (img_height == img_height2 && img_width == img_width2)
	{
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				int offset_rgba = i * img_width * 4 + j * 4;
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = img_data[offset_rgba + k] * img_data[offset_rgba + aa] / 255;
				img_data[offset_rgba + 3] = WHITE;
			}
		}
		mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
		renew();
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	// img_data out img_data2
	if (img_height == img_height2 && img_width == img_width2)
	{
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				int offset_rgba = i * img_width * 4 + j * 4;
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = img_data[offset_rgba + k] * (255 - img_data2[offset_rgba + aa]) / 255;
				img_data[offset_rgba + 3] = WHITE;
			}
		}
		mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
		renew();
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	// img_data atop img_data2
	if (img_height == img_height2 && img_width == img_width2)
	{
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				int offset_rgba = i * img_width * 4 + j * 4;
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = img_data[offset_rgba + k] * img_data2[offset_rgba + aa] / 255 + img_data2[offset_rgba + k] * (1 - img_data[offset_rgba + aa]) / 255;
				img_data[offset_rgba + 3] = WHITE;
			}
		}
		mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
		renew();
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	// img_data xor img_data2
	if (img_height == img_height2 && img_width == img_width2)
	{
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				int offset_rgba = i * img_width * 4 + j * 4;
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = img_data[offset_rgba + k] * (255 - img_data2[offset_rgba + aa]) / 255 + img_data2[offset_rgba + k] * (255 - img_data[offset_rgba + aa]) / 255;
				img_data[offset_rgba + 3] = WHITE;
			}
		}
		mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
		renew();
	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

// input rgb data; output rgb guassian data;
unsigned char *Application::getGaussianImgData(const unsigned char *sourceRGB,int n)
{
	// create Gaussian filter
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	for (int i = 0; i < img_width * img_height * 3; i++)
	{
		rgb[i] = sourceRGB[i];
	}
	vector<vector<double>> filter(n, vector<double>(n, 0));
	for (int i = 0; i < n; i++) {

		double tmpS = 1, tmpE = 1;
		for (int j = 0; j < i; j++)
			tmpS *= (n - 1 - j);
		for (int j = 1; j <= i; j++)
			tmpE *= j;

		filter[0][i] = tmpS / tmpE;
		filter[i][0] = tmpS / tmpE;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j != 0 && i != 0)
				filter[i][j] = filter[0][j] * filter[i][0];
		}
	}

	// apply filter
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			Stroke tmp = applyFilterRGB(rgb, img_width, img_height, j, i, filter);
			rgb[offset_rgb + rr] = tmp.r;
			rgb[offset_rgb + gg] = tmp.g;
			rgb[offset_rgb + bb] = tmp.b;
		}
	}
	return rgb;
}

unsigned char *Application::getGaussianImgData2(const unsigned char *sourceRGB, int n)
{
	double ** filter;
	filter = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
		filter[i] = (double *)malloc(n * sizeof(double));

	for (int i = 0; i < n; i++) {

		int tmpS = 1, tmpE = 1;
		for (int j = 0; j < i; j++)
			tmpS *= (n - 1 - j);
		for (int j = 1; j <= i; j++)
			tmpE *= j;

		filter[0][i] = tmpS / tmpE;
		filter[i][0] = tmpS / tmpE;
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j != 0 && i != 0)
				filter[i][j] = filter[0][j] * filter[i][0];
		}
	}

	bool edgeFlag = false;
	if (n == -5)
	{
		edgeFlag = true;
		n = 5;
	}
	unsigned char *rgb = this->To_RGB();
	unsigned char *outputData = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int startPixel = offset_rgb - 2 * img_width * 3 - 2 * 3;

			int sum = 0;

			if (!edgeFlag)
				for (int x = 0; x < n; x++)
				{
					for (int y = 0; y < n; y++)
					{
						sum += filter[x][y];
					}
				}
			else
				sum = 256;


			for (int k = 0; k < 3; k++)
			{
				int sttemp = startPixel + k;
				double colorSum = 0;
				for (int x = 0; x < n; x++)
				{
					for (int y = 0; y < n; y++)
					{
						int pixelAt = sttemp + x * img_width * 3 + y * 3;
						if (pixelAt < 0 || pixelAt>img_height*img_width * 3)
							colorSum += 0;
						else
							colorSum += rgb[pixelAt] * filter[x][y];
					}
				}
				if (colorSum < 0) {
					colorSum = 0;
				}
				outputData[offset_rgb + k] = colorSum / sum;
			}
		}
	}

	delete rgb;
	return outputData;
}


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	const unsigned char * sourceRGB = this->To_RGB();

	vector<int> radii = { 7, 3, 1 };  // 7 3 1
	for (int radius: radii)
	{
		unsigned char * reference__img = this->getGaussianImgData(sourceRGB, 2 * radius + 1);

		this->NPR_Paint_Layer(img_data, reference__img, radius);

		delete reference__img;
	}

	/*unsigned char * reference__img = this->getGaussianImgData(sourceRGB, 7 * 2 + 1);
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			img_data[offset_rgba + rr] = reference__img[offset_rgb + rr];
			img_data[offset_rgba + gg] = reference__img[offset_rgb + gg];
			img_data[offset_rgba + bb] = reference__img[offset_rgb + bb];
			img_data[offset_rgba + aa] = WHITE;
		}
	}*/

	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

// tCanvas RGBA; tReferenceImage Guassian RGB;
void Application::NPR_Paint_Layer(unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize)
{
	vector<Stroke> strokeList;
	float *distance = new float[img_width * img_height];

	// 參數
	float grid_2 = pow(tBrushSize * 2 + 1, 2.f);
	float threshold = 25.0f;

	int xStepSize = tBrushSize;
	int yStepSize = tBrushSize;

	// 計算與tCanvas 跟 tReferenceImage 的rgb距離
	ofstream ofs1, ofs2;
	ofs1.open(string("error") + to_string(tBrushSize) + ".csv") ;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset = i * img_width + j;
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			float dist = (float)sqrt(
				pow(tCanvas[offset_rgba + rr] - tReferenceImage[offset_rgb + rr], 2.f) +
				pow(tCanvas[offset_rgba + gg] - tReferenceImage[offset_rgb + gg], 2.f) +
				pow(tCanvas[offset_rgba + bb] - tReferenceImage[offset_rgb + bb], 2.f)
			);

			// debug
			distance[offset] = dist;
			ofs1 << dist << ',';
		}
		ofs1 << endl;
	}
	ofs1.close();


	ofs2.open(string("avgError") + to_string(tBrushSize) + ".csv");
	for (int i = tBrushSize; i < img_height; i += yStepSize)
	{
		for (int j = tBrushSize; j < img_width; j += xStepSize)
		{
			if (j + tBrushSize >= img_width || i + tBrushSize >= img_height)
			{
				ofs2 << 0.f << ',';
				continue;
			}
			else
			{
				float avg_error = 0.f;
				float max_error = 0.f;
				int x, y;
				float count = 0;
				// 取得在這個tBushSize的範圍內平均錯誤距離
				for (int new_i = i - tBrushSize; new_i <= i + tBrushSize; new_i++)
				{
					for (int new_j = j - tBrushSize; new_j <= j + tBrushSize; new_j++)
					{
						int new_offset = new_i * img_width + new_j;
						avg_error += distance[new_offset];
						if (distance[new_offset] >= max_error)
						{
							max_error = distance[new_offset];
							x = new_j; 
							y = new_i;
						}
						count++;
					}
				}
				assert(count == grid_2);
				avg_error =  avg_error / grid_2;

				ofs2 << avg_error << ',';
				// 如果平均錯誤距離高於門檻，就將當前範圍內error最大的座標當作stroke中心
				if (avg_error > threshold)
				{
					int offset_rgb = y * img_width * 3 + x * 3;
					strokeList.push_back(Stroke(tBrushSize, x, y, 
						tReferenceImage[offset_rgb + rr], 
						tReferenceImage[offset_rgb + gg], 
						tReferenceImage[offset_rgb + bb], 
						WHITE));
				}
			}
		}
		ofs2 << endl;
	}
	ofs2.close();

	// 打散畫筆順序
	random_shuffle(strokeList.begin(), strokeList.end());
	for (Stroke & stroke : strokeList)
	{
		this->Paint_Stroke(stroke);
	}


	delete distance;
}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke(const Stroke& s)
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++)
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++)
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height))
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared)
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				}
				else if (dist_squared == radius_squared + 1)
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}





///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
	radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}



