#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>

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
	unsigned char *rgb = this->To_RGB();



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

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N(unsigned int N)
{

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{

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

Stroke applyFilter(unsigned char *img_data_rgba, int img_width, int img_height, int x, int y, std::vector<std::vector<int>> &filter)
{
	int filterSizeY = (int)filter.size();
	int filterSizeX = (int)filter.front().size();
	Stroke result(0,0,0,0,0,0,0);
	int sum = 0;
	for (std::vector<int> &a : filter)
	{
		for (int &b : a)
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

void Application::Half_Size()
{
	unsigned char *rgb = this->To_RGB();

	int new_img_width = img_width  / 2;
	int new_img_height = img_height / 2;
	unsigned char *new_img_data = new unsigned char[new_img_width * new_img_height * 4];

	std::vector<std::vector<int>> filter = { {1,2,1},{2,4,2},{1,2,1} };
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
			/*int rrr = 0, ggg = 0, bbb = 0, aaa = 0;
			for (int k = -1; k < 2; k++)
			{
				for (int r = -1; r < 2; r++)
				{
					int origin_rgba_offset = (2 * i + r) * img_width * 4 + (j * 2 + k) * 4;

					if (i * 2 + r < 0 || i * 2 + r >= img_height || 2 * j + k < 0 || 2 * j + k >= img_width)
					{
						continue;
					}
					else
					{
						int coefficent = 1;
						if (k != 0 || r != 0)
						{
							if (abs(k + r) == 1)
							{
								coefficent = 2;
							}
						}
						else
						{
							coefficent = 4;
						}
						rrr += img_data[origin_rgba_offset + rr] * coefficent;
						ggg += img_data[origin_rgba_offset + gg] * coefficent;
						bbb += img_data[origin_rgba_offset + bb] * coefficent;
						aaa += img_data[origin_rgba_offset + aa] * coefficent;
					}
				}
			}
			new_img_data[new_offset_rgba + rr] = rrr / 16;
			new_img_data[new_offset_rgba + gg] = ggg / 16;
			new_img_data[new_offset_rgba + bb] = bbb / 16;
			new_img_data[new_offset_rgba + aa] = aaa / 16;*/
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

	std::vector<std::vector<int>> filter1 = { {1,2,1},{2,4,2}, {1,2,1} };
	std::vector<std::vector<int>> filter2 = { {1,3,3,1},{3,9,9,3},{3,9,9,3},{1,3,3,1} };
	std::vector<std::vector<int>> filter3 = { {1,2,1},{3,6,3},{3,6,3},{1,2,1} };
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


	for (int i = 0; i < new_img_height; i++)
	{
		for (int j = 0; j < new_img_width; j++)
		{
			int new_offset_rgba = i * new_img_width * 4 + j * 4;


			int rrr = 0, ggg = 0, bbb = 0, aaa = 0;
			for (int k = -1; k < 2; k++)
			{
				for (int r = -1; r < 2; r++)
				{
					int origin_rgba_offset = (i / scale + r) * img_width * 4 + (j / scale + k) * 4;

					if (i / scale + r < 0 || i / scale + r >= img_height || j / scale + k < 0 || j / scale + k >= img_width)
					{
						continue;
					}
					else
					{
						int coefficent = 1;
						if (k != 0 || r != 0)
						{
							if (abs(k + r) == 1)
							{
								coefficent = 2;
							}
						}
						else
						{
							coefficent = 4;
						}
						rrr += img_data[origin_rgba_offset + rr] * coefficent;
						ggg += img_data[origin_rgba_offset + gg] * coefficent;
						bbb += img_data[origin_rgba_offset + bb] * coefficent;
						aaa += img_data[origin_rgba_offset + aa] * coefficent;
					}
				}
			}

			new_img_data[new_offset_rgba + rr] = rrr / 16;
			new_img_data[new_offset_rgba + gg] = ggg / 16;
			new_img_data[new_offset_rgba + bb] = bbb / 16;
			new_img_data[new_offset_rgba + aa] = aaa / 16;

		}
	}

	delete[] rgb;
	delete[] img_data;
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
	if (img_height == img_height2 && img_width == img_width2)
	{

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
	if (img_height == img_height2 && img_width == img_width2)
	{

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
	if (img_height == img_height2 && img_width == img_width2)
	{

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
	if (img_height == img_height2 && img_width == img_width2)
	{

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
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

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
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::NPR_Paint_Layer(unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize)
{

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



