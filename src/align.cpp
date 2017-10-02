#include "align.h"
#include <string>

using std::string;
using std::cout;
using std::endl;
using std::get;

const int SHIFT = 15;

class convert
{
public:
	convert(Matrix<double> k):kernel(k){}

    std::tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;

        std::tuple<double, double, double> res;

        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                get<0>(res) += kernel(i, j)*get<0>(m(i, j));
                get<1>(res) += kernel(i, j)*get<1>(m(i, j));
                get<2>(res) += kernel(i, j)*get<2>(m(i, j));
            }
        }

        return res;
    }
    // Radius of neighbourhoud, which is passed to that operator
    static int radius;

    Matrix<double> kernel;
};

int convert::radius = 1;

Image one_dim_convert(Image src_image, Matrix<double> kernel, int radius, int dir){
	const int start_i = radius;
    const int n = (dir ? src_image.n_cols : src_image.n_rows);
    const int end_i = (dir ? src_image.n_rows : src_image.n_cols) - radius;

    Image dst_image(src_image.n_rows, src_image.n_cols);

    for(int j = 0; j < n; j++)
    	for(int i = start_i; i < end_i; i++){
    		std::tuple<double, double, double> sum = std::make_tuple(0, 0, 0);
    		for(int k = -radius; k <= radius; k++){
    			get<0>(sum) += (dir ? get<0>(src_image(i + k, j)) : get<0>(src_image(j, i + k)))*kernel(k + radius, 0);
    			get<1>(sum) += (dir ? get<1>(src_image(i + k, j)) : get<1>(src_image(j, i + k)))*kernel(k + radius, 0);
    			get<2>(sum) += (dir ? get<2>(src_image(i + k, j)) : get<2>(src_image(j, i + k)))*kernel(k + radius, 0);
    		}
    		(dir ? dst_image(i, j) : dst_image(j, i)) = sum;
    	}

    return dst_image;
}

std::tuple<int, int> get_dir(double theta){
	const double pi = 3.14;

	if(theta >= -pi/8 && theta < pi/8)
		return std::make_tuple(0, 1);

	if(theta >= pi/8 && theta < 3*pi/8)
		return std::make_tuple(1, 1);

	if(theta >= 3*pi/8 && theta < 5*pi/8)
		return std::make_tuple(1, 0);

	if(theta >= 5*pi/8 && theta < 7*pi/8)
		return std::make_tuple(1, -1);

	if(theta >= -3*pi/8 && theta < -pi/8)
		return std::make_tuple(-1, 1);

	if(theta >= -5*pi/8 && theta < -3*pi/8)
		return std::make_tuple(-1, 0);

	if(theta >= -7*pi/8 && theta < -5*pi/8)
		return std::make_tuple(-1, -1);

	if((theta >= 7*pi/8 && theta <= pi) || (theta >= -pi && theta < -7*pi/8))
		return std::make_tuple(0, -1);

	return std::make_tuple(0, 1);
}

std::tuple<int, int> MSE(Image image1, Image image2){

	int off_x = 0;
	int off_y = 0;
	uint min_sum = std::numeric_limits<uint>::max();

	Image im1;
	Image im2;

	for(int i = -SHIFT + 1; i < SHIFT; i++)
		for(int j = -SHIFT + 1; j < SHIFT; j++){
			int x = abs(i);
			int	y = abs(j);
			uint sum = 0;

			im1 = image1.submatrix(y*(j <= 0), x*(i >= 0), image1.n_rows - y, image1.n_cols - x);
			im2 = image2.submatrix(y*(j >= 0), x*(i <= 0), image1.n_rows - y, image1.n_cols - x);

			for(uint k1 = 0; k1 < im1.n_rows; k1++)
				for(uint k2 = 0; k2 < im1.n_cols; k2++)
					sum += (get<0>(im1(k1, k2)) - get<0>(im2(k1, k2)))*(get<0>(im1(k1, k2)) - get<0>(im2(k1, k2)));

			if(sum < min_sum){
				min_sum = sum;
				off_x = i;
				off_y = j;
			}
		}

	return std::make_tuple(off_x, off_y);
}

std::tuple<int, int> COR(Image image1, Image image2, int shift){

	int off_x = 0;
	int off_y = 0;
	uint max_sum = 0;

	Image im1;
	Image im2;

	for(int i = -SHIFT + 1; i < SHIFT; i++)
		for(int j = -SHIFT + 1; j < SHIFT; j++){
			int x = abs(i);
			int	y = abs(j);
			uint sum = 0;

			im1 = image1.submatrix(y*(j <= 0), x*(i >= 0), image1.n_rows - y, image1.n_cols - x);
			im2 = image2.submatrix(y*(j >= 0), x*(i <= 0), image1.n_rows - y, image1.n_cols - x);

			for(uint k1 = 0; k1 < im1.n_rows; k1++)
				for(uint k2 = 0; k2 < im1.n_cols; k2++)
					sum += get<0>(im1(k1, k2)) * get<0>(im2(k1, k2));

			if(sum > max_sum){
				max_sum = sum;
				off_x = i;
				off_y = j;
			}
		}

	return std::make_tuple(off_x, off_y);
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
	Image image_b = srcImage.submatrix(0, 0, srcImage.n_rows/3, srcImage.n_cols);
	Image image_g = srcImage.submatrix(srcImage.n_rows/3, 0, srcImage.n_rows/3, srcImage.n_cols);
	Image image_r = srcImage.submatrix(2*srcImage.n_rows/3, 0, srcImage.n_rows/3, srcImage.n_cols);

	Image dstImage(image_r.n_rows + 2*SHIFT, image_r.n_cols + 2*SHIFT);

	for(uint i = 0; i < image_r.n_rows; i++)
		for(uint j = 0; j < image_r.n_cols; j++)
			get<0>(dstImage(i + SHIFT, j + SHIFT)) = get<0>(image_r(i, j));

	std::tuple<int, int> offset = MSE(image_r, image_g);

	for(uint i = 0; i < image_g.n_rows; i++)
		for(uint j = 0; j < image_g.n_cols; j++)
			get<1>(dstImage(i + SHIFT - get<1>(offset), j + SHIFT + get<0>(offset))) = get<1>(image_g(i, j));

	offset = MSE(image_r, image_b);

	for(uint i = 0; i < image_b.n_rows; i++)
		for(uint j = 0; j < image_b.n_cols; j++)
			get<2>(dstImage(i + SHIFT - get<1>(offset), j + SHIFT + get<0>(offset))) = get<2>(image_b(i, j));

    return dstImage.submatrix(SHIFT, SHIFT, image_r.n_rows, image_r.n_cols);
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    convert::radius = 1;
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    convert::radius = 1;
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
	Matrix<double> kernel = {{ -1/6,  -2/3,  -1/6},
                             { -2/3,  13/3,  -2/3},
                             {-1/6, -2/3, -1/6}};
	convert::radius = 1;
    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {
	uint r, g, b, sum_r = 0, sum_g = 0, sum_b = 0; 

	for(uint i = 0; i < src_image.n_rows; i++)
		for(uint j = 0; j < src_image.n_cols; j++){
			 std::tie(r, g, b) = src_image(i, j);
             sum_r += r;
             sum_g += g;
             sum_b += b;
		}

	uint norm = src_image.n_rows * src_image.n_cols;

	sum_r /= norm;
	sum_g /= norm;
	sum_b /= norm;

	uint S = (sum_r + sum_g + sum_b)/3;

	Image dst_image = src_image;

	for(uint i = 0; i < dst_image.n_rows; i++)
		for(uint j = 0; j < dst_image.n_cols; j++){
			get<0>(dst_image(i, j)) *= S/sum_r;
			get<1>(dst_image(i, j)) *= S/sum_g;
			get<2>(dst_image(i, j)) *= S/sum_b;
		}

    return dst_image;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    Image dst_image = src_image.unary_map(convert(kernel));

    return dst_image;
}

Image autocontrast(Image src_image, double fraction) {
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
	int size = 2*radius + 1;
	Matrix<double> kernel(size, size);
	double sum = 0;

	for(int i = -size/2; i < size/2 + 1; i++)
		for(int j = -size/2; j < size/2 + 1; j++){
			kernel(i + size/2, j + size/2) = exp(-(i*i + j*j)/(2*sigma*sigma))/(2*3.14*sigma*sigma);
			sum += kernel(i + size/2, j + size/2);
		}

	for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++)
			kernel(i, j) /= sum;

	convert::radius = radius;
	Image dst_image = custom(src_image, kernel);

    return dst_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
	int size = 2*radius + 1;
	Matrix<double> kernel(size, 1);
	double sum = 0;

	for(int i = -size/2; i < size/2 + 1; i++){
		kernel(i + size/2, 0) = exp(-(i*i)/(2*sigma*sigma))/sqrt((2*3.14*sigma*sigma));
		sum += kernel(i + size/2, 0);
	}

	for(int i = 0; i < size; i++)
		kernel(i, 0) /= sum;

	Image dst_image = one_dim_convert(src_image, kernel, radius, 0);
	dst_image = one_dim_convert(dst_image, kernel, radius, 1);

    return dst_image;
}

Image median(Image src_image, int radius) {
    return src_image;
}

Image median_linear(Image src_image, int radius) {
    return src_image;
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

Image canny(Image src_image, int threshold1, int threshold2) {
	Image dst_image = gaussian_separable(src_image, 1.4, 2);

	Image imIx = sobel_x(dst_image);
	Image imIy = sobel_y(dst_image);

	Matrix<double> G(dst_image.n_rows, dst_image.n_cols);
	Matrix<double> theta(dst_image.n_rows, dst_image.n_cols);

	for(uint i = 0; i < dst_image.n_rows; i++)
		for(uint j = 0; j < dst_image.n_cols; j++){
			uint x = get<0>(imIx(i, j));
			uint y = get<0>(imIy(i, j));

			G(i, j) = sqrt(x*x + y*y);
			theta(i, j) = atan2(y, x);
		}

	Matrix<double> G_tres(G.n_rows, G.n_cols);
	Matrix<int> G_hyst(G.n_rows, G.n_cols);

	for(uint i = 1; i < G.n_rows - 1; i++)
		for(uint j = 1; j < G.n_cols - 1; j++){
			std::tuple<int, int> vec = get_dir(theta(i, j));

			double f_pix = G(i + get<0>(vec), j + get<1>(vec));
			double s_pix = G(i - get<0>(vec), j - get<1>(vec));

			if(G(i, j) >= f_pix && G(i, j) >= s_pix){
				if(G(i, j) > threshold2){
					G_tres(i, j) = G(i, j);
					G_hyst(i, j) = 2;
				}else if(G(i, j) > threshold1){
					G_hyst(i, j) = 1;
					G_tres(i, j) = 0;
				}
				else
					G_tres(i, j) = G_hyst(i, j) = 0;
			}
			else
				G_tres(i, j) = 0;
		}

	
	
    return dst_image;
}
