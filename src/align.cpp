#include "align.h"
#include <string>
#include <set>

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

        double res_r = 0, res_g = 0, res_b = 0;

        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                res_r += kernel(i, j)*get<0>(m(i, j));
                res_g += kernel(i, j)*get<1>(m(i, j));
                res_b += kernel(i, j)*get<2>(m(i, j));
            }
        }

        res_r = (res_r > 0 ? res_r : 0);
        res_r = (res_r < 255 ? res_r : 255);

        res_g = (res_g > 0 ? res_g : 0);
        res_g = (res_g < 255 ? res_g : 255);

        res_b = (res_b > 0 ? res_b : 0);
        res_b = (res_b < 255 ? res_b : 255);

        return std::make_tuple(res_r, res_g, res_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    static int radius;

    Matrix<double> kernel;
};

int convert::radius = 1;

//My functions
Image mirror(Image src_image){
	uint shift = convert::radius;

	Image dst_image(src_image.n_rows + 2*shift, src_image.n_cols + 2*shift);

	for(uint i = shift; i < src_image.n_rows + shift; i++)
		for(uint j = shift; j < src_image.n_cols + shift; j++)
			dst_image(i, j) = src_image(i - shift, j - shift);

	for(uint i = 0; i < shift; i++)
		for(uint j = 0; j < shift; j++)
			dst_image(i, j) = src_image(0, 0);

	for(uint i = 0; i < shift; i++)
		for(uint j = src_image.n_cols - shift; j < src_image.n_cols; j++)
			dst_image(i, j) = src_image(0, src_image.n_cols - 1);

	for(uint i = src_image.n_rows - shift; i < src_image.n_rows; i++)
		for(uint j = 0; j < shift; j++)
			dst_image(i, j) = src_image(src_image.n_rows - 1, 0);

	for(uint i = src_image.n_rows - shift; i < src_image.n_rows; i++)
		for(uint j = src_image.n_cols - shift; j < src_image.n_cols; j++)
			dst_image(i, j) = src_image(src_image.n_rows - 1, src_image.n_cols - 1);

	for(uint i = shift; i < src_image.n_cols + shift; i++)
		for(uint j = 0; j < shift; j++)
			dst_image(j, i) = src_image(j, i - shift);

	for(uint i = shift; i < src_image.n_cols + shift; i++)
		for(uint j = src_image.n_rows - shift; j < src_image.n_rows; j++)
			dst_image(j + shift, i) = src_image(j, i - shift);

	for(uint i = shift; i < src_image.n_rows + shift; i++)
		for(uint j = 0; j < shift; j++)
			dst_image(i, j) = src_image(i - shift, j);

	for(uint i = shift; i < src_image.n_rows + shift; i++)
		for(uint j = src_image.n_cols - shift; j < src_image.n_cols; j++)
			dst_image(i, j + shift) = src_image(i - shift, j);

	return dst_image;
}

//Median
std::tuple<uint, uint, uint> sort_neighb(Image n){
	int size = n.n_cols * n.n_rows;
	uint* mas_r = new uint[size];
	uint* mas_g = new uint[size];
	uint* mas_b = new uint[size];

	for(uint i = 0; i < n.n_rows; i++)
		for(uint j = 0; j < n.n_cols; j++){
			mas_r[i*n.n_cols + j] = get<0>(n(i, j));
			mas_g[i*n.n_cols + j] = get<1>(n(i, j));
			mas_b[i*n.n_cols + j] = get<2>(n(i, j));
		}

	for(int i = 0; i < size; i++)
		for(int j = 0; j < size - i - 1; j++){
			if(mas_r[j] > mas_r[j + 1]){
				uint tmp = mas_r[j];
				mas_r[j] = mas_r[j + 1];
				mas_r[j + 1] = tmp;
			}
			if(mas_g[j] > mas_g[j + 1]){
				uint tmp = mas_g[j];
				mas_g[j] = mas_g[j + 1];
				mas_g[j + 1] = tmp;
			}
			if(mas_b[j] > mas_b[j + 1]){
				uint tmp = mas_b[j];
				mas_b[j] = mas_b[j + 1];
				mas_b[j + 1] = tmp;
			}
		}

	uint median = size/2;

	return std::make_tuple(mas_r[median], mas_g[median], mas_b[median]);
}

//Canny
std::vector<std::set<int>> sets;

void find_strong_labels(std::vector<int>& strong_labels, std::set<int>& strong_pix){
	for(uint i = 0; i < sets.size(); i++){
		auto end = sets[i].end();
		for(uint j = 0; j < strong_labels.size(); j++)
			if(sets[i].find(strong_labels[j]) != end){
				strong_pix.insert(sets[i].begin(), end);
				break;
			}
		}
}

void add_eq(int A, int B, int C, int D){
	for(uint i = 0; i < sets.size(); i++){
		auto end = sets[i].end();
		if(sets[i].find(A) != end ||
		   sets[i].find(B) != end ||
		   sets[i].find(C) != end ||
		   sets[i].find(D) != end)
		{
			if(A) sets[i].insert(A);
			if(B) sets[i].insert(B);
			if(C) sets[i].insert(C);
			if(D) sets[i].insert(D);
			return;
		}
	}

	std::set<int> new_set;
	if(A) new_set.insert(A);
	if(B) new_set.insert(B);
	if(C) new_set.insert(C);
	if(D) new_set.insert(D);

	sets.push_back(new_set);
}

Matrix<int> get_connected(Matrix<int> mat, std::set<int>& strong_pix){
	Matrix<int> res(mat.n_rows, mat.n_cols);
	std::vector<int> strong_labels;
	uint label = 0;

	for(uint i = 0; i < mat.n_rows; i++)
		for(uint j = 0; j < mat.n_cols; j++){
			int Z = mat(i, j);
			
			if(!Z){
				res(i, j) = 0;
				continue;
			}

			int A = (j == 0? 0 : res(i, j - 1));
			int B = (j == 0 || i == 0? 0 : res(i - 1, j - 1));
			int C = (i == 0? 0 : res(i - 1, j));
			int D = (i == 0 || j == mat.n_cols - 1? 0 : res(i - 1, j + 1));

			if(A + B + C + D == 0)
				res(i, j) = ++label;
			else {
				add_eq(A, B, C, D);

				if(A)
					res(i, j) = A;
				else if(B)
					res(i, j) = B;
				else if(C)
					res(i, j) = C;
				else if(D)
					res(i, j) = D;
			}

			if(Z == 2)
				strong_labels.push_back(label);
		}

	for(uint i = 0; i < strong_labels.size(); i++)
		strong_pix.insert(strong_labels[i]);

	find_strong_labels(strong_labels, strong_pix);

	return res;
}

uint get_maxs_u_d(Matrix<int> G, bool ax){
	uint max1 = 0, max2 = 0;
	uint n_max1 = 0, n_max2 = 0;

	uint r = G.n_rows;
	uint c = G.n_cols;

	for(uint i = (ax?0:0.95*r); i < (ax?0.05*r:r); i++){
		uint sum = 0;
		for(uint j = 0; j < c; j++)
			if(G(i, j) == 3)
				sum++;

		if(sum > max1){
			max2 = max1;
			n_max2 = n_max1;

			max1 = sum;
			n_max1 = i;
		}else if(sum > max2){
			max2 = sum;
			n_max2 = i;
		}
	}

	uint n;

	if(ax)
		n = (n_max1 > n_max2 ? n_max1 : n_max2);
	else
		n = (n_max1 < n_max2 ? n_max1 : n_max2);

	return n;
}

uint get_maxs_l_r(Matrix<int> G, bool ax){
	uint max1 = 0, max2 = 0;
	uint n_max1 = 0, n_max2 = 0;

	uint r = G.n_rows;
	uint c = G.n_cols;

	for(uint i = (ax?0:0.95*c); i < (ax?0.05*r:c); i++){
		uint sum = 0;
		for(uint j = 0; j < r; j++)
			if(G(j, i) == 3)
				sum++;

		if(sum > max1){
			max2 = max1;
			n_max2 = n_max1;

			max1 = sum;
			n_max1 = i;
		}else if(sum > max2){
			max2 = sum;
			n_max2 = i;
		}
	}

	uint n;

	if(ax)
		n = (n_max1 > n_max2 ? n_max1 : n_max2);
	else
		n = (n_max1 < n_max2 ? n_max1 : n_max2);

	return n;
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

//Gause one dim
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

//Metrics
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

//Made before
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

	dstImage = dstImage.submatrix(SHIFT, SHIFT, image_r.n_rows, image_r.n_cols);

	if(isPostprocessing)
	{
		if(isMirror)
			dstImage = mirror(dstImage);

		if(postprocessingType == "--gray-world")
			dstImage = gray_world(dstImage);
		else if(postprocessingType == "--unsharp")
			dstImage = unsharp(dstImage);
		else if(postprocessingType == "--autocontrast")
			dstImage = autocontrast(dstImage, fraction);
	
		if(isMirror)
			dstImage = dstImage.submatrix(1, 1, image_r.n_rows, image_r.n_cols);
	}

    return dstImage;
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
	Matrix<double> kernel = {{ -1./6,  -2./3,  -1./6},
                             { -2./3,  13./3,  -2./3},
                             {-1./6, -2./3, -1./6}};
	convert::radius = 1;
    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {
	uint r_t, g_t, b_t, sum_r = 0, sum_g = 0, sum_b = 0; 

	for(uint i = 0; i < src_image.n_rows; i++)
		for(uint j = 0; j < src_image.n_cols; j++){
			 std::tie(r_t, g_t, b_t) = src_image(i, j);
             sum_r += r_t;
             sum_g += g_t;
             sum_b += b_t;
		}

	uint norm = src_image.n_rows * src_image.n_cols;

	sum_r /= norm;
	sum_g /= norm;
	sum_b /= norm;

	double S = (sum_r + sum_g + sum_b)/3;

	Image dst_image = src_image;

	for(uint i = 0; i < dst_image.n_rows; i++)
		for(uint j = 0; j < dst_image.n_cols; j++){
			uint r = get<0>(dst_image(i, j)) * S/sum_r;
			uint g = get<1>(dst_image(i, j)) * S/sum_g;
			uint b = get<2>(dst_image(i, j)) * S/sum_b;

			r = (r > 0 ? r : 0);
	        r = (r < 255 ? r : 255);

	        g = (g > 0 ? g : 0);
	        g = (g < 255 ? g : 255);

	        b = (b > 0 ? b : 0);
	        b = (b < 255 ? b : 255);
		
	        dst_image(i, j) = std::make_tuple(r, g, b);
		}

    return dst_image;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {

	Image dst_image = mirror(src_image);

    dst_image = dst_image.unary_map(convert(kernel));

    return dst_image.submatrix(convert::radius, convert::radius, src_image.n_rows, src_image.n_cols);
}

Image autocontrast(Image src_image, double fraction) {

	uint* hist = new uint[256];
	double min_v = 0, max_v = 255;
	bool min_flag = false, max_flag = false;
	uint min_n = fraction*src_image.n_rows*src_image.n_cols;
	uint max_n = (1 - fraction)*src_image.n_rows*src_image.n_cols;

	for(int i = 0; i < 256; i++)
		hist[i] = 0;

	for(uint i = 0; i < src_image.n_rows; i++)
		for(uint j = 0; j < src_image.n_cols; j++){
			uint intense = 0.2125*get<0>(src_image(i, j)) + 0.7154*get<1>(src_image(i, j)) + 0.0721*get<2>(src_image(i, j));
			hist[intense]++;
		}

	for(int i = 1; i < 256; i++){
		hist[i] += hist[i-1];
		
		if(hist[i] > min_n && !min_flag){
			min_v = i;
			min_flag = true;
		}

		if(hist[i] > max_n && !max_flag){
			max_v = i;
			max_flag = true;
		}
	}

	double b_c = 255.0/(1 - max_v/min_v);
	double a = -b_c/min_v;

	Image dst_image = src_image;

	for(uint i = 0; i < src_image.n_rows; i++)
		for(uint j = 0; j < src_image.n_cols; j++){
			int r = get<0>(src_image(i, j))*a + b_c;
			int g = get<1>(src_image(i, j))*a + b_c;
			int b = get<2>(src_image(i, j))*a + b_c;

			r = (r > 0 ? r : 0);
	        r = (r < 255 ? r : 255);

	        g = (g > 0 ? g : 0);
	        g = (g < 255 ? g : 255);

	        b = (b > 0 ? b : 0);
	        b = (b < 255 ? b : 255);

	        dst_image(i, j) = std::make_tuple(r, g, b);
		}

    return dst_image;
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

	Image dst_image = mirror(src_image);

	dst_image = one_dim_convert(dst_image, kernel, radius, 0);
	dst_image = one_dim_convert(dst_image, kernel, radius, 1);

    return dst_image.submatrix(radius, radius, src_image.n_rows, src_image.n_cols);
}

Image median(Image src_image, int radius) {
	uint size = 2*radius + 1;

	uint start_i = radius;
	uint start_j = radius;
	uint end_i = src_image.n_rows - radius;
	uint end_j = src_image.n_cols - radius;

	Image dst_image(src_image.n_rows, src_image.n_cols);

	for(uint i = start_i; i < end_i; i++)
		for(uint j = start_j; j < end_j; j++)
			dst_image(i, j) = sort_neighb(src_image.submatrix(i - radius, j - radius, size, size));

    return dst_image;
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


    std::set<int> strong_pix;

   	Matrix<int> labels = get_connected(G_hyst, strong_pix);

    auto end = strong_pix.end();

   	for(uint i = 0; i < G.n_rows; i++)
		for(uint j = 0; j < G.n_cols; j++)
			if(strong_pix.find(labels(i, j)) != end && G_hyst(i, j) == 1)
				G_hyst(i, j) = 3;

	uint n_u = get_maxs_u_d(G_hyst, 1);
	uint n_d = get_maxs_u_d(G_hyst, 0);
	uint n_l = get_maxs_l_r(G_hyst, 1);
	uint n_r = get_maxs_l_r(G_hyst, 0);

    return src_image.submatrix(n_u, n_l, n_d - n_u, n_r - n_l);
}
