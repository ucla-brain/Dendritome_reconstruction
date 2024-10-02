//
// Created by muyezhu on 3/3/17.
//

#ifndef MCP3D_LINE_FILTER_HPP
#define MCP3D_LINE_FILTER_HPP

#include <cstdlib>
#include <opencv2/core/core.hpp>
#include <vector>
#include <unordered_set>


namespace mcp3d {


const int DEFAULT_WINDOW_SIZE = 5;
const std::vector<float> DEFAULT_SIGMAS = {1.0, 3.0, 5.0};
/// more negative values favor curve lines with smaller diameter
const double GAUSSIAN_DERIV_NORMALIZER = -2;
const double ZERO_EIGEN_TOL = 1e-1;
const double LAMBDA2_HIGH = -1.0;
const double LAMBDA2_LOW = -20.0;
const double BETA = 0.5;
const double DEFAULT_FRANGI3D_RESCALE = 2.0;
const double W_LINE_3D = 1.0;
const double W_BLOB_3D = 1.0;
const double W_L2_3D = 2.0;

typedef std::pair<cv::Mat, cv::Mat> separate_filters_2d;

class Frangi2D
{
        public:
        /// construction of instance: input image, if input image is to be filtered in place,
        /// size of local window (must be odd if provided), scales to use (at scale t, t = sqrt(sigma))
        /// if k < 0, a window size is calculated for each scale
        Frangi2D(cv::Mat &input, int k = DEFAULT_WINDOW_SIZE, bool ipl = false,
        const std::vector<float> &_sigmas = DEFAULT_SIGMAS);

        /// matrix holding the lambda2 with the minimum value across scales
        /// the two eighen value l1 and l2, with l1 > l2, should satisfy l1 ~= 0, l2 << l1.
        /// for l1 not ~= 0, store 1 at position to indicate such. otherwise, store min(l2|scales)
        /// initialize to 1s before updating values
        cv::Mat minLambda2;

        cv::Mat maxVesselness;

        cv::Mat eigen_filtered;

        cv::Mat vesselness_filtered;

        /// at each scale, calls fill_level_hessians and updateImgMinLambda2
        void compute_hessian_eigen_through_scales();

        void compute_vesselness_through_scales();

        /// anisotropically filter the input image with eigen value response: I' = I * exp(max|minLambda2|)
        /// alternatively, I' = I * (1 / (1 + exp(-max|minLambda2|)))
        /// if ipl is true, the original input image is overwritten with the filtered image
        /// in DEV mode as defined in line_filter.cpp, the minLambda2 matrix is not modified by remapImgMinLambda2,
        /// instead, a new array is allocated to hold the remapped minLambda2 values.
        /// results from the filtering is adjusted to [0, 255], cast to uchar, and stored in eigen_viltered matrix
        void minLambda2_filter(double remap_factor = 0.15);

        void vesselness_filter(double remap_factor);

        std::vector<float> get_sigmas();

        std::pair<int, int> get_img_dimension();

        int get_window_size();

        std::tuple<cv::Mat, cv::Mat, cv::Mat> get_hessians();

        private:
        cv::Mat &input_img;
        int window_size;
        unsigned long numScales;
        int rows, cols;

        /// inplace: determines if the output filtered matrix will modify input_img matrix inplace
        bool inplace;

        bool dynamic_window_size;

        std::vector<float> sigmas;

        /* matrices holding image convolution with gaussian derivatives
        * the hessian matrix for each pixel at (i, j) is then [I_xx(i, j), I_xy(i, j); I_xy(i, j), I_yy(i, j)]
        * assuming I_xy = I_yx
        */
        cv::Mat I_xx;
        cv::Mat I_xy;
        cv::Mat I_yy;

        /// return the row and column filter for gaussian derivative to the order of dx and dy
        /// filters are normalized at different scales. first order derivative is normalized by
        /// sigma ^ GAUSSIAN_DERIV_NORMALIZER, 2nd order by sigma ^ (2 * GAUSSIAN_DERIV_NORMALIZER)
        separate_filters_2d get_level_kernels(int level, int dx, int dy);

        /// apply separate gaussian derivative filters to image
        cv::Mat apply_filters(separate_filters_2d filters);

        /// at a given scale/level, obtain image second order derivatives needed for a pixel hessian matrix
        void fill_level_hessians(int level);

        /// initialize minLambda2 as a matrix filled with 1s
        void initMinLambda2();

        /// initialize maxVesselness as a matrix filled with -100
        void initMaxVesselness();

        /// hessian l1 norm is approximated as sum (abs(elements in hessian))
        /// the greatest l1 norm acroos pixels is returned
        double get_max_hessian_l1_norm();

        /// update vesselness where the eigen value of greater magnitude is negative
        /// vesselness is computed for the pixel (i, j) and if greater than current entry in maxVesselness, replace it
        void updatePixelMaxVesselness(int i, int j, double max_hessian_l1, int level);

        void updateImgMaxVesselness(int level);

        /// compute vesselness as lambda2 / max_hessian_l1_norm
        /// where lambda2 is the eigenvalue with greater magnitude
        /// max_hessian_l1_norm is calculated for current level
        double computePixelVesselness(double lambda1, double lambda2, double max_hessian_l1);

        /// Update min eigen value across scales for pixel (i, j). the entry is only updated if at
        /// curren level, max(lambda1, lambda2) is smaller in magnitude than ZERO_EIGEN_TOL, and
        /// min(lambda1, lambda2) is smaller than current minLambda2(i, j)
        /// Because the hessian matrix is symmetrical assuming I_xy = I_yx, the 2x2 matrix is gauranteed
        /// to have all real solutions. use opencv eigen() to find eigenvalues
        /// \param i
        /// \param j
        void updatePixelMinLambda2(int i, int j);

        /// at each pixel, call updatePixelMinLambda2
        void updateImgMinLambda2();

        /// remap minLambda2: values between lambda2_low and lambda2_high are multiplied by remap_factor
        /// values lower than lambda2_low are replaced by remap_factor * lambda2_low (to prevent numerical underflow)
        /// values higher than lambda2_high are replaced by 10 (any positive integer should do)
        void remapImgMinLambda2(double remap_factor, cv::Mat &dest);

        /// values in maxVesselness multiplied by remap_factor but clamped above at 20
        void remapImgMaxVesselness(double remap_factor, cv::Mat &dest);
};

class Frangi3D
{
    public:
        /// the paths to all images contained in the whole 3D dataset
        std::vector<std::string> img_path_list;
        /// path to "anchor image", aka the image to filter
        std::string anchor_img_path;
        /// index of "anchor image", aka the image to filter, in the vector img_path_list
        int anchor_img_index;
        /// maximum vesselness at each pixel, computed through scale spaces
        cv::Mat max_vesselness;
        /// anisotropically filtered anchor image (by vesselness)
        cv::Mat filtered;

        bool filter_img(const std::string &img_path, double s = DEFAULT_FRANGI3D_RESCALE);

        bool filter_imgs(const std::vector<std::string>& img_paths, double s = DEFAULT_FRANGI3D_RESCALE);

        bool filter_all_imgs(double s);

        /// set values for anchor_img_path, anchor_img_index, prev_anchor_img_index. does not read imgs all alter
        /// img_stack / img_stack_vec
        bool init_stack(const std::string& img_path);

        /// calculates maximum vesselness at each pixel throughout scales space
        /// (1) call build_stack() to read necessary images and set anchor image
        /// (2) allocate the max_vesselness matrix and initilize all values to -10
        /// (3) at each scale space, call fill_hessians() to compute the image derivatives
        /// (4) at each scale space, call update_img_vesselness() to update the max_vesselness matrix based on
        ///     vesselness calculated at current scale
        void compute_vesselness();

        /// filter image by replacing pixel intensity I(i, j) as I(i, j) * exp(rescale(vesselness(i, j), s))
        /// normalize result to range [0, 255] and convert to uchar, store in filtered matrix
        void anisotropic_vesselness_filter(double s);

        /// check if each image path in img_list is valid.
        /// the sequence of images are loaded and processed sequentially.
        /// the initialization read the first window_size number of images, starting at the half_window_size th image
        /// in the sequence. after the image at the kernel anchor level is processed, the earliest image in sequence
        /// read in memory is released, and next image in sequence is read
        explicit Frangi3D(std::vector<std::string> l, int k = DEFAULT_WINDOW_SIZE, const std::vector<float> &_sigmas = DEFAULT_SIGMAS);

    private:
        std::unordered_set<std::string> bad_files;

        int window_size, half_window_size;
        std::vector<float> sigmas;
        unsigned long num_scales;
        int rows, cols;
        int prev_anchor_img_index;

        std::vector<cv::Mat> img_stack_vec;
        cv::Mat img_stack;
        cv::Mat anchor_img;
        cv::Mat I_xx, I_xy, I_xz, I_yy, I_yz, I_zz;

        double max_hessian_l1;

        cv::Mat get_gaussian1D_derivatives(int dx, int level);

        void fill_hessians(int level);

        void get_max_hessian_l1();

        void init_max_vesselness();

        void build_stack_from_empty();

        void build_stack();

        void update_pixel_vesselness(int i, int j, int level);

        void update_img_vesselness(int level);

        void rescale_max_vesselness(double s, cv::Mat& dst);

        void write_filtered_img(std::string dir = "");
};

}
#endif //MCP3D_LINE_FILTER_HPP
