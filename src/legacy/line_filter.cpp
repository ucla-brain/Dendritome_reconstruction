//
// Created by muyezhu on 3/1/17.
//
#include <iostream>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "line_filter.hpp"
#include "common/mcp3d_utility.hpp"

using namespace std;

#define DEV
#define PAD cout << "requesting image index out of range. padding processing stack with black images" << endl

mcp3d::Frangi2D::Frangi2D(cv::Mat& input, int k, bool ipl, const vector<float> &_sigmas) :
        input_img(input), window_size(k), inplace(ipl), sigmas(_sigmas)
{
    if (k > 3 && k % 2 == 0)
    {
        cout << "warning: even window size provided. k + 1 will be used instead." << endl;
        ++window_size;
    }
    dynamic_window_size = (k < 3);
    numScales = sigmas.size();
    rows = input_img.rows;
    cols = input_img.cols;
}

vector<float> mcp3d::Frangi2D::get_sigmas()
{
    vector<float> s = sigmas;
    return s;
}

pair<int, int> mcp3d::Frangi2D::get_img_dimension()
{
    return make_pair(rows, cols);
}

int mcp3d::Frangi2D::get_window_size()
{
    return window_size;
}

tuple<cv::Mat, cv::Mat, cv::Mat> mcp3d::Frangi2D::get_hessians()
{
    assert (I_xx.rows > 0 && I_xy.rows > 0 && I_yy.rows > 0);
    return make_tuple(I_xx, I_xy, I_yy);
}

/// Note: gaussian derivatives at different scales are normalized by either sigma (1st order derivative) or
/// sigma^2 (2nd order derivative)
mcp3d::separate_filters_2d mcp3d::Frangi2D::get_level_kernels(int level, int dx, int dy)
{
    double sigma = sigmas[level];
    if (dynamic_window_size)
    {
        if (sigma <= 0.8)
            window_size = DEFAULT_WINDOW_SIZE;
        else
        {
            //int m = (int) ceil((sigma - 0.35) / 0.15);
            //window_size = (m % 2 == 1 ? m : (m + 1));
            window_size = max(2 * (int)sigma * (int)sigma + 1, 3);
            window_size = min(window_size, min(rows, cols));
        }

    }
    cout << "window size = " << window_size << endl;
    // first get 1D gaussian kernel with windows size and sigma
    cv::Mat gk1D = cv::getGaussianKernel(window_size, sigma, CV_64F);
    cv::Mat ox = gk1D.clone();
    cv::Mat oy = gk1D.clone();
    // if no derivative in any direction, return the gaussian row and column filters
    if (dx + dy == 0)
        return make_pair(gk1D, gk1D);
    // derivatives are needed, calculate it based on f'(x) and f''(x) where f(x) = ce^(-x^2 / (2 sigma^2))
    // return separate filters for rows and columns
    assert(dx + dy == 2 && dx >= 0 && dy >= 0);
    if (dx == 2)
    {
        for (int i = 0; i < window_size; ++i)
        {
            double x = (double) (i - (window_size - 1) / 2);
            ox.at<double>(i) = (1 / pow(sigma, 2.0)) * gk1D.at<double>(i) * (x * x / pow(sigma, 2.0) - 1);
        }
        ox = ox * pow(sigma, 2 * GAUSSIAN_DERIV_NORMALIZER);

    }
    else if (dx == 1 && dy == 1)
    {
        for (int i = 0; i < window_size; ++i)
        {
            double x = (double) (i - (window_size - 1) / 2);
            ox.at<double>(i) = -gk1D.at<double>(i) * x / pow(sigma, 2.0);
        }
        ox = ox * pow(sigma, GAUSSIAN_DERIV_NORMALIZER);
        oy = ox;
    }
    else
    {
        for (int i = 0; i < window_size; ++i)
        {
            double y = (double) (i - (window_size - 1) / 2);
            oy.at<double>(i) = pow(sigma, -2.0) * gk1D.at<double>(i) * (y * y / pow(sigma, 2.0) - 1);
        }
        oy = oy * pow(sigma, 2 * GAUSSIAN_DERIV_NORMALIZER);
    }
    return make_pair(ox, oy);
}

cv::Mat mcp3d::Frangi2D::apply_filters(separate_filters_2d filters)
{
    cv::Mat filteredGaussianDeriv(rows, cols, CV_64F);
    cv::sepFilter2D(input_img, filteredGaussianDeriv, CV_64F, filters.first, filters.second);
    return filteredGaussianDeriv;
}

void mcp3d::Frangi2D::fill_level_hessians(int level)
{
    separate_filters_2d filters_gxx = get_level_kernels(level, 2, 0);
    I_xx = apply_filters(filters_gxx);
    separate_filters_2d filters_gxy = get_level_kernels(level, 1, 1);
    I_xy = apply_filters(filters_gxy);
    separate_filters_2d filters_gyy = get_level_kernels(level, 0, 2);
    I_yy = apply_filters(filters_gyy);
}

void mcp3d::Frangi2D::initMinLambda2()
{
    minLambda2 = cv::Mat(rows, cols, CV_64F, cv::Scalar(1.0));
}

void mcp3d::Frangi2D::updatePixelMinLambda2(int i, int j)
{
    /* real solutions are gauranteed: lambda * lambda - (ixx + iyy)lambda + (ixx * iyy - ixy * ixy) = 0
     * use opencv implementation to find eigenvalues
     */
    cv::Vec2d eigens;
    cv::Mat hessian(2, 2, CV_64F);
    hessian.at<double>(0, 0) = I_xx.at<double>(i, j);
    hessian.at<double>(0, 1) = I_xy.at<double>(i, j);
    hessian.at<double>(1, 0) = hessian.at<double>(0, 1);
    hessian.at<double>(1, 1) = I_yy.at<double>(i, j);
    bool b = cv::eigen(hessian, eigens);
    // update only if b is true
    //cout << eigens(0) << endl;
    //cout << eigens(1) << endl;
    if (b)
    {
        if (abs(eigens(0) > ZERO_EIGEN_TOL))
            return;
        if (eigens(1) < 0 && eigens(1) < minLambda2.at<double>(i, j))
            minLambda2.at<double>(i, j) = eigens(1);
    }
}

void mcp3d::Frangi2D::updateImgMinLambda2()
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
            updatePixelMinLambda2(i, j);
    }
}

void mcp3d::Frangi2D::compute_hessian_eigen_through_scales()
{
    initMinLambda2();
    for (int i = 0; i < (int)numScales; ++i)
    {
        fill_level_hessians(i);
        updateImgMinLambda2();
    }
}

void mcp3d::Frangi2D::remapImgMinLambda2(double remap_factor, cv::Mat& dest)
{
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (minLambda2.at<double>(i, j) < LAMBDA2_HIGH && minLambda2.at<double>(i, j) > LAMBDA2_LOW)
                dest.at<double>(i, j) = minLambda2.at<double>(i, j) * remap_factor;
            else if (minLambda2.at<double>(i, j) <= LAMBDA2_LOW)
                dest.at<double>(i, j) = LAMBDA2_LOW * remap_factor;
            else
                dest.at<double>(i, j) = 10;
}

void mcp3d::Frangi2D::remapImgMaxVesselness(double remap_factor, cv::Mat& dest)
{
    dest = cv::min(maxVesselness * remap_factor, 20.0);
}

/// I' = I * exp(max|minLambda2|)
void mcp3d::Frangi2D::minLambda2_filter(double remap_factor)
{
    eigen_filtered = cv::Mat(rows, cols, CV_64F);
    if (! inplace)
    {
        input_img.convertTo(eigen_filtered, CV_64F, 1.0 / 255);
    }
    else
    {
        input_img.convertTo(input_img, CV_64F, 1.0 / 255);
        eigen_filtered = input_img;
    }
    cv::Mat dest;
    #ifdef DEV
    {
        dest = cv::Mat(rows, cols, CV_64F);
    }
    #else
        cv::Mat dest = minLambda2;
    #endif
    remapImgMinLambda2(remap_factor, dest);
    cv::Mat S(rows, cols, CV_64F);
    cv::exp(dest * (-1), S);
    //cv::exp(dest, S);
    cv::Mat ones(rows, cols, CV_64F, cv::Scalar(1.0));
    //S = 1 / (S + ones);

    eigen_filtered = eigen_filtered.mul(S);
    double minVal, maxVal;
    cv::minMaxLoc(eigen_filtered, &minVal, &maxVal);
    cout << "maxval = " << maxVal << endl;
    eigen_filtered.convertTo(eigen_filtered, CV_8U, 255 / maxVal);
}

void mcp3d::Frangi2D::initMaxVesselness()
{
    maxVesselness = cv::Mat(rows, cols, CV_64F, cv::Scalar(-10.0));
}

double mcp3d::Frangi2D::computePixelVesselness(double lambda1, double lambda2, double max_hessian_l1)
{
    double l2 = abs(lambda1) <= abs(lambda2) ? lambda2 : lambda1;
    double vesselness = abs(l2) / sqrt(max_hessian_l1);
    return vesselness;
}

void mcp3d::Frangi2D::updatePixelMaxVesselness(int i, int j, double max_hessian_l1, int level)
{
    cv::Vec2d eigens;
    cv::Mat hessian(2, 2, CV_64F);
    hessian.at<double>(0, 0) = I_xx.at<double>(i, j);
    hessian.at<double>(0, 1) = I_xy.at<double>(i, j);
    hessian.at<double>(1, 0) = hessian.at<double>(0, 1);
    hessian.at<double>(1, 1) = I_yy.at<double>(i, j);
    bool b = cv::eigen(hessian, eigens);
    if (! b)
        return;
    if (((abs(eigens(0)) >= abs(eigens(1))) && eigens(0) >= 0) || ((abs(eigens(1)) >= abs(eigens(0))) && eigens(1) >= 0))
        return;
    double vesselness = computePixelVesselness(eigens(0), eigens(1), max_hessian_l1);
    maxVesselness.at<double>(i, j) = max(maxVesselness.at<double>(i, j), vesselness);
}

double mcp3d::Frangi2D::get_max_hessian_l1_norm()
{
    cv::Mat hessian_l1;
    hessian_l1 = cv::abs(I_xx) + cv::abs(I_xy) * 2 + cv::abs(I_yy);
    double max_hessian_l1, min_hessian_l2;
    cv::minMaxLoc(hessian_l1, &min_hessian_l2, &max_hessian_l1);
    return max_hessian_l1;
}

void mcp3d::Frangi2D::updateImgMaxVesselness(int level)
{
    double max_hessian_l1 = get_max_hessian_l1_norm();
    cout << "hessian l1 norm maimum = " << max_hessian_l1 << endl;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
        {
            updatePixelMaxVesselness(i, j, max_hessian_l1, level);
        }
}

void mcp3d::Frangi2D::compute_vesselness_through_scales()
{
    initMaxVesselness();
    for (int i = 0; i < (int)numScales; ++i)
    {
        fill_level_hessians(i);
        updateImgMaxVesselness(i);
    }
}

void mcp3d::Frangi2D::vesselness_filter(double remap_factor)
{
    vesselness_filtered = cv::Mat(rows, cols, CV_64F);
    if (! inplace)
    {
        input_img.convertTo(vesselness_filtered, CV_64F);
    }
    else
    {
        input_img.convertTo(input_img, CV_64F);
        vesselness_filtered = input_img;
    }
    cv::Mat dest;
    #ifdef DEV
    {
        dest = maxVesselness.clone();
    }
    #else
        cv::Mat dest = maxVesselness;
    #endif
    remapImgMaxVesselness(remap_factor, dest);
    cv::Mat S(rows, cols, CV_64F);
    cv::exp(dest, S);
    //cv::exp(dest, S);
    //cv::Mat ones(rows, cols, CV_64F, cv::Scalar(1.0));
    //S = 1 / (S + ones);

    vesselness_filtered = vesselness_filtered.mul(S);
    double minVal, maxVal;
    cv::minMaxLoc(vesselness_filtered, &minVal, &maxVal);
    cout << "maxval = " << maxVal << endl;
    vesselness_filtered.convertTo(vesselness_filtered, CV_8U, 255 / maxVal);
}

void mcp3d::Frangi3D::build_stack_from_empty()
{
    int start = max(anchor_img_index - half_window_size, 0);
    int end = min(anchor_img_index + half_window_size, (int)img_path_list.size() - 1);
    cv::Mat I;
    for (int i = anchor_img_index - half_window_size; i < anchor_img_index + half_window_size + 1; ++i)
    {
        cout << i << endl;
        if (i < start || i > end)
        {
            PAD;
            I = cv::Mat(rows, cols, CV_32F, cv::Scalar(0.0));
        }
        else
        {
            I = cv::imread(img_path_list[i], 0);
            cout << img_path_list[i] << endl;
            I.convertTo(I, CV_32F);
        }
        img_stack_vec.push_back(I);
    }
    rows = img_stack_vec[0].rows;
    cols = img_stack_vec[0].cols;
    assert(rows > 0 && cols > 0);
    img_stack = cv::Mat(rows, cols, CV_32FC(window_size));
    anchor_img = cv::Mat(rows, cols, CV_32F);
}

void mcp3d::Frangi3D::build_stack()
{
    cout << "building processing stack" << endl;
    if (img_stack.empty())
    {
        build_stack_from_empty();
    }
    else if (prev_anchor_img_index == anchor_img_index)
    {
        cout << "anchor image is same as previous. not rebuilding processing stack." << endl;
    }
    else if (abs(anchor_img_index - prev_anchor_img_index) <= half_window_size)
    {
        cout << "previous anchor image index = " << prev_anchor_img_index << endl;
        cout << "anchor image index = " << anchor_img_index << endl;
        int z_dif = anchor_img_index - prev_anchor_img_index;
        cv::Mat I;
        if (z_dif > 0)
        {
            for (int i = 0; i < z_dif; ++i)
            {
                img_stack_vec.erase(img_stack_vec.begin());
                int next = prev_anchor_img_index + half_window_size + 1 + i;
                cout << "index of next image to add to processing stack " << next << endl;

                if (next >= (int)img_path_list.size())
                {
                    PAD;
                    I = cv::Mat(rows, cols, CV_32F, cv::Scalar(0.0));
                }
                else
                {
                    I = cv::imread(img_path_list[next], 0);
                    I.convertTo(I, CV_32F);
                }
                img_stack_vec.push_back(I);
            }
        }
        else
        {
            for (int i = 0; i < abs(z_dif); ++i)
            {
                img_stack_vec.pop_back();
                int next = anchor_img_index - half_window_size - 1 - i;
                cv::Mat I;
                cout << "index of next image to add to processing stack " << next << endl;

                if (next < 0)
                {
                    PAD;
                    I = cv::Mat(rows, cols, CV_32F, cv::Scalar(0.0));
                }
                else
                {
                    I = cv::imread(img_path_list[next], 0);
                    I.convertTo(I, CV_32F);
                }
                img_stack_vec.insert(img_stack_vec.begin(), I);
            }
        }
    }
    cv::merge(img_stack_vec, img_stack);
    cv::extractChannel(img_stack, anchor_img, half_window_size);
    cout << "done" << endl;
}

cv::Mat mcp3d::Frangi3D::get_gaussian1D_derivatives(int dx, int level)
{
    assert (dx <= 3);
    float sigma = sigmas[level];
    cv::Mat gk = cv::getGaussianKernel(window_size, sigma, CV_32F);
    if (dx == 0)
        return gk;
    cv::Mat ox(gk.rows, gk.cols, CV_32F);
    if (dx == 1)
    {
        for (int i = 0; i < window_size; ++i)
        {
            float x = (float) (i - (window_size - 1) / 2);
            ox.at<float>(i) = -gk.at<float>(i) * x / (sigma * sigma);
        }
        ox = ox * pow(sigma, GAUSSIAN_DERIV_NORMALIZER);
        return ox;
    }
    else if (dx == 2)
    {
        for (int i = 0; i < window_size; ++i)
        {
            float x = (float) (i - (window_size - 1) / 2);
            ox.at<float>(i) = gk.at<float>(i) * (x * x * (float)pow(sigma, -4.0) - (float)pow(sigma, -2.0));
        }
        ox = ox * pow(sigma, 2 * GAUSSIAN_DERIV_NORMALIZER);
        return ox;
    }
    else
    {
        for (int i = 0; i < window_size; ++i)
        {
            float x = (float) (i - (window_size - 1) / 2);
            ox.at<float>(i) = gk.at<float>(i) * ( 3 * x / (float)pow(sigma, 4.0) - (float)pow(x, 3.0) / (float)pow(sigma, 6.0));
        }
        ox = ox * pow(sigma, 3 * GAUSSIAN_DERIV_NORMALIZER);
        return ox;
    }
}

void mcp3d::Frangi3D::fill_hessians(int level)
{
    // when possible, use separate filter 2D to have a 2D convolution
    // when only 1D result needed, for now use filter2D but pass the correct (k, 1) or (1, k) shape kernel
    // gkdn_x are row vectors of shape (1, k), gkdn_y are column vectors of shape (k, 1)

    // get gaussian 1d filters
    cv::Mat gkd0_x, gkd0_y, gkd1_x, gkd1_y, gkd2_x, gkd2_y;
    gkd0_y = get_gaussian1D_derivatives(0, level);
    gkd1_y = get_gaussian1D_derivatives(1, level);
    gkd2_y = get_gaussian1D_derivatives(2, level);
    cv::transpose(gkd0_y, gkd0_x);
    cv::transpose(gkd1_y, gkd1_x);
    cv::transpose(gkd2_y, gkd2_x);

    // prepare intermediate data matrices
    cv::Mat temp(rows, cols, CV_32F);
    cv::Mat temp_col(rows * cols, window_size, CV_32F);
    cv::Mat convolution(rows, cols, CV_32FC(window_size));
    cv::Mat convolution_reshape;
    vector<cv::Mat> convolution_vec, convolution_vec_temp;
    for (int i = 0; i < window_size; ++i)
    {
        convolution_vec.push_back(cv::Mat(rows, cols, CV_32F));
        convolution_vec_temp.push_back(cv::Mat(rows, cols, CV_32F));
    }

    // if uninitialized, initialize the II_ matrices
    I_xx.create(rows, cols, CV_32F);
    I_xy.create(rows, cols, CV_32F);
    I_xz.create(rows, cols, CV_32F);
    I_yy.create(rows, cols, CV_32F);
    I_yz.create(rows, cols, CV_32F);
    I_zz.create(rows, cols, CV_32F);


    // compute hessians with gkd0 in x direction: I_yy, I_yz, I_zz
    for (int i = 0; i < (int)img_stack_vec.size(); ++i)
    {
        // compute x direction gaussian convolution. convlution_vec_temp contain convolution results with gkd0_x
        cv::filter2D(img_stack_vec[i], convolution_vec_temp[i], CV_32F, gkd0_x);
    }
    // compute y direction 2nd order derivative of gaussian convolution and z direction gaussian convolution
    for (int i = 0; i < (int)convolution_vec_temp.size(); ++i)
        cv::filter2D(convolution_vec_temp[i], convolution_vec[i], CV_32F, gkd2_y);
    // construct the convolution stack from convolution vec
    cv::merge(convolution_vec, convolution);
    // reshape the convolution stack to oprate on z direction. in the reshaped matrix, channels are columns
    // therefore the convolution is along x axis in reshaped matrix
    convolution_reshape = convolution.reshape(1, rows * cols);
    cv::filter2D(convolution_reshape, temp_col, CV_32F, gkd0_x);
    convolution = temp_col.reshape(window_size, rows);
    cv::extractChannel(convolution, I_yy, half_window_size);
    // compute y and z direction 1st order derivative of gaussian convolution
    for (int i = 0; i < (int)convolution_vec_temp.size(); ++i)
        cv::filter2D(convolution_vec_temp[i], convolution_vec[i], CV_32F, gkd1_y);
    cv::merge(convolution_vec, convolution);
    convolution_reshape = convolution.reshape(1, rows * cols);
    cv::filter2D(convolution_reshape, temp_col, CV_32F, gkd1_x);
    convolution = temp_col.reshape(window_size, rows);
    cv::extractChannel(convolution, I_yz, half_window_size);
    // compute y direction gaussian convolution and z direction 2nd order derivative of gaussian convolution
    for (int i = 0; i < (int)convolution_vec_temp.size(); ++i)
        cv::filter2D(convolution_vec_temp[i], convolution_vec[i], CV_32F, gkd0_y);
    cv::merge(convolution_vec, convolution);
    convolution_reshape = convolution.reshape(1, rows * cols);
    cv::filter2D(convolution_reshape, temp_col, CV_32F, gkd2_x);
    convolution = temp_col.reshape(window_size, rows);
    cv::extractChannel(convolution, I_zz, half_window_size);


    // compute hessians with gkd0 in y direction: I_xx, I_xz
    for (int i = 0; i < (int)img_stack_vec.size(); ++i)
    {
        // compute y direction gaussian convolution. convlution_vec_temp contain convolution results with gk1d0_y
        cv::filter2D(img_stack_vec[i], convolution_vec_temp[i], CV_32F, gkd0_y);
    }
    // compute x direction 2nd order gaussian derivative convolution and z direction gaussian convolution
    for (int i = 0; i < (int)convolution_vec_temp.size(); ++i)
        cv::filter2D(convolution_vec_temp[i], convolution_vec[i], CV_32F, gkd2_x);
    cv::merge(convolution_vec, convolution);
    convolution_reshape = convolution.reshape(1, rows * cols);
    cv::filter2D(convolution_reshape, temp_col, CV_32F, gkd0_x);
    convolution = temp_col.reshape(window_size, rows);
    cv::extractChannel(convolution, I_xx, half_window_size);
    // compute x direction 1st order gaussian derivative convolution and z direction 1st order gaussian derivative
    // convolution
    for (int i = 0; i < (int)convolution_vec_temp.size(); ++i)
        cv::filter2D(convolution_vec_temp[i], convolution_vec[i], CV_32F, gkd1_x);
    cv::merge(convolution_vec, convolution);
    convolution_reshape = convolution.reshape(1, rows * cols);
    cv::filter2D(convolution_reshape, temp_col, CV_32F, gkd1_x);
    convolution = temp_col.reshape(window_size, rows);
    cv::extractChannel(convolution, I_xz, half_window_size);



    // compute hessians with gkd0 in z direction: I_xy
    for (int i = 0; i < (int)img_stack_vec.size(); ++i)
    {
        // compute x direction 1st order gaussian derivative convolution and y direction 1st order gaussian derivative
        // convolution with sepFilter2D
        cv::sepFilter2D(img_stack_vec[i], convolution_vec_temp[i], CV_32F, gkd1_x, gkd1_y);
    }
    cv::merge(convolution_vec_temp, convolution);
    convolution_reshape = convolution.reshape(1, rows * cols);
    cv::filter2D(convolution_reshape, temp_col, CV_32F, gkd0_x);
    convolution = temp_col.reshape(window_size, rows);
    cv::extractChannel(convolution, I_xy, half_window_size);

}

void mcp3d::Frangi3D::init_max_vesselness()
{
    max_vesselness = cv::Mat(rows, cols, CV_32F, cv::Scalar(-10.0));
}

void mcp3d::Frangi3D::rescale_max_vesselness(double s, cv::Mat& dst)
{
    dst = cv::min(max_vesselness * s, 20.0);
}

void mcp3d::Frangi3D::get_max_hessian_l1()
{
    cv::Mat diag, diag_reshape, non_diag, non_diag_reshape;
    diag.create(rows, cols, CV_32FC3);
    non_diag.create(rows, cols, CV_32FC3);
    cv::Mat col_diag, col_non_diag, L1;
    col_diag.create(rows * cols, 1, CV_32F);
    col_non_diag.create(rows * cols, 1, CV_32F);
    L1.create(rows * cols, 1, CV_32F);

    vector<cv::Mat> I = {I_xx, I_yy, I_zz};
    cv::merge(I, diag);
    diag_reshape = diag.reshape(1, rows * cols);
    cv::reduce(cv::abs(diag_reshape), col_diag, 1, CV_REDUCE_SUM, -1);

    vector<cv::Mat> J = {I_xy, I_xz, I_yz};
    cv::merge(J, non_diag);
    non_diag_reshape = non_diag.reshape(1, rows * cols);
    cv::reduce(cv::abs(non_diag_reshape), col_non_diag, 1, CV_REDUCE_SUM, -1);
    L1 = col_diag + col_non_diag + col_non_diag;

    double minVal, maxVal;
    cv::minMaxIdx(L1, &minVal, &maxVal);
    max_hessian_l1 = maxVal;
    cout << "max_hessian_l1 = " << max_hessian_l1 << endl;
}

void mcp3d::Frangi3D::update_pixel_vesselness(int i, int j, int level)
{
    // construct pixel hessian matrix and compute its eigen values
    cv::Mat H(3, 3, CV_32F);
    cv::Vec3d E;
    H.at<float>(0, 0) = I_xx.at<float>(i, j);
    H.at<float>(0, 1) = I_xy.at<float>(i, j);
    H.at<float>(0, 2) = I_xz.at<float>(i, j);
    H.at<float>(1, 0) = H.at<float>(0, 1);
    H.at<float>(1, 1) = I_yy.at<float>(i, j);
    H.at<float>(1, 2) = I_yz.at<float>(i, j);
    H.at<float>(2, 0) = H.at<float>(0, 2);
    H.at<float>(2, 1) = H.at<float>(1, 2);
    H.at<float>(2, 2) = I_zz.at<float>(i, j);
    //mcp3d::printMat<double>(H);
    bool b = cv::eigen(H, E);
    // compute vesselness based on eigen values
    if (! b)
        return;
    if (E(1) >= 0 || E(2) >= 0)
        return;
    //cout << "has negative lambda2 and lambda3" << endl;
    double blobness = abs(E(0)) * pow(E(1) * E(2), -0.5);
    double linelike = abs(E(1) / E(2));
    double norm = E(1) * E(1) + E(2) * E(2);
    double vesselness =  exp(-blobness * blobness / mcp3d::W_BLOB_3D) *
                        (1 - exp(-norm / mcp3d::W_L2_3D)) *
                        (1 - exp(-linelike * linelike / mcp3d::W_LINE_3D));
    //cout << vesselness << endl;

    vesselness /= max_hessian_l1;
    if (vesselness > max_vesselness.at<float>(i, j))
        max_vesselness.at<float>(i, j) = (float)vesselness;
}

void mcp3d::Frangi3D::update_img_vesselness(int level)
{
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            update_pixel_vesselness(i, j, level);
}

void mcp3d::Frangi3D::write_filtered_img(string dir)
{
    string out_dir, img_base;
    pair<string, string> splitted = mcp3d::SplitBaseName(anchor_img_path);
    if (dir.empty())
        out_dir = mcp3d::JoinPath({splitted.first, "filter"});
    else
        out_dir = dir;
    mcp3d::MakeDirectories(out_dir);
    img_base = mcp3d::RemoveFileExt(splitted.second) + "_filtered.tif";
    cout << "write to " << mcp3d::JoinPath({out_dir, img_base}) << endl;
    cv::imwrite(mcp3d::JoinPath({out_dir, img_base}), filtered);
}

////////////// public interface ///////////////////
mcp3d::Frangi3D::Frangi3D(vector<string> l, int k, const vector<float> &_sigmas): window_size(k), sigmas(_sigmas)
{
    num_scales = sigmas.size();
    for (const auto& path: l)
    {
        if (!mcp3d::IsFile(path))
        {
            cout << "WARNING:" << path << " is not a valid image path. skipping." << endl;
            bad_files.insert(path);
        }
    }
    for (auto it = l.begin(), end = l.end(); it != end; )
    {
        if (bad_files.find(*it) != bad_files.end())
        {
            it = l.erase(it);
            end = l.end();
        }
        else
            ++it;
    }
    img_path_list = l;
    if (k < 3)
    {
        cout << "window size must be odd and at least 3. setting window size to 3." << endl;
        window_size = 3;
    }
    else if (k % 2 == 0)
    {
        cout << "window size must be odd. setting window size to " + to_string(k + 1) << endl;
        ++window_size;
    }
    assert((int)img_path_list.size() >= window_size);
    half_window_size = (window_size - 1) / 2;
    prev_anchor_img_index = -1;
    cv::Mat I = cv::imread(img_path_list[0], 0);
    rows = I.rows;
    cols = I.cols;
}

bool mcp3d::Frangi3D::init_stack(const string& img_path)
{
    if (img_path != anchor_img_path)
    {
        // update previous and current anchor img index
        prev_anchor_img_index = anchor_img_index;
        for (int i = 0; i < (int)img_path_list.size(); ++i)
        {
            if (img_path_list[i] == img_path)
            {
                anchor_img_path = img_path;
                anchor_img_index = i;
                break;
            }
            if (i == (int)img_path_list.size() - 1)
            {
                cout << img_path << " is not in the list of images provided for the stack. do nothing." << endl;
                return false;
            }
        }
        cout << "anchor image index " << to_string(anchor_img_index) << endl;
        cout << anchor_img_path << endl;
    }
    else
    {
        // do nothing except for updating previous anchor img index
        prev_anchor_img_index = anchor_img_index;
        cout << "anchor image index " << to_string(anchor_img_index) << endl;
    }

    return true;
}

void mcp3d::Frangi3D::compute_vesselness()
{
    build_stack();
    init_max_vesselness();

    for (int i = 0; i < (int)num_scales; ++i)
    {
        fill_hessians(i);
        get_max_hessian_l1();
        update_img_vesselness(i);
    }
}

void mcp3d::Frangi3D::anisotropic_vesselness_filter(double s)
{
    filtered.create(rows, cols, CV_32F);
    //anchor_img.convertTo(filtered, CV_32F);
    cv::Mat dst;
    #ifdef DEV
    {
        dst = max_vesselness.clone();
        cv::Mat S(rows, cols, CV_32F);
        rescale_max_vesselness(s, dst);
        cv::exp(dst, S);
        filtered = anchor_img.mul(S);
    }
    #else
    {
        dst = max_vesselness;
        rescale_max_vesselness(s, dst);
        cv::exp(dst, dst);
        filtered = anchor_img.mul(exp(dst));
    }
    #endif
    double maxVal, minVal;
    cv::minMaxLoc(filtered, &minVal, &maxVal);
    filtered.convertTo(filtered, CV_8U, 255 / maxVal);
}

/// higher level public interface
bool mcp3d::Frangi3D::filter_img(const string& img_path, double s)
{
    bool b = init_stack(img_path);
    if (b)
    {
        compute_vesselness();
        anisotropic_vesselness_filter(s);
        return true;
    }
    else
        return false;
}

bool mcp3d::Frangi3D::filter_imgs(const vector<string> &img_paths, double s)
{
    bool b = true;
    for (const string &img_path: img_paths)
    {
        if (! filter_img(img_path, s))
        {
            b = false;
            cout << "processing of " << img_path << " is not successful" << endl;
        }
        else
            write_filtered_img();
    }
    return b;
}

bool mcp3d::Frangi3D::filter_all_imgs(double s)
{
    bool b = filter_imgs(img_path_list, s);
    return b;
}





















