/**
 * @file   mathfunctions.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup math
 * @ingroup    math
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef MATHFUNCTIONS_HPP
#define MATHFUNCTIONS_HPP

#define _USE_MATH_DEFINES
#include <cstddef>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "assert.hpp"
#include "triple.hpp"


namespace bip
{


/**
 * @var {variable-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
extern const float EPSILON;

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
bool fp_equal(float n1, float n2, float maxerror = EPSILON);
bool NONZEROF(float x);
/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float sqr(float n);

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
int sign(float n);

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float posmod(float n, float m);

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float rad2deg(float n);

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float deg2rad(float n);

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
triple<float> cart2sph(triple<float> n);

triple<float> cart2polar(triple<float> n);

triple<float> normalize(triple<float> n);

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
triple<float> sph2cart(triple<float> n);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
void fast_fourier_transform(fftwf_complex *M, triple<size_t> sizes, bool backward = false);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float hanning(triple<size_t> n, triple<size_t> sizes);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float gaussian(triple<float> n, triple<float> mu, float sigma);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
void eigens(double *M, double *eigenvalues, double *eigenvectors);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float median(float *array, size_t size);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
bool local_max(float *array, triple<size_t> pos, triple<size_t> sizes, size_t radius = 3);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
bool inside_region(triple<size_t> pos, triple<size_t> rstart, triple<size_t> rsizes);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float euclidean_dist(const std::vector<float> &v1, const std::vector<float> &v2);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
float chi_square_dist(const std::vector<float> &v1, const std::vector<float> &v2);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
void normalize_minmax(float *array, size_t size, float new_min, float new_max);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
void normalize_vnorm(std::vector<float> &v);

void normalize_descriptor(std::vector<float> &v);

/**
 * @fn {function-declaration}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @param[in]  {param-name} {param-description}
 * @param[out] {param-name} {param-description}
 *
 * @returns {return-value-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
triple<size_t> log_spherical_bin(triple<float> sph, triple<size_t> num_bins);

int spider_web_bin(triple<float> cart, float ri, float r, triple<size_t> nbins);

}

#endif
