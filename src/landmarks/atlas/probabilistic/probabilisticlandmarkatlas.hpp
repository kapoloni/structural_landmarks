/**
 * @file   probabilisticlandmarkatlas.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-probabilistic
 * @ingroup    landmark-atlas-probabilistic
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef PROBABILISTICLANDMARKATLAS_HPP
#define PROBABILISTICLANDMARKATLAS_HPP

#define BIP_VERBOSE_MODE
// #define BIP_DEBUG_MODE

#include <cstddef>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <itkImage.h>
#include <itkPoint.h>
#include <itkPointSet.h>
#include <itkManifoldParzenWindowsPointSetFunction.h>
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "mathfunctions.hpp"
#include "imageio.hpp"
#include "landmarkio.hpp"


namespace bip
{


/**
 * @class {class-name} {header-file}
 *
 * @brief {brief-description}
 * {detailed-description}
 * 
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
class probabilistic_landmark_atlas
{
    typedef itk::Image<float, 3>    itkImage;
    typedef itk::PointSet<float, 3> itkPointSet;

public:
    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @param[in]  {param-name} {param-description}
     * @param[out] {param-name} {param-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    probabilistic_landmark_atlas(std::string       filename_samples_prefix,
                                 triple<size_t>    sizes,
                                 size_t            num_samples,
                        typename itkImage::Pointer reference_img,
                                 float             kernel_sigma            = 1.0,
                                 float             regularization_sigma    = 1.0,
                                 size_t            evaluation_k            = 50,
                                 size_t            covariance_k            = 5,
                                 float             density_threshold       = 0.4,
                                 size_t            averaging_window_radius = 5);

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    ~probabilistic_landmark_atlas();

    /**
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
    void compute();

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    std::ostream& print(std::ostream &os) const;

public:
    /** @brief Getter for {attribute-name} */
    std::string get_filename_samples_prefix() const {
        return m_filename_samples_prefix;
    }

    /** @brief Getter for {attribute-name} */
    triple<size_t> get_sizes() const {
        return m_sizes;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_num_samples() const {
        return m_num_samples;
    }

    /** @brief Getter for {attribute-name} */
    itkImage::Pointer get_reference_img() const {
        return m_reference_img;
    }

    /** @brief Getter for {attribute-name} */
    float get_kernel_sigma() const {
        return m_kernel_sigma;
    }

    /** @brief Getter for {attribute-name} */
    float get_regularization_sigma() const {
        return m_regularization_sigma;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_evaluation_k() const {
        return m_evaluation_k;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_covariance_k() const {
        return m_covariance_k;
    }

    /** @brief Getter for {attribute-name} */
    float get_density_threshold() const {
        return m_density_threshold;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_averaging_window_radius() const {
        return m_averaging_window_radius;
    }

    /** @brief Getter for {attribute-name} */
    const std::vector<landmark>& get_average_landmarks() const {
        return m_average_landmarks;
    }

    /** @brief Getter for {attribute-name} */
    typename itkImage::Pointer get_density() const {
        return m_density;
    }

public:
    /** @brief Setter for {attribute-name} */
    void set_filename_samples_prefix(std::string filename_samples_prefix) {
        debug::assert2(!filename_samples_prefix.empty());
        m_filename_samples_prefix = filename_samples_prefix;
    }

    /** @brief Getter for {attribute-name} */
    void set_sizes(triple<size_t> sizes) {
        m_sizes = sizes;
    }

    /** @brief Getter for {attribute-name} */
    void set_num_samples(size_t num_samples) {
        m_num_samples = num_samples;
    }

    /** @brief Setter for {attribute-name} */
    void set_reference_img(typename itkImage::Pointer reference_img) {
        m_reference_img = reference_img;
    }

    /** @brief Setter for {attribute-name} */
    void set_kernel_sigma(float kernel_sigma) {
        debug::assert2(kernel_sigma > 0.0);
        m_kernel_sigma = kernel_sigma;
    }

    /** @brief Setter for {attribute-name} */
    void set_regularization_sigma(float regularization_sigma) {
        debug::assert2(regularization_sigma > 0.0);
        m_regularization_sigma = regularization_sigma;
    }

    /** @brief Setter for {attribute-name} */
    void set_evaluation_k(size_t evaluation_k) {
        debug::assert2(evaluation_k > 0);
        m_evaluation_k = evaluation_k;
    }

    /** @brief Setter for {attribute-name} */
    void set_covariance_k(size_t covariance_k) {
        debug::assert2(covariance_k > 0);
        m_covariance_k = covariance_k;
    }

    /** @brief Setter for {attribute-name} */
    void set_density_threshold(float density_threshold) {
        debug::assert2(density_threshold >= 0.0 && density_threshold <= 1.0);
        m_density_threshold = density_threshold;
    }

    /** @brief Setter for {attribute-name} */
    void set_averaging_window_radius(size_t averaging_window_radius) {
        debug::assert2(averaging_window_radius > 0);
        m_averaging_window_radius = averaging_window_radius;
    }

private:
    /**
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
    void compute_density();

    /**
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
    void compute_average_landmarks();

    /**
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
    std::vector<triple<size_t> > find_density_local_maxima();

private:
    /** @attention Purposely not implemented. */
    probabilistic_landmark_atlas(const probabilistic_landmark_atlas &other);

    /** @attention Purposely not implemented. */
    probabilistic_landmark_atlas& operator=(const probabilistic_landmark_atlas &rhs);

private:
    std::string           m_filename_samples_prefix; /**< {brief-description} */
    triple<size_t>        m_sizes;                   /**< {brief-description} */
    size_t                m_num_samples;             /**< {brief-description} */
    itkImage::Pointer     m_reference_img;           /**< {brief-description} */
    float                 m_kernel_sigma;            /**< {brief-description} */
    float                 m_regularization_sigma;    /**< {brief-description} */
    size_t                m_evaluation_k;            /**< {brief-description} */
    size_t                m_covariance_k;            /**< {brief-description} */
    float                 m_density_threshold;       /**< {brief-description} */
    size_t                m_averaging_window_radius; /**< {brief-description} */
    itkPointSet::Pointer  m_point_set;               /**< {brief-description} */
    itkImage::Pointer     m_density;                 /**< {brief-description} */
    std::vector<landmark> m_average_landmarks;       /**< {brief-description} */
};


// ------------------------------------------------------------------------------------------------

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
std::ostream& operator<<(std::ostream &os, const probabilistic_landmark_atlas &pla);


}

#endif
