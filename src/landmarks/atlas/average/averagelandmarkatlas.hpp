/**
 * @file   averagelandmarkatlas.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-average
 * @ingroup    landmark-atlas-average
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef AVERAGELANDMARKATLAS_HPP
#define AVERAGELANDMARKATLAS_HPP

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
#include <deque>
#include <itkImage.h>
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "mathfunctions.hpp"
#include "imageio.hpp"
#include "landmarkio.hpp"
#include "matchlandmarks.hpp"
#include "atlassample.hpp"


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
class average_landmark_atlas
{
    typedef itk::Image<float, 3> itkImage;

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
    average_landmark_atlas(std::string       filename_samples_prefix,
                           triple<size_t>    sizes,
                           size_t            num_samples,
                  typename itkImage::Pointer reference_img,
                           float             rel_frequency_threshold = 0.6,
                           size_t            max_countdown_value     = 50,
                           float             descriptors_tradeoff    = 0.5,
                           float             max_descriptor_dist     = FLT_MAX,
                           float             max_location_dist       = FLT_MAX);

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    ~average_landmark_atlas();

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
    float get_rel_frequency_threshold() const {
        return m_rel_frequency_threshold;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_max_countdown_value() const {
        return m_max_countdown_value;
    }

    /** @brief Getter for {attribute-name} */
    float get_descriptors_tradeoff() const {
        return m_descriptors_tradeoff;
    }

    /** @brief Getter for {attribute-name} */
    float get_max_descriptor_dist() const {
        return m_max_descriptor_dist;
    }

    /** @brief Getter for {attribute-name} */
    float get_max_location_dist() const {
        return m_max_location_dist;
    }

    /** @brief Getter for {attribute-name} */
    const atlas_sample* get_average_sample() const {
        return m_average_sample;
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
    void set_rel_frequency_threshold(float rel_frequency_threshold) {
        debug::assert2(rel_frequency_threshold >= 0.0 && rel_frequency_threshold <= 1.0);
        m_rel_frequency_threshold = rel_frequency_threshold;
    }

    /** @brief Setter for {attribute-name} */
    void set_max_countdown_value(size_t max_countdown_value) {
        debug::assert2(max_countdown_value > 0);
        m_max_countdown_value = max_countdown_value;
    }

    /** @brief Setter for {attribute-name} */
    void set_descriptors_tradeoff(float descriptors_tradeoff) {
        debug::assert2(descriptors_tradeoff >= 0.0 && descriptors_tradeoff <= 1.0);
        m_descriptors_tradeoff = descriptors_tradeoff;
    }

    /** @brief Setter for {attribute-name} */
    void set_max_descriptor_dist(float max_descriptor_dist) {
        m_max_descriptor_dist = (max_descriptor_dist > 0.0) ? max_descriptor_dist : FLT_MAX;
    }

    /** @brief Setter for {attribute-name} */
    void set_max_location_dist(float max_location_dist) {
        m_max_location_dist = (max_location_dist > 0.0) ? max_location_dist : FLT_MAX;
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
    void process_sample(const atlas_sample &sample);

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
    void discard_outliers();

private:
    /** @attention Purposely not implemented. */
    average_landmark_atlas(const average_landmark_atlas &other);

    /** @attention Purposely not implemented. */
    average_landmark_atlas& operator=(const average_landmark_atlas &rhs);

private:
    std::string       m_filename_samples_prefix; /**< {brief-description} */
    triple<size_t>    m_sizes;                   /**< {brief-description} */
    size_t            m_num_samples;             /**< {brief-description} */
    itkImage::Pointer m_reference_img;           /**< {brief-description} */
    float             m_rel_frequency_threshold; /**< {brief-description} */
    size_t            m_max_countdown_value;     /**< {brief-description} */
    float             m_descriptors_tradeoff;    /**< {brief-description} */
    float             m_max_descriptor_dist;     /**< {brief-description} */
    float             m_max_location_dist;       /**< {brief-description} */
    atlas_sample     *m_average_sample;          /**< {brief-description} */
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
std::ostream& operator<<(std::ostream &os, const average_landmark_atlas &ala);


}

#endif
