/**
 * @file   phasecongruency.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup phasecongruency
 * @ingroup    phasecongruency
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef PHASECONGRUENCY_HPP
#define PHASECONGRUENCY_HPP

#define BIP_VERBOSE_MODE
//#define CHECK_MAP
//#define DIRECTIONAL
//#define OTHERS
//#define BIP_DEBUG_MODE

#define DEFAULT_FFTW_NUMBER_OF_THREADS 8
#define _USE_MATH_DEFINES
#include <cstddef>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <new>
#include <fftw3.h>
#include <itkImage.h>
#include <itkVector.h>
#include "assert.hpp"
#include "triple.hpp"
#include "mathfunctions.hpp"
#include "imageio.hpp"
#include "loggaborfilterbank.hpp"


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
class phase_congruency
{
    typedef itk::Image<float, 3>  itkImage;
    typedef itk::Vector<float, 3> itkVector;

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
    phase_congruency(std::string           filename_prefix,
                     triple<size_t>        sizes,
                     loggabor_filter_bank *filter_bank,
                     float                *input_img,
                     unsigned char        *input_mask      = nullptr,
            typename itkImage::Pointer     reference_img   = nullptr,
                     float                 noise_threshold = -1.0,
                     float                 noise_std       = 3.0,
                     float                 sigmoid_gain    = 10.0,
                     float                 sigmoid_cutoff  = 0.5);

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    ~phase_congruency();

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

public:
    /** @brief Getter for {attribute-name} */
    std::string get_filename_prefix() const {
        return m_filename_prefix;
    }

    /** @brief Getter for {attribute-name} */
    loggabor_filter_bank* get_filter_bank() const {
        return m_filter_bank;
    }

    /** @brief Getter for {attribute-name} */
    float* get_input_img() const {
        return m_input_img;
    }

    /** @brief Getter for {attribute-name} */
    unsigned char* get_input_mask() const {
        return m_input_mask;
    }

    /** @brief Getter for {attribute-name} */
    typename itkImage::Pointer get_reference_img() const {
        return m_reference_img;
    }

    /** @brief Getter for {attribute-name} */
    triple<size_t> get_sizes() const {
        return m_sizes;
    }

    /** @brief Getter for {attribute-name} */
    float get_noise_threshold() const {
        return m_noise_threshold;
    }

    /** @brief Getter for {attribute-name} */
    float get_noise_std() const {
        return m_noise_std;
    }

    /** @brief Getter for {attribute-name} */
    float get_sigmoid_gain() const {
        return m_sigmoid_gain;
    }

    /** @brief Getter for {attribute-name} */
    float get_sigmoid_cutoff() const {
        return m_sigmoid_cutoff;
    }

private:
    /** @brief Setter for {attribute-name} */
    void set_filename_prefix(std::string filename_prefix) {
        debug::assert2(!filename_prefix.empty());
        m_filename_prefix = filename_prefix;
    }

    /** @brief Setter for {attribute-name} */
    void set_filter_bank(loggabor_filter_bank* filter_bank) {
        debug::assert2(filter_bank != nullptr);
        m_filter_bank = filter_bank;
    }

    /** @brief Setter for {attribute-name} */
    void set_input_img(float* input_img) {
        debug::assert2(input_img != nullptr);
        m_input_img = input_img;
    }

    /** @brief Setter for {attribute-name} */
    void set_input_mask(unsigned char* input_mask) {
        m_input_mask = input_mask;
    }

    /** @brief Setter for {attribute-name} */
    void set_reference_img(typename itkImage::Pointer reference_img) {
        debug::assert2(reference_img.IsNotNull());
        m_reference_img = reference_img;
    }

    /** @brief Setter for {attribute-name} */
    void set_sizes(triple<size_t> sizes) {
        m_sizes = sizes;
    }

    /** @brief Setter for {attribute-name} */
    void set_noise_threshold(float noise_threshold) {
        debug::assert2(noise_threshold > 0.0 || fp_equal(noise_threshold, -1.0));
        m_noise_threshold = noise_threshold;
    }

    /** @brief Setter for {attribute-name} */
    void set_noise_std(float noise_std) {
        debug::assert2(noise_std > 0.0);
        m_noise_std = noise_std;
    }

    /** @brief Setter for {attribute-name} */
    void set_sigmoid_gain(float sigmoid_gain) {
        m_sigmoid_gain = sigmoid_gain;
    }

    /** @brief Setter for {attribute-name} */
    void set_sigmoid_cutoff(float sigmoid_cutoff) {
        debug::assert2(sigmoid_cutoff >= 0.0 && sigmoid_cutoff <= 1.0);
        m_sigmoid_cutoff = sigmoid_cutoff;
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
    void write_scalar_map(const char *filename_suffix, float *map);

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
    void write_vectorial_map(const char *filename_suffix, itkVector *map);

private:
    /** @attention Purposely not implemented. */
    phase_congruency(const phase_congruency &other);

    /** @attention Purposely not implemented. */
    phase_congruency& operator=(const phase_congruency &rhs);

private:
    std::string           m_filename_prefix; /**< {brief-description} */
    loggabor_filter_bank *m_filter_bank;     /**< {brief-description} */
    triple<size_t>        m_sizes;           /**< {brief-description} */
    float                *m_input_img;       /**< {brief-description} */
    unsigned char        *m_input_mask;      /**< {brief-description} */
    itkImage::Pointer     m_reference_img;   /**< {brief-description} */
    float                 m_noise_threshold; /**< {brief-description} */
    float                 m_noise_std;       /**< {brief-description} */
    float                 m_sigmoid_gain;    /**< {brief-description} */
    float                 m_sigmoid_cutoff;  /**< {brief-description} */
};


// ------------------------------------------------------------------------------------------------

/**
 * @brief {brief-description}
 * {detailed-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
std::ostream& operator<<(std::ostream &os, const phase_congruency &pc);


}

#endif
