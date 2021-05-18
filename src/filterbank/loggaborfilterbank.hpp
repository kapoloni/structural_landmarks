/**
 * @file   loggaborfilterbank.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup filterbank
 * @ingroup    filterbank
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef BANKLOGGABOR_HPP
#define BANKLOGGABOR_HPP

#define BIP_VERBOSE_MODE
// #define BIP_DEBUG_MODE

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "assert.hpp"
#include "triple.hpp"
#include "mathfunctions.hpp"
#include "imageio.hpp"


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
class loggabor_filter_bank
{
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
    loggabor_filter_bank(std::string         filename_prefix,
                         bip::triple<size_t> sizes,
                         size_t              num_scales       = 4,
                         size_t              num_azimuths     = 6,
                         size_t              num_elevations   = 3,
                         float               max_frequency    = 0.3,
                         float               mult_factor      = 2.1,
                         float               frequency_ratio  = 0.55,
                         float               angular_ratio    = 1.2,
                         float               lowpass_order    = 15.0,
                         float               lowpass_cutoff   = 0.45,
                         bool                uniform_sampling = false);

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    ~loggabor_filter_bank();

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
    float* get_filter(size_t scale, size_t azimuth = 0, size_t elevation = 0);

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
    static loggabor_filter_bank* read_parameters(std::string filename);

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
    static void write_parameters(loggabor_filter_bank &bof);

public:
    /** @brief Getter for {member-name}. */
    std::string get_filename_prefix() const {
        return m_filename_prefix;
    }

    /** @brief Getter for {member-name}. */
    bip::triple<size_t> get_sizes() const {
        return m_sizes;
    }

    /** @brief Getter for {member-name}. */
    size_t get_num_scales() const {
        return m_num_scales;
    }

    /** @brief Getter for {member-name}. */
    size_t get_num_azimuths() const {
        return m_num_azimuths;
    }

    /** @brief Getter for {member-name}. */
    size_t get_num_azimuths_per_elevation(size_t elevation) const {
        bip::debug::assert2(elevation < m_num_elevations);
        return m_num_azimuths_per_elevation[elevation];
    }

    /** @brief Getter for {member-name}. */
    size_t get_num_elevations() const {
        return m_num_elevations;
    }

    /** @brief Getter for {member-name}. */
    float get_max_frequency() const {
        return m_max_frequency;
    }

    /** @brief Getter for {member-name}. */
    float get_mult_factor() const {
        return m_mult_factor;
    }

    /** @brief Getter for {member-name}. */
    float get_frequency_ratio() const {
        return m_frequency_ratio;
    }

    /** @brief Getter for {member-name}. */
    float get_angular_ratio() const {
        return m_angular_ratio;
    }

    /** @brief Getter for {member-name}. */
    float get_lowpass_order() const {
        return m_lowpass_order;
    }

    /** @brief Getter for {member-name}. */
    float get_lowpass_cutoff() const {
        return m_lowpass_cutoff;
    }

    /** @brief Getter for {member-name}. */
    bool get_uniform_sampling() const {
        return m_uniform_sampling;
    }

public:
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
    size_t get_num_orientations() const {
        size_t result = 0;
        for (size_t e = 0; e < m_num_elevations; ++e)
            result += m_num_azimuths_per_elevation[e];
        return result;
    }

private:
    /** @brief Setter for {member-name}. */
    void set_filename_prefix(std::string filename_prefix) {
        bip::debug::assert2(!filename_prefix.empty());
        m_filename_prefix = filename_prefix;
    }

    /** @brief Setter for {member-name}. */
    void set_sizes(bip::triple<size_t> sizes) {
        m_sizes = sizes;
    }

    /** @brief Setter for {member-name}. */
    void set_num_scales(size_t num_scales) {
        bip::debug::assert2(num_scales > 0);
        m_num_scales = num_scales;
    }

    /** @brief Setter for {member-name}. */
    void set_num_azimuths(size_t num_azimuths) {
        bip::debug::assert2(num_azimuths > 0);
        m_num_azimuths = num_azimuths;
    }

    /** @brief Setter for {member-name}. */
    void set_num_elevations(size_t num_elevations) {
        bip::debug::assert2(num_elevations > 0);
        m_num_elevations = num_elevations;
    }

    /** @brief Setter for {member-name}. */
    void set_max_frequency(float max_frequency) {
        bip::debug::assert2(max_frequency > 0.0 && max_frequency < 0.5);
        m_max_frequency = max_frequency;
    }

    /** @brief Setter for {member-name}. */
    void set_mult_factor(float mult_factor) {
        bip::debug::assert2(mult_factor > 0.0);
        m_mult_factor = mult_factor;
    }

    /** @brief Setter for {member-name}. */
    void set_frequency_ratio(float frequency_ratio) {
        bip::debug::assert2(frequency_ratio > 0.0);
        m_frequency_ratio = frequency_ratio;
    }

    /** @brief Setter for {member-name}. */
    void set_angular_ratio(float angular_ratio) {
        bip::debug::assert2(angular_ratio > 0.0);
        m_angular_ratio = angular_ratio;
    }

    /** @brief Setter for {member-name}. */
    void set_lowpass_order(float lowpass_order) {
        bip::debug::assert2(lowpass_order > 0.0);
        m_lowpass_order = lowpass_order;
    }

    /** @brief Setter for {member-name}. */
    void set_lowpass_cutoff(float lowpass_cutoff) {
        bip::debug::assert2(lowpass_cutoff > 0.0 && lowpass_cutoff <= 0.5);
        m_lowpass_cutoff = lowpass_cutoff;
    }

    /** @brief Setter for {member-name}. */
    void set_uniform_sampling(bool uniform_sampling) {
        m_uniform_sampling = uniform_sampling;
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
    float* create_filter(float freq0, float phi0, float theta0);

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
    float* read_filter(size_t scale, size_t azimuth = 0, size_t elevation = 0);

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
    void write_filter(float *filter, size_t scale, size_t azimuth = 0, size_t elevation = 0);

private:
    /** @attention Purposely not implemented. */
    loggabor_filter_bank(const loggabor_filter_bank &other);

    /** @attention Purposely not implemented. */
    loggabor_filter_bank& operator=(const loggabor_filter_bank &rhs);

private:
    std::string         m_filename_prefix;            /**< {brief-description} */
    bip::triple<size_t> m_sizes;                      /**< {brief-description} */
    size_t              m_num_scales;                 /**< {brief-description} */
    size_t              m_num_azimuths;               /**< {brief-description} */
    size_t             *m_num_azimuths_per_elevation; /**< {brief-description} */
    size_t              m_num_elevations;             /**< {brief-description} */
    float               m_max_frequency;              /**< {brief-description} */
    float               m_mult_factor;                /**< {brief-description} */
    float               m_frequency_ratio;            /**< {brief-description} */
    float               m_angular_ratio;              /**< {brief-description} */
    float               m_lowpass_order;              /**< {brief-description} */
    float               m_lowpass_cutoff;             /**< {brief-description} */
    bool                m_uniform_sampling;           /**< {brief-description} */
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
std::ostream& operator<<(std::ostream &os, const loggabor_filter_bank &lgab);


}

#endif
