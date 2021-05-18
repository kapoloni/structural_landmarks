/**
 * @file   landmarkdetector.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-detection
 * @ingroup    landmark-detection
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef LANDMARKDETECTOR_HPP
#define LANDMARKDETECTOR_HPP

#define BIP_VERBOSE_MODE
// #define BIP_DEBUG_MODE

#define _USE_MATH_DEFINES
#include <cstddef>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <itkImage.h>
#include <itkVector.h>
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
class landmark_detector
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
    landmark_detector(std::string       filename_prefix,
                      std::string map_name,
                      triple<size_t>    sizes,
             typename itkImage::Pointer reference_img,
                      unsigned char    *input_mask              = NULL,
                      float             saliency_threshold      = 0.5,
                      size_t            num_max_landmarks       = 10000,
                      size_t            local_descr_region_size = 16,
                      size_t            global_descr_max_radius = 64,
                      int                                    nr = 3,
                      int                                    na = 8,
                      unsigned char     options                 = OPTION_CALCULATE_ORIENTATION    +
                                                                  OPTION_CALCULATE_LANDMARK +
                                                                  OPTION_GET_GLOBAL_DESCRIPTORS +
                                                                  OPTION_GET_LOCAL_DESCRIPTORS);

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    ~landmark_detector();

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
    void compute(char *pointSetFilename);

public:
    /** @brief Getter for {attribute-name} */
    std::string get_filename_prefix() const {
        return m_filename_prefix;
    }

    /** @brief Getter for {attribute-name} */
    triple<size_t> get_sizes() const {
        return m_sizes;
    }

    /** @brief Getter for {attribute-name} */
    unsigned char* get_input_mask() const {
        return m_input_mask;
    }

    /** @brief Getter for {attribute-name} */
    std::string get_map_name() const
    {
        return m_map_name;
    }
    /** @brief Getter for {attribute-name} */
    typename itkImage::Pointer get_reference_img() const {
        return m_reference_img;
    }

    /** @brief Getter for {attribute-name} */
    float get_saliency_threshold() const {
        return m_saliency_threshold;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_num_max_landmarks() const {
        return m_num_max_landmarks;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_local_descr_region_size() const {
        return m_local_descr_region_size;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_global_descr_max_radius() const {
        return m_global_descr_max_radius;
    }

    /** @brief Getter for {attribute-name} */
    int get_nr() const {
        return m_nr;
    }

    /** @brief Getter for {attribute-name} */
    int get_na() const {
        return m_na;
    }

    /** @brief Getter for {attribute-name} */
    unsigned char get_options() const {
        return m_options;
    }

    /** @brief Getter for {attribute-name} */
    std::vector<landmark> get_landmarks() const {
        return m_landmarks;
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
    bool are_options_set(unsigned char options) const {
        return m_options & options;
    }

private:
    /** @brief Setter for {attribute-name} */
    void set_filename_prefix(std::string filename_prefix) {
        debug::assert2(!filename_prefix.empty());
        m_filename_prefix = filename_prefix;
    }

    /** @brief Setter for {attribute-name} */
    void set_sizes(triple<size_t> sizes) {
        m_sizes = sizes;
    }

    /** @brief Setter for {attribute-name} */
    void set_input_mask(unsigned char* input_mask) {
        m_input_mask = input_mask;
    }

    /** @brief Setter for {attribute-name} */
    void set_map_name(std::string map_name)
    {
        m_map_name = map_name;
    }
    /** @brief Setter for {attribute-name} */
    void set_reference_img(typename itkImage::Pointer reference_img) {
        debug::assert2(reference_img.IsNotNull());
        m_reference_img = reference_img;
    }

    /** @brief Setter for {attribute-name} */
    void set_saliency_threshold(float saliency_threshold) {
        debug::assert2(saliency_threshold >= 0.0 && saliency_threshold <= 1.0);
        m_saliency_threshold = saliency_threshold;
    }

    /** @brief Setter for {attribute-name} */
    void set_num_max_landmarks(size_t num_max_landmarks) {
        m_num_max_landmarks = num_max_landmarks;
    }

    /** @brief Setter for {attribute-name} */
    void set_local_descr_region_size(size_t local_descr_region_size) {
        m_local_descr_region_size = local_descr_region_size;
    }

    /** @brief Setter for {attribute-name} */
    void set_global_descr_max_radius(size_t global_descr_max_radius) {
        m_global_descr_max_radius = global_descr_max_radius;
    }

    /** @brief Setter for {attribute-name} */
    void set_nr(int nr) {
        m_nr = nr;
    }

    /** @brief Setter for {attribute-name} */
    void set_na(int na) {
        m_na = na;
    }

    /** @brief Setter for {attribute-name} */
    void set_options(unsigned char options) {
        m_options = options;
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
    void compute_landmark_locations();

    void read_landmark_locations(char * pointSetFilename);

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
    void adaptiveNonMaximalSuppression( const std::vector<landmark>& features, float c_robust );

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
    void compute_landmark_orientations();

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
    void compute_landmark_local_descriptors();

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
    void compute_landmark_global_descriptors();

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
    float *read_scalar_map(const char *filename_suffix);
    float *read_scalar_map_path(std::string filename);

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
    itkVector* read_vectorial_map(const char *filename_suffix);

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

private:
    /** @attention Purposely not implemented. */
    landmark_detector(const landmark_detector &other);

    /** @attention Purposely not implemented. */
    landmark_detector& operator=(const landmark_detector &rhs);

public:
    /** {brief-description} **/
    static const unsigned char OPTION_NONE = 0;

     /** {brief-description} **/
    static const unsigned char OPTION_CALCULATE_ORIENTATION = 1;

     /** {brief-description} **/
    static const unsigned char OPTION_CALCULATE_LANDMARK = 2;

     /** {brief-description} **/
    static const unsigned char OPTION_GET_GLOBAL_DESCRIPTORS = 4;

    /** {brief-description} **/
   static const unsigned char OPTION_GET_LOCAL_DESCRIPTORS = 8;

private:
    std::string           m_filename_prefix;         /**< {brief-description} **/
    triple<size_t> m_sizes;                          /**< {brief-description} **/
    std::string m_map_name;                          /**< {brief-description} **/
    unsigned char        *m_input_mask;              /**< {brief-description} **/
    itkImage::Pointer     m_reference_img;           /**< {brief-description} **/
    float                 m_saliency_threshold;      /**< {brief-description} **/
    size_t                m_num_max_landmarks;       /**< {brief-description} **/
    size_t                m_local_descr_region_size; /**< {brief-description} **/
    size_t                m_global_descr_max_radius; /**< {brief-description} **/
    int                                   m_nr;
    int                                    m_na;
    unsigned char         m_options;                 /**< {brief-description} **/
    std::vector<landmark> m_landmarks;               /**< {brief-description} **/
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
std::ostream& operator<<(std::ostream &os, const landmark_detector &ldet);


}

#endif
