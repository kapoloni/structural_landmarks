/**
 * @file   landmarkguidedregistration.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-registration
 * @ingroup    landmark-registration
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef LANDMARKGUIDEDREGISTRATION_HPP
#define LANDMARKGUIDEDREGISTRATION_HPP

#define BIP_VERBOSE_MODE
// #define BIP_DEBUG_MODE

#include <cstddef>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <itkImage.h>
#include <itkVector.h>
#include <itkMesh.h>
#include <itkBSplineTransform.h>
#include <itkLandmarkBasedTransformInitializer.h>
#include <itkImageRegionIterator.h>
#include <itkIterativeInverseDisplacementFieldImageFilter.h>
#include <itkWarpImageFilter.h>
#include <itkWarpMeshFilter.h>
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "imageio.hpp"
#include "meshio.hpp"
#include "matchlandmarks.hpp"


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
class landmark_guided_registration
{
    typedef itk::Image<float, 3>                                    itkImage;
    typedef itk::Vector<float, 3>                                   itkVector;
    typedef itk::Image<itkVector, 3>                                itkDisplField;
    typedef itk::Mesh<float, 3>                                     itkMesh;
    typedef std::vector<itkMesh::Pointer>                           itkMeshContainer;
    typedef itk::BSplineTransform<double, 3, 3>                     itkTransform;
    typedef itk::LandmarkBasedTransformInitializer<
                 itkTransform, itkImage, itkImage>                  itkTransformInitializer;
    typedef itk::IterativeInverseDisplacementFieldImageFilter<
                 itkDisplField, itkDisplField>                      itkInverseDisplFieldFilter;
    typedef itk::ImageRegionIterator<itkDisplField>                 itkDisplFieldIterator;
    typedef itk::WarpImageFilter<itkImage, itkImage, itkDisplField> itkWarpImageFilter;
    typedef itk::WarpMeshFilter<itkMesh, itkMesh, itkDisplField>    itkWarpMeshFilter;

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
    landmark_guided_registration(typename itkImage::Pointer      fixed_img,
                                 typename itkImage::Pointer      moving_img,
                                          itkMeshContainer       moving_meshes,
                                          std::vector<landmark> &fixed_img_landmarks,
                                          std::vector<landmark> &moving_img_landmarks,
                                          size_t                 num_control_points   = 7,
                                          float                  descriptors_tradeoff = 0.5,
                                          float                  max_descriptor_dist  = FLT_MAX,
                                          float                  max_location_dist    = FLT_MAX,
                                 typename itkImage::Pointer      landmark_weights_map = nullptr);

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    ~landmark_guided_registration();

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
    typename itkImage::Pointer get_fixed_img() const {
        return m_fixed_img;
    }

    /** @brief Getter for {attribute-name} */
    typename itkImage::Pointer get_moving_img() const {
        return m_moving_img;
    }

    /** @brief Getter for {attribute-name} */
    const itkMeshContainer& get_moving_meshes() const {
        return m_moving_meshes;
    }

    /** @brief Getter for {attribute-name} */
    const std::vector<landmark>& get_fixed_img_landmarks() const {
        return m_fixed_img_landmarks;
    }

    /** @brief Getter for {attribute-name} */
    const std::vector<landmark>& get_moving_img_landmarks() const {
        return m_moving_img_landmarks;
    }

    /** @brief Getter for {attribute-name} */
    size_t get_num_control_points() const {
        return m_num_control_points;
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
    typename itkImage::Pointer get_landmark_weights_map() const {
        return m_landmark_weights_map;
    }

    /** @brief Getter for {attribute-name} */
    typename itkTransform::Pointer get_transform() const {
        return m_transform;
    }

    /** @brief Getter for {attribute-name} */
    typename itkDisplField::Pointer get_displacement_field() const {
        return m_displacement_field;
    }

    /** @brief Getter for {attribute-name} */
    typename itkDisplField::Pointer get_inv_displacement_field() const {
        return m_inv_displacement_field;
    }

    /** @brief Getter for {attribute-name} */
    typename itkImage::Pointer get_output_img() const {
        return m_output_img;
    }

    /** @brief Getter for {attribute-name} */
    const itkMeshContainer& get_output_meshes() const {
        return m_output_meshes;
    }

private:
    /** @brief Setter for {attribute-name} */
    void set_fixed_img(typename itkImage::Pointer fixed_img) {
        debug::assert2(fixed_img.IsNotNull());
        m_fixed_img = fixed_img;
    }

    /** @brief Setter for {attribute-name} */
    void set_moving_img(typename itkImage::Pointer moving_img) {
        debug::assert2(moving_img.IsNotNull());
        m_moving_img = moving_img;
    }

    /** @brief Setter for {attribute-name} */
    void set_moving_meshes(const itkMeshContainer &moving_meshes) {
        m_moving_meshes = moving_meshes;
    }

    /** @brief Setter for {attribute-name} */
    void set_fixed_img_landmarks(const std::vector<landmark> &fixed_img_landmarks) {
        m_fixed_img_landmarks = fixed_img_landmarks;
    }

    /** @brief Setter for {attribute-name} */
    void set_moving_img_landmarks(const std::vector<landmark> &moving_img_landmarks) {
        m_moving_img_landmarks = moving_img_landmarks;
    }

    /** @brief Setter for {attribute-name} */
    void set_num_control_points(size_t num_control_points) {
        debug::assert2(num_control_points > 3);
        m_num_control_points = num_control_points;
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

    /** @brief Setter for {attribute-name} */
    void set_landmark_weights_map(typename itkImage::Pointer landmark_weights_map) {
        m_landmark_weights_map = landmark_weights_map;
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
    void initialize_transform();

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
    void generate_displacement_fields();

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
    void transform_image();

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
    void transform_meshes();

private:
    /** @attention Purposely not implemented. */
    landmark_guided_registration(const landmark_guided_registration &other);

    /** @attention Purposely not implemented. */
    landmark_guided_registration& operator=(const landmark_guided_registration &rhs);

private:
    itkImage::Pointer      m_fixed_img;              /**< {brief-description} */
    itkImage::Pointer      m_moving_img;             /**< {brief-description} */
    itkMeshContainer       m_moving_meshes;          /**< {brief-description} */
    std::vector<landmark>  m_fixed_img_landmarks;    /**< {brief-description} */
    std::vector<landmark>  m_moving_img_landmarks;   /**< {brief-description} */
    size_t                 m_num_control_points;     /**< {brief-description} */
    float                  m_descriptors_tradeoff;   /**< {brief-description} */
    float                  m_max_descriptor_dist;    /**< {brief-description} */
    float                  m_max_location_dist;      /**< {brief-description} */
    itkImage::Pointer      m_landmark_weights_map;   /**< {brief-description} */
    itkTransform::Pointer  m_transform;              /**< {brief-description} */
    itkDisplField::Pointer m_displacement_field;     /**< {brief-description} */
    itkDisplField::Pointer m_inv_displacement_field; /**< {brief-description} */
    itkImage::Pointer      m_output_img;             /**< {brief-description} */
    itkMeshContainer       m_output_meshes;          /**< {brief-description} */
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
std::ostream& operator<<(std::ostream &os, const landmark_guided_registration &reg);


}

#endif
