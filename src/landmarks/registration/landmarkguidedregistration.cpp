/**
 * @file   landmarkguidedregistration.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-registration
 * @ingroup    landmark-registration
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "landmarkguidedregistration.hpp"


namespace bip
{


landmark_guided_registration::
landmark_guided_registration(typename itkImage::Pointer      fixed_img,
                             typename itkImage::Pointer      moving_img,
                                      itkMeshContainer       moving_meshes,
                                      std::vector<landmark> &fixed_img_landmarks,
                                      std::vector<landmark> &moving_img_landmarks,
                                      size_t                 num_control_points,
                                      float                  descriptors_tradeoff,
                                      float                  max_descriptor_dist,
                                      float                  max_location_dist,
                             typename itkImage::Pointer      landmark_weights_map)
{
    // Initialize attributes.
    set_fixed_img(fixed_img);
    set_moving_img(moving_img);
    set_moving_meshes(moving_meshes);
    set_fixed_img_landmarks(fixed_img_landmarks);
    set_moving_img_landmarks(moving_img_landmarks);
    set_num_control_points(num_control_points);
    set_descriptors_tradeoff(descriptors_tradeoff);
    set_max_descriptor_dist(max_descriptor_dist);
    set_max_location_dist(max_location_dist);
    set_landmark_weights_map(landmark_weights_map);
}


landmark_guided_registration::
~landmark_guided_registration()
{
    // Nothing.
}


std::ostream&
landmark_guided_registration::
print(std::ostream &os) const
{
    os << "{"
       << "}";

    return os;
}


void
landmark_guided_registration::
compute()
{
    // Compute the pair of transforms.
    initialize_transform();
    generate_displacement_fields();

    // Do the registration.
    transform_image();
    if (m_moving_meshes.size() > 0)
        transform_meshes();
}


void
landmark_guided_registration::
initialize_transform()
{
    #ifdef BIP_VERBOSE_MODE
        puts("Computing the registration transform");
    #endif

    // Find matchings between the landmarks of the fixed image and the moving image.
    matches_vector landmark_matches = match_landmarks(
        m_fixed_img_landmarks, m_moving_img_landmarks, m_descriptors_tradeoff,
        m_max_descriptor_dist, m_max_location_dist);

    itkImage::Pointer normalized_landmark_weights_map = nullptr;

    // Normalize the landmark weights map (if any).
    if (m_landmark_weights_map.IsNotNull()) {
        itkImage::SizeType itkSize = m_landmark_weights_map->GetLargestPossibleRegion().GetSize();

        triple<size_t> sizes(itkSize[0], itkSize[1], itkSize[2]);
        size_t total_size = sizes[0] * sizes[1] * sizes[2];

        float *normalized_weights = image2array<float, 3>(m_landmark_weights_map);
        normalize_minmax(normalized_weights, total_size, 0.0, 1.0);

        normalized_landmark_weights_map = array2image<float, 3>(
            normalized_weights, sizes, m_landmark_weights_map);

        delete[] normalized_weights;
    }

    typename itkImage::IndexType                             fixed_index, moving_index;
    typename itkTransformInitializer::LandmarkPointType      fixed_point, moving_point;
    typename itkTransformInitializer::LandmarkPointContainer fixed_matched_landmark_points,
                                                             moving_matched_landmark_points;
    typename itkTransformInitializer::LandmarkWeightType     landmark_weights;

    // Add each pair of matched landmarks to the point containers.
    // Notice that we first need to convert the image coordinates (indices) to world coordinates
    // (physical points).
    // We also associate each landmark to a weight, which can be taken from the normalized
    // weights map.
    for (size_t i = 0; i < landmark_matches.size(); ++i) {
        size_t m1 = landmark_matches[i].first;
        size_t m2 = landmark_matches[i].second;

        for (size_t d = 0; d < 3; ++d) {
            fixed_index[d]  = m_fixed_img_landmarks[m1].get_location()[d];
            moving_index[d] = m_moving_img_landmarks[m2].get_location()[d];
        }

        if (normalized_landmark_weights_map.IsNotNull())
            landmark_weights.push_back(normalized_landmark_weights_map->GetPixel(moving_index));
        else
            landmark_weights.push_back(1.0);

        m_fixed_img->TransformIndexToPhysicalPoint(fixed_index, fixed_point);
        fixed_matched_landmark_points.push_back(fixed_point);

        m_moving_img->TransformIndexToPhysicalPoint(moving_index, moving_point);
        moving_matched_landmark_points.push_back(moving_point);

        #ifdef BIP_VERBOSE_MODE
            printf("  %lu. Transform point: (%g, %g, %g) <== (%g, %g, %g)\n",
                   i + 1,
                   fixed_point[0],  fixed_point[1],  fixed_point[2],
                   moving_point[0], moving_point[1], moving_point[2]);
        #endif
    }

    // Create the transform.
    m_transform = itkTransform::New();
    m_transform->SetIdentity();

    // Create the transform initializer and use it.
    typename itkTransformInitializer::Pointer initializer = itkTransformInitializer::New();
    initializer->SetFixedLandmarks(fixed_matched_landmark_points);
    initializer->SetMovingLandmarks(moving_matched_landmark_points);
    initializer->SetLandmarkWeight(landmark_weights);
    initializer->SetTransform(m_transform);
    initializer->SetBSplineNumberOfControlPoints(m_num_control_points);
    initializer->SetReferenceImage(m_fixed_img);

    initializer->InitializeTransform();
}


void
landmark_guided_registration::
generate_displacement_fields()
{
    // Create and configure the displacement field.
    m_displacement_field = itkDisplField::New();
    m_displacement_field->SetRegions(m_fixed_img->GetLargestPossibleRegion());
    m_displacement_field->SetOrigin(m_fixed_img->GetOrigin());
    m_displacement_field->SetSpacing(m_fixed_img->GetSpacing());
    m_displacement_field->SetDirection(m_fixed_img->GetDirection());
    m_displacement_field->Allocate();

    typename itkTransform::InputPointType  pin;
    typename itkTransform::OutputPointType pout;
    typename itkDisplField::IndexType      index;
    typename itkDisplField::PixelType      displacement;

    // Iterate through the whole displacement field space.
    itkDisplFieldIterator it(m_displacement_field,
                             m_displacement_field->GetLargestPossibleRegion());
    it.GoToBegin();
    while (!it.IsAtEnd())
    {
        index = it.GetIndex();

        // Set the difference between the transformed point and the original point as the
        // displacement vector for the current pixel.
        m_displacement_field->TransformIndexToPhysicalPoint(index, pin);
        pout = m_transform->TransformPoint(pin);
        displacement = pout - pin;
        it.Set(displacement);

        ++it;
    }

    // Now compute the inverse displacement field.
    itkInverseDisplFieldFilter::Pointer inverse_filter = itkInverseDisplFieldFilter::New();
    // inverse_filter->SetStopValue(...);
    // inverse_filter->SetNumberOfIterations(...);
    inverse_filter->SetInput(m_displacement_field);
    inverse_filter->Update();

    m_inv_displacement_field = inverse_filter->GetOutput();
}


void
landmark_guided_registration::
transform_image()
{
    debug::assert2(m_displacement_field.IsNotNull());

    #ifdef BIP_VERBOSE_MODE
        printf("Nonrigidly registering the moving image");
    #endif

    // The image is transformed using the inverse displacement field.
    typename itkWarpImageFilter::Pointer warp_filter = itkWarpImageFilter::New();
    warp_filter->SetInput(m_moving_img);
    warp_filter->SetDisplacementField(m_displacement_field);
    warp_filter->SetOutputStartIndex(m_fixed_img->GetLargestPossibleRegion().GetIndex());
    warp_filter->SetOutputOrigin(m_fixed_img->GetOrigin());
    warp_filter->SetOutputSpacing(m_fixed_img->GetSpacing());
    warp_filter->SetOutputDirection(m_fixed_img->GetDirection());
    warp_filter->Update();

    m_output_img = warp_filter->GetOutput();

    #ifdef BIP_VERBOSE_MODE
        puts(" - done");
    #endif
}


void
landmark_guided_registration::
transform_meshes()
{
    debug::assert2(m_inv_displacement_field.IsNotNull());

    #ifdef BIP_VERBOSE_MODE
        printf("Nonrigidly registering the moving meshes");
    #endif

    m_output_meshes.clear();

    // The meshes are transformed using the displacement field.
    for (size_t i = 0; i < m_moving_meshes.size(); ++i)
    {
        typename itkWarpMeshFilter::Pointer warp_filter = itkWarpMeshFilter::New();
        warp_filter->SetInput(m_moving_meshes[i]);
        warp_filter->SetDisplacementField(m_inv_displacement_field);
        warp_filter->Update();

        m_output_meshes.push_back(warp_filter->GetOutput());
    }

    #ifdef BIP_VERBOSE_MODE
        puts(" - done");
    #endif
}


// ------------------------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const landmark_guided_registration &reg)
{
    return reg.print(os);
}


}
