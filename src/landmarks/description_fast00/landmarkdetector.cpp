/**
 * @file   landmarkdetector.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-detection
 * @ingroup    landmark-detection
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "landmarkdetector.hpp"


namespace bip
{


landmark_detector::
landmark_detector(std::string       filename_prefix,
                  std::string       map_name,
                  triple<size_t>    sizes,
         typename itkImage::Pointer reference_img,
                  unsigned char    *input_mask,
                  float             saliency_threshold,
                  size_t            num_max_landmarks,
                  size_t            local_descr_region_size,
                  size_t            global_descr_max_radius,
                  unsigned char     options)
{
    // Initialization of members.
    set_filename_prefix(filename_prefix);
    set_map_name(map_name);
    set_input_mask(input_mask);
    set_sizes(sizes);
    set_reference_img(reference_img);
    set_saliency_threshold(saliency_threshold);
    set_num_max_landmarks(num_max_landmarks);
    set_local_descr_region_size(local_descr_region_size);
    set_global_descr_max_radius(global_descr_max_radius);
    set_options(options);
}


landmark_detector::
~landmark_detector()
{
    // Nothing.
}


std::ostream&
landmark_detector::
print(std::ostream &os) const
{
    os << "{"
       << "m_filename_prefix: "         << m_filename_prefix         << ", "
       << "m_map_name: "                << m_map_name                << ", "
       << "m_sizes: "                   << m_sizes                   << ", "
       << "m_input_mask: "              << m_input_mask              << ", "
       << "m_reference_img: "           << m_reference_img           << ", "
       << "m_saliency_threshold: "      << m_saliency_threshold      << ", "
       << "m_num_max_landmarks: "       << m_num_max_landmarks       << ", "
       << "m_local_descr_region_size: " << m_local_descr_region_size << ", "
       << "m_global_descr_max_radius: " << m_global_descr_max_radius << ", "
       << "}";

    return os;
}


void
landmark_detector::
compute()
{
    compute_landmark_locations();
    // compute_landmark_orientations();

    // if (are_options_set(OPTION_GET_LOCAL_DESCRIPTORS))
    //     compute_landmark_local_descriptors();

    if (are_options_set(OPTION_GET_GLOBAL_DESCRIPTORS))
        compute_landmark_global_descriptors();
}

void
landmark_detector::
compute_landmark_locations()
{
    /*
     * TODO:
     * Parallelize this method by dividing the computations below in multiple threads. This can be
     * done by dividing the saliencies matrix in regions and by assigning each region to a thread.
     * Obvioulsly, we must not make such regions smaller than the search neighborhood size. And
     * please notice that parallelizing the pixel access here would probably not be a good idea.
     *
     * PS: parallelizing this method might be more challenging than the others. Actually, maybe
     * it's not even worth it.
     */

    #ifdef BIP_VERBOSE_MODE
        std::cout << "Computing landmark locations:\n";
    #endif

    // Compute the saliencies map.
    // Actually, in this case the saliencies map is just the lower eigenvalues map. So we just
    // read it (and normalize it).

    const size_t total_size = m_sizes[0] * m_sizes[1] * m_sizes[2];
    const size_t num_dims   = (m_sizes[2] == 1) ? 2 : 3;

    char filename_suffix[512];
    sprintf(filename_suffix, "eigenvalues_%lu", num_dims - 1);

    float *saliencies_map = read_scalar_map(filename_suffix);
    normalize_minmax(saliencies_map, total_size, 0.0, 1.0);

    // This temporary vector is necessary to hold all detected landmarks.
    // Later we keep only the most prominent ones.
    std::vector<landmark> temp_landmarks;

    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i)
            {
                // Ignore locations outside the region of interest (mask).
                if (m_input_mask && !m_input_mask[i])
                    continue;

                // Ignore locations for which the saliency is below the threshold.
                if (saliencies_map[i] < m_saliency_threshold)
                    continue;

                triple<size_t> location(x, y, z);
                triple<float>  features(saliencies_map[i]);

                // The landmark locations are taken as the locations of saliency local maxima.
                // The only landmark feature that we have from now is the magnitude (rho), which
                // is the saliency value itself.
                if (local_max(saliencies_map, location, m_sizes)) {
                    temp_landmarks.push_back(landmark(location, features));

                    #ifdef BIP_VERBOSE_MODE
                        std::cout << "   [" << temp_landmarks.size()
                                  << "] landmark found at position "
                                  << location << "\n";
                    #endif
                }
            }

    // Write an image of the saliencies map (if required).
    if (are_options_set(OPTION_GET_SALIENCIES_MAP))
        write_scalar_map("saliencies", saliencies_map);

    delete[] saliencies_map;

    // Sort the temp_landmarks array in order to have the most prominent landmarks (the ones with
    // greater saliency values) at the end of the vector.
    std::sort(temp_landmarks.begin(), temp_landmarks.end(), landmark_comp());

    // Discard the excessive landmarks.
    size_t num_landmarks = std::min(temp_landmarks.size(), m_num_max_landmarks);
    for (long i = static_cast<long>(num_landmarks) - 1; i >= 0; --i)
        m_landmarks.push_back(temp_landmarks[i]);
}


void
landmark_detector::
compute_landmark_orientations()
{
    /*
     * TODO:
     * Parallelize this method by dividing the computations below in multiple threads. Each
     * landmark could be processed in a different thread, so the performance could be greatly
     * improved with a multithreaded approach.
     */

    #ifdef BIP_VERBOSE_MODE
        std::cout << "Computing landmark orientations:\n";
    #endif

    // Compute the map of feature magnitudes and orientations of each image element.
    // In this case, we use image gradients to get such information, since the phase congruency
    // orientations/moments alone can't give us all possible directions.
    // The gradients, however, are computed from the greatest eigenvalues map, which is already
    // a map of image edges.

    char filename_suffix[512];

    const size_t total_size = m_sizes[0] * m_sizes[1] * m_sizes[2];
    triple<float> *features_map = new triple<float>[total_size]();

    float *grad_mag_map   = new float[total_size]();
    float *grad_phi_map   = new float[total_size]();
    float *grad_theta_map = new float[total_size]();

    sprintf(filename_suffix, "PC");
    float *src_map = read_scalar_map(filename_suffix);

    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                const size_t x1 = (x == 0) ? 0 : x - 1;
                const size_t x2 = (x == m_sizes[0] - 1) ? m_sizes[0] - 1 : x + 1;
                const size_t y1 = (y == 0) ? 0 : y - 1;
                const size_t y2 = (y == m_sizes[1] - 1) ? m_sizes[1] - 1 : y + 1;
                const size_t z1 = (z == 0) ? 0 : z - 1;
                const size_t z2 = (z == m_sizes[2] - 1) ? m_sizes[2] - 1 : z + 1;

                const size_t ix1 = x1 + m_sizes[0] * (y  + m_sizes[1] * z);
                const size_t ix2 = x2 + m_sizes[0] * (y  + m_sizes[1] * z);
                const size_t iy1 = x  + m_sizes[0] * (y1 + m_sizes[1] * z);
                const size_t iy2 = x  + m_sizes[0] * (y2 + m_sizes[1] * z);
                const size_t iz1 = x  + m_sizes[0] * (y  + m_sizes[1] * z1);
                const size_t iz2 = x  + m_sizes[0] * (y  + m_sizes[1] * z2);

                features_map[i] = cart2sph(triple<float>(src_map[ix2] - src_map[ix1],
                                                         src_map[iy2] - src_map[iy1],
                                                         src_map[iz2] - src_map[iz1]));
                grad_mag_map[i]   = features_map[i][0];
                grad_phi_map[i]   = features_map[i][1];
                grad_theta_map[i] = features_map[i][2];
            }

    write_scalar_map("gradient_mag",   grad_mag_map);
    write_scalar_map("gradient_phi",   grad_phi_map);
    write_scalar_map("gradient_theta", grad_theta_map);

    delete[] grad_mag_map;
    delete[] grad_phi_map;
    delete[] grad_theta_map;

    // Image dimensions.
    const size_t num_dims = (m_sizes[2] == 1) ? 2 : 3;

    // Search radius.
    const size_t radius = 2;
    // const size_t radius = 0; // Use only the orientation at the landmarks's exact location.

    // Gaussian weight parameters.
    const float alpha = 2.1;
    const float sigma = (radius + EPSILON) / alpha;

    // Spherical space sampling parameters. We use 36 bins for azimuth angles and 18 bins for
    // elevation angles. And since we're partitioning the space only in orientations, the number
    // of radial bins is not important and can be set as 1.
    const size_t num_bins_azimuth   = 36;
    const size_t num_bins_elevation = (num_dims == 2) ? 1 : 18;
    const size_t num_bins_total     = num_bins_azimuth * num_bins_elevation;
    const triple<size_t> num_bins(1, num_bins_azimuth, num_bins_elevation);

    // Temporary array to include the replicated landmarks (if there is any).
    std::vector<landmark> new_landmarks;

    // Angular histogram. It accumulates the sampled orientation angles of each pixel around the
    // landmark and use such data to determine the dominant orientation of the landmark. It is
    // stored as a vector of triples because we need to know not only the accumulated values of
    // each bin but also their orientations.
    triple<float> histogram[num_bins_total];

    const float dtheta =   M_PI / num_bins[2];
    const float dphi   = 2*M_PI / num_bins[1];
    for (size_t b = 0, e = 0; e < num_bins[2]; ++e) {
        // const float theta = (num_dims == 2) ? 0.0 : -M_PI_2 + dtheta * e;
        const float theta = (num_dims == 2) ? 0.0 : -M_PI_2 + dtheta * (e + 0.5);

        for (size_t a = 0; a < num_bins[1]; ++a, ++b) {
            // const float phi = -M_PI + dphi * a;
            const float phi = -M_PI + dphi * (a + 0.5);

            histogram[b][1] = phi;
            histogram[b][2] = theta;
        }
    }

    for (size_t i = 0; i < m_landmarks.size(); ++i)
    {
        // Clear the histogram values accumulated in the previous iteration.
        for (size_t b = 0; b < num_bins_total; ++b)
            histogram[b][0] = 0.0;

        // Landmark location.
        const size_t x0 = m_landmarks[i].get_location()[0];
        const size_t y0 = m_landmarks[i].get_location()[1];
        const size_t z0 = m_landmarks[i].get_location()[2];

        const size_t zmin = (z0 < radius) ? 0 : z0 - radius;
        const size_t zmax = (z0 + radius >= m_sizes[2]) ? m_sizes[2]-1 : z0 + radius;
        const size_t ymin = (y0 < radius) ? 0 : y0 - radius;
        const size_t ymax = (y0 + radius >= m_sizes[1]) ? m_sizes[1]-1 : y0 + radius;
        const size_t xmin = (x0 < radius) ? 0 : x0 - radius;
        const size_t xmax = (x0 + radius >= m_sizes[0]) ? m_sizes[0]-1 : x0 + radius;

        for (size_t z = zmin; z <= zmax; ++z)
            for (size_t y = ymin; y <= ymax; ++y)
                for (size_t x = xmin; x <= xmax; ++x)
                {
                    // Get the index of the current neighbor's location.
                    const size_t j = x + m_sizes[0] * (y + m_sizes[1] * z);

                    // Ignore locations outside the region of interest (mask).
                    if (m_input_mask && !m_input_mask[j])
                        continue;

                    // The magnitude (rho) of the current neighbor is weighted by a 2D/3D Gaussian
                    // centered at the landmark location.
                    const triple<float> n(x, y, z);
                    const triple<float> mu(x0, y0, z0);
                    const float magnitude = features_map[j][0];
                    const float weight    = gaussian(n, mu, sigma);

                    // Find the bins of the orientation angles of the current neighbor.
                    // Notice that this function also gives us a radial bin index, but such value
                    // is not used in this context.
                    const triple<size_t> bins = log_spherical_bin(features_map[j], num_bins);

                    // Find the final bin index.
                    // PS: as stated above, bins[0] is not used and num_bins[0] == 1.
                    const size_t b = bins[1] + num_bins[1] * bins[2];

                    // Add the weighted magnitude to the angular histogram.
                    histogram[b][0] += magnitude * weight;
                }

        // Sort the angular bins in order to easily find the ones with greatest values.
        // We sort in order of ascending accumulated values (i.e., first element of the triples).
        std::sort(histogram, histogram + num_bins_total, triple_comp<float>(0));

        // Set the dominant orientation of the landmark.
        triple<size_t> location = m_landmarks[i].get_location();
        triple<float>  features = m_landmarks[i].get_features();
        size_t last = num_bins_total - 1;
        features[1] = histogram[last][1];
        features[2] = histogram[last][2];

        new_landmarks.push_back(landmark(location, features));

        #ifdef BIP_VERBOSE_MODE
            std::cout << "   [" << i + 1 << "/" << m_landmarks.size()
                      << "] orientation(s) at position " << location << " = {("
                      << rad2deg(features[1]) << ", "
                      << rad2deg(features[2]) << ")";
        #endif

        // Replicate landmarks which have more than one dominant orientation.
        // Each one receives a different orientation angle.
        for (long b = last - 1; b >= 0; --b) {
            if (histogram[b][0] >= 0.8 * histogram[last][0]) {
                features[1] = histogram[b][1];
                features[2] = histogram[b][2];
                new_landmarks.push_back(landmark(location, features));

                #ifdef BIP_VERBOSE_MODE
                    std::cout << ", ("
                              << rad2deg(features[1]) << ", "
                              << rad2deg(features[2]) << ")";
                #endif
            }
            else
                break;
        }
        #ifdef BIP_VERBOSE_MODE
            std::cout << "}\n";
        #endif
    }

    m_landmarks.clear();

    // Discard the excessive landmarks.
    size_t num_landmarks = std::min(new_landmarks.size(), m_num_max_landmarks);
    for (size_t i = 0; i < num_landmarks; ++i)
        m_landmarks.push_back(new_landmarks[i]);

    delete[] features_map;
}


void
landmark_detector::
compute_landmark_local_descriptors()
{
    // Not implemented.
}


void
landmark_detector::
compute_landmark_global_descriptors()
{
    /*
     * TODO:
     * Parallelize this method by dividing the computations below in multiple threads. Each
     * landmark could be processed in a different thread, so the performance could be greatly
     * improved with a multithreaded approach.
     */

    // The idea is to consider each pixel as a sample point and accumulate edge values in a single
    // log-spherical histogram centered at each landmark. The bins are chosen according to the
    // radial and angular positions of each neighbor pixel, relative to the landmark.
    // The greater the maximum radius, the more "global" is the resulting descriptor.
    // The map of edge values can be the PC map or the greatest eigenvalues map.

    #ifdef BIP_VERBOSE_MODE
        std::cout << "Computing global descriptors:\n";
    #endif

    // Image dimensions.
    const size_t num_dims = (m_sizes[2] == 1) ? 2 : 3;
    
    const size_t total_size = m_sizes[0] * m_sizes[1] * m_sizes[2];
    float *edge_map = read_scalar_map_path(m_map_name);
    normalize_minmax(edge_map, total_size, 0.0, 1.0);

    // char filename_suffix_em[128];
    // sprintf(filename_suffix_em, "edge_map.nii.gz");
    // write_scalar_map(filename_suffix_em, edge_map);

    // Log-spherical sampling parameters (number of bins).
    const size_t num_bins_radius    = 2;
    const size_t num_bins_azimuth   = 10;
    const size_t num_bins_elevation = (num_dims == 2) ? 1 : 5;
    const size_t num_bins_total     = num_bins_radius * num_bins_azimuth * num_bins_elevation;
    const triple<size_t> num_bins(num_bins_radius, num_bins_azimuth, num_bins_elevation);

    // Log-radial and angular histogram. It accumulates the sampled radial distances and
    // orientation angles of each pixel relative to the landmark and use such data to build
    // a descriptor vector.
    landmark::descriptor histogram(num_bins_total, 0.0);

    for (size_t i = 0; i < m_landmarks.size(); ++i)
    {
        #ifdef BIP_VERBOSE_MODE
            std::cout << "   [" << i + 1 << "/" << m_landmarks.size()
                      << "] global descriptor at position "
                      << m_landmarks[i].get_location();
        #endif

        // Clear the histogram values accumulated in the previous iteration.
        histogram.assign(num_bins_total, 0.0);

        // Landmark location and index.
        const size_t x0 = m_landmarks[i].get_location()[0];
        const size_t y0 = m_landmarks[i].get_location()[1];
        const size_t z0 = m_landmarks[i].get_location()[2];
        const size_t j0 = x0 + m_sizes[0] * (y0 + m_sizes[1] * z0);

        // Dominant orientation angles of the landmark.
        const double phi0   = m_landmarks[i].get_features()[1];
        const double theta0 = m_landmarks[i].get_features()[2];
        // std::cout << phi0 << " " << theta0 << "\n";

        const size_t radius = m_global_descr_max_radius;
        const size_t zmin   = (z0 < radius) ? 0 : z0 - radius;
        const size_t zmax   = (z0 + radius >= m_sizes[2]) ? m_sizes[2]-1 : z0 + radius;
        const size_t ymin   = (y0 < radius) ? 0 : y0 - radius;
        const size_t ymax   = (y0 + radius >= m_sizes[1]) ? m_sizes[1]-1 : y0 + radius;
        const size_t xmin   = (x0 < radius) ? 0 : x0 - radius;
        const size_t xmax   = (x0 + radius >= m_sizes[0]) ? m_sizes[0]-1 : x0 + radius;

        // Iterate through a squared region of size 2*radius in each direction.
        for (size_t z = zmin; z <= zmax; ++z)
            for (size_t y = ymin; y <= ymax; ++y)
                for (size_t x = xmin; x <= xmax; ++x)
                {
                    // Get the index of the current neighbor pixel.
                    size_t j = x + m_sizes[0] * (y + m_sizes[1] * z);

                    // Ignore the own landmark location.
                    if (j == j0)
                        continue;

                    // Ignore pixels with very low (< 5%) edge feature value.
                    if (edge_map[j] < 0.05)
                        continue;

                    // // Ignore locations outside the region of interest (mask).
                    // if (m_input_mask && !m_input_mask[j])
                    //     continue;

                    // Get the radial distance between the landmark and the current neighbor pixel.
                    // PS: we need to reverse the y and z axes in the image domain in order to make
                    // them equivalent to cartesian coordinates, and this makes the order of the
                    // subtracted terms change (it is the same as doing y = m_sizes[1] - y and
                    // analogously for y0).
                    const float dx = static_cast<float>(x)  - static_cast<float>(x0);
                    const float dy = static_cast<float>(y0) - static_cast<float>(y);
                    const float dz = static_cast<float>(z0) - static_cast<float>(z);

                    // Get the cartesian coordinates of the difference vector between the landmark
                    // and the current neighbor pixel. Then find its the azimuth and elevation
                    // angles.
                    triple<float> cart(dx, dy, dz);
                    triple<float> sph = cart2sph(cart);

                    // Ignore pixels farther (in rho distance) than the maximum radius.
                    if (sph[0] >= static_cast<float>(radius))
                        continue;

                    // The orientation of the current neighbor pixel is then rotated by the angle
                    // of the landmark orientation, in order to make the descriptor invariant to
                    // rotations. And the magnitude is set as the normalized radial distance.
                    sph[0] /= radius;
                    sph[2] -= theta0;
                    sph[1] -= phi0;

                    // Now we find the log-spherical bin of the current neighbor pixel.
                    const triple<size_t> bins = log_spherical_bin(sph, num_bins);
                    const size_t b = bins[0] + num_bins[0] * (bins[1] + num_bins[1] * bins[2]);

                    // And sum the weighted local PC value to the right bin.
                    // histogram[b] += edge_map[j] / bin_volumes[b];
                    histogram[b] += edge_map[j];
                }

        // Finally, the final histogram is normalized to have unit norm and is attached to the
        // landmark.
        normalize_vnorm(histogram);
        m_landmarks[i].set_global_descriptor(histogram);

        #ifdef BIP_VERBOSE_MODE
            std::cout << " - done\n";
        #endif
    }

    delete[] edge_map;
    // delete[] lambda0;
    // delete[] lambda1;
    // delete[] lambda2;
}


float*
landmark_detector::
read_scalar_map(const char *filename_suffix)
{
    debug::assert2(filename_suffix != NULL);

    char filename[512];
    sprintf(filename, "%s_%s.nii", m_filename_prefix.c_str(), filename_suffix);

    return image2array<float, 3>(read_image<float, 3>(std::string(filename)));
}

float*
landmark_detector::
read_scalar_map_path(std::string filename)
{
    // debug::assert2(filename != NULL);
    return image2array<float, 3>(read_image<float, 3>(filename));
}

landmark_detector::itkVector*
landmark_detector::
read_vectorial_map(const char *filename_suffix)
{
    debug::assert2(filename_suffix != NULL);

    char filename[512];
    sprintf(filename, "%s_%s.nii", m_filename_prefix.c_str(), filename_suffix);

    return image2array<itkVector, 3>(read_image<itkVector, 3>(std::string(filename)));
}


void
landmark_detector::
write_scalar_map(const char *filename_suffix, float *map)
{
    debug::assert2(filename_suffix != NULL);
    debug::assert2(map             != NULL);

    char filename[512];
    sprintf(filename, "%s_%s.nii", m_filename_prefix.c_str(), filename_suffix);

    write_image<float, 3>(std::string(filename),
        array2image<float, 3>(map, m_sizes, m_reference_img));
}


// ------------------------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const landmark_detector &ldet)
{
    return ldet.print(os);
}


}
