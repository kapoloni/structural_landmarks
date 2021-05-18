/**
 * @file   probabilisticlandmarkatlas.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-probabilistic
 * @ingroup    landmark-atlas-probabilistic
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "probabilisticlandmarkatlas.hpp"


namespace bip
{


probabilistic_landmark_atlas::
probabilistic_landmark_atlas(std::string       filename_samples_prefix,
                             triple<size_t>    sizes,
                             size_t            num_samples,
                    typename itkImage::Pointer reference_img,
                             float             kernel_sigma,
                             float             regularization_sigma,
                             size_t            evaluation_k,
                             size_t            covariance_k,
                             float             density_threshold,
                             size_t            averaging_window_radius)
{
    // Initialization of members.
    set_filename_samples_prefix(filename_samples_prefix);
    set_sizes(sizes);
    set_num_samples(num_samples);
    set_reference_img(reference_img);
    set_kernel_sigma(kernel_sigma);
    set_regularization_sigma(regularization_sigma);
    set_evaluation_k(evaluation_k);
    set_covariance_k(covariance_k);
    set_density_threshold(density_threshold);
    set_averaging_window_radius(averaging_window_radius);

    m_point_set = itkPointSet::New();
    m_density   = itkImage::New();

    typename itkImage::SizeType sizes2;
    sizes2[0] = sizes[0];
    sizes2[1] = sizes[1];
    sizes2[2] = sizes[2];

    typename itkImage::IndexType index;
    index[0] = index[1] = index[2] = 0;

    typename itkImage::RegionType region;
    region.SetSize(sizes2);
    region.SetIndex(index);

    m_density->SetRegions(region);
    m_density->Allocate();
    m_density->FillBuffer(0);
}


probabilistic_landmark_atlas::
~probabilistic_landmark_atlas()
{
    // Nothing.
}


void
probabilistic_landmark_atlas::
compute()
{
    compute_density();
    compute_average_landmarks();
}


std::ostream&
probabilistic_landmark_atlas::
print(std::ostream &os) const
{
    os << "{"
       << "}";

    return os;
}


void
probabilistic_landmark_atlas::
compute_density()
{
    #ifdef BIP_VERBOSE_MODE
        std::cout << "Estimating the p.d.f. of landmark locations";
    #endif

    // Repeatedly read samples (landmark files).
    for (size_t n = 0; n < m_num_samples; ++n) {
        char filename_sample[512];

        sprintf(filename_sample, "%s_%lu.txt", m_filename_samples_prefix.c_str(), n);
        std::vector<landmark> landmarks(read_landmarks(filename_sample));

        // Add each landmark of the current sample to the point set.
        for (size_t i = 0; i < landmarks.size(); ++i) {
            typename itkPointSet::PointType point;
            point[0] = landmarks[i].get_location()[0];
            point[1] = landmarks[i].get_location()[1];
            point[2] = landmarks[i].get_location()[2];

            m_point_set->GetPoints()->InsertElement(m_point_set->GetNumberOfPoints(), point);
        }
    }

    // There's no reason to continue if we have no point in the set, so check first.
    if (m_point_set->GetNumberOfPoints() > 0) {
        typedef itk::ManifoldParzenWindowsPointSetFunction<itkPointSet, float>
                itkManifoldParzenWindows;

        // Create and configure the manifold parzen windows object.
        typename itkManifoldParzenWindows::Pointer mpw = itkManifoldParzenWindows::New();
        mpw->SetInputPointSet(m_point_set);
        mpw->SetKernelSigma(m_kernel_sigma);
        mpw->SetRegularizationSigma(m_regularization_sigma);
        mpw->SetEvaluationKNeighborhood(m_evaluation_k);
        mpw->SetCovarianceKNeighborhood(m_covariance_k);
        mpw->NormalizeOn();
        mpw->UseAnisotropicCovariancesOn();

        // Set the value of the estimated density at each image location.
        // Such value is computed from the whole set of points which come from all added samples.
        // Depending on the Manifold Parzen Windows' parameters, this process can be very slow.
        for (size_t z = 0; z < m_sizes[2]; ++z)
            for (size_t y = 0; y < m_sizes[1]; ++y)
                for (size_t x = 0; x < m_sizes[0]; ++x)
                {
                    typename itkImage::IndexType    index;
                    typename itkPointSet::PointType point;

                    index[0] = point[0] = x;
                    index[1] = point[1] = y;
                    index[2] = point[2] = z;

                    m_density->SetPixel(index, mpw->Evaluate(point));
                }
    }

    // Adjust origin, spacing and direction accordiong to the reference image.
    m_density->SetOrigin(m_reference_img->GetOrigin());
    m_density->SetSpacing(m_reference_img->GetSpacing());
    m_density->SetDirection(m_reference_img->GetDirection());

    #ifdef BIP_VERBOSE_MODE
        std::cout << " - done\n";
    #endif
}


void
probabilistic_landmark_atlas::
compute_average_landmarks()
{
    #ifdef BIP_VERBOSE_MODE
        std::cout << "Computing average landmarks:\n";
    #endif

    // Get the local maxima in the densities map.
    std::vector<triple<size_t> > local_maxima = find_density_local_maxima();
    const size_t dims = (m_sizes[2] == 1) ? 2 : 3;

    // Number of possible dominant orientations for landmarks.
    // Here we may need to use a reduced set of angle ranges for the sake of performance.
    const size_t num_azimuths   = 8;                   // originally 36.
    const size_t num_elevations = (dims == 2) ? 1 : 4; // originally 1 or 18.
    const size_t num_ori_total  = num_azimuths * num_elevations;
    const triple<size_t> num_oris(1, num_azimuths, num_elevations);

    // // Parameters of the local descriptor.
    const size_t loc_ab = 10;                   // number of azimuth bins.
    const size_t loc_eb = 5;                   // number of elevation bins.
    const size_t loc_rb = 5;                   // number of radial bins.

    // Parameters of the global descriptor.
    // const size_t glob_rb = 4;                   // number of radial bins.
    // const size_t glob_ab = 8;                   // number of azimuth bins.
    // const size_t glob_eb = (dims == 2) ? 1 : 4; // number of elevation bins.
    const size_t radius = 32;
    const size_t glob_rb = 0;                   // number of radial bins.
    const size_t glob_ab = 0;                   // number of azimuth bins.
    const size_t glob_eb = 0; // number of elevation bins.

    //const size_t num_local_descr_bins  = pow(loc_sh, dims) * loc_ab * loc_eb;
    const size_t num_local_descr_bins = loc_rb * loc_ab * loc_eb;
    const size_t num_global_descr_bins = glob_rb * glob_ab * glob_eb;

    std::cout << "local " << num_local_descr_bins << "\n";
    // The average landmarks are computed separately for each discrete orientation. This is
    // necessary because, if we average landmarks with different dominant orientations, we would
    // lose the rotational invariance.
    // So we have, for each kind of descriptor, an MxO sized matrix of B-sized vectors, where M is
    // the number of local maxima points, O is the total number of possible landmark orientation
    // angles, and B is the total number of descriptor bins (so it is a matrix of descriptors).
    // The MxO matrix can be represented as a single array, but it must be dinamically allocated.
    landmark::descriptor *avg_local_descr =
        new landmark::descriptor[local_maxima.size() * num_ori_total]();
    landmark::descriptor *avg_global_descr =
        new landmark::descriptor[local_maxima.size() * num_ori_total]();

    // As we already know the length of the descriptors, we can initialize their vectors.
    for (size_t i = 0, m = 0; m < local_maxima.size(); ++m)
        for (size_t o = 0; o < num_ori_total; ++o, ++i) {
            avg_local_descr[o + num_ori_total * m].resize(num_local_descr_bins, 0.0);
            avg_global_descr[o + num_ori_total * m].resize(num_global_descr_bins, 0.0);
        }

    // We also need a similar structure to store the countings of landmarks added to each average
    // relative to orientation (so we can in fact average the results in the end).
    // In this case, a triple is used because later these countings will be partially sorted and
    // we need to keep track of the local maxima and orientation they represent. The counting
    // itself is stored in the first element of the triple.
    triple<size_t> *avg_countings = new triple<size_t>[local_maxima.size() * num_ori_total]();

    // Assign the local maxima and orientation indices to the counting triples.
    for (size_t i = 0, m = 0; m < local_maxima.size(); ++m)
        for (size_t o = 0; o < num_ori_total; ++o, ++i) {
            avg_countings[i][1] = o;
            avg_countings[i][2] = m;
        }

    // Repeatedly read samples (landmark files).
    for (size_t n = 0; n < m_num_samples; ++n) {
        #ifdef BIP_VERBOSE_MODE
            std::cout << "   Processing sample " << n + 1 << "/" << m_num_samples;
        #endif

        char filename_sample[512];
        std::cout << filename_sample << " ";

        sprintf(filename_sample, "%s_%lu.txt", m_filename_samples_prefix.c_str(), n);
        std::vector<landmark> landmarks(read_landmarks(filename_sample));

        // Check whether a landmark will take part in the average.
        for (size_t l = 0; l < landmarks.size(); ++l) {
            triple<size_t>       l_location     = landmarks[l].get_location();
            landmark::descriptor l_local_descr  = landmarks[l].get_local_descriptor();
            landmark::descriptor l_global_descr = landmarks[l].get_global_descriptor();
            //std::cout << filename_sample << " " << l_location << " ";
            //std::cout << filename_sample << " global " << l_local_descr[num_local_descr_bins-1];
            // We check circular areas around each local maxima. If a landmark is located
            // inside it, then it takes part in the average landmark located at such local
            // maximum point.
            for (size_t m = 0; m < local_maxima.size(); ++m) {
                //std::cout << "posic " << local_maxima[m][0] << " " << local_maxima[m][1] << " " << local_maxima[m][2] << " ";
                const float dist = sqrt(
                    sqr(static_cast<long>((l_location[0] - local_maxima[m][0]))) +
                    sqr(static_cast<long>((l_location[1] - local_maxima[m][1]))) +
                    sqr(static_cast<long>((l_location[2] - local_maxima[m][2]))));

                if (dist <= m_averaging_window_radius) {
                    // Get the dominant orientation of the landmark (in the original ranges).
                    triple<float> orientation = landmarks[l].get_features();
                    orientation[0] = EPSILON; // not REALLY necessary, but just in case...
                    //std::cout << "orientation " << orientation[0] << " " << orientation[1] << " " << orientation[2] << " ";

                    // Find the index of the dominant orientation of the landmark.
                    // const triple<size_t> oris = log_spherical_bin(orientation, num_oris);
                    // const size_t o = oris[0] + num_oris[0] * (oris[1] + num_oris[1] * oris[2]);

                    triple<float> cart = sph2cart(orientation);
                    const float lenght = cart[0];
                    const int o = spider_web_bin(cart, lenght, radius, num_oris);
                    // std::cout << "sw "<< o << "\n";

                    // Increment the counting of landmarks added to the average relative to the
                    // local maxima m and orientation o.
                    avg_countings[o + num_ori_total * m][0] += 1;

                    // Accumulate the landmark descriptors' data to the average relative to the
                    // local maxima m and orientation o.
                    for (size_t b = 0; b < num_local_descr_bins; ++b)
                        avg_local_descr[o + num_ori_total * m][b] += l_local_descr[b];
                    for (size_t b = 0; b < num_global_descr_bins; ++b)
                        avg_global_descr[o + num_ori_total * m][b] += l_global_descr[b];
                }
            }
        }
        #ifdef BIP_VERBOSE_MODE
            std::cout << " - done\n";
        #endif
    }

    // Now we need to normalize the accumulated values, dividing them by the countings.
    for (size_t i = 0, m = 0; m < local_maxima.size(); ++m)
        for (size_t o = 0; o < num_ori_total; ++o, ++i) {
            if (avg_countings[i][0] > 0) {
                for (size_t b = 0; b < num_local_descr_bins; ++b)
                    avg_local_descr[i][b] /= avg_countings[i][0];
                for (size_t b = 0; b < num_global_descr_bins; ++b)
                  avg_global_descr[i][b] /= avg_countings[i][0];
            }
        }

    // The countings vector of each local maximum point can also be seen as a histogram of the
    // landmark orientations. So we can use this histogram to determine which of all orientations
    // are going to be assigned to the averaged landmark, in the same way the landmark detector
    // does to determine the orientations of each detected landmark.
    // PS: notice that, for now, we DON'T assign any saliency or orientation to the averaged
    // landmarks. I think that the atlas will not need that...
    for (size_t m = 0; m < local_maxima.size(); ++m) {
        // Sort the countings for the current local maximum. Greatest values appear last.
        std::sort(&avg_countings[m     * num_ori_total],
                  &avg_countings[(m+1) * num_ori_total], triple_comp<size_t>(0));

        const size_t countmax = avg_countings[num_ori_total-1 + num_ori_total * m][0];
        const size_t omax     = avg_countings[num_ori_total-1 + num_ori_total * m][1];

        // This is VERY unlikely to happen, but just in case...
        if (countmax == 0)
            continue;

        // std::cout << avg_global_descr[omax + num_ori_total * m][93] << "\n";
        // Create the first averaged landmark. Its descriptors are those associated to the
        // orientation that has the greatest count.
        m_average_landmarks.push_back(landmark(local_maxima[m],
                                               triple<float>(),
                                               avg_local_descr[omax  + num_ori_total * m],
                                               avg_global_descr[omax + num_ori_total * m]));

        #ifdef BIP_VERBOSE_MODE
            std::cout << "   [" << m + 1 << "/" << local_maxima.size()
                      << "] average landmark(s) created at position "
                      << local_maxima[m];
        #endif

        size_t copies = 1;

        // Replicate landmarks that have high counts in other orientations too.
        for (long o = num_ori_total - 2; o >= 0; --o) {
            const size_t ocount   = avg_countings[o    + num_ori_total * m][0];
            const size_t countmax = avg_countings[omax + num_ori_total * m][0];

            if (ocount > 0 && ocount >= 0.8 * countmax) {
                // Its descriptors are associated to the orientation o.
                m_average_landmarks.push_back(landmark(local_maxima[m],
                                                       triple<float>(),
                                                       avg_local_descr[o  + num_ori_total * m],
                                                       avg_global_descr[o + num_ori_total * m]));
                ++copies;
            } else
                break;
        }

        #ifdef BIP_VERBOSE_MODE
            std::cout << " = " << copies << "\n";
        #endif
    }

    // Free memory.
    delete[] avg_countings;
    delete[] avg_local_descr;
    delete[] avg_global_descr;
}


std::vector<triple<size_t> >
probabilistic_landmark_atlas::
find_density_local_maxima()
{
    const size_t total_size = m_sizes[0] * m_sizes[1] * m_sizes[2];

    std::vector<triple<size_t> > local_maxima;

    // Put the density map in a normalized scale ([0, 1]). This makes the threshold easier.
    float *norm_density_map = image2array<float, 3>(m_density);
    normalize_minmax(norm_density_map, total_size, 0.0, 1.0);

    // Iterate through the whole image.
    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i)
            {
                // Ignore values which are below the relative threshold. They can't be considered
                // maximal points in this context.
                if (norm_density_map[i] < m_density_threshold)
                    continue;

                const triple<size_t> location(x, y, z);

                if (local_max(norm_density_map, location, m_sizes))
                    local_maxima.push_back(location);
            }

    delete[] norm_density_map;
    return local_maxima;
}


// ------------------------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const probabilistic_landmark_atlas &pla)
{
    return pla.print(os);
}


}
