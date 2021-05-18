/**
 * @file   phasecongruency.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup phasecongruency
 * @ingroup    phasecongruency
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "phasecongruency.hpp"


namespace bip
{


phase_congruency::
phase_congruency(std::string           filename_prefix,
                 triple<size_t>        sizes,
                 loggabor_filter_bank *filter_bank,
                 float                *input_img,
                 unsigned char        *input_mask,
        typename itkImage::Pointer     reference_img,
                 float                 noise_threshold,
                 float                 noise_std,
                 float                 sigmoid_gain,
                 float                 sigmoid_cutoff)
{
    // Initialization of members.
    set_filename_prefix(filename_prefix);
    set_filter_bank(filter_bank);
    set_sizes(sizes);
    set_input_img(input_img);
    set_input_mask(input_mask);
    set_reference_img(reference_img);
    set_noise_threshold(noise_threshold);
    set_noise_std(noise_std);
    set_sigmoid_gain(sigmoid_gain);
    set_sigmoid_cutoff(sigmoid_cutoff);
}


phase_congruency::
~phase_congruency()
{
    // Nothing.
}


std::ostream&
phase_congruency::
print(std::ostream &os) const
{
    os << "{"
       << "m_filename_prefix: " << m_filename_prefix << ", "
       << "m_filter_bank: "     << m_filter_bank     << ", "
       << "m_sizes: "           << m_sizes           << ", "
       << "m_input_img: "       << m_input_img       << ", "
       << "m_input_mask: "      << m_input_mask      << ", "
       << "m_reference_img: "   << m_reference_img   << ", "
       << "m_noise_threshold: " << m_noise_threshold << ", "
       << "m_noise_std: "       << m_noise_std       << ", "
       << "m_sigmoid_gain: "    << m_sigmoid_gain    << ", "
       << "m_sigmoid_cutoff: "  << m_sigmoid_cutoff  << ", "
       << "}";

    return os;

}


void
phase_congruency::
compute()
{
    /*
     * TODO:
     * Parallelize this method by dividing the computations below in multiple threads. There are
     * many places in which this can be done. For example, each log-Gabor filter is convolved by
     * the input image (in frequency domain) independently of the others. Also, most computations
     * are done pixelwise, so the pixel access could be parallelized too. The performance can be
     * severely improved in a multicore environment, but be wise in choosing what and where to
     * parallelize.
     */

    // --------------------------------------------------------------------------------------------
    // Step 1: compute the Fourier transform of the input image.
    // --------------------------------------------------------------------------------------------

    #ifdef BIP_VERBOSE_MODE
        printf("Computing FFT");
    #endif

    // Get total sizes in space and frequency domains.
    size_t total_size    = m_sizes[0] * m_sizes[1] * m_sizes[2];
    size_t ft_total_size = total_size * sizeof(fftwf_complex);

    // Image data for the Fourier transform of the input image.
    fftwf_complex *ft_input_img = (fftwf_complex*) fftwf_malloc(ft_total_size);

    // Array of images filtered in frequency domain.
    // It is one element per scale of the bank of filters.
    fftwf_complex *ft_filtered_imgs = (fftwf_complex*)
        fftwf_malloc(ft_total_size * m_filter_bank->get_num_scales());

    if (!ft_input_img || !ft_filtered_imgs)
        throw std::bad_alloc();

    // Shift the DC component to the center of image.
    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                ft_input_img[i][0] = m_input_img[i] * pow(-1.0, x + y + z);
                ft_input_img[i][1] = 0.0;
            }

    // Allow multithreading for computations of the FFT.
    fftwf_init_threads();
    fftwf_plan_with_nthreads(DEFAULT_FFTW_NUMBER_OF_THREADS);

    // Compute the forward FFT of the input image.
    fast_fourier_transform(ft_input_img, m_sizes);

    // (BIP_DEBUG_MODE) Save the Fourier spectrum.
    #ifdef BIP_DEBUG_MODE
    {
        float *ft_spectrum = new float[total_size]();

        for (size_t i = 0; i < total_size; ++i)
            ft_spectrum[i] = log(1.0 + sqrt(sqr(ft_input_img[i][0]) + sqr(ft_input_img[i][1])));

        normalize_minmax(ft_spectrum, total_size, 0.0, 1.0);

        write_scalar_map("spectrum", ft_spectrum);
        delete[] ft_spectrum;
    }
    #endif

    #ifdef BIP_VERBOSE_MODE
        puts(" - done");
    #endif

    // --------------------------------------------------------------------------------------------
    // Step 2: allocate memory for temporary and member arrays.
    // --------------------------------------------------------------------------------------------

    char filename_suffix[128];

    float     *sum_amplitude          = new float[total_size]();
    float     *max_amplitude          = new float[total_size]();
    float     *total_sum_amplitude    = new float[total_size]();
    float     *total_sum_energy       = new float[total_size]();
    float     *cov_xx                 = new float[total_size]();
    float     *cov_xy                 = new float[total_size]();
    float     *cov_xz                 = new float[total_size]();
    float     *cov_yy                 = new float[total_size]();
    float     *cov_yz                 = new float[total_size]();
    float     *cov_zz                 = new float[total_size]();
    // float     *mean_x                 = new float[total_size]();
    // float     *mean_y                 = new float[total_size]();
    // float     *mean_z                 = new float[total_size]();
    float     *pc_map                 = new float[total_size]();
    float     *directional_pc_map     = new float[total_size]();
    itkVector *directional_pc_max_map = new itkVector[total_size]();

    float *moments_eigenvalues_maps[3] = {
        new float[total_size](),
        new float[total_size](),
        new float[total_size]()
    };
    itkVector *moments_eigenvectors_maps[3] = {
        new itkVector[total_size](),
        new itkVector[total_size](),
        new itkVector[total_size]()
    };

    // --------------------------------------------------------------------------------------------
    // Step 3: compute the phase congruency.
    // --------------------------------------------------------------------------------------------

    #ifdef BIP_VERBOSE_MODE
        printf("Computing PC maps...\n");
    #endif

    // Current orientation's index.
    size_t o = 0;

    // Get the variation in elevation angles.
    float dtheta = (m_filter_bank->get_num_elevations() == 1) ?
                    0.0 :
                    M_PI_2 / (m_filter_bank->get_num_elevations() - 1);

    for (size_t e = 0; e < m_filter_bank->get_num_elevations(); ++e)
    {
        // Get the current elevation angle.
        float theta     = e * dtheta;
        float cos_theta = cos(theta);
        float sin_theta = sin(theta);

        // Get the variation in azimuth angles.
        float dphi = (m_filter_bank->get_num_azimuths() == 1) ?
                      0.0 :
                      (e == 0) ?
                          M_PI   / m_filter_bank->get_num_azimuths_per_elevation(0) :
                          M_PI*2 / m_filter_bank->get_num_azimuths_per_elevation(e);

        for (size_t a = 0; a < m_filter_bank->get_num_azimuths_per_elevation(e); ++a)
        {
            // Automatically calculated noise threshold constant.
            float T = 0.0;

            // Get the current azimuth angle.
            float phi     = a * dphi;
            float cos_phi = cos(phi);
            float sin_phi = sin(phi);

            // Initialize arrays used in the current orientation only.
            memset(sum_amplitude, 0, total_size * sizeof(float));
            memset(max_amplitude, 0, total_size * sizeof(float));

            for (size_t s = 0; s < m_filter_bank->get_num_scales(); ++s)
            {
                // Get a single log-Gabor filter for the current scale, azimuth and elevation.
                float *ft_filter = m_filter_bank->get_filter(s, a, e);

                #ifdef BIP_VERBOSE_MODE
                    printf("   Processing filter: s = %02lu, a = %02lu, e = %02lu", s, a, e);
                #endif

                // Pointer to the current scale's filtered image.
                fftwf_complex *p_ft_filtered_img = &ft_filtered_imgs[s * total_size];

                // Apply the log-Gabor filter to the input image in frequency domain. After that,
                // we don't need such filter anymore.
                for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
                    for (size_t y = 0; y < m_sizes[1]; ++y)
                        for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                            p_ft_filtered_img[i][0] = ft_filter[i] * ft_input_img[i][0];
                            p_ft_filtered_img[i][1] = ft_filter[i] * ft_input_img[i][1];
                        }

                delete[] ft_filter;

                // Use normalized backward FFT in order to get the filtered image in the spatial
                // domain.
                fast_fourier_transform(p_ft_filtered_img, m_sizes, true);

                // Shift the DC component back to the original location.
                for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
                    for (size_t y = 0; y < m_sizes[1]; ++y)
                        for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                            p_ft_filtered_img[i][0] *= pow(-1.0, x + y + z);
                            p_ft_filtered_img[i][1] *= pow(-1.0, x + y + z);
                        }

                // (BIP_DEBUG_MODE) Save the filter responses (even, odd, amplitude).
                #ifdef BIP_DEBUG_MODE
                {
                    float *ft_even      = new float[total_size]();
                    float *ft_odd       = new float[total_size]();
                    float *ft_amplitude = new float[total_size]();

                    for (size_t i = 0; i < total_size; ++i) {
                        ft_even[i]      = p_ft_filtered_img[i][0];
                        ft_odd[i]       = p_ft_filtered_img[i][1];
                        ft_amplitude[i] = sqrt(sqr(ft_even[i]) + sqr(ft_odd[i]));
                    }
                    normalize_minmax(ft_even,      total_size, -1.0, 1.0);
                    normalize_minmax(ft_odd,       total_size, -1.0, 1.0);
                    normalize_minmax(ft_amplitude, total_size,  0.0, 1.0);

                    sprintf(filename_suffix, "even_%02lu_%02lu_%02lu.nii", s, a, e);
                    write_scalar_map(filename_suffix, ft_even);

                    sprintf(filename_suffix, "odd_%02lu_%02lu_%02lu.nii", s, a, e);
                    write_scalar_map(filename_suffix, ft_odd);

                    sprintf(filename_suffix, "amplitude_%02lu_%02lu_%02lu.nii", s, a, e);
                    write_scalar_map(filename_suffix, ft_amplitude);

                    delete[] ft_even;
                    delete[] ft_odd;
                    delete[] ft_amplitude;
                }
                #endif

                // Accumulate amplitude responses along scales.
                for (size_t i = 0; i < total_size; ++i) {
                    // Ignore locations outside the region of interest (mask).
                    if (m_input_mask && !m_input_mask[i])
                      continue;


                    float even      = p_ft_filtered_img[i][0];
                    float odd       = p_ft_filtered_img[i][1];
                    float amplitude = sqrt(sqr(even) + sqr(odd));

                    sum_amplitude[i] += amplitude;
                    max_amplitude[i] = std::max(amplitude, max_amplitude[i]);
                }

                // Noise threshold estimation (if required).
                if (s == 0) {
                    if (m_noise_threshold < 0.0) {
                        float tau         = median(sum_amplitude, total_size) / sqrt(log(4.0));
                        float invmult     = 1.0 / m_filter_bank->get_mult_factor();
                        float nscales     = m_filter_bank->get_num_scales();
                        float total_tau   = tau * (1.0 - pow(invmult, nscales)) / (1.0 - invmult);
                        float noise_mean  = total_tau * sqrt(M_PI_2);
                        float noise_sigma = total_tau * sqrt((4.0 - M_PI) / 2.0);

                        T = noise_mean + m_noise_std * noise_sigma;

                        // (BIP_DEBUG_MODE) Print the estimated threshold.
                        #ifdef BIP_DEBUG_MODE
                            printf(" (estimated T = %g)", T);
                        #endif
                    }
                }

                #ifdef BIP_VERBOSE_MODE
                    puts(" - done");
                #endif
            }

            // We use the same block of memory for all directional PC maps, so we have to clean
            // the data from the previous orientation before we continue.
            memset(directional_pc_map, 0.0, sizeof(directional_pc_map));

            // Now the real thing begins...
            // We change the loops order to be able to save some memory space.
            for (size_t i = 0; i < total_size; ++i)
            {
                // Ignore locations outside the region of interest (mask).
                if (m_input_mask && !m_input_mask[i])
                    continue;

                float sum_even = 0.0;
                float sum_odd  = 0.0;

                // Compute the even and odd mean responses.
                for (size_t s = 0; s < m_filter_bank->get_num_scales(); ++s)
                {
                    fftwf_complex *p_ft_filtered_img = &ft_filtered_imgs[s * total_size];

                    sum_even += p_ft_filtered_img[i][0];
                    sum_odd  += p_ft_filtered_img[i][1];
                }
                float norm      = sqrt(sqr(sum_even) + sqr(sum_odd));
                float mean_even = sum_even / (norm + EPSILON);
                float mean_odd  = sum_odd / (norm + EPSILON);

                float local_energy = 0.0;

                // Compute the energy and phase congruency for the current orientation.
                for (size_t s = 0; s < m_filter_bank->get_num_scales(); ++s)
                {
                    fftwf_complex *p_ft_filtered_img = &ft_filtered_imgs[s * total_size];

                    float even = p_ft_filtered_img[i][0];
                    float odd  = p_ft_filtered_img[i][1];

                    local_energy += (even * mean_even + odd * mean_odd) -
                                fabs(even * mean_odd - odd * mean_even);
                }

                // Apply the energy threshold (the automatically computed T constant or a value
                // given by the user).
                if (m_noise_threshold < 0.0)
                    local_energy -= T;
                else
                    local_energy -= m_noise_threshold;

                // Apply the frequency distribution weighting on the energy value (if needed).
                if (local_energy > 0.0) {
                    // If there is only one non-zero component, width takes on a value of 0.
                    // If all components are equal, width is 1.
                    float width = (sum_amplitude[i] / (max_amplitude[i] + EPSILON) - 1.0) /
                                   (m_filter_bank->get_num_scales() - 1.0);

                    // Sigmoidal weighting function for this orientation.
                    float weight = 1.0 + exp(m_sigmoid_gain * (m_sigmoid_cutoff - width));

                    local_energy /= weight;
                }
                else
                    local_energy = 0.0;

                // Accumulate the total sums in amplitude and energy for the final phase congruency
                // map.
                total_sum_amplitude[i] += sum_amplitude[i];
                total_sum_energy[i]    += local_energy;

                // Compute the local phase congruency for the current location and orientation.
                float local_pc = local_energy / (sum_amplitude[i] + EPSILON);

                // Set the value of the current voxel for the directional PC map of the current
                // orientation.
                directional_pc_map[i] = local_pc;

                // Update the pixel at the directional PC maxima map.
                if (local_pc > directional_pc_max_map[i][0]) {
                    directional_pc_max_map[i][0] = local_pc;
                    directional_pc_max_map[i][1] = phi;
                    directional_pc_max_map[i][2] = theta;
                }

                // Projections in the cartesian space.
                float proj_x = local_pc * cos_theta * cos_phi;
                float proj_y = local_pc * cos_theta * sin_phi;
                float proj_z = local_pc * sin_theta;

                // Accumulate covariance values along the orientations.
                // mean_x[i] += proj_x;
                // mean_y[i] += proj_y;
                // mean_z[i] += proj_z;
                cov_xx[i] += sqr(proj_x);
                cov_yy[i] += sqr(proj_y);
                cov_zz[i] += sqr(proj_z);
                cov_xy[i] += proj_x * proj_y;
                cov_xz[i] += proj_x * proj_z;
                cov_yz[i] += proj_y * proj_z;
            }

            // Write the directional PC map.
            // sprintf(filename_suffix, "directional_PC_%lu", o);
            // write_scalar_map(filename_suffix, directional_pc_map);

            ++o; // Next orientation.
        }
    }

    // Write the directional PC maxima map.
    // write_vectorial_map("directional_PC_max", directional_pc_max_map);

    // Compute and write the final phase congruency map.
    for (size_t i = 0; i < total_size; ++i) {
        pc_map[i] = total_sum_energy[i] / (total_sum_amplitude[i] + EPSILON);
    }
    write_scalar_map("PC", pc_map);

    // Free memory.
    delete[] sum_amplitude;
    delete[] max_amplitude;
    delete[] total_sum_amplitude;
    delete[] total_sum_energy;
    delete[] pc_map;
    delete[] directional_pc_map;
    delete[] directional_pc_max_map;

    fftwf_free(ft_input_img);
    fftwf_free(ft_filtered_imgs);
    fftwf_cleanup_threads();

    // --------------------------------------------------------------------------------------------
    // Step 4: compute the maps of moments eigenvalues and eigenvectors.
    // --------------------------------------------------------------------------------------------

    #ifdef BIP_VERBOSE_MODE
        printf("Computing PC moments, eigenvalues and eigenvectors");
    #endif

    // Covariance normalization factor.
    float orientations      = static_cast<float>(m_filter_bank->get_num_orientations());
    float half_orientations = 0.5 * orientations;
    // float norm = 3.0 * orientations;

    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i)
            {
                // Ignore locations outside the region of interest (mask).
                if (m_input_mask && !m_input_mask[i])
                    continue;

                // Set the coefficients of the real covariance matrix. These coefficients are
                // calculated from the image raw moments, which we already have.
                // cov_xx[i] = cov_xx[i] / norm - sqr(mean_x[i] / norm);
                // cov_xy[i] = (cov_xy[i] - mean_x[i] * mean_y[i]) / norm;
                // cov_xz[i] = (cov_xz[i] - mean_x[i] * mean_z[i]) / norm;
                // cov_yy[i] = cov_yy[i] / norm - sqr(mean_y[i] / norm);
                // cov_yz[i] = (cov_yz[i] - mean_y[i] * mean_z[i]) / norm;
                // cov_zz[i] = cov_zz[i] / norm - sqr(mean_z[i] / norm);

                cov_xx[i] /= half_orientations;
                cov_xy[i] /= orientations;
                cov_xz[i] /= orientations;
                cov_yy[i] /= half_orientations;
                cov_yz[i] /= orientations;
                cov_zz[i] /= half_orientations;

                // Compute the eigenvalues and eigenvectors.
                // When we have 1D or 2D data, the covariances related to Y and/or Z directions
                // are going to be 0, so it will be the same as computing the eigenvalue and
                // eigenvector of an 1D or 2D matrix.
                double M[9] = {
                    cov_xx[i], cov_xy[i], cov_xz[i],
                    cov_xy[i], cov_yy[i], cov_yz[i],
                    cov_xz[i], cov_yz[i], cov_zz[i]
                };
                double eigenvalues[3];
                double eigenvectors[9];
                eigens(M, eigenvalues, eigenvectors);

                // Set the values at the current pixel for the eigenvalue and eigenvector maps.
                for (size_t d = 0; d < 3; ++d) {
                    moments_eigenvalues_maps[d][i]     = eigenvalues[d];
                    // moments_eigenvectors_maps[d][i][0] = eigenvalues[d] * eigenvectors[3*d];
                    // moments_eigenvectors_maps[d][i][1] = eigenvalues[d] * eigenvectors[3*d + 1];
                    // moments_eigenvectors_maps[d][i][2] = eigenvalues[d] * eigenvectors[3*d + 2];
                    moments_eigenvectors_maps[d][i][0] = eigenvectors[3*d];
                    moments_eigenvectors_maps[d][i][1] = eigenvectors[3*d + 1];
                    moments_eigenvectors_maps[d][i][2] = eigenvectors[3*d + 2];
                }
            }

    // Write all the eigenvalues and eigenvectors maps.
    for (size_t d = 0; d < 3; ++d) {
        sprintf(filename_suffix, "eigenvalues_%lu", d);
        write_scalar_map(filename_suffix, moments_eigenvalues_maps[d]);

        // sprintf(filename_suffix, "eigenvectors_%lu", d);
        // write_vectorial_map(filename_suffix, moments_eigenvectors_maps[d]);
    }

    // Free memory.
    // delete[] mean_x;
    // delete[] mean_y;
    // delete[] mean_z;
    delete[] cov_xx;
    delete[] cov_xy;
    delete[] cov_xz;
    delete[] cov_yy;
    delete[] cov_yz;
    delete[] cov_zz;
    for (size_t d = 0; d < 3; ++d) {
        delete[] moments_eigenvalues_maps[d];
        delete[] moments_eigenvectors_maps[d];
    }

    #ifdef BIP_VERBOSE_MODE
        puts(" - done");
    #endif
}


void
phase_congruency::
write_scalar_map(const char *filename_suffix, float *map)
{
    debug::assert2(filename_suffix != NULL);
    debug::assert2(map             != NULL);

    char filename[512];
    sprintf(filename, "%s_%s.nii.gz", m_filename_prefix.c_str(), filename_suffix);

    write_image<float, 3>(std::string(filename),
        array2image<float, 3>(map, m_sizes, m_reference_img));
}


void
phase_congruency::
write_vectorial_map(const char *filename_suffix, phase_congruency::itkVector *map)
{
    debug::assert2(filename_suffix != NULL);
    debug::assert2(map             != NULL);

    char filename[512];
    sprintf(filename, "%s_%s.nii.gz", m_filename_prefix.c_str(), filename_suffix);

    write_image<itkVector, 3>(std::string(filename),
        array2image<itkVector, 3>(map, m_sizes));
}


// ------------------------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const phase_congruency &pc)
{
    return pc.print(os);
}


}
