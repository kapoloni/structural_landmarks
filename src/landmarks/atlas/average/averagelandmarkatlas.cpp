/**
 * @file   averagelandmarkatlas.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-average
 * @ingroup    landmark-atlas-average
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "averagelandmarkatlas.hpp"


namespace bip
{


average_landmark_atlas::
average_landmark_atlas(std::string       filename_samples_prefix,
                       triple<size_t>    sizes,
                       size_t            num_samples,
              typename itkImage::Pointer reference_img,
                       float             rel_frequency_threshold,
                       size_t            max_countdown_value,
                       float             descriptors_tradeoff,
                       float             max_descriptor_dist,
                       float             max_location_dist) :
    m_average_sample(NULL)
{
    // Initialization of members.
    set_filename_samples_prefix(filename_samples_prefix);
    set_sizes(sizes);
    set_num_samples(num_samples);
    set_reference_img(reference_img);
    set_rel_frequency_threshold(rel_frequency_threshold);
    set_max_countdown_value(max_countdown_value);
    set_descriptors_tradeoff(descriptors_tradeoff);
    set_max_descriptor_dist(max_descriptor_dist);
    set_max_location_dist(max_location_dist);
}


average_landmark_atlas::
~average_landmark_atlas()
{
    if (m_average_sample != NULL)
        delete m_average_sample;
}


void
average_landmark_atlas::
compute()
{
    #ifdef BIP_VERBOSE_MODE
        std::cout << "Averaging matching landmarks from samples\n";
    #endif

    for (size_t n = 0; n < m_num_samples; ++n) {
        char filename_sample[512];
        sprintf(filename_sample, "%s_%lu.txt", m_filename_samples_prefix.c_str(), n);

        // Read landmarks and put them into a deque.
        std::vector<landmark> temp(read_landmarks(filename_sample));
        std::deque<landmark> landmarks(temp.begin(), temp.end());

        // Create a new atlas sample containing the landmarks.
        atlas_sample sample(landmarks, m_max_countdown_value);

        #ifdef BIP_VERBOSE_MODE
            std::cout << "   Processing sample " << n+1 << "/" << m_num_samples
                      << " (" << sample.get_landmarks().size() << " landmarks)";
        #endif

        process_sample(sample);

        #ifdef BIP_VERBOSE_MODE
            std::cout << " - done\n";
        #endif
    }

    discard_outliers();
}


void
average_landmark_atlas::
process_sample(const atlas_sample &sample)
{
    if (m_average_sample == NULL) {
        // First sample.
        m_average_sample = new atlas_sample(sample);
    } else {
        // Convert from deque to vector.
        std::vector<landmark> landmarks1(m_average_sample->get_landmarks().begin(),
                                         m_average_sample->get_landmarks().end());
        std::vector<landmark> landmarks2(sample.get_landmarks().begin(),
                                         sample.get_landmarks().end());

        // Find matchings between the landmarks of the current average sample and the new sample.
        matches_vector landmark_matches = match_landmarks(landmarks1, landmarks2,
            m_descriptors_tradeoff, m_max_descriptor_dist, m_max_location_dist);

        // Merge the current average sample to the new sample. This updates the set of landmarks,
        // adds new landmarks from the new sample, and possibly eliminates some old non-matched
        // landmarks (outliers).
        m_average_sample->merge_data(sample, landmark_matches, m_max_countdown_value);
    }
}


std::ostream&
average_landmark_atlas::
print(std::ostream &os) const
{
    os << "{"
       << "}";

    return os;
}


void
average_landmark_atlas::
discard_outliers()
{
    #ifdef BIP_VERBOSE_MODE
        std::cout << "Discarding outliers (least frequent landmarks)";
    #endif

    // Nothing to do if no sample was computed.
    if (m_average_sample == NULL)
        return;

    // Eliminate "outliers": landmarks whose relative frequency are lower than the frequency
    // threshold. So only the most relevant landmarks detected in the population are kept.
    size_t i = 0;
    while (i < m_average_sample->get_landmarks().size()) {
        float rel_frequency = (float) m_average_sample->get_frequencies()[i] / m_num_samples;

        if (rel_frequency < m_rel_frequency_threshold)
            m_average_sample->erase_data(i);
        else
            ++i;
    }

    #ifdef BIP_VERBOSE_MODE
        std::cout << " - done\n";
    #endif
}


// ------------------------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const average_landmark_atlas &ala)
{
    return ala.print(os);
}


}
