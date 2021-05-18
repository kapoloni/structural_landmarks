/**
 * @file   atlassample.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-average
 * @ingroup    landmark-atlas-average
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "atlassample.hpp"


namespace bip
{


atlas_sample::
atlas_sample(std::deque<landmark> &landmarks, size_t countdown_max) :
    m_landmarks(landmarks),
    m_frequencies(std::deque<size_t>(landmarks.size(), 1)),
    m_countdowns(std::deque<size_t>(landmarks.size(), countdown_max))
{
    // Nothing.
}


atlas_sample::
atlas_sample(const atlas_sample &other) :
    m_landmarks(other.m_landmarks),
    m_frequencies(other.m_frequencies),
    m_countdowns(other.m_countdowns)
{
    // Nothing.
}


atlas_sample::
~atlas_sample()
{
    // Nothing.
}


atlas_sample&
atlas_sample::
operator=(const atlas_sample &rhs)
{
    if (this != &rhs) {
        m_landmarks   = rhs.m_landmarks;
        m_frequencies = rhs.m_frequencies;
        m_countdowns  = rhs.m_countdowns;
    }
    return *this;
}


std::ostream&
atlas_sample::
print(std::ostream &os) const
{
    os << "{"
       << "m_landmarks: [";
    for (size_t i = 0; i < m_landmarks.size(); ++i) {
        os << m_landmarks[i] << ", ";
    }
    os << "], "
       << "m_frequencies: [";
    for (size_t i = 0; i < m_frequencies.size(); ++i) {
        os << m_frequencies[i] << ", ";
    }
    os << "], "
       << "m_countdowns: [";
    for (size_t i = 0; i < m_countdowns.size(); ++i) {
        os << m_countdowns[i] << ", ";
    }
    os << "]}";
    
    return os;
}


void
atlas_sample::
merge_data(const atlas_sample   &other,
           const matches_vector &landmark_matches,
           size_t               countdown_reset_value)
{
    const size_t num_landmarks_this  = m_landmarks.size();
    const size_t num_landmarks_other = other.m_landmarks.size();

    // These arrays are initialized with 0 (false).
    bool *matched_this  = new bool[num_landmarks_this]();
    bool *matched_other = new bool[num_landmarks_other]();

    // Modify the matched landmarks of this sample, averaging their data (location, descriptors)
    // with the data of the corresponding matched landmarks of the other sample.
    // For these matched landmarks, the frequency is the sum of their individual frequencies, and
    // the countdown is reset to its initial (max) value.
    for (size_t i = 0; i < landmark_matches.size(); ++i) {
        size_t m1 = landmark_matches[i].first;
        size_t m2 = landmark_matches[i].second;

        matched_this[m1] = matched_other[m2] = true;

        landmark       *lm_this  = &m_landmarks[m1];
        const landmark *lm_other = &other.m_landmarks[m2];

        // Average landmark locations.
        lm_this->set_location(triple<size_t>(
                                 (lm_this->get_location()[0] + lm_other->get_location()[0]) / 2,
                                 (lm_this->get_location()[1] + lm_other->get_location()[1]) / 2,
                                 (lm_this->get_location()[2] + lm_other->get_location()[2]) / 2));
        
        // Average landmark features.
        lm_this->set_features(triple<float>(
                                 (lm_this->get_features()[0] + lm_other->get_features()[0]) / 2,
                                 (lm_this->get_features()[1] + lm_other->get_features()[1]) / 2,
                                 (lm_this->get_features()[2] + lm_other->get_features()[2]) / 2));

        landmark::descriptor local_descriptor_this(lm_this->get_local_descriptor());
        landmark::descriptor local_descriptor_other(lm_other->get_local_descriptor());
        landmark::descriptor global_descriptor_this(lm_this->get_global_descriptor());
        landmark::descriptor global_descriptor_other(lm_other->get_global_descriptor());

        // Average local descriptors.
        for (size_t b = 0; b < local_descriptor_this.size(); ++b) {
            local_descriptor_this[b] = 0.5 * (local_descriptor_this[b] +
                                              local_descriptor_other[b]);
        }
        lm_this->set_local_descriptor(local_descriptor_this);

        // Average global descriptors.
        for (size_t b = 0; b < global_descriptor_this.size(); ++b) {
            global_descriptor_this[b] = 0.5 * (global_descriptor_this[b] +
                                               global_descriptor_other[b]);
        }
        lm_this->set_global_descriptor(global_descriptor_this);

        // Update the frequency and reset the countdown.
        m_frequencies[m1] += other.m_frequencies[m2];
        m_countdowns[m1] = countdown_reset_value;
    }

    // Unmatched landmarks from this sample have their countdowns decremented.
    // If such countdown is about to reach 0, then the landmark data must be eliminated.
    size_t i = 0, k = 0;
    while (i < m_landmarks.size()) {
        if (matched_this[k] || m_countdowns[i] > 1) {
            m_countdowns[i] -= !matched_this[k];
            ++i;
        }
        else
            erase_data(i);
        ++k;
    }

    // Unmatched landmarks from the other sample are all added to this one, unless their data
    // are about to be eliminated.
    for (size_t j = 0; j < other.m_landmarks.size(); ++j) {
        if (!matched_other[j] && other.m_countdowns[j] > 1) {
            m_landmarks.push_back(other.m_landmarks[j]);
            m_frequencies.push_back(other.m_frequencies[j]);
            m_countdowns.push_back(other.m_countdowns[j] - 1);
        }
    }

    delete[] matched_this;
    delete[] matched_other;
}


void
atlas_sample::
erase_data(size_t index)
{
    debug::assert2(index < m_landmarks.size());

    m_landmarks.erase(m_landmarks.begin() + index);
    m_frequencies.erase(m_frequencies.begin() + index);
    m_countdowns.erase(m_countdowns.begin() + index);
}


}
