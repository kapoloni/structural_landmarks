/**
 * @file   landmark.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup types
 * @ingroup    types
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "landmark.hpp"


namespace bip
{


landmark::
landmark(triple<size_t>             location,
         triple<float>              features,
         const landmark::descriptor &local_descriptor,
         const landmark::descriptor &global_descriptor) :
    m_location(location),
    m_features(features),
    m_local_descriptor(local_descriptor),
    m_global_descriptor(global_descriptor)
{
    // Nothing.
}


landmark::
landmark(const landmark &other) :
    m_location(other.m_location),
    m_features(other.m_features),
    m_local_descriptor(other.m_local_descriptor),
    m_global_descriptor(other.m_global_descriptor)
{
    // Nothing.
}


landmark::
~landmark()
{
    // Nothing.
}


landmark&
landmark::
operator=(const landmark &rhs)
{
    if (this != &rhs) {
        m_location          = rhs.m_location;
        m_features          = rhs.m_features;
        m_local_descriptor  = rhs.m_local_descriptor;
        m_global_descriptor = rhs.m_global_descriptor;
    }
    return *this;
}


std::ostream&
landmark::
print(std::ostream &os) const
{
    size_t i;

    os << "{"
       << "m_location: " << m_location << ", "
       << "m_features: " << m_features << ", "
       << "m_local_descriptor: [";

    if (!m_local_descriptor.empty()) {
        for (i = 0; i < m_local_descriptor.size() - 1; ++i) {
           os << m_local_descriptor[i] << ", ";
        }
        os << m_local_descriptor[i];
    }
    os << "], "
       << "m_global_descriptor: [";
    if (!m_global_descriptor.empty()) {
        for (i = 0; i < m_global_descriptor.size() - 1; ++i) {
           os << m_global_descriptor[i] << ", ";
        }
        os << m_global_descriptor[i];
    }
    os << "]}";

    return os;
}


// ------------------------------------------------------------------------------------------------

bool
landmark_comp::
operator()(const landmark &l0, const landmark &l1) const
{
    return l0.m_features[0] < l1.m_features[0];
}


// ------------------------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const landmark &l)
{
    return l.print(os);
}


}
