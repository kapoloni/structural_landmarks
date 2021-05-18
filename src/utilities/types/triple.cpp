/**
 * @file   triple.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup types
 * @ingroup    types
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

namespace bip
{


template <class T>
triple<T>::
triple(T a0, T a1, T a2)
{
    m_data[0] = a0;
    m_data[1] = a1;
    m_data[2] = a2;
}


template <class T>
triple<T>::
triple(const triple<T> &other)
{
    m_data[0] = other.m_data[0];
    m_data[1] = other.m_data[1];
    m_data[2] = other.m_data[2];
}


template <class T>
triple<T>::
~triple()
{
    // Nothing.
}


template <class T>
triple<T>&
triple<T>::
operator=(const triple<T> &rhs)
{
    if (this != &rhs) {
        m_data[0] = rhs.m_data[0];
        m_data[1] = rhs.m_data[1];
        m_data[2] = rhs.m_data[2];
    }
    return *this;
}


template <class T>
T&
triple<T>::
operator[](size_t idx)
{
    bip::debug::assert2(idx < 3);
    return m_data[idx];
}


template <class T>
const T&
triple<T>::
operator[](size_t idx) const
{
    bip::debug::assert2(idx < 3);
    return m_data[idx];
}


template <class T>
std::ostream&
triple<T>::
print(std::ostream &os) const
{
    return os << "(" << m_data[0] << ", " << m_data[1] << ", " << m_data[2] << ")";
}


// ------------------------------------------------------------------------------------------------

template <class T>
triple_comp<T>::
triple_comp(size_t idx)
{
    bip::debug::assert2(idx < 3);
    m_idx = idx;
}


template <class T>
bool
triple_comp<T>::
operator()(const triple<T> &t0, const triple<T> &t1) const
{
    return t0[m_idx] < t1[m_idx];
}


// ------------------------------------------------------------------------------------------------

template <class T>
std::ostream&
operator<<(std::ostream &os, const triple<T> &t)
{
    return t.print(os);
}


}
