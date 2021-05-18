/**
 * @file   triple.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup types
 * @ingroup    types
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef TRIPLE_HPP
#define TRIPLE_HPP

#include <cstddef>
#include <cstdarg>
#include <iostream>
#include "assert.hpp"


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
template <class T>
class triple
{
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
    explicit triple(T a0 = 0, T a1 = 0, T a2 = 0);

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    triple(const triple<T> &other);

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    virtual ~triple();

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    triple<T>& operator=(const triple<T> &rhs);

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
    T& operator[](size_t idx);

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
    const T& operator[](size_t idx) const;

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


    triple<T> operator + ( triple<T> const & p) const
    { 
        return triple<T>( m_data[0]+p.m_data[0], m_data[1]+p.m_data[1], m_data[2]+p.m_data[2] );
    }

    triple<T> operator - ( triple<T> const & p) const
    {
        return triple<T>( m_data[0]-p.m_data[0], m_data[1]-p.m_data[1], m_data[2]-p.m_data[2] );
    }

    triple<T> operator * ( const float s ) const
    {
        return triple<T>( m_data[0]*s, m_data[1]*s, m_data[2]*s );
    }

    triple<T> operator / ( const float s ) const
    {
        return triple<T>( m_data[0]/s, m_data[1]/s, m_data[2]/s );
    }

    // dot product
    triple<T> operator * ( triple<T> const & p ) const
    {
        return ( m_data[0]*p.m_data[0] + m_data[1]*p.m_data[1] + m_data[2]*p.m_data[2] );
    }  

    // cross product
    triple<T> operator ^ ( triple<T> const & p ) const
    {
        return triple<T>
        (
            m_data[1]*p.m_data[2] - m_data[2]*p.m_data[1],
            m_data[2]*p.m_data[0] - m_data[0]*p.m_data[2],
            m_data[0]*p.m_data[1] - m_data[1]*p.m_data[0]
        );
    }

    triple<T> & operator += ( triple<T> const & p)
    {
        m_data[0] += p.m_data[0];
        m_data[1] += p.m_data[1];
        m_data[2] += p.m_data[2];
        return *this;
    }

    triple<T> & operator -= ( triple<T> const & p)
    {
        m_data[0] -= p.m_data[0];
        m_data[1] -= p.m_data[1];
        m_data[2] -= p.m_data[2];
        return *this;
    }

    triple<T> & operator *= ( const float s )
    {
        m_data[0] *= s;
        m_data[1] *= s;
        m_data[2] *= s;
        return *this;
    }

    triple<T> & operator /= ( const float s )
    {
        m_data[0] /= s;
        m_data[1] /= s;
        m_data[2] /= s;
        return *this;
    }
    
    inline float Norm() const
    {
        return sqrt( m_data[0]*m_data[0] + m_data[1]*m_data[1] + m_data[2]*m_data[2] );
    }
    
    inline float SquaredNorm() const
    {
        return ( m_data[0]*m_data[0] + m_data[1]*m_data[1] + m_data[2]*m_data[2] );
    }

    triple<T> Normalize() const
    {
        float norm = sqrt( m_data[0]*m_data[0] + m_data[1]*m_data[1] + m_data[2]*m_data[2]);
        return triple<T>( m_data[0]/norm, m_data[1]/norm, m_data[2]/norm );
    }    

private:
    T m_data[3]; /**< {brief-description} **/
};


// ------------------------------------------------------------------------------------------------

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
template <class T>
class triple_comp
{
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
    triple_comp(size_t idx);

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
    bool operator()(const triple<T> &t0, const triple<T> &t1) const;

private:
    size_t m_idx; /**< {brief-description} **/
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
template <class T>
std::ostream& operator<<(std::ostream &os, const triple<T> &t);


}

#include "triple.cpp"

#endif
