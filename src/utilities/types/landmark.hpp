/**
 * @file   landmark.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup types
 * @ingroup    types
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include <cstddef>
#include <iostream>
#include <vector>
#include "assert.hpp"
#include "triple.hpp"


namespace bip
{


class landmark_comp;


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
class landmark
{
    friend class landmark_comp;

public:
    typedef std::vector<float> descriptor;

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
    explicit landmark(triple<size_t>   location           = triple<size_t>(),
                      triple<float>    features           = triple<float>(),
                      const descriptor &local_descriptor  = descriptor(),
                      const descriptor &global_descriptor = descriptor()
    );

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    landmark(const landmark &other);

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    virtual ~landmark();

    /**
     * @brief {brief-description}
     * {detailed-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    landmark& operator=(const landmark &rhs);

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

public:
    /** Getter for {attribute-name} **/
    triple<size_t> get_location() const {
        return m_location;
    }

    /** Getter for {attribute-name} **/
    triple<float> get_features() const {
        return m_features;
    }

    /** Getter for {attribute-name} **/
    descriptor get_local_descriptor() const {
        return m_local_descriptor;
    }

    /** Getter for {attribute-name} **/
    descriptor get_global_descriptor() const {
        return m_global_descriptor;
    }

    /** Setter for {attribute-name} **/
    void set_location(triple<size_t> location) {
        m_location = location;
    }

    /** Setter for {attribute-name} **/
    void set_features(triple<float> features) {
        m_features = features;
    }

    /** Setter for {attribute-name} **/
    void set_local_descriptor(const descriptor &local_descriptor) {
        m_local_descriptor = local_descriptor;
    }

    /** Setter for {attribute-name} **/
    void set_global_descriptor(const descriptor &global_descriptor) {
        m_global_descriptor = global_descriptor;
    }

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
    size_t get_num_local_descriptor_bins() const {
        return m_local_descriptor.size();
    }

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
    size_t get_num_global_descriptor_bins() const {
        return m_global_descriptor.size();
    }

private:
    triple<size_t> m_location;          /**< {brief-description} **/
    triple<float>  m_features;          /**< {brief-description} **/
    descriptor     m_local_descriptor;  /**< {brief-description} **/
    descriptor     m_global_descriptor; /**< {brief-description} **/
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
class landmark_comp
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
    bool operator()(const landmark &l0, const landmark &l1) const;
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
std::ostream& operator<<(std::ostream &os, const landmark &l);


}

#endif
