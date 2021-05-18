/**
 * @file   atlassample.hpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-atlas-average
 * @ingroup    landmark-atlas-average
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#ifndef ATLASSAMPLE_HPP
#define ATLASSAMPLE_HPP

#define BIP_VERBOSE_MODE
// #define BIP_DEBUG_MODE

#include <cstddef>
#include <iostream>
#include <deque>
#include <utility>
#include "assert.hpp"
#include "triple.hpp"
#include "landmark.hpp"
#include "matchlandmarks.hpp"


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
class atlas_sample
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
    atlas_sample(std::deque<landmark> &landmarks, size_t countdown_max);

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    atlas_sample(const atlas_sample &other);

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    ~atlas_sample();

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    atlas_sample& operator=(const atlas_sample &rhs);

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
    /** @brief Getter for {attribute-name} */
    const std::deque<landmark>& get_landmarks() const {
        return m_landmarks;
    }

    /** @brief Getter for {attribute-name} */
    const std::deque<size_t>& get_frequencies() const {
        return m_frequencies;
    }

    /** @brief Getter for {attribute-name} */
    const std::deque<size_t>& get_countdowns() const {
        return m_countdowns;
    }

public:
    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @param[in]  {param-name} {param-description}
     * @param[out] {param-name} {param-description}
     *
     * @returns {return-value-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    size_t get_size() const {
        return m_landmarks.size();
    }

public:
    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @param[in]  {param-name} {param-description}
     * @param[out] {param-name} {param-description}
     *
     * @returns {return-value-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    void merge_data(const atlas_sample   &other,
                    const matches_vector &landmark_matches,
                    size_t               countdown_reset_value);

    /**
     * @brief {brief-description}
     * {detailed-description}
     * 
     * @param[in]  {param-name} {param-description}
     * @param[out] {param-name} {param-description}
     *
     * @returns {return-value-description}
     *
     * @warning {warning-text}
     * @attention {attention-text}
     *
     * @see {references}
     */
    void erase_data(size_t index);

private:
    std::deque<landmark> m_landmarks;   /**< {brief-description} **/
    std::deque<size_t>   m_frequencies; /**< {brief-description} **/
    std::deque<size_t>   m_countdowns;  /**< {brief-description} **/
};


// ------------------------------------------------------------------------------------------------

/**
 * @fn {definition-name}
 *
 * @brief {brief-description}
 * {detailed-description}
 *
 * @warning {warning-text}
 * @attention {attention-text}
 *
 * @see {references}
 */
std::ostream& operator<<(std::ostream &os, const atlas_sample &as);


}

#endif
