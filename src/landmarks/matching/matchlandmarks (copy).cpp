/**
 * @file   matchlandmarks.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-matching
 * @ingroup    landmark-matching
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "matchlandmarks.hpp"


namespace bip
{


matches_vector
match_landmarks(const std::vector<landmark> &landmarks1,
                const std::vector<landmark> &landmarks2,
                float                       descriptors_tradeoff,
                float                       max_descriptor_dist,
                float                       max_location_dist)
{
    /*
     * TODO:
     * Parallelize this method by dividing the computations below in multiple threads. There are
     * many places in which this can be done. For example, in the computation of the distance
     * matrices.
     *
     * Also, it would be very nice to optimize the computing time and storage requirements of this
     * function by using a less naive strategy. If you can find a more efficient algorithm that
     * computes the nearest neighbors of two sets of landmarks without violating the distance
     * constraints and being UNBIASED (same result no matter the order of the landmarks lists),
     * your changes to this code would be appreciated.
     */

    debug::assert2(descriptors_tradeoff >= 0.0 && descriptors_tradeoff <= 1.0);

    // Consider non-positive maximum distances as "infinity".
    if (max_descriptor_dist <= 0.0)
        max_descriptor_dist = FLT_MAX;
    if (max_location_dist <= 0.0)
        max_location_dist = FLT_MAX;

    // Data related to the (first and second) nearest neighbors are stored in a triple whose first
    // element is the neighbor's index (stored as a floating point number, but not problem), the
    // second one is the location distance and the third one is the final (weighted) descriptor
    // distance (see below).
    typedef triple<float> nn_data;

    #ifdef BIP_VERBOSE_MODE
        // std::cout << "Searching for matching landmarks:\n";
    #endif

    matches_vector matches;
    const size_t num_max_matches = std::min(landmarks1.size(), landmarks2.size());
    double dx, dy, dz;
    double vx, vy, vz;
    double x1, x2, y1, y2, z1, z2;

    // std::cout << "Number of max matches = " << num_max_matches << "\n";
    
    // Make sure that none of the landmark lists is empry.
    if (num_max_matches > 0)
    {
        matches.reserve(num_max_matches);

        // Lists with the nearest neighbors data of each landmark.
        std::pair<nn_data, nn_data> *landmarks1_nn =
                                    new std::pair<nn_data, nn_data>[landmarks1.size()]();
        std::pair<nn_data, nn_data> *landmarks2_nn =
                                    new std::pair<nn_data, nn_data>[landmarks2.size()]();

        #ifdef BIP_DEBUG_MODE
            std::ofstream ofs("distances.txt");
        #endif

        // We need to compute the nearest neighbors of each landmark from both the lists landmarks1
        // and landmarks2. In order to save some lines of code (without creating another function),
        // we use a 2-iterations loop and a set of pointers which either point to data related to
        // landmarks1 or landmarks2. In the end of an iteration, we swap these pointers. So the
        // first iteration will compute the nearest neighbors of each element from landmarks1, and
        // the second iteration will compute the nearest neighbors of each element from landmarks2.
        const std::vector<landmark>  *la    = &landmarks1;
        const std::vector<landmark>  *lb    = &landmarks2;
        std::pair<nn_data, nn_data> **la_nn = &landmarks1_nn;
        std::pair<nn_data, nn_data> **lb_nn = &landmarks2_nn;

        for (size_t k = 0; k < 2; ++k) {
            // Temporarily store the nearest neighbors data for the current landmark.
            std::pair<nn_data, nn_data> nn;

            // Compute distances and find the nearest neighbors of each landmark from both lists.
            for (size_t i = 0; i < (*la).size(); ++i) {
                nn.first[1]  = FLT_MAX; // Location distance.
                nn.first[2]  = FLT_MAX; // Descriptor distance.
                nn.second[1] = FLT_MAX; // Location distance.
                nn.second[2] = FLT_MAX; // Descriptor distance.

                for (size_t j = 0; j < (*lb).size(); ++j) {
                    // Get the euclidean distance between the local descriptors.
                    float local_descr_dist = euclidean_dist((*la)[i].get_local_descriptor(),
                                                            (*lb)[j].get_local_descriptor());
                    // std::cout<<"Local: "<<local_descr_dist<<"\n";
                    
                    // Get the chi-square distance between the global descriptors.
                    float global_descr_dist = chi_square_dist((*la)[i].get_global_descriptor(),
                                                              (*lb)[j].get_global_descriptor());
                    // std::cout<<"Global: "<<global_descr_dist<<"\n";

                    // The final weighted distance is a linear combination of the two distances
                    // above, weighted by the tradeoff parameter.
                    float final_descr_dist = descriptors_tradeoff * local_descr_dist +
                                             (1 - descriptors_tradeoff) * global_descr_dist;

                    // Get the distance of landmarks in terms of location.
                    triple<size_t> xa = (*la)[i].get_location();
                    triple<size_t> xb = (*lb)[j].get_location();
            
                    
                    float location_dist = 0;
                    if(((sqr(static_cast<float>(xa[0]) - static_cast<float>(xb[0]))< max_location_dist)
                        &&  (sqr(static_cast<float>(xa[1]) - static_cast<float>(xb[1]))< max_location_dist)) &&
                         (sqr(static_cast<float>(xa[2]) - static_cast<float>(xb[2]))< max_location_dist)){
                        location_dist = max_location_dist;
                    }else{
                        location_dist = max_location_dist + 10;
                     }

                    // float location_dist = sqrt(
                    //     sqr(static_cast<float>(xa[0]) - static_cast<float>(xb[0])) +
                    //     sqr(static_cast<float>(xa[1]) - static_cast<float>(xb[1])) +
                    //     sqr(static_cast<float>(xa[2]) - static_cast<float>(xb[2]))
                    // );
                    

                    // Update the nearest neighbor data of (*la)'s i-th landmark.
                    if (final_descr_dist < nn.first[2]) {
                        nn.first[0] = static_cast<float>(j); // 1st nn's index.
                        nn.first[1] = location_dist;         // 1st nn's location distance.
                        nn.first[2] = final_descr_dist;      // 1st nn's descriptor distance.
                    }
                    else if (final_descr_dist < nn.second[2]) {
                        nn.second[0] = static_cast<float>(j); // 2nd nn's index.
                        nn.second[1] = location_dist;         // 2nd nn's location distance.
                        nn.second[2] = final_descr_dist;      // 2nd nn's descriptor distance.
                    }

                    #ifdef BIP_DEBUG_MODE
                        // Since only the distances matter here, we only need to do this in the
                        // first iteration.
                        if (k == 0) {
                            ofs << location_dist     << " "
                                << local_descr_dist  << " "
                                << global_descr_dist << " "
                                << final_descr_dist  << "\n";
                        }
                    #endif
                }
                #ifdef BIP_DEBUG_MODE
                    ofs << "\n";
                #endif

                // Set the nearest neighbors data for the current landmark of the current list.
                (*la_nn)[i] = nn;
            }

            // Swap the pointers.
            std::swap<const std::vector<landmark>*>(la, lb);
            std::swap<std::pair<nn_data, nn_data>**>(la_nn, lb_nn);
        }

        #ifdef BIP_DEBUG_MODE
            ofs.close();
        #endif

        // Now find the matchings.
        for (size_t i = 0; i < landmarks1.size(); ++i) {
            // If the first nearest neighbor is too distant in terms of location, we can't take
            // it as a valid match.
            if (landmarks1_nn[i].first[1] > max_location_dist) // || landmarks1_nn[i].first[1] < 5)
                continue;

            // If the first nearest neighbor is too distant in terms of descriptor, we can't take
            // it as a valid match.
            if (landmarks1_nn[i].first[2] > max_descriptor_dist)
                continue;

            // Indices of first nearest neighbors.
            const size_t j  = static_cast<size_t>(landmarks1_nn[i].first[0]);
            const size_t ii = static_cast<size_t>(landmarks2_nn[j].first[0]);

            // So the i-th landmark of landmarks1 has the j-th landmark of landmarks2 as its
            // first nearest neighbor. If there is a reciprocity in that (i.e., the j-th landmark
            // of landmarks2 also has the the i-th landmark of landmarks1 as its first nearest
            // neighbor), then we have a valid match.
            if (i == ii) {
                matches.push_back(std::pair<size_t, size_t>(i, j));
                // std::cout<<landmarks1_nn[i].first[2]<<"\n";     
                // std::cout<<"Global landmarks Distances: "<<landmarks1_nn[i].first[2]<<"\n";                	
                // std::cout<<"Location Distances: "<<landmarks1_nn[i].first[1]<<"\n";
                #ifdef BIP_VERBOSE_MODE
                    std::cout 
                              << landmarks1[i].get_location()[0]<<","
                              << landmarks1[i].get_location()[1]<<","
                              << landmarks1[i].get_location()[2]<<","
                              << landmarks2[j].get_location()[0]<<","
                              << landmarks2[j].get_location()[1]<<","
                              << landmarks2[j].get_location()[2]<<"\n";
                #endif
            } else {
                #ifdef BIP_DEBUG_MODE
                    std::cout << "  " << matches.size()
                              << ". Match fail: "
                              << landmarks1[i].get_location()
                              << " <==> "
                              << landmarks2[j].get_location()
                              << " (location_dist = "
                              << landmarks1_nn[i].first[1]
                              << "; descriptor_dist = "
                              << landmarks1_nn[i].first[2]
                              << ")\n";
                #endif
            }

        }

        // Free memory.
        delete[] landmarks2_nn;
        delete[] landmarks1_nn;
    }

    #ifdef BIP_VERBOSE_MODE
        // std::cout << "Number of matches: " << matches.size() << "\n";
    #endif

    return matches;
}


}
