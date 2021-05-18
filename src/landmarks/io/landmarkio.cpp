/**
 * @file   landmarkio.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup landmark-io
 * @ingroup    landmark-io
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "landmarkio.hpp"


namespace bip
{


std::vector<landmark>
read_landmarks(std::string filename)
{
    size_t num_landmarks              = 0;
    size_t num_local_descriptor_bins  = 0;
    size_t num_global_descriptor_bins = 0;

    // The landmark data are read from a simple text file.
    std::ifstream ifs(filename.c_str());

    // Read the number of landmarks, and local/global descriptor bins from first line of the file.
    ifs >> num_landmarks
        >> num_local_descriptor_bins
        >> num_global_descriptor_bins;

    // std::cout << "read landmarks" <<num_landmarks << " " << num_local_descriptor_bins << " " << num_global_descriptor_bins << "\n";

    if (ifs.fail()) throw "read_landmarks(): Failed to open the file";

    std::vector<landmark> landmarks;

    for (size_t i = 0; i < num_landmarks; ++i) {
        triple<size_t>       location;
        triple<float>        features;
        landmark::descriptor local_descriptor;
        landmark::descriptor global_descriptor;

        // Read the location, saliency and orientation id of the current landmark.
        ifs >> location[0] >> location[1] >> location[2]
            >> features[0] >> features[1] >> features[2];

        if (ifs.fail())
            throw "read_landmarks(): Failed to read from file";

        if (num_local_descriptor_bins > 0)
            local_descriptor.resize(num_local_descriptor_bins, 0.0);

        // Read the local descriptor of the current landmark.
        for (size_t b = 0; b < num_local_descriptor_bins; ++b) {
            ifs >> local_descriptor[b];
            if (ifs.fail())
                throw "read_landmarks(): Failed to read from file";
        }

        if (num_global_descriptor_bins > 0)
            global_descriptor.resize(num_global_descriptor_bins, 0.0);

        // Read the global descriptor of the current landmark.
        for (size_t b = 0; b < num_global_descriptor_bins; ++b) {
            ifs >> global_descriptor[b];
            if (ifs.fail())
                throw "read_landmarks(): Failed to read from file";
        }

        // Create the landmark object and add it to the results.
        landmarks.push_back(landmark(location, features, local_descriptor, global_descriptor));
    }
    ifs.close();

    return landmarks;
}


void
write_landmarks(std::string filename, const std::vector<landmark> &landmarks)
{
    const size_t num_local_descriptor_bins =
                (landmarks.size() == 0) ? 0 : landmarks[0].get_num_local_descriptor_bins();
    const size_t num_global_descriptor_bins =
                (landmarks.size() == 0) ? 0 : landmarks[0].get_num_global_descriptor_bins();

    // The landmark data are written into a simple text file.
    std::ofstream ofs(filename.c_str());

    // The first line of the data file contains three numbers: the total count of landmarks and
    // the number of local and global descriptor bins per landmark.
    ofs << landmarks.size()           << " "
        << num_local_descriptor_bins  << " "
        << num_global_descriptor_bins << "\n";

    if (ofs.fail())
        throw "write_landmarks(): Failed to write in file";

    // The following lines contain one landmark data each. Each landmark data is stored as a
    // simple sequence of numbers, with no special format (harder for humans to read, but easier
    // to parse later).
    for (size_t i = 0; i < landmarks.size(); ++i)
    {
        const triple<size_t> location = landmarks[i].get_location();
        const triple<float>  features = landmarks[i].get_features();

        // First write the landmark location (x,y,z) and features (saliency, orientations).
        ofs << location[0] << " " << location[1] << " " << location[2] << " "
            << features[0] << " " << features[1] << " " << features[2] << " ";

        if (ofs.fail())
            throw "write_landmarks(): Failed to write in file";

        // Then the local descriptor bins.
        landmark::descriptor local_descriptor(landmarks[i].get_local_descriptor());
        for (size_t b = 0; b < num_local_descriptor_bins; ++b) {
            ofs << local_descriptor[b] << " ";

            if (ofs.fail())
                throw "write_landmarks(): Failed to write in file";
        }

        // And then the global descriptor bins.
        landmark::descriptor global_descriptor(landmarks[i].get_global_descriptor());
        for (size_t b = 0; b < num_global_descriptor_bins; ++b) {
            ofs << global_descriptor[b] << " ";

            if (ofs.fail())
                throw "write_landmarks(): Failed to write in file";
        }
        ofs << "\n";
    }
    ofs.close();
}


}
