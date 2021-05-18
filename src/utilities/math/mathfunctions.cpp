/**
 * @file   mathfunctions.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup math
 * @ingroup    math
 *
 * @copyright {copyright-information}
 *
 * {SOFTWARE-LICENSE-GOES-HERE}
 */

#include "mathfunctions.hpp"


namespace bip
{


const float EPSILON = 1E-5;


bool
fp_equal(float n1, float n2, float maxerror)
{
    // First try a simple equality test.
    if (n1 == n2)
        return true;
    // Then try a "fuzzy" comparison using the absolute error.
    return fabs(n1 - n2) < maxerror;
}


float
sqr(float n)
{
    return n * n;
}


int
sign(float n)
{
    return (n > 0.0) - (n < 0.0);
}


float
posmod(float n, float m)
{
    float quotient = n / m;
    return (quotient - floor(quotient)) * m;
}


float
rad2deg(float n)
{
    return n / M_PI * 180;
}


float
deg2rad(float n)
{
    return n / 180 * M_PI;
}


triple<float>
cart2sph(triple<float> n)
{
    const float x     = n[0];
    const float y     = n[1];
    const float z     = n[2];
    const float rho   = sqrt(sqr(x) + sqr(y) + sqr(z));
    const float phi   = atan2(y, x);
    const float theta = atan2(z, sqrt(sqr(x) + sqr(y)));

    return triple<float>(rho, phi, theta);
}


bool NONZEROF(float x){
  if (fabs(x) > 1.0e-4f)
    return true;
  else
    return false;
}

triple<float>
normalize(triple<float> n)
{
  const float x     = n[0];
  const float y     = n[1];
  const float z     = n[2];
  const float lenght = sqrt(sqr(x) + sqr(y) + sqr(z));
  return triple<float>(x/lenght, y/lenght, z/lenght);
}


triple<float>
cart2polar(triple<float> n)
{
    const float x     = n[0];
    const float y     = n[1];
    const float z     = n[2];
    const float lenght = sqrt(sqr(x) + sqr(y) + sqr(z));
    float a     = acos(z/lenght);
    float b = 0;
    const float phi   = atan2(y, x);
    if (NONZEROF(x)) // x <> 0
		{
			if (x<0)        // x<0 -> b += 180°; b in [-90°, 90°) -> b in [90°, 270°)
				b = M_PI + atan(y/x);
			else
				if (y<0)    // x>0 && y<0 -> b += 360°; b in [-90°, 0°) -> b in [270°, 360°)
					b = 2*M_PI + atan(y/x);
				else          // x> 0 && y>= 0; b in [0°, 90°)
					b = atan(y/x);
		}
		else // x == 0
		{
			if (NONZEROF(y))
				if (y < 0) b = M_PI;  // y<0  -> b=180°
				else b = M_PI_2;           // y>=0 -> b=90°
			else // x==y==0; z <> 0
				if (z<0)
				{
					a = 0;
					b = M_PI;  // x=y=0; z<0  -> a=0, b=180°
				}
				else
					a = b = 0;          // x=y=0; z>=0 -> a = b = 0°
		}

    return triple<float>(a, b, 0);
}


triple<float>
sph2cart(triple<float> n)
{
    const float rho   = n[0];
    const float phi   = n[1];
    const float theta = n[2];
    const float x     = rho * cos(theta) * cos(phi);
    const float y     = rho * cos(theta) * sin(phi);
    const float z     = rho * sin(theta);

    return triple<float>(x, y, z);
}


void
fast_fourier_transform(fftwf_complex *M, triple<size_t> sizes, bool backward)
{
    const size_t dimensions = 3;

    const char WISDOM_FILENAME_FORWARD[32]  = "wisdom_fftwf_forward.txt";
    const char WISDOM_FILENAME_BACKWARD[32] = "wisdom_fftwf_backward.txt";

    fftwf_plan plan;
    FILE *wisdom_file = NULL;

    // Open the correct wisdom file according to the desired transform.
    if (backward)
        wisdom_file = fopen(WISDOM_FILENAME_BACKWARD, "r");
    else
        wisdom_file = fopen(WISDOM_FILENAME_FORWARD, "r");

    // Import FFTW plan settings from wisdom file if possible.
    // This will save a lot of runtime in this FFT computation.
    if (wisdom_file) {
        fftwf_import_wisdom_from_file(wisdom_file);
        fclose(wisdom_file);
    }

    // Reveser the order of the sizes array (so the size in Z axis is the first element).
    int ft_sizes[3];
    for (size_t d = 0; d < dimensions; ++d)
        ft_sizes[d] = static_cast<int>(sizes[dimensions-d-1]);

    // Create the correct FFTW plan according to the desired transform.
    if (backward)
        plan = fftwf_plan_dft(dimensions, ft_sizes, M, M, FFTW_BACKWARD, FFTW_ESTIMATE);
    else
        plan = fftwf_plan_dft(dimensions, ft_sizes, M, M, FFTW_FORWARD, FFTW_ESTIMATE);

    // Export FFTW plan settings to wisdom file.
    // This will save a lot of runtime in future FFT computations.
    if (wisdom_file) {
        fftwf_export_wisdom_to_file(wisdom_file);
        fclose(wisdom_file);
    }

    // Compute the transform and clean up.
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // The normalization step is applied only in backward FFT.
    if (backward) {
        size_t total_size = 1;
        for (size_t d = 0; d < dimensions; ++d)
            total_size *= sizes[d];

        for (size_t i = 0; i < total_size; ++i) {
            M[i][0] /= total_size;
            M[i][1] /= total_size;
        }
    }
}


float
hanning(triple<size_t> n, triple<size_t> sizes)
{
    const float hx = (sizes[0] == 1) ? 1 : 0.5 * (1 - cos(2*M_PI * n[0] / (sizes[0]-1)));
    const float hy = (sizes[1] == 1) ? 1 : 0.5 * (1 - cos(2*M_PI * n[1] / (sizes[1]-1)));
    const float hz = (sizes[2] == 1) ? 1 : 0.5 * (1 - cos(2*M_PI * n[2] / (sizes[2]-1)));

    return hx * hy * hz;
}


float
gaussian(triple<float> n, triple<float> mu, float sigma)
{
    return exp(-(0.5 * sqr((n[0] - mu[0]) / sigma) +
                 0.5 * sqr((n[1] - mu[1]) / sigma) +
                 0.5 * sqr((n[2] - mu[2]) / sigma)));
}


void
eigens(double *M, double *eigenvalues, double *eigenvectors)
{
    const size_t dimensions = 3;

    // Allocate gsl data.
    gsl_matrix_view m = gsl_matrix_view_array(M, dimensions, dimensions);
    gsl_vector *evals = gsl_vector_alloc(dimensions);
    gsl_matrix *evecs = gsl_matrix_alloc(dimensions, dimensions);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(dimensions);

    // Calculate the eigenvalues and eigenvectors of M.
    gsl_eigen_symmv(&m.matrix, evals, evecs, w);
    gsl_eigen_symmv_free(w);

    // Sort in descending order by signed eigenvalues.
    gsl_eigen_symmv_sort(evals, evecs, GSL_EIGEN_SORT_VAL_DESC);

    // Each eigenvector is in a column of the resulting matrix.
    for (size_t i = 0; i < dimensions; ++i) {
        eigenvalues[i] = gsl_vector_get(evals, i);

        for (size_t j = 0; j < dimensions; ++j)
            eigenvectors[i*dimensions + j] = gsl_matrix_get(evecs, j, i);
    }
}


float
median(float *array, size_t size)
{
    debug::assert2(array != NULL);

    // Copy data to an auxiliary vector.
    std::vector<float> aux;
    aux.assign(array, array + size);

    // Find the element in the middle of the data array.
    // Notice that this function rearranges the vector. That's why we need to copy the original
    // array first.
    size_t n = aux.size() / 2;
    std::nth_element(aux.begin(), aux.begin() + n, aux.end());

    // If the data size is odd, we have the median element.
    if(aux.size() % 2)
        return aux[n];

    // Otherwise, we compute the average of the nth and (n-1)th elements.
    std::nth_element(aux.begin(), aux.begin() + n-1, aux.end());
    return 0.5 * (aux[n] + aux[n-1]);
}


bool
local_max(float *array, triple<size_t> pos, triple<size_t> sizes, size_t radius)
{
    // Get the coordinates and index of the central position.
    const size_t x0 = pos[0];
    const size_t y0 = pos[1];
    const size_t z0 = pos[2];
    const size_t i0 = x0 + sizes[0] * (y0 + sizes[1] * z0);

    const size_t zmin = (z0 < radius) ? 0 : z0 - radius;
    const size_t zmax = (z0 + radius >= sizes[2]) ? sizes[2]-1 : z0 + radius;
    const size_t ymin = (y0 < radius) ? 0 : y0 - radius;
    const size_t ymax = (y0 + radius >= sizes[1]) ? sizes[1]-1 : y0 + radius;
    const size_t xmin = (x0 < radius) ? 0 : x0 - radius;
    const size_t xmax = (x0 + radius >= sizes[0]) ? sizes[0]-1 : x0 + radius;

    for (size_t z = zmin; z <= zmax; ++z)
        for (size_t y = ymin; y <= ymax; ++y)
            for (size_t x = xmin; x <= xmax; ++x)
            {
                // Ignore the central position.
                if (x == x0 && y == y0 && z == z0)
                    continue;

                // Get the index of the current neighbor.
                size_t i = x + sizes[0] * (y + sizes[1] * z);

                // No neighbor may have a value greater than or equal to the central point's.
                if (array[i] >= array[i0])
                    return false;
            }

    return true;
}


bool
inside_region(triple<size_t> pos, triple<size_t> rstart, triple<size_t> rsizes)
{
    return pos[0] >= rstart[0] && pos[0] < rstart[0] + rsizes[0] &&
           pos[1] >= rstart[1] && pos[1] < rstart[1] + rsizes[1] &&
           pos[2] >= rstart[2] && pos[2] < rstart[2] + rsizes[2];
}


float
euclidean_dist(const std::vector<float> &v1, const std::vector<float> &v2)
{
    debug::assert2(v1.size() == v2.size());

    float result = 0.0;
    for (size_t i = 0; i < v1.size(); ++i)
        result += sqr(v1[i] - v2[i]);

    return sqrt(result);
}


float
chi_square_dist(const std::vector<float> &v1, const std::vector<float> &v2)
{
    debug::assert2(v1.size() == v2.size());

    float result = 0.0;
    for (size_t i = 0; i < v1.size(); ++i)
        result += sqr(v1[i] - v2[i]) / (v1[i] + v2[i] + EPSILON);

    return 0.5 * result;
}


void
normalize_minmax(float *array, size_t size, float new_min, float new_max)
{
    float old_min =  FLT_MAX;
    float old_max = -FLT_MAX;

    // Find minimum and maximum.
    for (size_t i = 0; i < size; ++i) {
        if (old_min > array[i])
            old_min = array[i];
        if (old_max < array[i])
            old_max = array[i];
    }

    float old_range = old_max - old_min;
    float new_range = new_max - new_min;

    // Put data in the range [new_min, new_max].
    for (size_t i = 0; i < size; ++i)
        array[i] = (array[i] - old_min) / old_range;
    for (size_t i = 0; i < size; ++i)
        array[i] = array[i] * new_range + new_min;
}


void
normalize_vnorm(std::vector<float> &v)
{
    float vnorm = 0.0;

    // Calculate the vector norm.
    for (size_t i = 0; i < v.size(); ++i)
        vnorm += sqr(v[i]);
    vnorm = sqrt(vnorm);

    // Normalize to make the norm become unitary.
    if (vnorm > 0.0) {
        for (size_t i = 0; i < v.size(); ++i)
            v[i] /= vnorm;
    }
}

void
normalize_descriptor(std::vector<float> &v)
{
    float lenght = 0.0;

    // Calculate the vector norm.
    for (size_t i = 0; i < v.size(); ++i){
        lenght += v[i];
    }

    // Normalize to make the sum unitary.
    if (lenght > 0.0) {
        for (size_t i = 0; i < v.size(); ++i)
            v[i] /= lenght;
    }
}

int
spider_web_bin(triple<float> cart, float ri, float r, triple<size_t> nbins)
{
  int rlogbase = 2;
  // num_bins = [num_bins_radius, num_bins_azimuth, num_bins_elevation]
  const triple<float> polar = cart2polar(cart);
	float a = polar[0];
  float b = polar[1];

	int secx = (int)( a*M_1_PI * nbins[2] );
  int secy = (int)( b*M_1_PI*0.5 * nbins[1] );
	int shell;

	if (rlogbase>1)  // shellindex ~ 1- r ; inner shells bigger than outer shells
		shell = (int)( nbins[0] * pow(rlogbase, nbins[0]*((ri/r)-1) ) );
	// for log-base=1 it is linear -> shells are equi-distant
	else
		shell = (int)( (ri/r) * nbins[0] );

	return (int)(secx + (secy * nbins[2]) + shell * (nbins[2]*nbins[1]));
}


triple<size_t>
log_spherical_bin(triple<float> sph, triple<size_t> nbins)
{
    float rho   = sph[0];
    float phi   = sph[1];
    float theta = sph[2];

    // Put the azimuth angles in the range [-pi, pi).
    while (phi < -M_PI)
        phi += 2*M_PI;
    while (phi >= M_PI)
        phi -= 2*M_PI;

    // Put the elevation angles in the range [-pi/2, pi/2].
    while (theta < -M_PI)
        theta += 2*M_PI;
    while (theta >= M_PI)
        theta -= 2*M_PI;
    if (theta >= M_PI_2 || theta < -M_PI_2) {
        theta = sign(theta) * M_PI - theta;
        phi  += (phi < 0.0) ? M_PI : -M_PI;
    }

    const float dphi   = 2*M_PI / nbins[1];
    const float dtheta =   M_PI / nbins[2];

    // Compute the radial and angular bin indices in the log-spherical space.
    const size_t bin_rho = (fp_equal(rho, 0.0)) ?
                            0.0 :
                            static_cast<size_t>(std::max(0.0, nbins[0] + floor(log(rho)/log(2.0))));
    const size_t bin_phi = posmod(floor(nbins[1] * (M_PI + phi) / (2*M_PI)), nbins[1]);
    const size_t bin_theta = (fp_equal(theta, M_PI_2)) ?
                              nbins[2] - 1 :
                              posmod(floor(nbins[2] * (M_PI_2 + theta) / M_PI), nbins[2]);

    return triple<size_t>(bin_rho, bin_phi, bin_theta);
}


}
