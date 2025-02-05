#include "Pgrm_setup.hpp"

fftw_complex *set::fft_temp_x1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NX2);
fftw_complex *set::fft_temp_x2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NX2);
fftw_complex *set::fft_temp_y1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NY2);
fftw_complex *set::fft_temp_y2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NY2);
fftw_complex *set::fft_temp_z1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NZ2);
fftw_complex *set::fft_temp_z2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NZ2);

fftw_plan set::plan_xf = fftw_plan_dft_1d(NX2, fft_temp_x1, fft_temp_x2, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan set::plan_yf = fftw_plan_dft_1d(NY2, fft_temp_y1, fft_temp_y2, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan set::plan_zf = fftw_plan_dft_1d(NZ2, fft_temp_z1, fft_temp_z2, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan set::plan_xb = fftw_plan_dft_1d(NX2, fft_temp_x1, fft_temp_x2, FFTW_BACKWARD, FFTW_ESTIMATE);
fftw_plan set::plan_yb = fftw_plan_dft_1d(NY2, fft_temp_y1, fft_temp_y2, FFTW_BACKWARD, FFTW_ESTIMATE);
fftw_plan set::plan_zb = fftw_plan_dft_1d(NZ2, fft_temp_z1, fft_temp_z2, FFTW_BACKWARD, FFTW_ESTIMATE);

double set::kx_local[NXHPX] = {};
double set::ky_local[NY] = {};
double set::kz_local[NZPZ] = {};

Complex set::u[NXHPX][NY][NZPZ] = {};
Complex set::v[NXHPX][NY][NZPZ] = {};
Complex set::w[NXHPX][NY][NZPZ] = {};

Complex set::vorx[NXHPX][NY][NZPZ] = {};
Complex set::vory[NXHPX][NY][NZPZ] = {};
Complex set::vorz[NXHPX][NY][NZPZ] = {};
Complex set::vorx0[NXHPX][NY][NZPZ] = {};
Complex set::vory0[NXHPX][NY][NZPZ] = {};
Complex set::vorz0[NXHPX][NY][NZPZ] = {};
Complex set::vorx1[NXHPX][NY][NZPZ] = {};
Complex set::vory1[NXHPX][NY][NZPZ] = {};
Complex set::vorz1[NXHPX][NY][NZPZ] = {};
Complex set::vorx2[NXHPX][NY][NZPZ] = {};
Complex set::vory2[NXHPX][NY][NZPZ] = {};
Complex set::vorz2[NXHPX][NY][NZPZ] = {};

Complex set::jacx[NXHPX][NY][NZPZ] = {};
Complex set::jacy[NXHPX][NY][NZPZ] = {};
Complex set::jacz[NXHPX][NY][NZPZ] = {};
Complex set::jacx0[NXHPX][NY][NZPZ] = {};
Complex set::jacy0[NXHPX][NY][NZPZ] = {};
Complex set::jacz0[NXHPX][NY][NZPZ] = {};
Complex set::jacx1[NXHPX][NY][NZPZ] = {};
Complex set::jacy1[NXHPX][NY][NZPZ] = {};
Complex set::jacz1[NXHPX][NY][NZPZ] = {};
Complex set::jacx2[NXHPX][NY][NZPZ] = {};
Complex set::jacy2[NXHPX][NY][NZPZ] = {};
Complex set::jacz2[NXHPX][NY][NZPZ] = {};

Complex set::P12[NXHPX][NY][NZPZ] = {};
Complex set::P13[NXHPX][NY][NZPZ] = {};
Complex set::P23[NXHPX][NY][NZPZ] = {};