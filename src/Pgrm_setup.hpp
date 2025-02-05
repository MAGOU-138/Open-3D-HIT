#ifndef SIMULATION_CONFIG_HPP
#define SIMULATION_CONFIG_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <hdf5.h>
#include <complex>
#include <fftw3.h>
#include <mpi.h>
#include <iostream>
#include <filesystem>
using Complex = std::complex<double>;
namespace fs = std::filesystem;

namespace set
{
    // These are what you can modify
    constexpr int NX = 48;
    constexpr int NY = 64;
    constexpr int NZ = 64;
    constexpr int PX = 4;
    constexpr int PZ = 2;
    constexpr double dt = 1.e-3;
    constexpr int Nsave = 100;
    constexpr int NT = 100000;
    constexpr double nu = 1.e-2;

    // The rest can not be modified
    constexpr int NXH = NX / 2;
    constexpr int NYH = NY / 2;
    constexpr int NZH = NZ / 2;
    constexpr int NX2 = NX * 3 / 2;
    constexpr int NY2 = NY * 3 / 2;
    constexpr int NZ2 = NZ * 3 / 2;
    // constexpr int NX2NPX = NX2 / NPX;
    constexpr int NXHPX = NXH / PX;
    constexpr int NZPZ = NZ / PZ;
    constexpr int NX2PX = NX2 / PX;
    constexpr int NY2PX = NY2 / PX;
    constexpr int NY2PZ = NY2 / PZ;
    constexpr int NZ2PX = NZ2 / PX;
    constexpr int NZ2PZ = NZ2 / PZ;

    constexpr Complex imag_unit = Complex(0.0, 1.0);

    extern double kx_local[NXHPX];
    extern double ky_local[NY];
    extern double kz_local[NZPZ];
    extern Complex vorx[NXHPX][NY][NZPZ], vorx0[NXHPX][NY][NZPZ], vorx1[NXHPX][NY][NZPZ], vorx2[NXHPX][NY][NZPZ];
    extern Complex vory[NXHPX][NY][NZPZ], vory0[NXHPX][NY][NZPZ], vory1[NXHPX][NY][NZPZ], vory2[NXHPX][NY][NZPZ];
    extern Complex vorz[NXHPX][NY][NZPZ], vorz0[NXHPX][NY][NZPZ], vorz1[NXHPX][NY][NZPZ], vorz2[NXHPX][NY][NZPZ];
    extern Complex jacx[NXHPX][NY][NZPZ], jacx0[NXHPX][NY][NZPZ], jacx1[NXHPX][NY][NZPZ], jacx2[NXHPX][NY][NZPZ];
    extern Complex jacy[NXHPX][NY][NZPZ], jacy0[NXHPX][NY][NZPZ], jacy1[NXHPX][NY][NZPZ], jacy2[NXHPX][NY][NZPZ];
    extern Complex jacz[NXHPX][NY][NZPZ], jacz0[NXHPX][NY][NZPZ], jacz1[NXHPX][NY][NZPZ], jacz2[NXHPX][NY][NZPZ];
    extern Complex u[NXHPX][NY][NZPZ];
    extern Complex v[NXHPX][NY][NZPZ];
    extern Complex w[NXHPX][NY][NZPZ];
    extern Complex P12[NXHPX][NY][NZPZ];
    extern Complex P13[NXHPX][NY][NZPZ];
    extern Complex P23[NXHPX][NY][NZPZ];

    extern int id_step;

    extern fftw_plan plan_xf, plan_xb, plan_yf, plan_yb, plan_zf, plan_zb;
    extern fftw_complex *fft_temp_x1;
    extern fftw_complex *fft_temp_x2;
    extern fftw_complex *fft_temp_y1;
    extern fftw_complex *fft_temp_y2;
    extern fftw_complex *fft_temp_z1;
    extern fftw_complex *fft_temp_z2;

    extern int my_rank;
    extern int numProcs;
    extern int PZID;
    extern int PXID;
    extern MPI_Comm MPI_COMM_X;
    extern MPI_Comm MPI_COMM_Z;

    constexpr double pi = 3.141592653589793238462643383279502884197;

    // Init wave numbers
}

using namespace set;
// program

// void checkCase();
void initCase();
void initField();
void get_nonlinear_term();
void timeIntegration();

void Fourier2Physical_dealiased(const Complex (&var_)[NXHPX][NY][NZPZ], double (&var_p)[NX2][NY2PZ][NZ2PX]);
void Physical2Fourier_dealiased(const double (&var_p)[NX2][NY2PZ][NZ2PX], Complex (&var_)[NXHPX][NY][NZPZ]);
void evaluate(const Complex (&var_send)[NXHPX][NY][NZPZ], Complex (&var_recv)[NXHPX][NY][NZPZ]);

void EXCHANGE_Y2Z(const Complex (&var_y)[NXHPX][NY2][NZPZ], Complex (&var_z)[NXHPX][NY2PZ][NZ]);
void EXCHANGE_Z2Y(const Complex (&var_z)[NXHPX][NY2PZ][NZ], Complex (&var_y)[NXHPX][NY2][NZPZ]);
void EXCHANGE_Z2X(const Complex (&var_z)[NXHPX][NY2PZ][NZ2], Complex (&var_x)[NXH][NY2PZ][NZ2PX]);
void EXCHANGE_X2Z(const Complex (&var_x)[NXH][NY2PZ][NZ2PX], Complex (&var_z)[NXHPX][NY2PZ][NZ2]);

void multiplication();

// hid_t FILE_CREATE(const std::string &IO_FILENAME);
// void WRITE_2D(hid_t file_id, const std::string &dataset_name, const double (&data)[NXHNP][NY]);
void get_velocity_from_vorticity();
void get_vorticity_from_velocity();
void save_Instantaneous(const int id_step_);

#endif