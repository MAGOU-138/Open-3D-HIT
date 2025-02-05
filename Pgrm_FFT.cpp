#include "Pgrm_setup.hpp"
#include "algorithm"
#include "cstring"
using namespace set;

void Fourier2Physical_dealiased(const Complex (&var_)[NXHPX][NY][NZPZ],
                                double (&var_p)[NX2][NY2PZ][NZ2PX])
{
    // transfrom the variable from Fourier space to physical space using 3/2 rule
    // This is the most difficult part of the code, as it involves a lot of MPI communication and I/FFT

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // IFFT along y direction
    Complex var_x1y2z1_py[NXHPX][NY2][NZPZ];
    for (size_t k = 0; k < NZPZ; k++)
    {
        for (size_t i = 0; i < NXHPX; i++)
        {
            // zero-pad
            std::memset(fft_temp_y1, 0, sizeof(Complex) * NY2);
            for (size_t j = 0; j < NYH; j++)
            {
                fft_temp_y1[j][0] = var_[i][j][k].real();
                fft_temp_y1[j][1] = var_[i][j][k].imag();

                fft_temp_y1[j + NY][0] = var_[i][j + NYH][k].real();
                fft_temp_y1[j + NY][1] = var_[i][j + NYH][k].imag();
            }
            fft_temp_y1[NY][0] = 0.0; // set Nyquist wave zero
            fft_temp_y1[NY][1] = 0.0; // set Nyquist wave zero

            fftw_execute(plan_yb); // IFFT in y direction
            for (size_t j = 0; j < NY2; j++)
            {
                var_x1y2z1_py[i][j][k] = Complex(fft_temp_y2[j][0], fft_temp_y2[j][1]);
            }
        }
    }

    // change parallel direction
    Complex var_x1y2z1_pz[NXHPX][NY2PZ][NZ];
    EXCHANGE_Y2Z(var_x1y2z1_py, var_x1y2z1_pz);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // IFFT along z direction
    Complex var_x1y2z2_pz[NXHPX][NY2PZ][NZ2];
    // zero-pad and IFFT along z direction
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            // zero-pad
            std::memset(fft_temp_z1, 0, sizeof(Complex) * NZ2);
            for (size_t k = 0; k < NZH; k++)
            {
                fft_temp_z1[k][0] = var_x1y2z1_pz[i][j][k].real();
                fft_temp_z1[k][1] = var_x1y2z1_pz[i][j][k].imag();

                fft_temp_z1[k + NZ][0] = var_x1y2z1_pz[i][j][k + NZH].real();
                fft_temp_z1[k + NZ][1] = var_x1y2z1_pz[i][j][k + NZH].imag();
            }
            fft_temp_z1[NZ][0] = 0.0; // set Nyquist wave zero
            fft_temp_z1[NZ][1] = 0.0; // set Nyquist wave zero

            fftw_execute(plan_zb); // IFFT in z direction
            for (size_t k = 0; k < NZ2; k++)
            {
                var_x1y2z2_pz[i][j][k] = Complex(fft_temp_z2[k][0], fft_temp_z2[k][1]);
            }
        }
    }

    // change paralle direction
    Complex var_x1y2z2_px[NXH][NY2PZ][NZ2PX];
    EXCHANGE_Z2X(var_x1y2z2_pz, var_x1y2z2_px);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // zero-pad and IFFT along x direction
    for (size_t j = 0; j < NY2PZ; j++)
    {
        for (size_t k = 0; k < NZ2PX; k++)
        {
            // zero-pad
            std::memset(fft_temp_x1, 0, sizeof(Complex) * NX2);
            for (size_t i = 0; i < NXH; i++)
            {
                fft_temp_x1[i][0] = var_x1y2z2_px[i][j][k].real();
                fft_temp_x1[i][1] = var_x1y2z2_px[i][j][k].imag();

                fft_temp_x1[i + NX][0] = var_x1y2z2_px[NXH - i][j][k].real();
                fft_temp_x1[i + NX][1] = -var_x1y2z2_px[NXH - i][j][k].imag(); // conjugate
            }
            fft_temp_x1[NX][0] = 0.0; // set Nyquist wave zero
            fft_temp_x1[NX][1] = 0.0; // set Nyquist wave zero

            fftw_execute(plan_xb); // IFFT in x direction
            for (size_t i = 0; i < NX2; i++)
            {
                var_p[i][j][k] = fft_temp_x2[i][0];
            }
        }
    }
}

void Physical2Fourier_dealiased(const double (&var_p)[NX2][NY2PZ][NZ2PX],
                                Complex (&var_)[NXHPX][NY][NZPZ])
{
    // transfrom the variable from physical space to Fourier space using 3/2 rule
    // This is the most difficult part of the code, as it involves a lot of MPI communication and I/FFT

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // FFT along x direction and zero-pad
    Complex var_x1y2z2_px[NXH][NY2PZ][NZ2PX];
    for (size_t k = 0; k < NZ2PX; k++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            // FFT x direction
            for (size_t i = 0; i < NX2; i++)
            {
                fft_temp_x1[i][0] = var_p[i][j][k];
                fft_temp_x1[i][1] = 0.0;
            }
            fftw_execute(plan_xf);
            // zero-pad
            for (size_t i = 0; i < NXH; i++)
            {
                var_x1y2z2_px[i][j][k] = Complex(fft_temp_x2[i][0], fft_temp_x2[i][1]) / (double)NX2;
            }
        }
    }

    // change paralle direction
    Complex var_x1y2z2_pz[NXHPX][NY2PZ][NZ2];
    EXCHANGE_X2Z(var_x1y2z2_px, var_x1y2z2_pz);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // FFT along z direction
    Complex var_x1y2z1_pz[NXHPX][NY2PZ][NZ];
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            // FFT z direction
            for (size_t k = 0; k < NZ2; k++)
            {
                fft_temp_z1[k][0] = var_x1y2z2_pz[i][j][k].real();
                fft_temp_z1[k][1] = var_x1y2z2_pz[i][j][k].imag();
            }
            fftw_execute(plan_zf);
            // zero-pad
            for (size_t k = 0; k < NZH; k++)
            {
                var_x1y2z1_pz[i][j][k] = Complex(fft_temp_z2[k][0], fft_temp_z2[k][1]) / (double)NZ2;
                var_x1y2z1_pz[i][j][k + NZH] = Complex(fft_temp_z2[k + NZ][0], fft_temp_z2[k + NZ][1]) / (double)NZ2;
            }
        }
    }

    // change paralle direction
    Complex var_x1y2z1_py[NXHPX][NY2][NZPZ];
    EXCHANGE_Z2Y(var_x1y2z1_pz, var_x1y2z1_py);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // FFT along y direction
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t k = 0; k < NZPZ; k++)
        {
            // FFT y direction
            for (size_t j = 0; j < NY2; j++)
            {
                fft_temp_y1[j][0] = var_x1y2z1_py[i][j][k].real();
                fft_temp_y1[j][1] = var_x1y2z1_py[i][j][k].imag();
            }
            fftw_execute(plan_yf);
            // zero-pad
            for (size_t j = 0; j < NYH; j++)
            {
                var_[i][j][k] = Complex(fft_temp_y2[j][0], fft_temp_y2[j][1]) / (double)NY2;
                var_[i][j + NYH][k] = Complex(fft_temp_y2[j + NY][0], fft_temp_y2[j + NY][1]) / (double)NY2;
            }
        }
    }
}