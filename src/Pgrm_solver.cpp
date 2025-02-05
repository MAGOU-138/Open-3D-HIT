#include "Pgrm_setup.hpp"
using namespace std;
using namespace set;

void initCase()
{
    if (my_rank == 0)
    {
        if (!fs::exists("DATA"))
        {
            fs::create_directory("DATA");
        }
    }
}

void get_vorticity_from_velocity()
{
    // Calculate the vorticity from the velocity field
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                double kx = kx_local[i];
                double ky = ky_local[j];
                double kz = kz_local[k];

                vorx[i][j][k] = (w[i][j][k] * ky - v[i][j][k] * kz) * imag_unit;
                vory[i][j][k] = (u[i][j][k] * kz - w[i][j][k] * kx) * imag_unit;
                vorz[i][j][k] = (v[i][j][k] * kx - u[i][j][k] * ky) * imag_unit;
            }
        }
    }
}

void get_velocity_from_vorticity()
{
    // Calculate the velocity from the vorticity field
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                double kx = kx_local[i];
                double ky = ky_local[j];
                double kz = kz_local[k];
                double k2 = kx * kx + ky * ky + kz * kz;
                if (k2 == 0.0)
                {
                    u[i][j][k] = Complex(0.0, 0.0);
                    v[i][j][k] = Complex(0.0, 0.0);
                    w[i][j][k] = Complex(0.0, 0.0);
                }
                else
                {
                    Complex phi_1 = vorx[i][j][k] / k2;
                    Complex phi_2 = vory[i][j][k] / k2;
                    Complex phi_3 = vorz[i][j][k] / k2;

                    u[i][j][k] = (phi_3 * ky - phi_2 * kz) * imag_unit;
                    v[i][j][k] = (phi_1 * kz - phi_3 * kx) * imag_unit;
                    w[i][j][k] = (phi_2 * kx - phi_1 * ky) * imag_unit;
                }
            }
        }
    }
}

void evaluate(const Complex (&var_send)[NXHPX][NY][NZPZ], Complex (&var_recv)[NXHPX][NY][NZPZ])
{
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                var_recv[i][j][k] = var_send[i][j][k];
            }
        }
    }
}

void multiplication()
{
    // Calculate the multiplication of velocity and vorticity

    /////////////////////////////////////////////////////////
    // Transform the velocity and vorticity from Fourier space to physical space
    double up[NX2][NY2PZ][NZ2PX];
    double vp[NX2][NY2PZ][NZ2PX];
    double wp[NX2][NY2PZ][NZ2PX];
    double vorxp[NX2][NY2PZ][NZ2PX];
    double voryp[NX2][NY2PZ][NZ2PX];
    double vorzp[NX2][NY2PZ][NZ2PX];
    double P12p[NX2][NY2PZ][NZ2PX];
    double P13p[NX2][NY2PZ][NZ2PX];
    double P23p[NX2][NY2PZ][NZ2PX];

    // Transform the velocity field
    Fourier2Physical_dealiased(u, up);
    Fourier2Physical_dealiased(v, vp);
    Fourier2Physical_dealiased(w, wp);

    // Transform the vorticity field
    Fourier2Physical_dealiased(vorx, vorxp);
    Fourier2Physical_dealiased(vory, voryp);
    Fourier2Physical_dealiased(vorz, vorzp);

    // Multiply the velocity and vorticity fields
    for (size_t i = 0; i < NX2; i++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            for (size_t k = 0; k < NZ2PX; k++)
            {
                P12p[i][j][k] = up[i][j][k] * voryp[i][j][k] - vp[i][j][k] * vorxp[i][j][k];
                P13p[i][j][k] = up[i][j][k] * vorzp[i][j][k] - wp[i][j][k] * vorxp[i][j][k];
                P23p[i][j][k] = vp[i][j][k] * vorzp[i][j][k] - wp[i][j][k] * voryp[i][j][k];
            }
        }
    }

    // Transform the result back to Fourier space using 3/2 rule
    Physical2Fourier_dealiased(P12p, P12);
    Physical2Fourier_dealiased(P13p, P13);
    Physical2Fourier_dealiased(P23p, P23);
}

void get_nonlinear_term()
{
    /** This function calculates the nonlinear term of the equation of motion. */
    multiplication();

    // Calculate the nonlinear term
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                double kx = kx_local[i];
                double ky = ky_local[j];
                double kz = kz_local[k];
                jacx[i][j][k] = (+P12[i][j][k] * ky + P13[i][j][k] * kz) * imag_unit;
                jacy[i][j][k] = (-P12[i][j][k] * kx + P23[i][j][k] * kz) * imag_unit;
                jacz[i][j][k] = (-P13[i][j][k] * kx - P23[i][j][k] * ky) * imag_unit;
            }
        }
    }
}

void initField()
{
    // The wavenumbers
    for (size_t i = 0; i < NXHPX; i++)
    {
        kx_local[i] = (double)(PXID * NXHPX + i);
    }
    for (size_t j = 0; j < NYH; j++)
    {
        ky_local[j] = (double)j;
        ky_local[j + NYH] = (double)j - (double)NYH;
    }
    for (size_t k = 0; k < NZPZ; k++)
    {
        double kz_temp = (double)(NZPZ * PZID + k);
        if (kz_temp >= NZH)
        {
            kz_temp -= NZ;
        }
        kz_local[k] = kz_temp;
    }

    /**** Check wavenumbers for each rank***/
    // if (rank == 7)
    // {
    //     std::cout << "kx_local: ";
    //     for (size_t i = 0; i < NXHPX; i++)
    //     {
    //         std::cout << kx_local[i] << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "ky_local: ";
    //     for (size_t j = 0; j < NY; j++)
    //     {
    //         std::cout << ky_local[j] << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "kz_local: ";
    //     for (size_t k = 0; k < NZPZ; k++)
    //     {
    //         std::cout << kz_local[k] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // Initialize the velocity field by given kinetic energy spectra
    if (my_rank == 0)
    {
        std::cout << "Initializing the velocity field..." << std::endl;
    }

    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                double kx = kx_local[i];
                double ky = ky_local[j];
                double kz = kz_local[k];
                // vor x
                if (kx == 0.0 && ky == 0.0 && kz == 2.0)
                {
                    vorx[i][j][k] = Complex(0.0, -1.0);
                }
                if (kx == 0.0 && ky == 0.0 && kz == -2.0)
                {
                    vorx[i][j][k] = Complex(0.0, 1.0);
                }
                if (kx == 0.0 && ky == 3.0 && kz == 0.0)
                {
                    vorx[i][j][k] = Complex(0.0, 1.5);
                }
                if (kx == 0.0 && ky == -3.0 && kz == 0.0)
                {
                    vorx[i][j][k] = Complex(0.0, -1.5);
                }
                // vor y
                if (kx == 0.0 && ky == 0.0 && kz == 1.0)
                {
                    vory[i][j][k] = Complex(0.0, 0.5);
                }
                if (kx == 0.0 && ky == 0.0 && kz == -1.0)
                {
                    vory[i][j][k] = Complex(0.0, -0.5);
                }
                if (kx == 3.0 && ky == 0.0 && kz == 0.0)
                {
                    vory[i][j][k] = Complex(-1.5, 0.0);
                }
                // vor z
                if (kx == 0.0 && ky == 1.0 && kz == 0.0)
                {
                    vorz[i][j][k] = Complex(-0.5, 0.0);
                }
                if (kx == 0.0 && ky == -1.0 && kz == 0.0)
                {
                    vorz[i][j][k] = Complex(-0.5, 0.0);
                }
                if (kx == 2.0 && ky == 0.0 && kz == 0.0)
                {
                    vorz[i][j][k] = Complex(1.0, 0.0);
                }
            }
        }
    }

    // for (size_t i = 0; i < NXHPX; i++)
    // {
    //     for (size_t j = 0; j < NY; j++)
    //     {
    //         for (size_t k = 0; k < NZPZ; k++)
    //         {
    //             double kx = kx_local[i];
    //             double ky = ky_local[j];
    //             double kz = kz_local[k];
    //             double k2 = kx * kx + ky * ky + kz * kz;
    //             double E = 0.0;
    //             double phi_u = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * 2.0 * M_PI;
    //             double phi_v = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * 2.0 * M_PI;
    //             double phi_w = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * 2.0 * M_PI;

    //             bool valid_wave = true;
    //             if (k2 == 0.0)
    //             {
    //                 valid_wave = false;
    //             }
    //             if (kx == -NXH * 1.0 or ky == -NYH * 1.0 or kz == -NZH * 1.0)
    //             {
    //                 std::cout << "rank: " << rank << " Nyquist wave: " << kx << " " << ky << " " << kz << std::endl;
    //                 valid_wave = false;
    //             }

    //             if (valid_wave)
    //             {
    //                 double kabs = sqrt(k2);
    //                 double P11 = -kx * kx / k2 + 1.0;
    //                 double P12 = -kx * ky / k2;
    //                 double P13 = -kx * kz / k2;
    //                 double P21 = -ky * kx / k2;
    //                 double P22 = -ky * ky / k2 + 1.0;
    //                 double P23 = -ky * kz / k2;
    //                 double P31 = -kz * kx / k2;
    //                 double P32 = -kz * ky / k2;
    //                 double P33 = -kz * kz / k2 + 1.0;
    //                 E = 1.0 / (kabs * (1.0 + powf(kabs, 4) / powf(6.0, 4)));
    //                 double theta = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * M_PI;
    //                 double phi = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * 2.0 * M_PI;
    //                 double Ex = sqrt(2.0 * E) * sin(theta) * cos(phi);
    //                 double Ey = sqrt(2.0 * E) * sin(theta) * sin(phi);
    //                 double Ez = sqrt(2.0 * E) * cos(theta);

    //                 u[i][j][k] = (Ex * P11 + Ey * P12 + Ez * P13) * Complex(cos(phi_u), sin(phi_u));
    //                 v[i][j][k] = (Ex * P21 + Ey * P22 + Ez * P23) * Complex(cos(phi_v), sin(phi_v));
    //                 w[i][j][k] = (Ex * P31 + Ey * P32 + Ez * P33) * Complex(cos(phi_w), sin(phi_w));

    //                 if (rank == 0 && i == 0 && k == 0)
    //                 {
    //                     std::cout << "kx: " << kx << " ky: " << ky << " kz: " << kz << " E: " << E << " Ex: " << Ex << " Ey: " << Ey << " Ez: " << Ez << "P11: " << P11 << " v: " << v[i][j][k] << std::endl;
    //                 }
    //             }
    //         }
    //     }
    // }

    get_velocity_from_vorticity();

    // get_vorticity_from_velocity();
    get_nonlinear_term();

    // initialize the multi-steps
    evaluate(vorx, vorx0);
    evaluate(vorx, vorx1);
    evaluate(vorx, vorx2);

    evaluate(vory, vory0);
    evaluate(vory, vory1);
    evaluate(vory, vory2);

    evaluate(vorz, vorz0);
    evaluate(vorz, vorz1);
    evaluate(vorz, vorz2);

    evaluate(jacx, jacx0);
    evaluate(jacx, jacx1);
    evaluate(jacx, jacx2);

    evaluate(jacy, jacy0);
    evaluate(jacy, jacy1);
    evaluate(jacy, jacy2);

    evaluate(jacz, jacz0);
    evaluate(jacz, jacz1);
    evaluate(jacz, jacz2);
}

void timeIntegration()
{
    // Time integration using numerical splitting method
    /***
     * Karniadakis G E, Israeli M, Orszag S A. High-order splitting methods for the incompressible Navier-Stokes equations[J].
     * Journal of computational physics, 1991, 97(2): 414-443.
     ***/

    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                Complex residual_x = 3.0 * vorx0[i][j][k] - 1.5 * vorx1[i][j][k] + 1.0 / 3.0 * vorx2[i][j][k] - dt * (3.0 * jacx0[i][j][k] - 3.0 * jacx1[i][j][k] + jacx2[i][j][k]);
                Complex residual_y = 3.0 * vory0[i][j][k] - 1.5 * vory1[i][j][k] + 1.0 / 3.0 * vory2[i][j][k] - dt * (3.0 * jacy0[i][j][k] - 3.0 * jacy1[i][j][k] + jacy2[i][j][k]);
                Complex residual_z = 3.0 * vorz0[i][j][k] - 1.5 * vorz1[i][j][k] + 1.0 / 3.0 * vorz2[i][j][k] - dt * (3.0 * jacz0[i][j][k] - 3.0 * jacz1[i][j][k] + jacz2[i][j][k]);
                double k2 = kx_local[i] * kx_local[i] + ky_local[j] * ky_local[j] + kz_local[k] * kz_local[k];
                vorx[i][j][k] = residual_x / (11.0 / 6.0 + dt * (nu * k2));
                vory[i][j][k] = residual_y / (11.0 / 6.0 + dt * (nu * k2));
                vorz[i][j][k] = residual_z / (11.0 / 6.0 + dt * (nu * k2));
            }
        }
    }

    get_velocity_from_vorticity();
    get_nonlinear_term();

    evaluate(vorx1, vorx2);
    evaluate(vorx0, vorx1);
    evaluate(vorx, vorx0);

    evaluate(vory1, vory2);
    evaluate(vory0, vory1);
    evaluate(vory, vory0);

    evaluate(vorz1, vorz2);
    evaluate(vorz0, vorz1);
    evaluate(vorz, vorz0);

    evaluate(jacx1, jacx2);
    evaluate(jacx0, jacx1);
    evaluate(jacx, jacx0);

    evaluate(jacy1, jacy2);
    evaluate(jacy0, jacy1);
    evaluate(jacy, jacy0);

    evaluate(jacz1, jacz2);
    evaluate(jacz0, jacz1);
    evaluate(jacz, jacz0);
}