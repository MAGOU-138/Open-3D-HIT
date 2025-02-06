

#include <iostream>
#include <stdexcept>
#include <string>
#include "Pgrm_setup.hpp"

/*
3d turbulence
programmed by Chutian Wu
*/

namespace set
{
    int numProcs = 1;
    int my_rank = 0;
    int id_step = 0;
    int PZID = 0;
    int PXID = 0;
    MPI_Comm MPI_COMM_X = 0;
    MPI_Comm MPI_COMM_Z = 0;
}

int main(int argc, char **argv)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &set::numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &set::my_rank);

    set::PZID = set::my_rank / set::PX;
    set::PXID = set::my_rank - set::PZID * set::PX;
    int new_rank, new_size;

    MPI_Comm_split(MPI_COMM_WORLD, set::PXID, set::PZID, &set::MPI_COMM_X);
    MPI_Comm_rank(set::MPI_COMM_X, &new_rank);
    MPI_Comm_size(set::MPI_COMM_X, &new_size);

    MPI_Comm_split(MPI_COMM_WORLD, set::PZID, set::PXID, &set::MPI_COMM_Z);
    MPI_Comm_rank(set::MPI_COMM_Z, &new_rank);
    MPI_Comm_size(set::MPI_COMM_Z, &new_size);
    // std::cout << set::rank << ' ' << set::PXID << ' ' << set::PZID << std::endl;
    ///////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////
    // Init case
    // checkCase();
    initCase();

    // init field
    initField();
    ///////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////
    // time integration
    for (size_t i = 0; i < set::NT; i++)
    {
        if (set::my_rank == 0)
        {
            std::cout << "time step " << set::id_step << std::endl;
        }

        if (set::id_step % set::Nsave == 0)
        {
            save_Instantaneous(set::id_step);
        }
        timeIntegration();
        set::id_step += 1;
    }
    /////////////////////////////////////////////////////////////////////////////////////////

    // Free FFT
    fftw_destroy_plan(set::plan_xf);
    fftw_destroy_plan(set::plan_xb);
    fftw_destroy_plan(set::plan_yf);
    fftw_destroy_plan(set::plan_yb);
    fftw_free(set::fft_temp_x1);
    fftw_free(set::fft_temp_x2);
    fftw_free(set::fft_temp_y1);
    fftw_free(set::fft_temp_y2);
    // Finalize MPI
    MPI_Comm_free(&set::MPI_COMM_X);
    MPI_Comm_free(&set::MPI_COMM_Z);
    MPI_Finalize();

    return 0;
}
