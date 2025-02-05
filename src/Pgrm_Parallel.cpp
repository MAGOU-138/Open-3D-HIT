#include "Pgrm_setup.hpp"
using namespace set;

void EXCHANGE_Y2Z(const Complex (&var_y)[NXHPX][NY2][NZPZ],
                  Complex (&var_z)[NXHPX][NY2PZ][NZ])
{
    // This function changes parallel direction from Y to Z

    Complex send[PZ][NZPZ][NY2PZ][NXHPX];
    Complex recv[PZ][NZPZ][NY2PZ][NXHPX];
    constexpr int all2all_size = NXHPX * NY2PZ * NZPZ;

    // Copy data from var_y to send buffer
    for (int i = 0; i < NXHPX; i++)
    {
        for (int j = 0; j < NY2PZ; j++)
        {
            for (int k = 0; k < NZPZ; k++)
            {
                for (int p = 0; p < PZ; p++)
                {
                    // global y index
                    int y_idx = p * NY2PZ + j;
                    send[p][k][j][i] = var_y[i][y_idx][k];
                }
            }
        }
    }

    // MPI send/recv
    MPI_Alltoall(send, all2all_size, MPI_DOUBLE_COMPLEX,
                 recv, all2all_size, MPI_DOUBLE_COMPLEX,
                 MPI_COMM_X);

    // recover data from recv buffer to var_z
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                for (size_t p = 0; p < PZ; p++)
                {
                    // global z index
                    int z_idx = p * NZPZ + k;
                    var_z[i][j][z_idx] = recv[p][k][j][i];
                }
            }
        }
    }
}

void EXCHANGE_Z2Y(const Complex (&var_z)[NXHPX][NY2PZ][NZ], Complex (&var_y)[NXHPX][NY2][NZPZ])
{

    // This function changes parallel direction from Z to Y

    Complex send[PZ][NZPZ][NY2PZ][NXHPX];
    Complex recv[PZ][NZPZ][NY2PZ][NXHPX];
    constexpr int all2all_size = NXHPX * NY2PZ * NZPZ;

    // recover data from recv buffer to var_z
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                for (size_t p = 0; p < PZ; p++)
                {
                    // global z index
                    int z_idx = p * NZPZ + k;
                    send[p][k][j][i] = var_z[i][j][z_idx];
                }
            }
        }
    }

    // MPI send/recv
    MPI_Alltoall(send, all2all_size, MPI_DOUBLE_COMPLEX,
                 recv, all2all_size, MPI_DOUBLE_COMPLEX,
                 MPI_COMM_X);

    // Copy data from var_y to send buffer
    for (int i = 0; i < NXHPX; i++)
    {
        for (int j = 0; j < NY2PZ; j++)
        {
            for (int k = 0; k < NZPZ; k++)
            {
                for (int p = 0; p < PZ; p++)
                {
                    // global y index
                    int y_idx = p * NY2PZ + j;
                    var_y[i][y_idx][k] = recv[p][k][j][i];
                }
            }
        }
    }
}

void EXCHANGE_Z2X(const Complex (&var_z)[NXHPX][NY2PZ][NZ2],
                  Complex (&var_x)[NXH][NY2PZ][NZ2PX])
{
    // This function changes parallel direction from Z to X
    Complex send[PX][NZ2PX][NY2PZ][NXHPX];
    Complex recv[PX][NZ2PX][NY2PZ][NXHPX];
    constexpr int all2all_size = NZ2PX * NY2PZ * NXHPX;

    // Copy data from var_z to send buffer
    for (int i = 0; i < NXHPX; i++)
    {
        for (int j = 0; j < NY2PZ; j++)
        {
            for (int k = 0; k < NZ2PX; k++)
            {
                for (int p = 0; p < PX; p++)
                {
                    // global z index
                    int z_idx = p * NZ2PX + k;
                    send[p][k][j][i] = var_z[i][j][z_idx];
                }
            }
        }
    }

    // MPI send/recv
    MPI_Alltoall(send, all2all_size, MPI_DOUBLE_COMPLEX,
                 recv, all2all_size, MPI_DOUBLE_COMPLEX,
                 MPI_COMM_Z);

    // recover data from recv buffer to var_x
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            for (size_t k = 0; k < NZ2PX; k++)
            {
                for (size_t p = 0; p < PX; p++)
                {
                    // global x index
                    int x_idx = p * NXHPX + i;
                    var_x[x_idx][j][k] = recv[p][k][j][i];
                }
            }
        }
    }
}

void EXCHANGE_X2Z(const Complex (&var_x)[NXH][NY2PZ][NZ2PX], Complex (&var_z)[NXHPX][NY2PZ][NZ2])
{
    // This function changes parallel direction from X  to Z
    Complex send[PX][NZ2PX][NY2PZ][NXHPX];
    Complex recv[PX][NZ2PX][NY2PZ][NXHPX];
    constexpr int all2all_size = NZ2PX * NY2PZ * NXHPX;

    // recover data from recv buffer to var_x
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY2PZ; j++)
        {
            for (size_t k = 0; k < NZ2PX; k++)
            {
                for (size_t p = 0; p < PX; p++)
                {
                    // global x index
                    int x_idx = p * NXHPX + i;
                    send[p][k][j][i] = var_x[x_idx][j][k];
                }
            }
        }
    }

    // MPI send/recv
    MPI_Alltoall(send, all2all_size, MPI_DOUBLE_COMPLEX,
                 recv, all2all_size, MPI_DOUBLE_COMPLEX,
                 MPI_COMM_Z);

    // Copy data from var_z to send buffer
    for (int i = 0; i < NXHPX; i++)
    {
        for (int j = 0; j < NY2PZ; j++)
        {
            for (int k = 0; k < NZ2PX; k++)
            {
                for (int p = 0; p < PX; p++)
                {
                    // global z index
                    int z_idx = p * NZ2PX + k;
                    var_z[i][j][z_idx] = recv[p][k][j][i];
                }
            }
        }
    }
}