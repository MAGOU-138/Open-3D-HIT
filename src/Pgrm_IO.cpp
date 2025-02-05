#include <iomanip>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include "Pgrm_setup.hpp"
using namespace set;

// evaluate array to a temp array
void evaluate_array(const Complex (&in)[NXHPX][NY][NZPZ],
                    double (&out_real)[NXHPX][NY][NZPZ],
                    double (&out_imag)[NXHPX][NY][NZPZ])
{
    for (size_t i = 0; i < NXHPX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZPZ; k++)
            {
                out_real[i][j][k] = in[i][j][k].real();
                out_imag[i][j][k] = in[i][j][k].imag();
            }
        }
    }
}

hid_t FILE_CREATE(const std::string &io_filename)
{
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if (plist_id < 0)
    {
        return -1;
    }

    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    hid_t file_id = H5Fcreate(io_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    H5Pclose(plist_id);

    return file_id;
}

void WRITE_3D(hid_t file_id, const std::string &dataset_name, const double (&data)[NXHPX][NY][NZPZ])
{
    int dim_rank = 3;
    hsize_t dims[3] = {(hsize_t)NXH, (hsize_t)NY, (hsize_t)NZ};
    hsize_t count[3] = {(hsize_t)NXHPX,
                        (hsize_t)NY,
                        (hsize_t)NZPZ};
    hsize_t offset[3] = {(hsize_t)NXHPX * PXID, 0, (hsize_t)NZPZ * PZID};

    hid_t file_space_id = H5Screate_simple(dim_rank, dims, NULL);
    hid_t dset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE, file_space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(file_space_id);

    hid_t mem_space_id = H5Screate_simple(dim_rank, count, NULL);
    file_space_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, plist_id, data);

    H5Pclose(plist_id);
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);
    H5Dclose(dset_id);
}

void save_Instantaneous(const int id_step_)
{
    std::stringstream filename_stream;
    filename_stream << "./DATA/F-" << std::setw(8) << std::setfill('0') << id_step_ << ".H5";
    std::string filename = filename_stream.str();
    if (my_rank == 0)
    {
        std::cout << "step " << id_step << ", saving data to file:" << filename << std::endl;
    }

    // IO init
    H5open();

    hid_t file_id = FILE_CREATE(filename);
    if (file_id < 0)
    {
        return;
    }

    double var_real[NXHPX][NY][NZPZ];
    double var_imag[NXHPX][NY][NZPZ];

    // save velocity u
    evaluate_array(u, var_real, var_imag);
    WRITE_3D(file_id, "u_real", var_real);
    WRITE_3D(file_id, "u_imag", var_imag);

    // save velocity v
    evaluate_array(v, var_real, var_imag);
    WRITE_3D(file_id, "v_real", var_real);
    WRITE_3D(file_id, "v_imag", var_imag);

    // save velocity w
    evaluate_array(w, var_real, var_imag);
    WRITE_3D(file_id, "w_real", var_real);
    WRITE_3D(file_id, "w_imag", var_imag);

    // save vorticity omega_x
    evaluate_array(vorx, var_real, var_imag);
    WRITE_3D(file_id, "vorx_real", var_real);
    WRITE_3D(file_id, "vorx_imag", var_imag);

    // save vorticity omega_y
    evaluate_array(vory, var_real, var_imag);
    WRITE_3D(file_id, "vory_real", var_real);
    WRITE_3D(file_id, "vory_imag", var_imag);

    // save vorticity omega_z
    evaluate_array(vorz, var_real, var_imag);
    WRITE_3D(file_id, "vorz_real", var_real);
    WRITE_3D(file_id, "vorz_imag", var_imag);

    // save nonlinear jacx
    evaluate_array(jacx, var_real, var_imag);
    WRITE_3D(file_id, "jacx_real", var_real);
    WRITE_3D(file_id, "jacx_imag", var_imag);

    // save nonlinear jacy
    evaluate_array(jacy, var_real, var_imag);
    WRITE_3D(file_id, "jacy_real", var_real);
    WRITE_3D(file_id, "jacy_imag", var_imag);

    // save nonlinear jacz
    evaluate_array(jacz, var_real, var_imag);
    WRITE_3D(file_id, "jacz_real", var_real);
    WRITE_3D(file_id, "jacz_imag", var_imag);

    // IO finish
    H5Fclose(file_id);
    H5close();
}