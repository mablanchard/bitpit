/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2020 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_volcartesian.hpp"

#include "unit_test_volcartesian_common.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing patch constructors.
*/
void subtest_001()
{
    {
        std::unique_ptr<VolCartesian> patch = std::unique_ptr<VolCartesian>(new VolCartesian());
        if (!patch) {
            throw std::runtime_error("Function 'VolCartesian()' failed unit test.");
        }
    }

    {
        int dimension = DIMENSION;
        std::array<double, 3> origin = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};
        std::array<double, 3> lengths = {{LENGTH_X, LENGTH_Y, LENGTH_Z}};
        std::array<int, 3> nCells = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};
        std::unique_ptr<VolCartesian> patch = std::unique_ptr<VolCartesian>(new VolCartesian(dimension, origin, lengths, nCells));
        if (!patch) {
            throw std::runtime_error("Function 'VolCartesian(int dimension, const std::array<double, 3> &origin, const std::array<double, 3> &lengths, const std::array<int, 3> &nCells)' failed unit test.");
        }
    }

    {
        int id = 0;
        int dimension = DIMENSION;
        std::array<double, 3> origin = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};
        std::array<double, 3> lengths = {{LENGTH_X, LENGTH_Y, LENGTH_Z}};
        std::array<int, 3> nCells = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};
        std::unique_ptr<VolCartesian> patch = std::unique_ptr<VolCartesian>(new VolCartesian(id, dimension, origin, lengths, nCells));
        if (!patch) {
            throw std::runtime_error("Function 'VolCartesian(int id, int dimension, const std::array<double, 3> &origin, const std::array<double, 3> &lengths, const std::array<int, 3> &nCells)' failed unit test.");
        }
    }

    {
        int dimension = DIMENSION;
        std::array<double, 3> origin = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};
        double length = LENGTH_X;
        int nCells1D = N_CELLS_X;
        std::unique_ptr<VolCartesian> patch = std::unique_ptr<VolCartesian>(new VolCartesian(dimension, origin, length, nCells1D));
        if (!patch) {
            throw std::runtime_error("Function 'VolCartesian(int dimension, const std::array<double, 3> &origin, double length, int nCells1D)' failed unit test.");
        }
    }

    {
        int id = 0;
        int dimension = DIMENSION;
        std::array<double, 3> origin = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};
        double length = LENGTH_X;
        int nCells1D = N_CELLS_X;
        std::unique_ptr<VolCartesian> patch = std::unique_ptr<VolCartesian>(new VolCartesian(id, dimension, origin, length, nCells1D));
        if (!patch) {
            throw std::runtime_error("Function 'VolCartesian(int id, int dimension, const std::array<double, 3> &origin, double length, int nCells1D)' failed unit test.");
        }
    }

    {
        int dimension = DIMENSION;
        std::array<double, 3> origin = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};
        double length = LENGTH_X;
        double dh = SPACING_X;
        std::unique_ptr<VolCartesian> patch = std::unique_ptr<VolCartesian>(new VolCartesian(dimension, origin, length, dh));
        if (!patch) {
            throw std::runtime_error("Function 'VolCartesian(int dimension, const std::array<double, 3> &origin, double length, double dh)' failed unit test.");
        }
    }

    {
        int id = 0;
        int dimension = DIMENSION;
        std::array<double, 3> origin = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};
        double length = LENGTH_X;
        double dh = SPACING_X;
        std::unique_ptr<VolCartesian> patch = std::unique_ptr<VolCartesian>(new VolCartesian(id, dimension, origin, length, dh));
        if (!patch) {
            throw std::runtime_error("Function 'VolCartesian(int id, int dimension, const std::array<double, 3> &origin, double length, double dh)' failed unit test.");
        }
    }

    {
        std::unique_ptr<VolCartesian> initialPatch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        int archiveVersion = 1;
        std::string archiveHeader = "VolCartesian - Unit Test";
        OBinaryArchive binaryWriter("unit_test", archiveVersion, archiveHeader);
        initialPatch->dump(binaryWriter.getStream());
        binaryWriter.close();

        IBinaryArchive binaryReader("unit_test");
        std::unique_ptr<VolCartesian> restoredPatch = std::unique_ptr<VolCartesian>(new VolCartesian(binaryReader.getStream()));
        binaryReader.close();
        if (!restoredPatch) {
            throw std::runtime_error("Function 'VolCartesian(std::istream &stream)' failed unit test.");
        }
    }
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    // Initialize the logger
    log::manager().initialize(log::COMBINED);

    // Run the subtests
    log::cout() << "Unit tests for VolCartesian constructors" << std::endl;

    try {
        subtest_001();
    } catch (const std::exception &exception) {
        log::cout() << exception.what() << std::endl;
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
