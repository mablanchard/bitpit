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

#include "bitpit_common.hpp"
#include "bitpit_volcartesian.hpp"

#include "unit_test_volcartesian_common.hpp"

using namespace bitpit;

/*!
 * Subtest 001
 *
 * Testing cell-related methods - Normal memory mode.
 */
void subtest_001()
{
    {
        const long EXPECTED_RESULT = N_CELLS;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->getCellCount();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellCount() const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        for (int d = 0; d < 3; ++d) {
            int result = patch->getCellCount(d);
            if (result != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getCellCount(int direction) const' failed unit test.");
            }
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::VOXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        ElementType result = patch->getCellType();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellType() const' failed unit test.");
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::VOXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        ElementType result = patch->getCellType(0);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellType(long id) const' failed unit test.");
        }
    }

    {
        const double EXPECTED_RESULT = CELL_VOLUME;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        double result = patch->evalCellVolume(0);
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'evalCellVolume(long id) const' failed unit test.");
        }
    }

    {
        const double EXPECTED_RESULT = CELL_SIZE;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        double result = patch->evalCellSize(0);
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'evalCellSize(long id) const' failed unit test.");
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<double, 3> result = patch->evalCellCentroid(0);
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'evalCellCentroid(long id) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        for (int d = 0; d < 3; ++d) {
            double result = patch->getCellCentroids(d)[0];
            if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getCellCentroids(int direction) const' failed unit test.");
            }
        }
    }

    {
        const long EXPECTED_RESULT = 101;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->getCellLinearId(1, 2, 3);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellLinearId(int i, int j, int k) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 101;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long id = patch->getCellLinearId({{1, 2, 3}});
        if (id != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<int, 3> cartesianId = patch->getCellCartesianId(101);
        for (int d = 0; d < 3; ++d) {
            if (cartesianId[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getCellCartesianId(long idx) const' failed unit test.");
            }
        }
    }

    {
        const bool EXPECTED_RESULT = true;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        bool result = patch->isCellCartesianIdValid({{1, 2, 3}});
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'isCellCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 102;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->locateClosestCell({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'locateClosestCell(std::array<double, 3> const &point) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{2, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<int, 3> result = patch->locateClosestCellCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'locateClosestCellCartesian(std::array<double, 3> const &point) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<long> result = patch->extractCellSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<long> result = patch->extractCellSubSet(101, 107);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractCellSubSet(int idxMin, int idxMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 2> EXPECTED_RESULT = {{102, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
        const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
        std::vector<long> result = patch->extractCellSubSet(pointMin, pointMax);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
            }
        }
    }
}

/*!
 * Subtest 002
 *
 * Testing cell-related methods - Light memory mode.
 */
void subtest_002()
{
    {
        const long EXPECTED_RESULT = N_CELLS;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long result = patch->getCellCount();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellCount() const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        for (int d = 0; d < 3; ++d) {
            int result = patch->getCellCount(d);
            if (result != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getCellCount(int direction) const' failed unit test.");
            }
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::VOXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        ElementType result = patch->getCellType();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellType() const' failed unit test.");
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::VOXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        ElementType result = patch->getCellType(0);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellType(long id) const' failed unit test.");
        }
    }

    {
        const double EXPECTED_RESULT = CELL_VOLUME;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        double result = patch->evalCellVolume(0);
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'evalCellVolume(long id) const' failed unit test.");
        }
    }

    {
        const double EXPECTED_RESULT = CELL_SIZE;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        double result = patch->evalCellSize(0);
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'evalCellSize(long id) const' failed unit test.");
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<double, 3> result = patch->evalCellCentroid(0);
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'evalCellCentroid(long id) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{0.5 * SPACING_X, 0.5 * SPACING_Y, 0.5 * SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        for (int d = 0; d < 3; ++d) {
            double result = patch->getCellCentroids(d)[0];
            if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getCellCentroids(int direction) const' failed unit test.");
            }
        }
    }

    {
        const long EXPECTED_RESULT = 101;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long result = patch->getCellLinearId(1, 2, 3);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellLinearId(int i, int j, int k) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 101;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long id = patch->getCellLinearId({{1, 2, 3}});
        if (id != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getCellLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<int, 3> cartesianId = patch->getCellCartesianId(101);
        for (int d = 0; d < 3; ++d) {
            if (cartesianId[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getCellCartesianId(long idx) const' failed unit test.");
            }
        }
    }

    {
        const bool EXPECTED_RESULT = true;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        bool result = patch->isCellCartesianIdValid({{1, 2, 3}});
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'isCellCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 102;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long result = patch->locateClosestCell({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'locateClosestCell(std::array<double, 3> const &point) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{2, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<int, 3> result = patch->locateClosestCellCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'locateClosestCellCartesian(std::array<double, 3> const &point) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::vector<long> result = patch->extractCellSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractCellSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 4> EXPECTED_RESULT = {{101, 102, 106, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::vector<long> result = patch->extractCellSubSet(101, 107);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractCellSubSet(int idxMin, int idxMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 2> EXPECTED_RESULT = {{102, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
        const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
        std::vector<long> result = patch->extractCellSubSet(pointMin, pointMax);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractCellSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
            }
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
    log::cout() << "Unit tests for VolCartesian cell-related methods" << std::endl;

    try {
        subtest_001();
        subtest_002();
    } catch (const std::exception &exception) {
        log::cout() << exception.what() << std::endl;
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
