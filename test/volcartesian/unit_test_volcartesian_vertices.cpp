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
 * Testing vertex-related methods - Normal memory mode.
 */
void subtest_001()
{
    {
        const long EXPECTED_RESULT = N_VERTICES;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->getVertexCount();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getVertexCount() const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{N_VERTICES_X, N_VERTICES_Y, N_VERTICES_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        for (int d = 0; d < 3; ++d) {
            int result = patch->getVertexCount(d);
            if (result != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCount(int direction) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{1 * SPACING_X, 2 * SPACING_Y, 3 * SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<double, 3> result = patch->evalVertexCoords(139);
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'evalVertexCoords(long id) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        for (int d = 0; d < 3; ++d) {
            double result = patch->getVertexCoords(d)[1];
            if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getVertexCoords(int direction) const' failed unit test.");
            }
        }
    }

    {
        const long EXPECTED_RESULT = 139;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->getVertexLinearId(1, 2, 3);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getVertexLinearId(int i, int j, int k) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 139;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->getVertexLinearId({{1, 2, 3}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getVertexLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<int, 3> result = patch->getVertexCartesianId(139);
        for (int d = 0; d < 3; ++d) {
            if (result[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCartesianId(long idx) const' failed unit test.");
            }
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<int, 3> result = patch->getVertexCartesianId(101, 0);
        for (int d = 0; d < 3; ++d) {
            if (result[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCartesianId(long cellIdx, int vertex) const' failed unit test.");
            }
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<int, 3> result = patch->getVertexCartesianId({{1, 2, 3}}, 0);
        for (int d = 0; d < 3; ++d) {
            if (result[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCartesianId(const std::array<int, 3> &cellIjk, int vertex) const' failed unit test.");
            }
        }
    }

    {
        const bool EXPECTED_RESULT = true;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        bool result = patch->isVertexCartesianIdValid({{1, 2, 3}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'isVertexCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 146;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->locateClosestVertex({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'locateClosestVertex(std::array<double, 3> const &point) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{2, 3, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<int, 3> result = patch->locateClosestVertexCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'locateClosestVertexCartesian(std::array<double, 3> const &point) const' failed unit test.");
            }
        }
    }


    {
        const std::array<long, 4> EXPECTED_RESULT = {{139, 140, 145, 146}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<long> result = patch->extractVertexSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 2> EXPECTED_RESULT = {{101, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<long> result = patch->extractVertexSubSet(101, 107);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractVertexSubSet(int idxMin, int idxMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 2> EXPECTED_RESULT = {{140, 146}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
        const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
        std::vector<long> result = patch->extractVertexSubSet(pointMin, pointMax);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
            }
        }
    }
}

/*!
 * Subtest 002
 *
 * Testing vertex-related methods - Light memory mode.
 */
void subtest_002()
{
    {
        const long EXPECTED_RESULT = N_VERTICES;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long result = patch->getVertexCount();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getVertexCount() const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{N_VERTICES_X, N_VERTICES_Y, N_VERTICES_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        for (int d = 0; d < 3; ++d) {
            int result = patch->getVertexCount(d);
            if (result != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCount(int direction) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{1 * SPACING_X, 2 * SPACING_Y, 3 * SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<double, 3> result = patch->evalVertexCoords(139);
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'evalVertexCoords(long id) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        for (int d = 0; d < 3; ++d) {
            double result = patch->getVertexCoords(d)[1];
            if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getVertexCoords(int direction) const' failed unit test.");
            }
        }
    }

    {
        const long EXPECTED_RESULT = 139;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long result = patch->getVertexLinearId(1, 2, 3);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getVertexLinearId(int i, int j, int k) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 139;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long result = patch->getVertexLinearId({{1, 2, 3}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getVertexLinearId(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<int, 3> result = patch->getVertexCartesianId(139);
        for (int d = 0; d < 3; ++d) {
            if (result[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCartesianId(long idx) const' failed unit test.");
            }
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<int, 3> result = patch->getVertexCartesianId(101, 0);
        for (int d = 0; d < 3; ++d) {
            if (result[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCartesianId(long cellIdx, int vertex) const' failed unit test.");
            }
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{1, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<int, 3> result = patch->getVertexCartesianId({{1, 2, 3}}, 0);
        for (int d = 0; d < 3; ++d) {
            if (result[d] != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'getVertexCartesianId(const std::array<int, 3> &cellIjk, int vertex) const' failed unit test.");
            }
        }
    }

    {
        const bool EXPECTED_RESULT = true;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        bool result = patch->isVertexCartesianIdValid({{1, 2, 3}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'isVertexCartesianIdValid(const std::array<int, 3> &ijk) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 146;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        long result = patch->locateClosestVertex({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'locateClosestVertex(std::array<double, 3> const &point) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{2, 3, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::array<int, 3> result = patch->locateClosestVertexCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'locateClosestVertexCartesian(std::array<double, 3> const &point) const' failed unit test.");
            }
        }
    }


    {
        const std::array<long, 4> EXPECTED_RESULT = {{139, 140, 145, 146}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::vector<long> result = patch->extractVertexSubSet(std::array<int, 3>{{1, 2, 3}}, std::array<int, 3>{{2, 3, 3}});
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractVertexSubSet(std::array<int, 3> const &ijkMin, std::array<int, 3> const &ijkMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 2> EXPECTED_RESULT = {{101, 107}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        std::vector<long> result = patch->extractVertexSubSet(101, 107);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractVertexSubSet(int idxMin, int idxMax) const' failed unit test.");
            }
        }
    }

    {
        const std::array<long, 2> EXPECTED_RESULT = {{140, 146}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        const std::array<double, 3> pointMin = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
        const std::array<double, 3> pointMax = {{ORIGIN_X + 0.45 * LENGTH_X + 1., ORIGIN_Y + 0.45 * LENGTH_Y + 1., ORIGIN_Z + 0.45 * LENGTH_Z + 1.}};
        std::vector<long> result = patch->extractVertexSubSet(pointMin, pointMax);
        for (std::size_t n = 0; n < EXPECTED_RESULT.size(); ++n) {
            if (!utils::DoubleFloatingEqual()(result[n], EXPECTED_RESULT[n])) {
                throw std::runtime_error("Function 'extractVertexSubSet(std::array<double, 3> const &pointMin, std::array<double, 3> const &pointMax) const' failed unit test.");
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
    log::cout() << "Unit tests for VolCartesian vertex-related methods" << std::endl;

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
