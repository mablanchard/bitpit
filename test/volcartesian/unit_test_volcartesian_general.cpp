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
#include <vector>
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
 * Testing patch general methods.
 */
void subtest_001()
{
    {
        std::unique_ptr<VolCartesian> initialPatch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::unique_ptr<PatchKernel> clonedPatch = initialPatch->clone();
        if (!clonedPatch) {
            throw std::runtime_error("Function 'clone()' failed unit test.");
        }
    }

    {
        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->reset();
        if (patch->getCellCount() != 0) {
            throw std::runtime_error("Function 'reset()' failed unit test.");
        } else if (patch->getVertexCount() != 0) {
            throw std::runtime_error("Function 'reset()' failed unit test.");
        } else if (patch->getInterfaceCount() != 0) {
            throw std::runtime_error("Function 'reset()' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{N_CELLS_X, N_CELLS_Y, N_CELLS_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->setDiscretization(EXPECTED_RESULT);
        for (int d = 0; d < 3; ++d) {
            if (patch->getCellCount(d) != EXPECTED_RESULT[d]) {
                throw std::runtime_error("Function 'setDiscretization(const std::array<int, 3> &nCells)' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<double, 3> result = patch->getSpacing();
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getSpacing() const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{SPACING_X, SPACING_Y, SPACING_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        for (int d = 0; d < 3; ++d) {
            double result = patch->getSpacing(d);
            if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getSpacing(int direction) const' failed unit test.");
            }
        }
    }

    {
        const VolCartesian::MemoryMode EXPECTED_RESULT = VolCartesian::MEMORY_LIGHT;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->switchMemoryMode(VolCartesian::MEMORY_LIGHT);
        VolCartesian::MemoryMode result = patch->getMemoryMode();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'switchMemoryMode(MemoryMode mode)' failed unit test.");
        }
    }

    {
        const VolCartesian::MemoryMode EXPECTED_RESULT = VolCartesian::MEMORY_NORMAL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        VolCartesian::MemoryMode result = patch->getMemoryMode();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getMemoryMode() const' failed unit test.");
        }
    }

    {
        const bool EXPECTED_RESULT = true;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        bool result = patch->isPointInside({{ORIGIN_X + 0.5 * LENGTH_X, ORIGIN_Y + 0.5 * LENGTH_Y, ORIGIN_Z + 0.5 * LENGTH_Z}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'isPointInside(const std::array<double, 3> &point) const' failed unit test.");
        }
    }

    {
        const bool EXPECTED_RESULT = true;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        bool result = patch->isPointInside(102, {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'isPointInside(const std::array<double, 3> &point) const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 102;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        long result = patch->locatePoint({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'locatePoint(const std::array<double, 3> &point) const' failed unit test.");
        }
    }

    {
        const std::array<int, 3> EXPECTED_RESULT = {{2, 2, 3}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<int, 3> result = patch->locatePointCartesian({{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}});
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'locatePointCartesian(const std::array<double, 3> &point) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{ORIGIN_X, ORIGIN_Y, ORIGIN_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<double, 3> result = patch->getOrigin();
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getOrigin() const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{10., 10., 10.}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->setOrigin(EXPECTED_RESULT);
        std::array<double, 3> result = patch->getOrigin();
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'setOrigin(const std::array<double, 3> &origin)' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{ORIGIN_X + 10., ORIGIN_Y + 10., ORIGIN_Z + 10.}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->translate({{10., 10., 10.}});
        std::array<double, 3> result = patch->getOrigin();
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'translate(const std::array<double, 3> &translation)' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{LENGTH_X, LENGTH_Y, LENGTH_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::array<double, 3> result = patch->getLengths();
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'getLengths() const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{LENGTH_X, LENGTH_Y, LENGTH_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->setLengths(EXPECTED_RESULT);
        std::array<double, 3> result = patch->getLengths();
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'setLengths(const std::array<double, 3> &lengths)' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{10. * LENGTH_X, 10. * LENGTH_Y, 10. * LENGTH_Z}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->scale({{10., 10., 10.}}, {{10., 10., 10.}});
        std::array<double, 3> result = patch->getLengths();
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'scale(const std::array<double, 3> &scaling, const std::array<double, 3> &center)' failed unit test.");
            }
        }
    }

    {
        const double EXPECTED_RESULT = 56.5;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<double> cellData(patch->getCellCount());
        for (long n = 0; n < patch->getCellCount(); ++n) {
            cellData[n] = n;
        }

        std::vector<double> result = patch->convertToVertexData(cellData);
        if (!utils::DoubleFloatingEqual()(result[101], EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'convertToVertexData(const std::vector<double> &cellData) const' failed unit test.");
        }
    }

    {
        const double EXPECTED_RESULT = 163.5;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<double> vertexData(patch->getVertexCount());
        for (long n = 0; n < patch->getVertexCount(); ++n) {
            vertexData[n] = n;
        }

        std::vector<double> result = patch->convertToCellData(vertexData);
        if (!utils::DoubleFloatingEqual()(result[101], EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'convertToCellData(const std::vector<double> &vertexData) const' failed unit test.");
        }
    }

    {
        const std::array<double, 8> EXPECTED_RESULT = {{0.07, 0.21, 0.0175, 0.0525, 0.13, 0.39, 0.0325, 0.0975}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<double> vertexData(patch->getVertexCount());
        for (long n = 0; n < patch->getVertexCount(); ++n) {
            vertexData[n] = n;
        }

        std::vector<int> stencil;
        std::vector<double> weights;
        std::array<double, 3> point = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
        int stencilSize = patch->linearCellInterpolation(point, &stencil, &weights);
        for (int i = 0; i < stencilSize; ++i) {
            double result = weights[i];
            if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[i])) {
                throw std::runtime_error("Function 'linearCellInterpolation(std::array<double, 3> &point, std::vector<int> &stencil, std::vector<double> &weights) const' failed unit test.");
            }
        }
    }

    {
        const std::array<double, 8> EXPECTED_RESULT = {{0.19125, 0.06375, 0.44625, 0.14875, 0.03375, 0.01125, 0.07875, 0.02625}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        std::vector<double> vertexData(patch->getVertexCount());
        for (long n = 0; n < patch->getVertexCount(); ++n) {
            vertexData[n] = n;
        }

        std::vector<int> stencil;
        std::vector<double> weights;
        std::array<double, 3> point = {{ORIGIN_X + 0.45 * LENGTH_X, ORIGIN_Y + 0.45 * LENGTH_Y, ORIGIN_Z + 0.45 * LENGTH_Z}};
        int stencilSize = patch->linearVertexInterpolation(point, &stencil, &weights);
        for (int i = 0; i < stencilSize; ++i) {
            double result = weights[i];
            if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT[i])) {
                throw std::runtime_error("Function 'linearVertexInterpolation(std::array<double, 3> &point, std::vector<int> &stencil, std::vector<double> &weights) const' failed unit test.");
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
    log::cout() << "Unit tests for VolCartesian general methods" << std::endl;

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
