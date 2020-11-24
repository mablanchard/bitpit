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
 * Testing interface-related methods - Normal memory mode.
 */
void subtest_001()
{
    {
        const long EXPECTED_RESULT = N_INTERFACES;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->buildInterfaces();
        long result = patch->getInterfaceCount();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getInterfaceCount() const' failed unit test.");
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::PIXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        ElementType result = patch->getInterfaceType();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::PIXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        ElementType result = patch->getInterfaceType(0);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = 0;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->buildInterfaces();
        patch->resetInterfaces();
        long result = patch->getInterfaces().size();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'resetInterfaces()' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = N_INTERFACES;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->resetInterfaces();
        patch->buildInterfaces();
        long result = patch->getInterfaces().size();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'buildInterfaces()' failed unit test.");
        }
    }

    {
        const long EXPECTED_RESULT = N_INTERFACES;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->resetInterfaces();
        patch->updateInterfaces(patch->getCells().getIds());
        long result = patch->getInterfaces().size();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'updateInterfaces(const std::vector<long> &cellIds)' failed unit test.");
        }
    }

    {
        const double EXPECTED_RESULT = INTERFACE_AREA_X;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->buildInterfaces();
        double result = patch->evalInterfaceArea(0);
        if (!utils::DoubleFloatingEqual()(result, EXPECTED_RESULT)) {
            throw std::runtime_error("Function 'evalInterfaceArea(long id) const' failed unit test.");
        }
    }

    {
        const std::array<double, 3> EXPECTED_RESULT = {{-1., 0., 0.}};

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_NORMAL);

        patch->buildInterfaces();
        std::array<double, 3> result = patch->evalInterfaceNormal(10);
        for (int d = 0; d < 3; ++d) {
            if (!utils::DoubleFloatingEqual()(result[d], EXPECTED_RESULT[d])) {
                throw std::runtime_error("Function 'evalInterfaceNormal(long id) const' failed unit test.");
            }
        }
    }
}

/*!
 * Subtest 002
 *
 * Testing interface-related methods - Light memory mode.
 */
void subtest_002()
{
    {
        const long EXPECTED_RESULT = N_INTERFACES;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        patch->buildInterfaces();
        long result = patch->getInterfaceCount();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getInterfaceCount() const' failed unit test.");
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::PIXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        ElementType result = patch->getInterfaceType();
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
        }
    }

    {
        const ElementType EXPECTED_RESULT = ElementType::PIXEL;

        std::unique_ptr<VolCartesian> patch = createTestPatch(VolCartesian::MEMORY_LIGHT);

        ElementType result = patch->getInterfaceType(0);
        if (result != EXPECTED_RESULT) {
            throw std::runtime_error("Function 'getInterfaceType() const' failed unit test.");
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
    log::cout() << "Unit tests for VolCartesian interfaces-related methods" << std::endl;

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
