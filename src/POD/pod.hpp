/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#ifndef __BITPIT_POD_HPP__
#define __BITPIT_POD_HPP__

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "pod_kernel.hpp"
#include "pod_voloctree.hpp"

namespace bitpit {

class PODKernel;

class POD : public VTKBaseStreamer {

public:

    enum MemoryMode {
        MEMORY_NORMAL,
        MEMORY_LIGHT
    };

    enum RunMode {
        RESTORE,
        COMPUTE
    };

    enum WriteMode {
        DUMP,
        DEBUG,
        NONE
    };

    enum ReconstructionMode {
        PROJECTION,
        MINIMIZATION
    };

    enum MeshType {
        UNDEFINED,
        VOLOCTREE
    };

public:
# if BITPIT_ENABLE_MPI
    POD(MPI_Comm comm = MPI_COMM_WORLD);
# else
    POD();
# endif

    ~POD();

    POD(POD&& other) = default;

    void clear();

    void setDirectory(const std::string &directory);
    const std::string & getDirectory();
    void setName(const std::string &name);
    const std::string & getName();
    void addSnapshot(const std::string &directory, const std::string &name);
    void addSnapshot(const pod::SnapshotFile &file);
    void setSnapshots(const std::vector<pod::SnapshotFile> &database);
    void addReconstructionSnapshot(const std::string &directory, const std::string &name);
    void addReconstructionSnapshot(const pod::SnapshotFile &file);
    void setModeCount(std::size_t nmodes);
    std::size_t getModeCount();
    void setEnergyLevel(double energy);
    double getEnergyLevel();
    void setMeshType(MeshType type);
    MeshType getMeshType();
    void setStaticMesh(bool flag);

    void setMemoryMode(MemoryMode mode);
    MemoryMode getMemoryMode();
    void setRunMode(RunMode mode);
    RunMode getRunMode();
    void setWriteMode(WriteMode mode);
    WriteMode getWriteMode();
    void setReconstructionMode(ReconstructionMode mode);
    ReconstructionMode getReconstructionMode();

    void setSensorMask(const PiercedStorage<bool> &mask);

    std::size_t getSnapshotCount();
    std::vector<std::string> getScalarNames();
    std::vector<std::array<std::string,3>> getVectorNames();
    std::vector<std::string> getFieldsNames();

    const VolumeKernel* getMesh();
    const pod::PODMode & getMean();
    const std::vector<pod::PODMode> & getModes();
    std::vector<std::vector<double> > getReconstructionCoeffs();
    const std::vector<long int> & getListID();
    std::size_t getListIDInternalCount();

    void run();
    void dump();
    void restore();

    void evalMeanMesh();
    void fillListID(const PiercedStorage<bool> &bfield);
    void evalCorrelation();
    void evalModes();
    void evalEigen();
    void evalReconstruction();

    void reconstructFields(pod::PODField &field, pod::PODField &recon);
    void dumpField(const std::string &name, const pod::PODField &field) const;

    void reconstructFields(PiercedStorage<double> &fields, const VolumeKernel *mesh,
            const std::map<std::string, std::size_t> &targetScalars,
            const std::map<std::array<std::string, 3>, std::array<std::size_t, 3>> &targetVector,
            const std::unordered_set<long> *targetCells);
    void reconstructFields(PiercedStorage<double> &fields, const VolumeKernel *mesh,
            std::map<std::string, std::size_t> targetFields,
            const std::unordered_set<long> *targetCells);

protected:
    std::unique_ptr<PODKernel>              m_podkernel;                /**< POD computational kernel */
    MeshType                                m_meshType;                 /**< Type of POD mesh*/
    bool                                    m_staticMesh;               /**< If true the mesh is unique and the same for each snapshot and for POD modes [it is read one time together with the first snapshot].*/
    std::string                             m_directory;                /**< Input/output directory.*/
    std::string                             m_name;                     /**< POD session name.*/
    std::vector<pod::SnapshotFile>          m_database;                 /**< Vector of snapshots (directory and file name structure) */
    std::vector<pod::SnapshotFile>          m_databaseReconstruction;   /**< Vector of snapshots to be reconstructed (directory and file name structure) */
    std::size_t                             m_nSnapshots;               /**< Number of snapshots*/
    std::size_t                             m_nReconstructionSnapshots; /**< Number of snapshots to be reconstructed*/
    std::size_t                             m_nScalarFields;            /**< Number of scalar fields (note. first fields in dumped file)*/
    std::size_t                             m_nVectorFields;            /**< Number of vector fields (note. last fields in dumped file)*/
    std::size_t                             m_nFields;                  /**< Number of total fields*/
    std::vector<std::string>                m_nameScalarFields;         /**< Names of scalar fields. */
    std::vector<std::array<std::string,3>>  m_nameVectorFields;         /**< Names of vector fields. */
    bool                                    m_toUpdate;                 /**< If true the pod structures need to be updated.*/
    PiercedStorage<bool>                    m_filter;                   /**< Filter field (!=0 fluid cell, ==0 solid cell) used to compute POD modes (no POD on solid cells).*/
    PiercedStorage<bool>                    m_sensorMask;               /**< Sensor mask field (!=0 solve cell, ==0 no-solve cell) used to project (orthogonally and non-orthogonally) on POD modes.*/
    pod::PODMode                            m_mean;                     /**< Mean field of the snapshots database.*/
    std::vector<pod::PODMode>               m_modes;                    /**< POD Modes*/
    std::size_t                             m_nModes;                   /**< Number of retained POD modes*/
    double                                  m_energyLevel;              /**< Level of percentage energy of the retained POD modes*/

    std::vector<std::vector<double>>               m_correlationMatrices;   /**< Correlation matrices (internal use)*/
    std::vector<std::vector<double>>               m_minimizationMatrices;  /**< Least-squares minimization matrices (internal use)*/
    std::vector<std::vector<double>>               m_lambda;                /**< Eigenvalue of correlation matrix. */
    std::vector<std::vector<std::vector<double>>>  m_podCoeffs;             /**< Eigenvectors of correlation matrix (i.e. pod coefficients of database snapshots).*/
    std::vector<std::vector<double>>               m_reconstructionCoeffs;  /**< Pod coefficients of last reconstructed snapshot.*/

    std::vector<long int>                          m_listActiveIDs;           /**<List of ID of active cells [to be updated when filter/mask change]. */
    std::size_t                                    m_sizeInternal;            /**<Number of internal cells in the list of ID of active cells [the internal cells are placed first in the list of active IDs].*/

#if BITPIT_ENABLE_MPI
    MPI_Comm            m_communicator; /**< MPI communicator */
#endif
    int                 m_rank;         /**< Local rank of process. */
    int                 m_nProcs;       /**< Number of processes. */

    //pod options
    MemoryMode          m_memoryMode;           /**<Memory mode: MEMORY_NORMAL - pod modes always in memory, MEMORY_LIGHT - pod modes read from file. */
    RunMode             m_runMode;              /**<Restore or compute pod modes, mean field and pod mesh. */
    WriteMode           m_writeMode;            /**<Write mode: dump write pod info, modes, mean field and pod mesh on dump files only, DEBUG write even vtu files and NONE to dump/write nothing. [Default = DUMP] */
    ReconstructionMode  m_reconstructionMode;   /**<Evaluate reconstruction by PROJECTION or by MINIMIZATION. [Default = MINIMIZATION] */

    std::vector<std::size_t>    _m_nr;  /**<Temporary number of modes to track the energy level of retained number of modes for different fields.*/

    const static int    ARCHIVE_VERSION = 0;

    void evalMeanStaticMesh();
    void checkModeCount(double *alambda, std::size_t ifield);
    void evalModesStaticMesh();
    void initCorrelation();
    void evalCorrelationTerm(int i, pod::PODField &snapi, int j, pod::PODField &snapj);
    void evalReconstructionCoeffs(pod::PODField &snapi);
    void evalReconstructionCoeffsStaticMesh(pod::PODField &snapi);
    void buildFields(pod::PODField &recon);
    void evalMinimizationMatrices();
    void initMinimization();
    void solveMinimization(std::vector<std::vector<double>> &rhs);

    void dumpMode(std::size_t ir);

    void readSnapshot(pod::SnapshotFile snap, VolumeKernel *mesh, pod::PODField &fieldr);
    void readSnapshot(pod::SnapshotFile snap, pod::PODField &fieldr);
    void readMode(std::size_t ir);

    void diff(pod::PODField &a, const pod::PODMode &b);
    void sum(pod::PODField &a, const pod::PODMode &b);

#if BITPIT_ENABLE_MPI
    void initializeCommunicator(MPI_Comm communicator);
    MPI_Comm getCommunicator() const;
    bool isCommunicatorSet() const;
    void freeCommunicator();
#endif

    void evalReconstructionCoeffs(PiercedStorage<double> &fields, const VolumeKernel *mesh,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::size_t> &podscalarIds,
            const std::vector<std::array<std::size_t, 3>> &vectorIds, const std::vector<std::size_t> &podvectorIds);
    void evalReconstructionCoeffsStaticMesh(PiercedStorage<double> &fields,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::size_t> &podscalarIds,
            const std::vector<std::array<std::size_t, 3>> &vectorIds, const std::vector<std::size_t> &podvectorIds);
    void buildFields(PiercedStorage<double> &fields,
            const std::vector<std::size_t> &scalarIds, const std::vector<std::size_t> &podscalarIds,
            const std::vector<std::array<std::size_t, 3>> &vectorIds, const std::vector<std::size_t> &podvectorIds,
            const std::unordered_set<long> *targetCells = nullptr);

};

}

#endif