#include "FEMManager.hpp"

bool DEBUGFEM = true;
namespace gf {

    FEMManager::FEMManager(const getfem::mesh& mesh)
        : M_mfU1(mesh, mesh.dim()), M_mfU2(mesh, mesh.dim()),
        M_mfStress1(mesh), M_mfStress2(mesh),
        M_mfRhs(mesh) {
        }

    void FEMManager::setMeshFem(const Numerics &n, const getfem::mesh& mesh){
        
        if (getfem::MPI_IS_MASTER())
            std::cout << "Setting Finite Element... ";

        // /** \DEBUG: */
        // std::cout << "FEMtypeDisp: " << n.FEMTypeDisplacement << std::endl;
        // std::cout << "FEMtypeStress: " << n.FEMTypeStress << std::endl;
        // std::cout << "FEMtypeRhs: " << n.FEMTypeRhs << std::endl;
        // /** \end debug */


        getfem::pfem pfU = getfem::fem_descriptor(n.FEMTypeDisplacement);
        getfem::pfem pfStress = getfem::fem_descriptor(n.FEMTypeStress);
        getfem::pfem pfRhs = getfem::fem_descriptor(n.FEMTypeRhs);

        M_mfU1.set_finite_element(pfU);
        M_mfU2.set_finite_element(pfU);
        M_mfStress1.set_qdim(3,3);
        M_mfStress2.set_qdim(3,3);
        M_mfRhs.set_finite_element(pfRhs);

        if (getfem::MPI_IS_MASTER())
            std::cout << "done." << std::endl;

    }

} // namespace gf