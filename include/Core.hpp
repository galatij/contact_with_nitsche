#ifndef _CORE_HPP_
#define _CORE_HPP_
#define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG

// Level Set and Xfem stuff:
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_config.h>
#include <getfem/getfem_assembling.h> // import assembly methods (and comp. of 
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>     // export functions (save the solution in a file)
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_MUMPS_interface.h>
#include <gmm/gmm_superlu_interface.h>

#include <getfem/bgeot_mesh.h>

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <memory>

#include <getfem/getfem_mesh_fem_product.h>
#include <getfem/getfem_mesh_fem_global_function.h>
#include <getfem/getfem_models.h>
#include <getfem/getfem_model_solvers.h>
#include <getfem/getfem_plasticity.h>

#include "GetPot"
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <iomanip>
#include <functional>
#include <mpi.h>



namespace gf {
    
    /* some Getfem++ types that we will be using */
    using bgeot::scalar_type; ///< = double
    using bgeot::base_small_vector; ///< special class for small (dim < 16) vectors
    using bgeot::base_vector;
    using bgeot::base_node; ///< geometrical nodes (derived from base_small_vector)
    using bgeot::size_type; ///< basically std::vector<scalar_type>::size_type
    using bgeot::dim_type;

    using ScalarFunctionType = std::function<scalar_type(base_node,scalar_type)>; ///< f(\vector{x},t)
    using VectorFunctionType = std::function<base_small_vector(base_node, scalar_type)>; ///< \vector{f}(\vector{x},t)
    
    /* definition of some matrix/vector types. These ones are built
    using the predefined types in Gmm++ */
    using sparse_vector = getfem::modeling_standard_sparse_vector;
    using sparse_matrix = getfem::modeling_standard_sparse_matrix;
    using plain_vector = getfem::modeling_standard_plain_vector;

    using varnamelist = getfem::model::varnamelist;


    enum BCType {
        Dirichlet,
        Neumann,
        Mixed
    };

    enum RegionType {
        // BottomLeft = 1,
        // YmaxLeft = 2,
        // TopLeft = 3,
        // YminLeft = 4,
        // BottomRight = 5,
        // YmaxRight = 6,
        // TopRight = 7,
        // YminRight = 8,
        // Xmin = 9,
        // Xmax = 10,
        BulkLeft = 11,
        BulkRight = 12,
        Fault = 13
    };

} // namespace gf

#endif // _CORE_HPP_

