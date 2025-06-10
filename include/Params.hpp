#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include "Core.hpp"
#include "GetPot"

namespace gf {

    struct Domain {
        size_type dim;
        scalar_type Lx;
        scalar_type Ly;
        scalar_type Lz;
        size_type Nx;
        size_type Ny;
        size_type Nz;
        scalar_type h;
        scalar_type angle;
        std::string meshType;
    };
    struct Physics {
        scalar_type M_E0;
        scalar_type M_nu;
        scalar_type M_lambda;
        scalar_type M_mu;
        base_small_vector M_gravity;
        scalar_type M_mu_friction;
    };
    struct It {
        scalar_type maxIt;
        scalar_type rtol;
        scalar_type atol;
        /** \todo:
         * add other Newton Params, if needed
         */
    };
    struct Nitsche {
        scalar_type theta;
        scalar_type gamma0;
    };
    struct Numerics {
        std::string integration;
        std::string FEMTypeDisplacement;
        std::string FEMTypeStress;
        std::string FEMTypeRhs;
    }; 
    struct Time {
        scalar_type t0;
        scalar_type tend;
        scalar_type dt;
        /** \todo:
         * add backtracking parameters, if needed
         * */
    };
    
    struct BCStrings {
        std::vector<size_type> regionsDirID;
        std::vector<size_type> regionsNeuID;
        std::vector<size_type> regionsMixID;
        std::vector<std::string> stringsDir;
        std::vector<std::string> stringsNeu;
        std::vector<std::string> stringsMix;
    };

    struct Params {
        GetPot datafile;
        Domain domain;
        Physics physics;
        It it;
        Nitsche nitsche;
        Time time;
        Numerics numerics;
        BCStrings bc;
        bool verbose = false;
        bool gmsh = false;
        
        Params(int argc, char* argv[]);

        void broadcast(MPI_Comm comm);
        
    };

    
    std::ostream& operator<<(std::ostream&, const Params&);

} // namespace gf


#endif // _PARAMS_HPP_