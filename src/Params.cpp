#include "Params.hpp"
#include "Utils.hpp"
#include "muParserXInterface.hpp"

namespace gf {

    Params::Params(int argc, char* argv[])
    /* : datafile("data.pot")*/
    {
        if (getfem::MPI_IS_MASTER()){
            GetPot command_line(argc, argv);
            const std::string dataFileName = command_line.follow("data.pot", 2, "-f",
                    "--file");
            
            verbose = command_line.search("-v");
            gmsh = command_line.search("-m");
            
            std::ifstream check(dataFileName);
            if(!check.is_open())
                throw std::runtime_error("Could not open the datafile!");
            check.close();
            GetPot datafile("data.pot");
            
            // GetPot datafile {filename.c_str()};
            
            domain.dim = datafile("domain/dim", 3);
            domain.Lx = datafile("domain/Lx", 1.0);
            domain.Ly = datafile("domain/Ly", 1.0);
            domain.Lz = datafile("domain/Lz", 1.0);
            domain.h = datafile("domain/h", 0.2);
            domain.angle = datafile("domain/angle", 0.);
            domain.Nx = static_cast<int>(domain.Lx/domain.h);
            domain.Ny = static_cast<int>(domain.Ly/domain.h / (domain.Ly/domain.Lx));
            domain.Nz = static_cast<int>(domain.Lz/domain.h / (domain.Lz/domain.Lx));
            domain.meshType = datafile("domain/meshType", "GT_QK(3,1)");
            
            physics.M_E0 = datafile("physics/E", 0.0);
            physics.M_nu = datafile("physics/nu", 0.0);
            physics.M_mu_friction = datafile("physics/mu_friction", 0.5);
            physics.M_lambda = (physics.M_E0 * physics.M_nu) /
                                ((1 + physics.M_nu) * (1 - 2 * physics.M_nu));
            physics.M_mu = physics.M_E0 / (2 * (1 + physics.M_nu));
            // Load gravity vector
            std::string gravity_str = datafile("physics/bulkLoad", "[0., 0., 0.]"); // Note: expecting [ ... ]
            std::vector<std::string> gVecStr = splitString(gravity_str);
            for (size_type i{}; i < 3; ++i)
                physics.M_gravity.push_back(std::stod(gVecStr[i]));

            if (physics.M_gravity.size() != 3)
                throw std::runtime_error("Expected 3 components in physics/bulkLoad");

            it.maxIt = datafile("it/maxit", 30);
            it.rtol = datafile("it/tol", 1e-6);
            it.atol = datafile("it/atol", 0.); // atol: not in file, default to 0
            
            nitsche.theta = datafile("it/theta", 0.0);
            nitsche.gamma0 = datafile("it/gamma", 10.0);

            time.t0 = datafile("time/t0", 0.0);
            time.tend = datafile("time/tend", 1.0);
            time.dt = datafile("time/dt", 0.1);

            if (domain.meshType == "GT_PK(3,1)")
            {
                numerics.integration = "IM_TETRAHEDRON(5)";
                numerics.FEMTypeDisplacement = "FEM_PK(3,1)";
                numerics.FEMTypeRhs = "FEM_PK(3,1)";
                numerics.FEMTypeStress = "FEM_PK(3,1)";
            }
            else if (domain.meshType == "GT_QK(3,1)")
            {
                numerics.integration = "IM_HEXAHEDRON(5)";
                numerics.FEMTypeDisplacement = "FEM_QK(3,1)";
                numerics.FEMTypeRhs = "FEM_QK(3,1)";
                numerics.FEMTypeStress = "FEM_QK(3,1)";
                // std::cout << "Inside Params constructor: numerics.FEMTypeRhs = " << numerics.FEMTypeRhs << std::endl;
            }
            else 
                throw std::runtime_error("Select either GT_PK(3,1) or GT_QK(3,1)");

            std::string regionsDir;
            std::string regionsNeu;
            std::string regionsMix;
            regionsDir = datafile("physics/regionDisp", "");
            regionsNeu = datafile("physics/regionLoad", "");
            regionsMix = datafile("physics/regionDispNormal", "");

            bc.regionsDirID = gf::toVec(regionsDir);
            bc.regionsNeuID = gf::toVec(regionsNeu);
            bc.regionsMixID = gf::toVec(regionsMix);

            for (size_t i = 0; i < bc.regionsDirID.size(); ++i) {
                std::ostringstream varname;
                varname << "physics/bdDisp" << (i + 1); // bdDisp1, bdDisp2, ...
                std::string stringValue = datafile(varname.str().c_str(), "");
                if (stringValue.empty())
                    throw std::runtime_error("Dirichlet String Value undetected!");
                bc.stringsDir.emplace_back(stringValue);
            }

            for (size_t i = 0; i < bc.regionsNeuID.size(); ++i) {
                std::ostringstream varname;
                varname << "physics/bdLoad" << (i + 1); // bdLoad1, bdLoad2, ...
                std::string stringValue = datafile(varname.str().c_str(), "");
                if (stringValue.empty())
                    throw std::runtime_error("Neumann String Value undetected!");
                bc.stringsNeu.emplace_back(stringValue);
            }

            for (size_t i = 0; i < bc.regionsMixID.size(); ++i) {
                std::ostringstream varname;
                varname << "physics/bdDispN" << (i + 1);
                std::string stringValue = datafile(varname.str().c_str(), "");
                if (stringValue.empty())
                    throw std::runtime_error("Mixed String Value undetected!");
                bc.stringsMix.emplace_back(stringValue);
            }

            if (verbose) std::cout << *this << std::endl;
        }

        // Broadcast
        broadcast(MPI_COMM_WORLD);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 1)
            std::cout << "From rank 1: " << *this << std::endl;

    }

    std::ostream& operator<<(std::ostream& os, const Params& p){
        os << "================ PARAMETERS ================\n";
        os << "DOMAIN:\n";
        os << "-- [Lx, Ly, Lz] = ["
           << p.domain.Lx << ", "
           << p.domain.Ly << ", "
           << p.domain.Lz << "]\n";
        os << "-- [Nx, Ny, Nz] = ["
           << p.domain.Nx << ", "
           << p.domain.Ny << ", "
           << p.domain.Nz << "]\n";
        os << "PHYSICS:\n";
        os << "-- E0 = " << p.physics.M_E0 << "\n";
        os << "-- nu = " << p.physics.M_nu << "\n";
        os << "-- g = [" 
           << p.physics.M_gravity[0] << ", "
           << p.physics.M_gravity[1] << ", "
           << p.physics.M_gravity[2] << "]\n";
        os << "-- mu_friction = " << p.physics.M_mu_friction << "\n";
        os << "IT:\n";
        os << "--maxit = " << p.it.maxIt << "\n";
        os << "--rtol = " << p.it.rtol << "\n";
        os << "--atol = " << p.it.atol << "\n";
        os << "NITSCHE:\n";
        os << "--theta = " << p.nitsche.theta << "\n";
        os << "--gamma0 = " << p.nitsche.gamma0 << "\n";
        os << "TIME:\n";
        os << "--t0 = " << p.time.t0 << "\n";
        os << "--tend = " << p.time.tend << "\n";
        os << "--dt = " << p.time.dt << "\n";
        os << "NUMERICS:\n";
        os << "--integration = " << p.numerics.integration << "\n";
        os << "--FEMTypeDisplacement = " << p.numerics.FEMTypeDisplacement << "\n";
        os << "--FEMTypeStress = " << p.numerics.FEMTypeStress << "\n";
        os << "--FEMTypeRhs = " << p.numerics.FEMTypeRhs << "\n";
        os << "BC STRINGS:\n";
        os << "-- regionsDirID: ";
        for (const auto& id : p.bc.regionsDirID) {
            os << id << " ";
        }
        os << "\n-- regionsNeuID: ";
        for (const auto& id : p.bc.regionsNeuID) {
            os << id << " ";
        }
        os << "\n-- regionsMixID: ";
        for (const auto& id : p.bc.regionsMixID) {
            os << id << " ";
        }
        os << "\n-- stringsDir: ";
        for (const auto& str : p.bc.stringsDir) {
            os << str << " ";
        }
        os << "\n-- stringsNeu: ";
        for (const auto& str : p.bc.stringsNeu) {
            os << str << " ";
        }
        os << "\n-- stringsMix: ";
        for (const auto& str : p.bc.stringsMix) {
            os << str << " ";
        }
        os << "\nVERBOSE: " << (p.verbose ? "true" : "false") << "\n";
        os << "GMSH: " << (p.gmsh ? "true" : "false") << "\n";
        os << "============================================\n";
        /** \todo: print other parameters */
    
        return os;
    };

    void Params::broadcast(MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);

        if (rank == 0) {
            // Pack all data into a std::vector<char>
            std::vector<char> buf;

            // Helper lambda to append raw data to buf
            auto append = [&](auto& x) {
                const char* p = reinterpret_cast<const char*>(&x);
                buf.insert(buf.end(), p, p + sizeof(x));
            };

            append(domain.dim);
            append(domain.Lx);
            append(domain.Ly);
            append(domain.Lz);
            append(domain.Nx);
            append(domain.Ny);
            append(domain.Nz);
            append(domain.h);
            append(domain.angle);
            serializeString(domain.meshType, buf);

            append(physics.M_E0);
            append(physics.M_nu);
            append(physics.M_lambda);
            append(physics.M_mu);

            // physics.M_gravity is std::vector<scalar_type>, size = 3
            size_t grav_size = physics.M_gravity.size();
            append(grav_size);
            for (auto& val : physics.M_gravity) {
                append(val);
            }

            append(physics.M_mu_friction);

            append(it.maxIt);
            append(it.rtol);
            append(it.atol);

            append(nitsche.theta);
            append(nitsche.gamma0);

            append(time.t0);
            append(time.tend);
            append(time.dt);

            serializeString(numerics.integration, buf);
            serializeString(numerics.FEMTypeDisplacement, buf);
            serializeString(numerics.FEMTypeStress, buf);
            serializeString(numerics.FEMTypeRhs, buf);

            append(verbose);
            append(gmsh);

            // Now broadcast buffer size
            size_t buf_size = buf.size();
            MPI_Bcast(&buf_size, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

            // Broadcast buffer data
            MPI_Bcast(buf.data(), buf_size, MPI_CHAR, 0, comm);

        } else {
            // Receive buffer size
            size_t buf_size;
            MPI_Bcast(&buf_size, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

            std::vector<char> buf(buf_size);
            MPI_Bcast(buf.data(), buf_size, MPI_CHAR, 0, comm);

            // Unpack from buf
            const char* ptr = buf.data();

            auto extract = [&](auto& x) {
                std::memcpy(&x, ptr, sizeof(x));
                ptr += sizeof(x);
            };

            extract(domain.dim);
            extract(domain.Lx);
            extract(domain.Ly);
            extract(domain.Lz);
            extract(domain.Nx);
            extract(domain.Ny);
            extract(domain.Nz);
            extract(domain.h);
            extract(domain.angle);

            deserializeString(domain.meshType, ptr);

            extract(physics.M_E0);
            extract(physics.M_nu);
            extract(physics.M_lambda);
            extract(physics.M_mu);

            size_t grav_size;
            extract(grav_size);
            physics.M_gravity.resize(grav_size);
            for (size_t i=0; i<grav_size; ++i) {
                extract(physics.M_gravity[i]);
            }

            extract(physics.M_mu_friction);

            extract(it.maxIt);
            extract(it.rtol);
            extract(it.atol);

            extract(nitsche.theta);
            extract(nitsche.gamma0);

            extract(time.t0);
            extract(time.tend);
            extract(time.dt);

            deserializeString(numerics.integration, ptr);
            deserializeString(numerics.FEMTypeDisplacement, ptr);
            deserializeString(numerics.FEMTypeStress, ptr);
            deserializeString(numerics.FEMTypeRhs, ptr);

            extract(verbose);
            extract(gmsh);
        }

        int root = 0;
        
        broadcastSizetVector(bc.regionsDirID, root, comm);
        broadcastSizetVector(bc.regionsNeuID, root, comm);
        broadcastSizetVector(bc.regionsMixID, root, comm);

        broadcastStringVector(bc.stringsDir, root, comm);
        broadcastStringVector(bc.stringsNeu, root, comm);
        broadcastStringVector(bc.stringsMix, root, comm);
    }

} // namespace gf