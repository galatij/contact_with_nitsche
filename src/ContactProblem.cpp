#include "ContactProblem.hpp"

namespace gf{

    bool contact_flag = false;

    ContactProblem::ContactProblem(const Mesh& m, const Params& p)
    : M_mesh(m),
    M_params(p),
    M_BC(m.get()),
    M_FEM(m.get()),
    M_integrationMethod(m.get()),
    M_integrationMethodSurface(m.get()),
    M_imData(M_integrationMethod)
    {
    }

    void
    ContactProblem::init()
    {
        // std::cout << "Inside ConstactProblem::init()" << std::endl;
        // std::cout << M_params.numerics.FEMTypeRhs << std::endl;

        M_BC.readBC(M_params.bc);

        M_FEM.setMeshFem(M_params.numerics, M_mesh.get());

        getfem::pintegration_method ppi = getfem::int_method_descriptor(M_params.numerics.integration);
        size_type N = M_params.domain.dim;
        M_integrationMethod.set_integration_method(M_mesh.get().convex_index(), ppi);
        // M_imData.set_tensor_size(bgeot::multi_index(N,N)); /** \todo */

    }


    void
    ContactProblem::assemble()
    {
        using getfem::MPI_IS_MASTER;

        if (MPI_IS_MASTER()) std::cout << "Preparing the assembly phase:\n";
        
        // gmm::set_traces_level(1);
        
        dim_type dim = M_mesh.get().dim();
        size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();

        // Main unknown of the problem (displacement)
        if (MPI_IS_MASTER()) std::cout << "  Defining variables...";
        
        M_model.add_fem_variable("uL", M_FEM.mf_u1());
        M_model.add_fem_variable("uR", M_FEM.mf_u2());
        
        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Add scalar data to the model
        if (MPI_IS_MASTER()) std::cout << "  Initializing scalar data...";
        
        M_model.add_initialized_scalar_data("lambda", M_params.physics.M_lambda);
        M_model.add_initialized_scalar_data("mu", M_params.physics.M_mu);
        M_model.add_initialized_scalar_data("gammaN", 10*M_params.physics.M_E0/M_params.domain.h);
        M_model.add_initialized_scalar_data("theta", M_params.nitsche.theta);  // symmetric variant
        M_model.add_initialized_scalar_data("mu_fric", M_params.physics.M_mu_friction);
        
        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Define some useful macro for Nitsche contact integrals
        if (MPI_IS_MASTER()) std::cout << "  Defining macros...";

        M_model.add_initialized_scalar_data("eps", 1.e-20);
        M_model.add_macro("n", "Normal"); // use normal on contact face
        M_model.add_macro("eps_dir", "1e-4"); // to avoid numerical degeneracies
        M_model.add_macro("a", "[0, 1, 0]"); /** \todo Modify to be more robust... */

        // Project a onto plane orthogonal to n
        M_model.add_macro("t1_raw", "a - (a . n) * n");
        M_model.add_macro("t1", "Normalized(t1_raw)");

        // Cross to get the second tangent
        M_model.add_macro("t2", "Cross_product(n, t1)");

        // Define macros for jump
        M_model.add_macro("u_jump", "(uL - Interpolate(uR,neighbor_element))");
        M_model.add_macro("v_jump", "(Test_uL - Interpolate(Test_uR,neighbor_element))");
        M_model.add_macro("un_jump", "u_jump . n"); // normal trial jump
        M_model.add_macro("ut1_jump", "u_jump . t1");
        M_model.add_macro("ut2_jump", "u_jump . t2");
        M_model.add_macro("vn_jump", "v_jump . n"); // normal test jump
        M_model.add_macro("vt1_jump", "v_jump . t1");
        M_model.add_macro("vt2_jump", "v_jump . t2");

        M_model.add_macro("traction_u", "((lambda*Trace(Grad_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_uL)) * n)");
        M_model.add_macro("traction_v", "((lambda*Trace(Grad_Test_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_Test_uL)) * n)");
        M_model.add_macro("sig_u_nL", "traction_u . n");
        M_model.add_macro("sig_u_t1", "traction_u . t1");
        M_model.add_macro("sig_u_t2", "traction_u . t2");
        M_model.add_macro("sig_v_nL", "traction_v . n");
        M_model.add_macro("sig_v_t1", "traction_v . t1");
        M_model.add_macro("sig_v_t2", "traction_v . t2");

        // Normal gap and stress
        M_model.add_macro("Pn_u", "(gammaN * un_jump - sig_u_nL)");
        M_model.add_macro("Pt1_u", "(gammaN * ut1_jump - sig_u_t1)");
        M_model.add_macro("Pt2_u", "(gammaN * ut2_jump - sig_u_t2)");
        M_model.add_macro("Pn_v_theta", "(gammaN * vn_jump - theta*sig_v_nL)");
        M_model.add_macro("Pt1_v_theta", "(gammaN * vt1_jump - theta*sig_v_t1)");
        M_model.add_macro("Pt2_v_theta", "(gammaN * vt2_jump - theta*sig_v_t2)");
        
        // Friction threshold
        M_model.add_macro("Sh", "mu_fric * pos_part(Pn_u)");
        M_model.add_macro("norm_Pt", "sqrt(Pt1_u*Pt1_u + Pt2_u*Pt2_u)");
        M_model.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
        M_model.add_macro("proj_Pt2_u", "Pt2_u * min(1, Sh / (norm_Pt + eps))");

        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Add isotropic elasticity bricks
        if (MPI_IS_MASTER()) std::cout << "  Adding elasticity bricks...";
        
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uL", "lambda", "mu", RegionType::BulkLeft);
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uR", "lambda", "mu", RegionType::BulkRight);
        
        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Add linear stress brick
        if (MPI_IS_MASTER()) std::cout << "  Adding linear stress brick...";

        getfem::add_linear_term(
            M_model,
            M_integrationMethod,
            "- theta/gammaN * sig_u_nL * sig_v_nL", /** expression */
            Fault, /** region */
            false, /** symmetric */
            false, /** coercive */
            "linear_stress",
            false /** check */
        );

        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Add KKT condition brick
        if (MPI_IS_MASTER()) std::cout << "  Adding KKT condition brick...";

        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "1/gammaN * pos_part(Pn_u) * Pn_v_theta",
            Fault,
            false,
            false,
            "KKTbrick"
        );

        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Add Coulomb condition brick
        if (MPI_IS_MASTER()) std::cout << "  Adding Coulomb friction brick...";
        
        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "(1/gammaN) * (proj_Pt1_u * Pt1_v_theta + proj_Pt2_u * Pt2_v_theta)",
            Fault,
            false,
            false,
            "CoulombBrick"
        );
        
        if (MPI_IS_MASTER()) std::cout << "done.\n";


        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "eps * uL . Test_uL + eps * uR . Test_uR",
            BulkLeft,
            false,
            false,
            "stabilizer1"
        );
        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "eps * uL . Test_uL + eps * uR . Test_uR",
            BulkRight,
            false,
            false,
            "stabilizer2"
        );

        
        // Volumic source term (gravity)
        if (MPI_IS_MASTER()) std::cout << "  Adding volumic source term brick...";
        
        plain_vector G(M_FEM.mf_rhs().nb_dof()*dim);

        for (size_type i = 0; i < nb_dof_rhs; ++i)
            gmm::copy(M_params.physics.M_gravity, gmm::sub_vector(G, gmm::sub_interval(i*dim, dim)));

        M_model.add_initialized_fem_data("VolumicData", M_FEM.mf_rhs(), G);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uL", "VolumicData", BulkLeft);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uR", "VolumicData", BulkRight);
        
        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Neumann conditions
        if (MPI_IS_MASTER()) std::cout << "  Adding Neumann condition bricks...";
        
        const auto & NeumannBCs = M_BC.Neumann();
        plain_vector F(nb_dof_rhs*dim);
        for (const auto& bc: NeumannBCs){
            const auto& rg = bc->getRegion();
            auto& f = bc->f();
            // just for debug
            auto f0 = [&f](const base_node&x){ return f(x, 0); };
            getfem::interpolation_function(M_FEM.mf_rhs(),F,f0, rg);
            M_model.add_initialized_fem_data(bc->name(), M_FEM.mf_rhs(),F);

            if (bc->isLeftBoundary())
                getfem::add_source_term_brick(M_model, M_integrationMethod,"uL",bc->name(),bc->ID());
            else
                getfem::add_source_term_brick(M_model, M_integrationMethod,"uR",bc->name(),bc->ID());
        }
        if (MPI_IS_MASTER()) std::cout << "done.\n";

        // Dirichlet conditions
        if (MPI_IS_MASTER()) std::cout << "  Adding Dirichlet condition bricks...";
        const auto & DirichletBCs = M_BC.Dirichlet();
        plain_vector D(nb_dof_rhs*dim);

        for (const auto& bc: DirichletBCs){
            const auto& rg = bc->getRegion();
            auto& f = bc->f();
            // take time 0
            auto f0 = [&f](const base_node&x){ return f(x, 0); };

            /// @DEBUG:            
            getfem::interpolation_function(M_FEM.mf_rhs(),D,f0, rg);
            M_model.add_initialized_fem_data(bc->name(), M_FEM.mf_rhs(),D);

            if (bc->isLeftBoundary())
                getfem::add_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uL", M_FEM.mf_u1(), bc->ID(), 
                     bc->name());
            else
                getfem::add_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uR", M_FEM.mf_u2(), bc->ID(), 
                     bc->name());
        }

        if (MPI_IS_MASTER()) std::cout << "done.\n";


        // Mixed conditions (normal Dirichlet)
        if (MPI_IS_MASTER()) std::cout << "  Adding normal Dirichlet condition bricks...";

        const auto & MixedBCs = M_BC.Mixed();
        plain_vector M(nb_dof_rhs);

        for (const auto& bc: MixedBCs){
            const auto& rg = bc->getRegion();
            auto& f = bc->f();

            auto f0 = [&f](const base_node&x){ return f(x, 0); };
            getfem::interpolation_function(M_FEM.mf_rhs(), M, f0, rg);

            M_model.add_initialized_fem_data(bc->name(), M_FEM.mf_rhs(),M);

            if (bc->isLeftBoundary())
                // getfem::add_normal_Dirichlet_condition_with_penalization
                //     (M_model, M_integrationMethod, "uL", 1e-5, bc->ID(), bc->name());
                getfem::add_normal_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uL", M_FEM.mf_rhs(), bc->ID(), bc->name());
            else
                // getfem::add_normal_Dirichlet_condition_with_penalization
                //     (M_model, M_integrationMethod, "uR", 1e-5, bc->ID(), bc->name());
                getfem::add_normal_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uR", M_FEM.mf_rhs(), bc->ID(), bc->name());  
        }

        if (MPI_IS_MASTER()) std::cout << "done.\n";

    }

    void
    ContactProblem::solve() {

        using getfem::MPI_IS_MASTER;

        if (MPI_IS_MASTER()) std::cout << "Solving the problem..." << std::endl;

        const auto & NeumannBCs = M_BC.Neumann();
        gmm::set_traces_level(0);
        
        dim_type dim = M_mesh.get().dim();
        size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();
        plain_vector F(nb_dof_rhs*dim);

        // Time loop
        size_type n_timesteps = static_cast<size_type>(M_params.time.tend - M_params.time.t0)/M_params.time.dt;
        
        scalar_type t {};

        for (size_type i{}; i < n_timesteps; ++i)
        {
            if (MPI_IS_MASTER()) std::cout << "t = " << t << std::endl;
            // Neumann conditions
            for (const auto& bc: NeumannBCs){
                const auto& rg = bc->getRegion();
                auto& f = bc->f();
                // just for debug
                auto ft = [&f,t](const base_node&x){ return f(x, t); };
                getfem::interpolation_function(M_FEM.mf_rhs(),F,ft, rg);
                gmm::copy(F, M_model.set_real_variable(bc->name()));
            }

            /** \todo Customize the choice of the solver */
            // Solve the problem
            gmm::iteration iter(M_params.it.rtol, 1, M_params.it.maxIt);
            getfem::standard_solve(M_model,iter);

            // Export results
            exportVtk(i);

            // Advance in time
            t += M_params.time.dt;
        }

    }

    
    void ContactProblem::exportVtk(size_type i) {

        dim_type dim = M_mesh.get().dim();

        // const plain_vector& uL = M_model.real_variable("uL");
        // const plain_vector& uR = M_model.real_variable("uR");

        // dal::bit_vector dofsL = M_FEM.mf_u1().basic_dof_on_region(BulkLeft);
        // dal::bit_vector dofsR = M_FEM.mf_u2().basic_dof_on_region(BulkRight);
        // dal::bit_vector dofsF = M_FEM.mf_u1().basic_dof_on_region(Fault);

        // size_type nb_dof_L = dofsL.card();
        // size_type nb_dof_R = dofsR.card();
        // size_type nb_dof_F = dofsF.card();
        // size_type nb_dof = M_FEM.mf_u1().nb_dof();
        // size_type nb_dof_tot = nb_dof_L+nb_dof_R; // == nb_dof + nb_dof_L

        // // for (dal::bv_visitor bv(dofsL); !bv.finished();++bv) std::cout << bv << " ";
        // // std::cout << std::endl;
        
        // // Compute the actual intersection of Left and Right
        // dal::bit_vector intersection = dofsL & dofsR;

        // // std::cout << "\nDofsF: " << std::endl;
        // // for (dal::bv_visitor bv(dofsF); !bv.finished();++bv) std::cout << bv << " ";


        // plain_vector U;
        // U.reserve(nb_dof_tot);
        // std::vector<base_node> dof_coords;
        // dof_coords.reserve(nb_dof_tot);

        // // Map original dof index to new index in U
        // std::map<size_type, size_type> dof_map;
        // std::map<size_type, size_type> dof_to_point;

        // // Insert uL values
        // for (dal::bv_visitor ii(dofsL); !ii.finished(); ++ii) {
        //     size_type dof = ii;
        //     if (dof % dim == 0) {  // only once per node
        //         base_node pt = M_FEM.mf_u1().point_of_basic_dof(dof);
        //         pt[0] -= 1e-16;
        //         dof_coords.push_back(pt);
        //         size_type pt_idx = dof_coords.size() - 1;
        //         dof_to_point[dof] = pt_idx;
        //     }
        //     U.push_back(uL[dof]);
        // }

        // // Insert uR values
        // for (dal::bv_visitor jj(dofsR); !jj.finished(); ++jj) {
        //     size_type dof = jj;
        //     if (dof % dim == 0) {  // only once per node
        //         base_node pt = M_FEM.mf_u2().point_of_basic_dof(dof);
        //         pt[0] += 1e-16;
        //         dof_coords.push_back(pt);
        //         size_type pt_idx = dof_coords.size() - 1;
        //         dof_to_point[dof] = pt_idx;
        //     }
        //     U.push_back(uL[dof]);
        // }


        // std::vector<base_small_vector> U_vector;
        // U_vector.reserve(U.size() / dim);

        // for (size_type i = 0; i < U.size(); i += dim) {
        //     base_small_vector vec(dim);
        //     for (size_type d = 0; d < dim; ++d)
        //         vec[d] = U[i + d];
        //     U_vector.emplace_back(vec);
        // }

        // assert(U_vector.size() == dof_coords.size());

        // // Build connectivity matrix
        // std::vector<std::vector<size_type>> connectivity;

        // for (dal::bv_visitor cv(M_mesh.get().convex_index()); !cv.finished(); ++cv) {
        //     auto dofs = M_FEM.mf_u1().ind_basic_dof_of_element(cv);
        //     std::vector<size_type> cell;

        //     for (size_type d = 0; d < dofs.size(); d += dim) {
        //         size_type dof_id = dofs[d]; // one per point
        //         auto it = dof_to_point.find(dof_id);
        //         if (it != dof_to_point.end()) {
        //             cell.push_back(it->second);
        //         } else {
        //             std::cerr << "DOF " << dof_id << " not found in dof_to_point!\n";
        //         }
        //     }

        //     if ((dim == 3 && (cell.size() == 4 || cell.size() == 8)) ||
        //         (dim == 2 && (cell.size() == 3 || cell.size() == 4))) {
        //         connectivity.push_back(cell);
        //     }
        // }

        // // Export to vtk (custom)
        // std::string filename = "result_" + std::to_string(i) + ".vtk";

        // std::ofstream out(filename);
        // if (!out) {
        //     std::cerr << "Cannot open file " << filename << " for writing.\n";
        //     return;
        // }

        // out << "# vtk DataFile Version 3.0\n";
        // out << "Custom GetFEM Export\n";
        // out << "ASCII\n";
        // out << "DATASET UNSTRUCTURED_GRID\n";

        // out << "POINTS " << dof_coords.size() << " float\n";
        // for (const auto& pt : dof_coords) {
        //     for (size_t d = 0; d < dim; ++d)
        //         out << pt[d] << " ";
        //     out << "\n";
        // }


        // size_t num_cells = connectivity.size();
        // size_t total_ints = 0;
        // for (const auto& c : connectivity)
        //     total_ints += (1 + c.size()); // size prefix + indices

        // out << "CELLS " << num_cells << " " << total_ints << "\n";
        // for (const auto& c : connectivity) {
        //     out << c.size();
        //     for (auto idx : c)
        //         out << " " << idx;
        //     out << "\n";
        // }

        // out << "CELL_TYPES " << num_cells << "\n";
        // for (const auto& c : connectivity) {
        //     int vtk_type = -1;
        //     if (c.size() == 3) vtk_type = 5;        // triangle
        //     else if (c.size() == 4) vtk_type = 10;  // tetrahedron
        //     else if (c.size() == 8) vtk_type = 12;  // hexahedron
        //     else if (c.size() == 6) vtk_type = 13;  // wedge (prism), optional
        //     else if (c.size() == 5) vtk_type = 14;  // pyramid, optional
        //     else throw std::runtime_error("Unsupported element type with " + std::to_string(c.size()) + " points.");
        //     out << vtk_type << "\n";
        // }

        // out << "POINT_DATA " << U_vector.size() << "\n";
        // out << "VECTORS U float\n";
        // for (const auto& v : U_vector) {
        //     for (size_t d = 0; d < dim; ++d)
        //         out << v[d] << " ";
        //     out << "\n";
        // }

        // std::cout << "VTK file written to: " << filename << "\n";

        // getfem::vtk_export exp("result_" + std::to_string(i) + ".vtk");
        // exp.exporting(M_FEM.mf_u1());
        // exp.write_mesh();
        // exp.write_point_data(M_FEM.mf_u1(), U, "u");

        getfem::vtk_export exp("result_" + std::to_string(i) + ".vtk");
        exp.exporting(M_FEM.mf_u1());
        exp.write_mesh();
        exp.write_point_data(M_FEM.mf_u1(), M_model.real_variable("uL"), "uL");
        exp.write_point_data(M_FEM.mf_u2(), M_model.real_variable("uR"), "uR");

    }

} // namespace gf
