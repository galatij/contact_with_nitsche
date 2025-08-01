#include "ContactEnforcementStrategy.hpp"


namespace gf {

    void
    NitscheContactEnforcement::enforce(getfem::model& md, const getfem::mesh_im& im, bool verbose) const
    {
        md.add_initialized_scalar_data("gammaN", M_gammaN);
        md.add_initialized_scalar_data("theta", M_theta);  // symmetric variant
        
        // Normal gap and stress
        md.add_macro("Pn_u", "(gammaN * un_jump - sig_u_nL)");
        md.add_macro("Pt1_u", "(gammaN * ut1_jump - sig_u_t1)");
        md.add_macro("Pt2_u", "(gammaN * ut2_jump - sig_u_t2)");
        md.add_macro("Pn_v_theta", "(gammaN * vn_jump - theta*sig_v_nL)");
        md.add_macro("Pt1_v_theta", "(gammaN * vt1_jump - theta*sig_v_t1)");
        md.add_macro("Pt2_v_theta", "(gammaN * vt2_jump - theta*sig_v_t2)");

        // Friction threshold
        md.add_macro("Sh", "mu_fric * pos_part(Pn_u)");
        md.add_macro("norm_Pt", "sqrt(Pt1_u*Pt1_u + Pt2_u*Pt2_u)");
        md.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
        md.add_macro("proj_Pt2_u", "Pt2_u * min(1, Sh / (norm_Pt + eps))");
        if (verbose) std::cout << "done.\n";
        
        // Add linear stress brick
        if (verbose) std::cout << "  Adding linear stress brick...";
        getfem::add_linear_term(
            md,
            im,
            "- theta/gammaN * sig_u_nL * sig_v_nL", /** expression */
            Fault, /** region */
            false, /** symmetric */
            false, /** coercive */
            "linear_stress",
            false /** check */
        );
        if (verbose) std::cout << "done.\n";

        // Add KKT condition brick
        if (verbose) std::cout << "  Adding KKT condition brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "1/gammaN * pos_part(Pn_u) * Pn_v_theta",
            Fault,
            false,
            false,
            "KKTbrick"
        );
        if (verbose) std::cout << "done.\n";

        // Add Coulomb condition brick
        if (verbose) std::cout << "  Adding Coulomb friction brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "(1/gammaN) * (proj_Pt1_u * Pt1_v_theta + proj_Pt2_u * Pt2_v_theta)",
            Fault,
            false,
            false,
            "CoulombBrick"
        );
        if (verbose) std::cout << "done.\n";

    }


    void
    PenaltyContactEnforcement::enforce(getfem::model& md, const getfem::mesh_im& im, bool verbose) const
    {   
        md.add_initialized_scalar_data("gammaP", M_gammaP);
        md.add_macro("Sh", "mu_fric * pos_part(un_jump)");
        md.add_macro("norm_ut_jump", "sqrt(ut1_jump*ut1_jump + ut2_jump*ut2_jump)");
        md.add_macro("proj_ut1_jump", "ut1_jump * min(1, Sh / (norm_ut_jump + eps))");
        md.add_macro("proj_ut2_jump", "ut2_jump* min(1, Sh / (norm_ut_jump + eps))");
        if (verbose) std::cout << "done.\n";
        
        // Add KKT condition brick
        if (verbose) std::cout << "  Adding KKT condition brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "gammaP * pos_part(un_jump) * vn_jump",
            Fault,
            false,
            false,
            "KKTbrick"
        );
        if (verbose) std::cout << "done.\n";

        // Add Coulomb condition brick
        if (verbose) std::cout << "  Adding Coulomb friction brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "gammaP*(proj_ut1_jump * vt1_jump + proj_ut2_jump * vt2_jump)",
            Fault,
            false,
            false,
            "CoulombBrick"
        );
        if (verbose) std::cout << "done.\n";

    }


    void
    AugmentedLagrangianContactEnforcement::enforce(getfem::model& md, const getfem::mesh_im& im, bool verbose) const
    {    
        md.add_filtered_fem_variable("mult", M_mfLM, Fault);

        // Initialize the Lagrange multiplier variable with a small value to avoid singularity issues
        plain_vector mult0 (M_mfLM.dof_on_region(Fault).card(), 1.e-3);
        md.set_real_variable("mult") = mult0;

        md.add_initialized_scalar_data("gammaL", M_gammaL);
        
        // Normal gap and stress
        md.add_macro("lambdan", "mult . n");
        md.add_macro("lambdat1", "mult . t1");
        md.add_macro("lambdat2", "mult . t2");
        md.add_macro("mun", "Test_mult . n");
        md.add_macro("mut1", "Test_mult . t1");
        md.add_macro("mut2", "Test_mult . t2");

        md.add_macro("Pn_u", "(gammaL * un_jump - lambdan)");
        md.add_macro("Pt1_u", "(gammaL * ut1_jump - lambdat1)");
        md.add_macro("Pt2_u", "(gammaL * ut2_jump - lambdat2)");

        // Friction threshold
        md.add_macro("Sh", "mu_fric * pos_part(Pn_u)");
        md.add_macro("norm_Pt", "sqrt(Pt1_u*Pt1_u + Pt2_u*Pt2_u)");
        md.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
        md.add_macro("proj_Pt2_u", "Pt2_u * min(1, Sh / (norm_Pt + eps))");
        if (verbose) std::cout << "done.\n";
        
        // Add KKT condition brick
        if (verbose) std::cout << "  Adding KKT condition brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "pos_part(Pn_u) * vn_jump",
            Fault,
            false,
            false,
            "KKTbrick"
        );
        if (verbose) std::cout << "done.\n";

        // Add Coulomb condition brick
        if (verbose) std::cout << "  Adding Coulomb friction brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "proj_Pt1_u * vt1_jump + proj_Pt2_u * vt2_jump",
            Fault,
            false,
            false,
            "CoulombBrick"
        );
        if (verbose) std::cout << "done.\n"; 
        
        // Add LM term for normal gap
        if (verbose) std::cout << "  Adding Lagrange multiplier term for normal gap...";
        getfem::add_nonlinear_term(
            md,
            im,
            " - 1/gammaL * (lambdan + pos_part(Pn_u)) * mun",
            Fault,
            false,
            false,
            "LM_NormalGapBrick"
        );
        if (verbose) std::cout << "done.\n";
        
        // Add LM term for tangential gap
        if (verbose) std::cout << "  Adding Lagrange multiplier term for tangential gap...";
        getfem::add_nonlinear_term(
            md,
            im,
            " -1/gammaL * (lambdat1 + proj_Pt1_u) * mut1 + (lambdat2 + proj_Pt2_u) * mut2",
            Fault,
            false,
            false,
            "LM_TangentialGapBrick"
        );
        if (verbose) std::cout << "done.\n";

        // Add stabilization terms for Lagrange multipliers (negligible, eps = 1.e-20)
        getfem::add_nonlinear_term(
            md, im,
            "eps * lambdan * mun",
            Fault,
            false, false,
            "LM_NormalStab"
        );

        getfem::add_nonlinear_term(
            md, im,
            "eps * lambdat1 . mut1 + eps * lambdat2 . mut2",
            Fault,
            false, false,
            "LM_TangentialStab"
        );

        if (verbose) std::cout << "done.\n";          
    }


} // namespace gf