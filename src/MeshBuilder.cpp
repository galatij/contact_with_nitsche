#include "MeshBuilder.hpp"

bool DEBUGIR = true;
namespace gf {

    void
    MeshBuilderStrategy::construct(getfem::mesh& mesh) const
    {
        buildMesh(mesh);
        initRegions(mesh);
    }


    // BuiltInBuilder:: IMPLEMENTATION
    void
    BuiltInBuilder::buildMesh(getfem::mesh& mesh) const
    {
        MPI_Barrier(MPI_COMM_WORLD);

        if (getfem::MPI_IS_MASTER())
            std::cout << "Building the mesh internally... ";
        size_type N = M_domain.dim;
        std::string meshType = M_domain.meshType; // element type for meshing, linear by defaults

        auto pgt = bgeot::geometric_trans_descriptor(meshType);

        if (true){
            auto pbDim = pgt->dim();
            assert(N==pbDim);
        }

        std::vector<size_type> nsubdiv(N);
        std::vector<scalar_type> lengths(N);
        nsubdiv[0] = M_domain.Nx;
        nsubdiv[1] = M_domain.Ny;
        nsubdiv[2] = M_domain.Nz;
        lengths[0] = M_domain.Lx;
        lengths[1] = M_domain.Ly;
        lengths[2] = M_domain.Lz;
    
        // create the mesh
        getfem::regular_unit_mesh(mesh, nsubdiv, pgt);

        // transform the mesh (scale the unit mesh to [M_domain.Lx, M_domain.Ly, M_domain.Lz])
        //!\todo: maybe I also need to translate it, I want the fault to be positioned in x=0
        bgeot::base_matrix M(N,N);
        for (size_type i = 0; i < N; ++i)
            M(i,i) = lengths[i];
        mesh.transformation(M);

        base_small_vector v {-0.5*M_domain.Lx,0.,0.};
        mesh.translation(v);
        if (getfem::MPI_IS_MASTER())
            std::cout << "done.\n";

    }


    void
    BuiltInBuilder::initRegions(getfem::mesh& mesh) const
    {
        if (getfem::MPI_IS_MASTER())
            std::cout << "Initializing regions... ";

        // Create the internal regions BulkLeft, BulkRight, Fault
        dal::bit_vector leftConvexesList;
        dal::bit_vector rightConvexesList;
        dal::bit_vector centeredConvexesList;
        
        // loop over all convexes in mesh
        for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
            
            bool insertedLeft = false, insertedRight = false, insertedCentered = false;

            // get the coords of the points of convex i
            auto local_points = mesh.points_of_convex(i);
            
            if (DEBUGIR){
                // std::clog << "Element: " << i << std::endl;
                // std::clog << "--Points\n";
            }

            // loop over the points of a convex
            for (const auto &pt: local_points){
                if (DEBUGIR)
                    // std::clog <<"["<<pt[0]<<","<<pt[1]<<","<<pt[2]<<"]\n";
                
                if (pt[0] < -1e-6 && !insertedLeft){ // if one point has x < 0 --> left bulk region
                    leftConvexesList.add(i);
                    insertedLeft = true;
                }
                else if (pt[0] > 1.e-6 && !insertedRight){ // "   "   "   x > 0 --> right   "   "
                    rightConvexesList.add(i);
                    insertedRight = true;
                }
                if (std::abs(pt[0]) < 1.e-6 && !insertedCentered){
                    if (DEBUGIR) {
                        // std::clog << "Adding an element in Fault convexList: "
                        // << "pt[0] = " << pt[0]
                        // << ", abs(pt[0]) = " << std::abs(pt[0]) << " < 1.e-6? "
                        // << std::boolalpha << (std::abs(pt[0]) < 1.e-6) << std::endl;
                    }
                    centeredConvexesList.add(i);
                    insertedCentered = true;
                }

            }
        }

        mesh.region(RegionType::BulkLeft).add(leftConvexesList);

        mesh.region(RegionType::BulkRight).add(rightConvexesList);

        for (dal::bv_visitor i(centeredConvexesList); !i.finished(); ++i){

            auto psc = mesh.structure_of_convex(i);

            // check the centroid (x < 0 || x > 0)      
            auto pts = mesh.points_of_convex(i);
            base_node centroid(pts[0].size());

            for (const auto& pt : pts)
                centroid += pt;
            centroid /= static_cast<double>(pts.size());

            // loop over the faces of a convex
            for (size_type f = 0; f < psc->nb_faces(); ++f){

                // check the normal
                auto pgt = mesh.trans_of_convex(i);
                auto face_pts = pgt->convex_ref()->points_of_face(f);

                base_small_vector n = mesh.mean_normal_of_face_of_convex(i, f);

                if (std::abs(n[0] - 1.0) < 1.e-6 && std::abs(n[1]) < 1.e-6 &&
                    std::abs(n[2]) < 1.e-6 && centroid[0] < -1.e-6){

                    // add the face to the fault's meshRegion
                    mesh.region(RegionType::Fault).add(i,f);
                }
            }
        }

        
        // Create the boundary regions
        getfem::mesh_region borderFaces;
        getfem::outer_faces_of_mesh(mesh, borderFaces);

        double eps = 1.e-2;

        auto make_box = [eps](const base_node& a, const base_node& b) {
            base_node pt_min(a.size()), pt_max(a.size());
            for (size_type i = 0; i < a.size(); ++i) {
                pt_min[i] = std::min(a[i], b[i]) - eps;
                pt_max[i] = std::max(a[i], b[i]) + eps;
            }
            return std::make_pair(pt_min, pt_max);
        };
    {
        auto [min1, max1] = make_box(
            base_node{-M_domain.Lx/2., 0., 0.},
            base_node{0., M_domain.Ly, 0.}
        );
        mesh.region(1) = getfem::select_faces_in_box(mesh, borderFaces, min1, max1);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{0., M_domain.Ly, 0.},
            base_node{-M_domain.Lx, M_domain.Ly, M_domain.Lz});
        mesh.region(2) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{-M_domain.Lx/2, M_domain.Ly, M_domain.Lz},
            base_node{0., 0., M_domain.Lz});
        mesh.region(3) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{0., 0., M_domain.Lz},
            base_node{-M_domain.Lx/2., 0., 0.});
        mesh.region(4) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }

    {       
        auto [pt_min, pt_max] = make_box(
            base_node{0., 0., 0.},
            base_node{M_domain.Lx/2., M_domain.Ly, 0.});
        mesh.region(5) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{M_domain.Lx/2., M_domain.Ly, 0.},
            base_node{0., M_domain.Ly, M_domain.Lz});
        mesh.region(6) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{0., M_domain.Ly, M_domain.Lz},
            base_node{M_domain.Lx/2., 0., M_domain.Lz});
        mesh.region(7) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {    
        auto [pt_min, pt_max] = make_box(
            base_node{M_domain.Lx/2., 0., M_domain.Lz},
            base_node{0., 0., 0.});
        mesh.region(8) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{-M_domain.Lx/2., 0., 0.},
            base_node{-M_domain.Lx/2., M_domain.Ly, M_domain.Lz});
        mesh.region(9) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{M_domain.Lx/2., 0., 0.},
            base_node{M_domain.Lx/2., M_domain.Ly, M_domain.Lz});
        mesh.region(10) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }

        // // alternatively:
        // // bds(9) = getfem::select_faces_of_normal(mesh, border_faces,
        // //    base_small_vector(-1.,0.,0.), 0.01);
        // // bds(10) = getfem::select_faces_of_normal(mesh, border_faces,
        // //    base_small_vector(1.,0.,0.), 0.01);
        if (getfem::MPI_IS_MASTER())
            std::cout << "done." << std::endl;

    }


    // GmshBuilder:: IMPLEMENTATION
    void
    GmshBuilder::buildMesh(getfem::mesh& mesh) const
    {

        if (getfem::MPI_IS_MASTER()){
            std::cout << "Generating mesh file... ";
            generateMeshFile();
            std::cout << "done." << std::endl;
        }

        // Synchronize all processes: wait for the file to be ready
        MPI_Barrier(MPI_COMM_WORLD);

        if (getfem::MPI_IS_MASTER()) {
            std::cout << "Importing the mesh file... ";
            std::ifstream test("fractured_mesh.msh");
            if (!test.is_open()) {
                throw std::runtime_error("Cannot open fractured_mesh.msh");
            }
        }

        using RegMap = std::map<std::string, size_type>;
        RegMap regmap;

        // Sequentially import the mesh in each rank
        /** \todo: may be done better with a custom serialized buffer broadcast */
        int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk);
        int sz; MPI_Comm_size(MPI_COMM_WORLD, &sz);
        for (int r{}; r < sz; ++r)
            if (r == rk)
                getfem::import_mesh_gmsh("fractured_mesh.msh", mesh, regmap);

        if (getfem::MPI_IS_MASTER()){
            // std::cout << regmap.size() << "\n";
            // for (RegMap::iterator i = regmap.begin(); i != regmap.end(); i++) {
            //     std::cout << i->first << " " << i->second << "\n";
            // }
            std::cout << "done." << std::endl;
        }
    }

    void
    GmshBuilder::initRegions(getfem::mesh& mesh) const
    {

        if (getfem::MPI_IS_MASTER())
            std::cout << "Initializing regions... ";

        getfem::mesh_region &fault_region = mesh.region(Fault);
        const getfem::mesh_region &bulk_region = mesh.region(BulkRight);

        // Loop through (convex, face) pairs
        // std::clog << "Faces in region Fault (convex, face): ";
        for (getfem::mr_visitor it(fault_region); !it.finished(); ++it) {
            if (bulk_region.is_in(it.cv())){
                // std::cout << "Removing (" << it.cv() << "," << it.f() << ")...";
                fault_region.sup(it.cv(),it.f());
            }
            // std::cout << "done." << std::endl;
        }

        if (getfem::MPI_IS_MASTER())
            std::cout << "done." << std::endl;

    }
    
    
    void
    GmshBuilder::generateMeshFile() const
    {
        // Extract parameters from domain (already read from .pot)
        scalar_type Lx = M_domain.Lx, Ly = M_domain.Ly, Lz = M_domain.Lz, h = M_domain.h;
        scalar_type faultX0 = Lx/2;
        scalar_type faultX1 = Lx/2 + Lz*std::tan(M_domain.angle*M_PI/180.0);

        std::ofstream meshfile("fractured_mesh.geo");
        if (!meshfile.is_open())
            throw std::runtime_error("Could not open fractured_mesh.geo for writing.");

        // Write .geo file using domain parameters
        meshfile << "Mesh.MshFileVersion = 2.2;\n"
            << "Mesh.Format = 1;\n\n"
            << "// ------------ Parameters ------------\n"
            << "Lx = " << Lx << ";\n"
            << "Ly = " << Ly << ";\n"
            << "Lz = " << Lz << ";\n"
            << "faultX0 = " << faultX0 << ";\n"
            << "faultX1 = " << faultX1 << ";\n"
            << "h = " << h << ";\n\n"
            << "// ------------ Points ------------\n"
            << "Point(1) = {0, 0, 0, h};\n"
            << "Point(2) = {faultX0, 0, 0, h};\n"
            << "Point(3) = {Lx, 0, 0, h};\n"
            << "Point(4) = {0, Ly, 0, h};\n"
            << "Point(5) = {faultX0, Ly, 0, h};\n"
            << "Point(6) = {Lx, Ly, 0, h};\n"
            << "Point(7) = {0, 0, Lz, h};\n"
            << "Point(8) = {faultX1, 0, Lz, h};\n"
            << "Point(9) = {Lx, 0, Lz, h};\n"
            << "Point(10) = {0, Ly, Lz, h};\n"
            << "Point(11) = {faultX1, Ly, Lz, h};\n"
            << "Point(12) = {Lx, Ly, Lz, h};\n\n"

            << "// ------------ Fault Surface ------------\n"
            << "Line(1) = {1, 2};\n"
            << "Line(2) = {2, 3};\n"
            << "Line(3) = {3, 6};\n"
            << "Line(4) = {6, 5};\n"
            << "Line(5) = {5, 4};\n"
            << "Line(6) = {4,1};\n"
            << "Line(7) = {7,8};\n"
            << "Line(8) = {8,9};\n"
            << "Line(9) = {9,12};\n"
            << "Line(10) = {12,11};\n"
            << "Line(11) = {11,10};\n"
            << "Line(12) = {10,7};\n"
            << "Line(13) = {1,7};\n"
            << "Line(14) = {3,9};\n"
            << "Line(15) = {6,12};\n"
            << "Line(16) = {4,10};\n"
            << "Line(17) = {2,5};\n"
            << "Line(18) = {5,11};\n"
            << "Line(19) = {11,8};\n"
            << "Line(20) = {8,2};\n\n"

            << "// ------------ Left Block Faces ------------\n"
            << "Line Loop(1) = {-6,-5,-17,-1};\n"
            << "Plane Surface(1) = {1};\n"
            << "Line Loop(2) = {5,16,-11,-18};\n"
            << "Plane Surface(2) = {2};\n"
            << "Line Loop(3) = {7,-19,11,12};\n"
            << "Plane Surface(3) = {3};\n"
            << "Line Loop(4) = {1,-20,-7,-13};\n"
            << "Plane Surface(4) = {4};\n"
            << "Line Loop(9) = {6,13,-12,-16};\n"
            << "Plane Surface(9) = {9};\n\n"

            << "// ------------ Right Block Faces ------------\n"
            << "Line Loop(5) = {-2,17,-4,-3};\n"
            << "Plane Surface(5) = {5};\n"
            << "Line Loop(6) = {4,18,-10,-15};\n"
            << "Plane Surface(6) = {6};\n"
            << "Line Loop(7) = {8,9,10,19};\n"
            << "Plane Surface(7) = {7};\n"
            << "Line Loop(8) = {2,14,-8,20};\n"
            << "Plane Surface(8) = {8};\n"
            << "Line Loop(10) = {3,15,-9,-14};\n"
            << "Plane Surface(10) = {10};\n\n"

            << "// ------------ Fault Surface ------------\n"
            << "Line Loop(11) = {17,18,19,20};\n"
            << "Plane Surface(1000) = {11};\n\n"

            << "// ------------ Volume Definitions ------------\n"
            << "Surface Loop(1001) = {1,2,3,4,9,1000};\n"
            << "Volume(101) = {1001};\n"
            << "Surface Loop(1002) = {5,6,7,8,10,-1000};\n"
            << "Volume(102) = {1002};\n\n"

            << "// ------------ Physical Regions ------------\n"
            << "Physical Surface(\"BottomLeft\") = {1};\n"
            << "Physical Surface(\"YmaxLeft\") = {2};\n"
            << "Physical Surface(\"TopLeft\") = {3};\n"
            << "Physical Surface(\"YminLeft\") = {4};\n"
            << "Physical Surface(\"BottomRight\") = {5};\n"
            << "Physical Surface(\"YmaxRight\") = {6};\n"
            << "Physical Surface(\"TopRight\") = {7};\n"
            << "Physical Surface(\"YminRight\") = {8};\n"
            << "Physical Surface(\"Xmin\") = {9};\n"
            << "Physical Surface(\"Xmax\") = {10};\n"
            << "Physical Volume(\"BulkLeft\") = {101};\n"
            << "Physical Volume(\"BulkRight\") = {102};\n"
            << "Physical Surface(\"Fault\") = {1000};\n\n";

        if (M_domain.meshType == "GT_QK(3,1)")         
            meshfile << "// ------------ Mesh Settings ------------\n"
            << "Transfinite Line {1:20} = 10 Using Progression 1;\n"
            << "Transfinite Surface {1,2,3,4,5,6,7,8,9,10,11,1000};\n"
            << "Transfinite Volume {101, 102};\n"
            << "Recombine Surface {1,2,3,4,5,6,7,8,9,10,11,1000};\n"
            << "Recombine Volume {101, 102};\n";

        meshfile.close();

        // Run Gmsh to generate the mesh
        std::string gmsh_cmd = "gmsh fractured_mesh.geo -3 -format msh2 -o fractured_mesh.msh > gmsh.log 2>&1";
        int gmsh_status = std::system(gmsh_cmd.c_str());
        // if (gmsh_status != 0)
        if (gmsh_status == -1 || !WIFEXITED(gmsh_status) || WEXITSTATUS(gmsh_status) != 0)
            throw std::runtime_error("Gmsh mesh generation failed.");

    }


} // namespace gf