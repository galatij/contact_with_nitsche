#include "Core.hpp"
#include "Mesh.hpp"
#include "ContactProblem.hpp"

int main(int argc, char * argv[]){
    using namespace getfem;
    using namespace gf;

    GETFEM_MPI_INIT(argc, argv);

    Params p(argc, argv);

    // int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
    // if (r == 2)
    //     std::cout << "\nr = "<< r <<"\nVerbosity: " << p.verbose 
    //         << "\nFEMDisp = "<< p.numerics.FEMTypeDisplacement
    //         << "\nbdLoad2 = "<< p.bc.stringsNeu[1]
    //         << " on region bd = " << p.bc.regionsNeuID[1]
    //         << std::endl;


    Mesh mesh(p);
    ContactProblem pb(mesh, p);
    
    pb.init();
    pb.assemble();
    pb.solve();

    GETFEM_MPI_FINALIZE;

    return 0;

}