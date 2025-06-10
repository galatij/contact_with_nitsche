#include "BCHandler.hpp"
#include "Utils.hpp"
#include <sstream>
#include <algorithm>

bool DEBUGBC = false;

namespace gf {

    BCHandler::BCHandler(const getfem::mesh& m)
    : M_mesh(m) {
    }

    void BCHandler::readBC(const BCStrings& bcStr) {
        read<BCType::Dirichlet>(bcStr);
        read<BCType::Neumann>(bcStr);
        read<BCType::Mixed>(bcStr);
    }


    template <BCType T>
    void
    BCHandler::read(const BCStrings& bcStr)
    {
        if constexpr (T == BCType::Dirichlet)
        {
            for (size_t i = 0; i < bcStr.regionsDirID.size(); ++i) {

                // Use muparser
                M_parser.set_expression(bcStr.stringsDir[i]);
                // M_BCStrings[T].emplace_back(stringValue);

                getfem::mr_visitor it(M_mesh.region(bcStr.regionsDirID[i]));
                base_small_vector n = M_mesh.mean_normal_of_face_of_convex(it.cv(), it.f());

                // Build the BCDir object bc
                auto bc = std::make_unique<BCDir>(M_mesh.region(bcStr.regionsDirID[i]), bcStr.regionsDirID[i], M_parser, T, n);
                
                // Add to map
                M_BCList[T].emplace_back(std::move(bc));

                if (getfem::MPI_IS_MASTER() && DEBUGBC)
                    std::cout << "DEBUG: evaluating bdDisp" << i <<" at ((1,2,3),100): ("
                    << M_BCList[BCType::Dirichlet][i]->eval({1,2,3},100)[0] << ", "
                    << M_BCList[BCType::Dirichlet][i]->eval({1,2,3},100)[1] << ", "
                    << M_BCList[BCType::Dirichlet][i]->eval({1,2,3},100)[2]
                    << ")" << std::endl;
            }
        }
        else if constexpr (T == BCType::Neumann)
        {
            for (size_t i = 0; i < bcStr.regionsNeuID.size(); ++i) {

                // Use muparser
                M_parser.set_expression(bcStr.stringsNeu[i]);
                // M_BCStrings[T].emplace_back(stringValue);

                getfem::mr_visitor it(M_mesh.region(bcStr.regionsNeuID[i]));
                base_small_vector n = M_mesh.mean_normal_of_face_of_convex(it.cv(), it.f());

                // Build the BCNeu object bc
                auto bc = std::make_unique<BCNeu>(M_mesh.region(bcStr.regionsNeuID[i]), bcStr.regionsNeuID[i], M_parser, T, n);
                // Add to map
                M_BCList[T].emplace_back(std::move(bc));
            }
        }
        else if constexpr (T == BCType::Mixed)
        {
            for (size_t i = 0; i < bcStr.regionsMixID.size(); ++i) {

                // Use muparser
                M_parser.set_expression(bcStr.stringsMix[i]);
                // M_BCStrings[T].emplace_back(stringValue);

                getfem::mr_visitor it(M_mesh.region(bcStr.regionsMixID[i]));
                base_small_vector n = M_mesh.mean_normal_of_face_of_convex(it.cv(), it.f());

                // Build the BCMixed object bc
                auto bc = std::make_unique<BCMix>(M_mesh.region(bcStr.regionsMixID[i]), bcStr.regionsMixID[i], M_parser, T, n);
                // Add to map
                M_BCList[BCType::Mixed].emplace_back(std::move(bc));
            }
        }

    }



    const std::vector<std::unique_ptr<BC>> & 
    BCHandler::Neumann() const {
        static const std::vector<std::unique_ptr<BC>> empty;
        auto it = M_BCList.find(BCType::Neumann);
        return (it != M_BCList.end()) ? it->second : empty;
        // return M_BCList.at(BCType::Neumann);
    }

    const std::vector<std::unique_ptr<BC>> &
    BCHandler::Dirichlet() const {
        static const std::vector<std::unique_ptr<BC>> empty;
        auto it = M_BCList.find(BCType::Dirichlet);
        return (it != M_BCList.end()) ? it->second : empty;
        // return M_BCList.at(BCType::Dirichlet);
    }

    const std::vector<std::unique_ptr<BC>> & 
    BCHandler::Mixed() const {
        static const std::vector<std::unique_ptr<BC>> empty;
        auto it = M_BCList.find(BCType::Mixed);
        return (it != M_BCList.end()) ? it->second : empty;
    }


} // namespace gf