#ifndef _BC_HANDLER_HPP_
#define _BC_HANDLER_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "Mesh.hpp"
#include "BC.hpp"
#include "muParserXInterface.hpp"

#include <unordered_map>

namespace gf {

    class BCHandler {

        using BCListType = std::unordered_map<BCType, std::vector<std::unique_ptr<BC>>>;
        using BCstringsType = std::unordered_map<BCType, std::vector<std::string>>;

        const getfem::mesh& M_mesh; ///< The mesh object
        BCListType M_BCList; ///< The BCobjects
        BCstringsType M_BCStrings; ///< The strings of the prescribed BC functions (just for output)
        
        muParserXInterface M_parser;

        // VectorFunctionType buildBCFunctionFromExpressions(const std::vector<std::string>&);

        template <BCType T>
        void read(const BCStrings&);

    public:
        BCHandler(const getfem::mesh&);

        void readBC(const BCStrings&);

        const std::vector<std::unique_ptr<BC>> & Neumann() const;

        const std::vector<std::unique_ptr<BC>> & Dirichlet() const;

        const std::vector<std::unique_ptr<BC>> & Mixed() const;

    };

} // namespace gf



#endif // _BC_HANDLER_HPP_