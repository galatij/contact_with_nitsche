#include "Utils.hpp"

namespace gf {

    // Remove brackets
    std::vector<std::string>
    splitString(std::string stringValue) {
        stringValue.erase(std::remove(stringValue.begin(), stringValue.end(), '['), stringValue.end());
        stringValue.erase(std::remove(stringValue.begin(), stringValue.end(), ']'), stringValue.end());

        // Split by commas
        std::vector<std::string> components;
        std::istringstream ss(stringValue);
        std::string item;

        while (std::getline(ss, item, ',')) {
            components.push_back(item);
        }

        return components;
    }

    std::vector<size_type>
    toVec(std::string stringValue) {
        auto components = splitString(stringValue);
        std::vector<size_t> result; 
        for (const auto& c:components)
            result.emplace_back(std::stod(c));
        return result;
    }


    void serializeString(const std::string& str, std::vector<char>& buffer) {
        size_t len = str.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + sizeof(size_t));
        buffer.insert(buffer.end(), str.begin(), str.end());
    }


    void deserializeString(std::string& str, const char* &buffer_ptr) {
        size_t len;
        std::memcpy(&len, buffer_ptr, sizeof(size_t));
        buffer_ptr += sizeof(size_t);
        str.assign(buffer_ptr, len);
        buffer_ptr += len;
    }


    void broadcastSizetVector(std::vector<size_type>& vec, int root, MPI_Comm comm) {
        size_type size = vec.size();
        MPI_Bcast(&size, 1, MPI_SIZE_T, root, comm);
        if (!getfem::MPI_IS_MASTER()) vec.resize(size);
        MPI_Bcast(vec.data(), size, MPI_SIZE_T, root, comm);
    }


    void broadcastStringVector(std::vector<std::string>& vec, int root, MPI_Comm comm) {
        int size = static_cast<int>(vec.size());
        MPI_Bcast(&size, 1, MPI_INT, root, comm);
        if (!getfem::MPI_IS_MASTER()) vec.resize(size);
        for (int i = 0; i < size; ++i) {
            int len = 0;
            if (getfem::MPI_IS_MASTER()) len = vec[i].size();
            MPI_Bcast(&len, 1, MPI_INT, root, comm);
            if (!getfem::MPI_IS_MASTER()) vec[i].resize(len);
            MPI_Bcast(&vec[i][0], len, MPI_CHAR, root, comm);
        }
    }


} // namespace gf


    // VectorFunctionType
    // BCHandler::buildBCFunctionFromExpressions(const std::vector<std::string>& components)
    // {
    //     // RMK: i assume that each call to buildBCFunctionFromExpressions resets the parser
        
    //     // Build a function for each component and return the combined VectorFunctionType
    //     VectorFunctionType parsedVectorFunction;
    //     std::vector<ScalarFunctionType> parsedFunctionsVec;

    //     for (size_type k {}; k < 3; ++k) {
    //         // Reset the expression of the parser for the new component
    //         M_parser.set_expression(components[k]);

    //         // Define a lambda function that binds to the parsed expression
    //         auto func = [this](base_node x, scalar_type t) -> scalar_type {
    //             // evaluate the expression for x[0], x[1], x[2], t and return the result as a base_small_vector
    //             std::array<double, 4> inputs = { x[0], x[1], x[2], t }; // assuming base_node has x[0], x[1], x[2]
    //             scalar_type result = M_parser(inputs);  // evaluating the function at the input values
    //             return result;  // return result as scalar_type
    //         };

    //         parsedFunctionsVec[k] = std::move(func);

    //     }

    //     // Combine and return the result
    //     return [parsedFunctionsVec](base_node node, scalar_type t) -> base_small_vector {
    //         base_small_vector result(3);
    //         for (size_type i {}; i < 3 ; ++i) {
    //             result[i] = parsedFunctionsVec[i](node, t);
    //         }
    //         return result;
    //     };

    //     return parsedVectorFunction;
    // }

