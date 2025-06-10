#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include "Core.hpp"

namespace gf {

    std::vector<std::string> splitString(std::string stringValue);

    std::vector<size_type> toVec(std::string stringValue);

    void serializeString(const std::string& str, std::vector<char>& buffer);

    void deserializeString(std::string& str, const char* &buffer_ptr);

    void broadcastSizetVector(std::vector<size_type>& vec, int root, MPI_Comm comm);
    
    void broadcastStringVector(std::vector<std::string>& vec, int root, MPI_Comm comm);
    

} // namespace gf

#endif // _UTILS_HPP_