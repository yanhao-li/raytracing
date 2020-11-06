#include "types.h"

// JSON parser library (https://github.com/nlohmann/json)
#include "json.hpp"
using json = nlohmann::json;

void load_off(const std::string &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
Scene load_scene(const std::string &filename);