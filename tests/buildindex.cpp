#include "construction.h"

std::unordered_map<std::string, std::string> paths;

int M;
int ef_construction;
int threads;

int main(int argc, char **argv)
{
    for (int i = 0; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--data_path")
            paths["data_vector"] = argv[i + 1];
        if (arg == "--index_file")
            paths["index_save"] = argv[i + 1];
        if (arg == "--M")
            M = std::stoi(argv[i + 1]);
        if (arg == "--ef_construction")
            ef_construction = std::stoi(argv[i + 1]);
        if (arg == "--threads")
            threads = std::stoi(argv[i + 1]);
    }

    if (paths["data_vector"] == "")
        throw Exception("data path is empty");
    if (paths["index_save"] == "")
        throw Exception("index path is empty");
    if (M <= 0)
        throw Exception("M should be a positive integer");
    if (ef_construction <= 0)
        throw Exception("ef_construction should be a positive integer");
    if (threads <= 0)
        throw Exception("threads should be a positive integer");

    iRangeGraph::DataLoader storage;
    storage.LoadData(paths["data_vector"]);
    iRangeGraph::iRangeGraph_Build<float> index(&storage, M, ef_construction);
    index.max_threads = threads;
    index.buildandsave(paths["index_save"]);
}