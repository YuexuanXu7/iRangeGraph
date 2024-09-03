#include "iRG_search_multi.h"

std::unordered_map<std::string, std::string> paths;

const int query_K = 10;
int M;

void Generate(iRangeGraph_multi::DataLoader &storage)
{
    storage.synthesize_2Dranges(paths["range_prefix"]);
    storage.LoadRanges(paths["range_prefix"]);
    storage.Generate_Groundtruth(paths["groundtruth_prefix"]);
}

int main(int argc, char **argv)
{
    for (int i = 0; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--data_path")
            paths["data_vector"] = argv[i + 1];
        if (arg == "--query_path")
            paths["query_vector"] = argv[i + 1];
        if (arg == "--range_saveprefix")
            paths["range_prefix"] = argv[i + 1];
        if (arg == "--groundtruth_saveprefix")
            paths["groundtruth_prefix"] = argv[i + 1];
        if (arg == "--index_file")
            paths["index"] = argv[i + 1];
        if (arg == "--result_saveprefix")
            paths["result_saveprefix"] = argv[i + 1];
        if (arg == "--attribute1")
            paths["attribute1"] = argv[i + 1];
        if (arg == "--attribute2")
            paths["attribute2"] = argv[i + 1];
        if (arg == "--M")
            M = std::stoi(argv[i + 1]);
    }

    if (argc != 19)
        throw Exception("please check input parameters");

    iRangeGraph_multi::DataLoader storage;
    storage.query_K = query_K;
    storage.LoadQuery(paths["query_vector"]);
    storage.LoadData(paths["data_vector"]);
    // the order of load in attribute values should not be switched
    // the first attribute should be first loaded
    storage.LoadAttribute(paths["attribute1"]);
    storage.LoadAttribute(paths["attribute2"]);

    // Generate should be called when running for the first time; otherwise, it can be skipped.
    Generate(storage);

    storage.LoadRanges(paths["range_prefix"]);
    storage.LoadGroundtruth(paths["groundtruth_prefix"]);
    storage.Sort_by_Attr(0);

    iRangeGraph_multi::iRangeGraph_Search_Multi<float> index(paths["index"], &storage, M);
    index.setprob();
    std::vector<int>
        SearchEF = {1400, 700, 400, 300, 250, 200, 180, 160, 140, 120, 100, 90, 80, 70, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10};
    index.search(SearchEF, paths["result_saveprefix"], M);
}