#pragma once

#include "space_l2.h"
#include <filesystem>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <sys/time.h>
#include <map>

class Exception : public std::runtime_error
{
public:
    Exception(const std::string &msg) : std::runtime_error(msg) {}
};

void CheckPath(std::string filename)
{
    std::filesystem::path pathObj(filename);
    std::filesystem::path dirPath = pathObj.parent_path();
    if (!std::filesystem::exists(dirPath))
    {
        try
        {
            if (std::filesystem::create_directories(dirPath))
            {
                std::cout << "Directory created: " << dirPath << std::endl;
            }
            else
            {
                std::cerr << "Failed to create directory: " << dirPath << std::endl;
            }
        }
        catch (std::filesystem::filesystem_error &e)
        {
            throw Exception(e.what());
        }
    }
}

float GetTime(timeval &begin, timeval &end)
{
    return end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;
}

namespace iRangeGraph
{
    typedef std::pair<float, int> PFI;
    typedef unsigned int tableint;
    typedef unsigned int linklistsizeint;

    class DataLoader
    {
    public:
        int Dim, query_nb, query_K;
        std::vector<std::vector<float>> query_points;
        int data_nb;
        std::vector<std::vector<float>> data_points;
        std::unordered_map<int, std::vector<std::pair<int, int>>> query_range;
        std::unordered_map<int, std::vector<std::vector<int>>> groundtruth;

        DataLoader() {}
        ~DataLoader() {}

        // query vector filename format: 4 bytes: query number; 4 bytes: dimension; query_nb*Dim vectors
        void LoadQuery(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
                throw Exception("cannot open " + filename);
            infile.read((char *)&query_nb, sizeof(int));
            infile.read((char *)&Dim, sizeof(int));
            query_points.resize(query_nb);
            for (int i = 0; i < query_nb; i++)
            {
                query_points[i].resize(Dim);
                infile.read((char *)query_points[i].data(), Dim * sizeof(float));
            }
            infile.close();
        }

        // Used only when computing groundtruth and constructing index. Do not use this to load data for search process
        void LoadData(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
                throw Exception("cannot open " + filename);
            infile.read((char *)&data_nb, sizeof(int));
            infile.read((char *)&Dim, sizeof(int));
            data_points.resize(data_nb);
            for (int i = 0; i < data_nb; i++)
            {
                data_points[i].resize(Dim);
                infile.read((char *)data_points[i].data(), Dim * sizeof(float));
            }
            infile.close();
        }

        // By default generation, 0.bin~9.bin denotes 2^0~2^-9 range fractions, 17.bin denotes mixed range fraction.
        // Before reading the query ranges, make sure query vectors have been read.
        void LoadQueryRange(std::string fileprefix)
        {
            std::vector<int> s;
            for (int i = 0; i < 10; i++)
                s.emplace_back(i);
            s.emplace_back(17);
            for (auto suffix : s)
            {
                std::string filename = fileprefix + std::to_string(suffix) + ".bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                for (int i = 0; i < query_nb; i++)
                {
                    int ql, qr;
                    infile.read((char *)&ql, sizeof(int));
                    infile.read((char *)&qr, sizeof(int));
                    query_range[suffix].emplace_back(ql, qr);
                }
                infile.close();
            }
        }

        // 0.bin~9.bin correspond to groundtruth for 2^0~2^-9 range fractions, 17.bin for mixed fraction
        void LoadGroundtruth(std::string fileprefix)
        {
            for (auto t : query_range)
            {
                int suffix = t.first;
                std::string filename = fileprefix + std::to_string(suffix) + ".bin";
                std::ifstream infile(filename, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                    throw Exception("cannot open " + filename);
                groundtruth[suffix].resize(query_nb);
                for (int i = 0; i < query_nb; i++)
                {
                    groundtruth[suffix][i].resize(query_K);
                    infile.read((char *)groundtruth[suffix][i].data(), query_K * sizeof(int));
                }
                infile.close();
            }
        }
    };

    class QueryGenerator
    {
    public:
        int data_nb, query_nb;
        hnswlib::L2Space *space;

        QueryGenerator(int data_num, int query_num) : data_nb(data_num), query_nb(query_num) {}
        ~QueryGenerator() {}

        void GenerateRange(std::string saveprefix)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::vector<std::pair<int, int>> rs;
            int current_len = data_nb;
            for (int i = 0; i < 10; i++)
            {
                if (current_len < 10)
                    throw Exception("dataset size is too small, increase the amount of data objects!");
                rs.emplace_back(current_len, i);
                current_len /= 2;
            }
            for (auto t : rs)
            {
                int len = t.first, suffix = t.second;
                std::string savepath = saveprefix + std::to_string(suffix) + ".bin";
                CheckPath(savepath);
                std::cout << "save query range to" << savepath << std::endl;
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);
                std::uniform_int_distribution<int> u_start(0, data_nb - len);
                for (int i = 0; i < query_nb; i++)
                {
                    int ql = u_start(e);
                    int qr = ql + len - 1;
                    if (ql >= data_nb || qr >= data_nb)
                        throw Exception("Query range out of bound");
                    outfile.write((char *)&ql, sizeof(int));
                    outfile.write((char *)&qr, sizeof(int));
                }
                outfile.close();
            }

            rs.clear();
            current_len = data_nb;
            for (int i = 0; i < 10; i++)
            {
                rs.emplace_back(current_len, i);
                current_len /= 2;
            }
            std::string savepath = saveprefix + "17.bin";
            CheckPath(savepath);
            std::cout << "save query range to" << savepath << std::endl;
            std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
            if (!outfile.is_open())
                throw Exception("cannot open " + savepath);

            for (auto t : rs)
            {
                int len = t.first;
                std::uniform_int_distribution<int> u_start(0, data_nb - len);

                for (int i = 0; i < query_nb / 10; i++)
                {
                    int ql = u_start(e);
                    int qr = ql + len - 1;
                    if (ql >= data_nb || qr >= data_nb)
                        throw Exception("Query range out of bound");
                    outfile.write((char *)&ql, sizeof(int));
                    outfile.write((char *)&qr, sizeof(int));
                }
            }
            outfile.close();
        }

        float dis_compute(std::vector<float> &v1, std::vector<float> &v2)
        {
            hnswlib::DISTFUNC<float> fstdistfunc_ = space->get_dist_func();
            float dis = fstdistfunc_((char *)v1.data(), (char *)v2.data(), space->get_dist_func_param());
            return dis;
        }

        void GenerateGroundtruth(std::string saveprefix, DataLoader &storage)
        {
            space = new hnswlib::L2Space(storage.Dim);
            for (auto t : storage.query_range)
            {
                int suffix = t.first;
                std::string savepath = saveprefix + std::to_string(suffix) + ".bin";
                CheckPath(savepath);
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                    throw Exception("cannot open " + savepath);
                std::cout << "generating for " << t.first << std::endl;
                for (int i = 0; i < query_nb; i++)
                {
                    auto rp = t.second[i];
                    int ql = rp.first, qr = rp.second;
                    std::priority_queue<std::pair<float, int>> ans;
                    for (int j = ql; j <= qr; j++)
                    {
                        float dis = dis_compute(storage.query_points[i], storage.data_points[j]);
                        ans.emplace(dis, j);
                        if (ans.size() > storage.query_K)
                            ans.pop();
                    }
                    while (ans.size())
                    {
                        auto id = ans.top().second;
                        ans.pop();
                        outfile.write((char *)&id, sizeof(int));
                    }
                }
                outfile.close();
            }
        }
    };

    class TreeNode
    {
    public:
        int node_id;
        int lbound, rbound;
        int depth;
        std::vector<TreeNode *> childs;
        TreeNode(int l, int r, int d) : lbound(l), rbound(r), depth(d) {}
    };

    class SegmentTree
    {
    public:
        int ways_ = 2;
        TreeNode *root{nullptr};
        int max_depth{-1};
        std::vector<TreeNode *> treenodes;

        SegmentTree(int data_nb)
        {
            root = new TreeNode(0, data_nb - 1, 0);
        }

        void BuildTree(TreeNode *u)
        {
            if (u == nullptr)
                throw Exception("Tree node is a nullptr");
            treenodes.emplace_back(u);
            max_depth = std::max(max_depth, u->depth);
            int L = u->lbound, R = u->rbound;
            size_t Len = R - L + 1;
            if (L == R)
                return;
            int gap = (R - L + 1) / ways_;
            int res = (R - L + 1) % ways_;

            for (int l = L; l <= R;)
            {
                int r = l + gap - 1;
                if (res > 0)
                {
                    r++;
                    res--;
                }
                r = std::min(r, R);
                TreeNode *childnode = new TreeNode(l, r, u->depth + 1);
                u->childs.emplace_back(childnode);
                BuildTree(childnode);
                l = r + 1;
            }
        }

        std::vector<TreeNode *> range_filter(TreeNode *u, int ql, int qr)
        {
            if (u->lbound >= ql && u->rbound <= qr)
                return {u};
            std::vector<TreeNode *> res;
            if (u->lbound > qr)
                return res;
            if (u->rbound < ql)
                return res;
            for (auto child : u->childs)
            {
                auto t = range_filter(child, ql, qr);
                while (t.size())
                {
                    res.emplace_back(t.back());
                    t.pop_back();
                }
            }
            return res;
        }
    };
}