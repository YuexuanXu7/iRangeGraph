#include "utils.h"

namespace iRangeGraph_multi
{
    typedef std::pair<float, int> PFI;

    struct TwoRangeQuery
    {
        int l1, r1, l2, r2;
    };

    class DataLoader
    {
    public:
        int Dim, query_nb, query_K;
        std::vector<std::vector<float>> query_points;
        int data_nb;
        std::vector<std::vector<float>> data_points;
        std::vector<int> original_id;

        int attr_nb{0};
        std::vector<std::vector<int>> attributes;

        hnswlib::L2Space *space;

        struct Attr_Constraint
        {
            std::vector<std::pair<int, int>> attr_constraints;
        };

        std::unordered_map<std::string, std::vector<Attr_Constraint>> query_range;
        std::unordered_map<std::string, std::vector<std::vector<int>>> ground_truth;

        std::unordered_map<std::string, std::vector<std::pair<int, int>>> mapped_queryrange;

        DataLoader() {}
        ~DataLoader() {}

        float dis_compute(std::vector<float> &v1, std::vector<float> &v2)
        {
            hnswlib::DISTFUNC<float> fstdistfunc_ = space->get_dist_func();
            float dis = fstdistfunc_((char *)v1.data(), (char *)v2.data(), space->get_dist_func_param());
            return dis;
        }

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
            space = new hnswlib::L2Space(Dim);
            infile.close();
        }

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
            attributes.resize(data_nb);
            infile.close();
        }

        void LoadAttribute(std::string filename)
        {
            std::ifstream infile(filename, std::ios::in | std::ios::binary);
            if (!infile.is_open())
            {
                throw Exception("cannot open " + filename);
            }
            for (int i = 0; i < data_nb; i++)
            {
                int val;
                infile.read((char *)&val, sizeof(int));
                attributes[i].emplace_back(val);
            }
            infile.close();
            attr_nb++;
        }

        bool check_amount(std::map<std::pair<std::string, std::string>, std::vector<TwoRangeQuery>> &mp)
        {
            for (auto t : mp)
            {
                if (t.second.size() < query_nb)
                    return false;
            }
            return true;
        }

        void synthesize_2Dranges(std::string saveprefix)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::uniform_int_distribution<int> u_start(0, data_nb - 1);

            std::string savepath = saveprefix + "mixed.bin";
            CheckPath(savepath);
            std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
            if (!outfile.is_open())
            {
                throw Exception("cannot open " + savepath);
            }
            for (int i = 0; i < query_nb; i++)
            {
                int l1 = u_start(e);
                int r1 = u_start(e);
                int l2 = u_start(e);
                int r2 = u_start(e);
                if (l1 > r1)
                    std::swap(l1, r1);
                if (l2 > r2)
                    std::swap(l2, r2);

                outfile.write((char *)&l1, sizeof(int));
                outfile.write((char *)&r1, sizeof(int));
                outfile.write((char *)&l2, sizeof(int));
                outfile.write((char *)&r2, sizeof(int));
            }
            outfile.close();
        }

        void LoadRanges(std::string saveprefix)
        {
            std::string savepath = saveprefix + "mixed.bin";
            std::ifstream infile(savepath, std::ios::in | std::ios::binary);
            if (!infile.is_open())
            {
                throw Exception("cannot open " + savepath);
            }
            for (int k = 0; k < query_nb; k++)
            {
                int l1, r1, l2, r2;
                infile.read((char *)&l1, sizeof(int));
                infile.read((char *)&r1, sizeof(int));
                infile.read((char *)&l2, sizeof(int));
                infile.read((char *)&r2, sizeof(int));
                Attr_Constraint t;
                t.attr_constraints.emplace_back(l1, r1);
                t.attr_constraints.emplace_back(l2, r2);
                query_range["mixed"].emplace_back(t);
            }
            infile.close();
        }

        void Generate_Groundtruth(std::string saveprefix)
        {
            for (auto t : query_range)
            {
                std::string domain = t.first;
                std::string savepath = saveprefix + domain + ".bin";
                CheckPath(savepath);
                std::ofstream outfile(savepath, std::ios::out | std::ios::binary);
                if (!outfile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }
                for (int i = 0; i < query_nb; i++)
                {
                    auto constraints = t.second[i];
                    std::priority_queue<std::pair<float, int>> ans;
                    for (int pid = 0; pid < data_nb; pid++)
                    {
                        bool flag = true;
                        for (int j = 0; j < attr_nb; j++)
                        {
                            int ql = constraints.attr_constraints[j].first, qr = constraints.attr_constraints[j].second;
                            if (attributes[pid][j] < ql || attributes[pid][j] > qr)
                                flag = false;
                        }
                        if (!flag)
                            continue;
                        float dis = dis_compute(query_points[i], data_points[pid]);
                        ans.emplace(dis, pid);
                        if (ans.size() > query_K)
                            ans.pop();
                    }
                    while (ans.size() < query_K)
                        ans.emplace(0, -1);
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

        void LoadGroundtruth(std::string saveprefix)
        {
            for (auto t : query_range)
            {
                std::string domain = t.first;
                std::string savepath = saveprefix + domain + ".bin";
                std::ifstream infile(savepath, std::ios::in | std::ios::binary);
                if (!infile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }
                ground_truth[domain].resize(query_nb);
                for (int qid = 0; qid < query_nb; qid++)
                {
                    ground_truth[domain][qid].resize(query_K);
                    infile.read((char *)ground_truth[domain][qid].data(), query_K * sizeof(int));
                }
                infile.close();
            }
        }

        void Sort_by_Attr(int aid)
        {
            if (aid >= attr_nb)
                throw Exception("out of attribute number limit");
            std::vector<std::pair<int, int>> p;
            for (int i = 0; i < data_nb; i++)
            {
                p.push_back({attributes[i][aid], i});
            }
            sort(p.begin(), p.end());
            std::vector<std::vector<float>> vec_tmp;
            vec_tmp.resize(data_nb);
            original_id.resize(data_nb);

            for (int i = 0; i < data_nb; i++)
            {
                int pid = p[i].second;
                original_id[i] = pid;
                vec_tmp[i] = data_points[pid];
            }
            std::swap(data_points, vec_tmp);
            for (auto t : query_range)
            {
                std::string domain = t.first;
                for (int qid = 0; qid < query_nb; qid++)
                {
                    int ql = t.second[qid].attr_constraints[aid].first;
                    int qr = t.second[qid].attr_constraints[aid].second;
                    int l_bound = std::lower_bound(p.begin(), p.end(), std::make_pair(ql, -1)) - p.begin();
                    int r_bound = std::upper_bound(p.begin(), p.end(), std::make_pair(qr, data_nb)) - p.begin() - 1;
                    mapped_queryrange[domain].emplace_back(l_bound, r_bound);
                }
            }
            std::cout << "sorted data points by " << aid << "th attribute" << std::endl;
        }
    };
}
