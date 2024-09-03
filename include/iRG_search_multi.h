#include "utils_multi.h"

namespace iRangeGraph_multi
{
    typedef unsigned int tableint;
    typedef unsigned int linklistsizeint;

    typedef std::pair<float, int> PFI;
    typedef std::pair<float, std::pair<int, int>> PFII;

    template <typename dist_t>
    class iRangeGraph_Search_Multi
    {
    public:
        int MaxStep{20};
        std::vector<double> probability;
        iRangeGraph::SegmentTree *tree;
        DataLoader *storage;
        size_t max_elements_{0};
        size_t dim_{0};
        size_t M_out{0};
        size_t ef_construction{0};

        size_t size_data_per_element_{0};
        size_t size_links_per_element_{0};
        size_t data_size_{0};

        size_t size_links_per_layer_{0};
        size_t offsetData_{0};

        char *data_memory_{nullptr};

        hnswlib::L2Space *space;
        hnswlib::DISTFUNC<dist_t> fstdistfunc_;
        void *dist_func_param_{nullptr};

        size_t metric_distance_computations{0};
        size_t metric_hops{0};

        std::vector<int> visitedpool;
        size_t visited_tag{0};

        // purepost = True -> p=1   purepost =  False -> 0<=p<=1
        bool purepost{true};

        iRangeGraph_Search_Multi(std::string edgefilename, DataLoader *store, int M) : storage(store)
        {
            std::ifstream edgefile(edgefilename, std::ios::in | std::ios::binary);
            if (!edgefile.is_open())
                throw Exception("cannot open " + edgefilename);

            max_elements_ = storage->data_nb;
            dim_ = storage->Dim;
            tree = new iRangeGraph::SegmentTree(max_elements_);
            tree->BuildTree(tree->root);

            space = new hnswlib::L2Space(dim_);
            fstdistfunc_ = space->get_dist_func();
            dist_func_param_ = space->get_dist_func_param();
            M_out = M;
            visitedpool.resize(max_elements_);

            data_size_ = dim_ * sizeof(float);
            size_links_per_layer_ = M_out * sizeof(tableint) + sizeof(linklistsizeint);
            size_links_per_element_ = size_links_per_layer_ * (tree->max_depth + 1);
            size_data_per_element_ = size_links_per_element_ + data_size_;
            offsetData_ = size_links_per_element_;

            data_memory_ = (char *)malloc(max_elements_ * size_data_per_element_);
            if (data_memory_ == nullptr)
                throw std::runtime_error("Not enough memory");

            for (int pid = 0; pid < max_elements_; pid++)
            {
                for (int layer = 0; layer <= tree->max_depth; layer++)
                {
                    linklistsizeint *data = get_linklist(pid, layer);
                    edgefile.read((char *)data, sizeof(tableint));
                    int size = getListCount(data);
                    if (size > M_out)
                        throw Exception("real linklist size is bigger than defined M_out");
                    for (int i = 0; i < size; i++)
                    {
                        char *current_neighbor_ = (char *)(data + 1 + i);
                        edgefile.read(current_neighbor_, sizeof(tableint));
                    }
                }

                size_t size_in_bytes = dim_ * sizeof(float);
                char *data = getDataByInternalId(pid);
                std::memcpy(data, reinterpret_cast<char *>(storage->data_points[pid].data()), size_in_bytes);
            }
            edgefile.close();
        }

        ~iRangeGraph_Search_Multi()
        {
            free(data_memory_);
            data_memory_ = nullptr;
        }

        void setprob()
        {
            probability.resize(MaxStep);

            for (int x = 0; x < MaxStep; x++)
            {
                probability[x] = 1 / (1 + exp(x));
            }
        }

        int ProbFunc(int x)
        {
            if (purepost)
                return 1;
            if (x >= MaxStep)
                return 0;
            double randNum = (double)rand() / RAND_MAX;
            return randNum < probability[x] ? 1 : 0;
        }

        inline char *getDataByInternalId(tableint internal_id) const
        {
            return (data_memory_ + internal_id * size_data_per_element_ + offsetData_);
        }

        linklistsizeint *get_linklist(tableint internal_id, int layer) const
        {
            return (linklistsizeint *)(data_memory_ + internal_id * size_data_per_element_ + layer * size_links_per_layer_);
        }

        int getListCount(linklistsizeint *ptr) const
        {
            return *((int *)ptr);
        }

        int GetOverLap(int l, int r, int ql, int qr)
        {
            int L = std::max(l, ql);
            int R = std::min(r, qr);
            return R - L + 1;
        }

        inline bool CheckInQueryRange(int pid, std::vector<std::pair<int, int>> &queryrange)
        {
            int originalId = storage->original_id[pid];
            for (int i = 0; i < storage->attr_nb; i++)
            {
                int val = storage->attributes[originalId][i];
                if (val < queryrange[i].first || val > queryrange[i].second)
                    return false;
            }
            return true;
        }

        std::vector<std::pair<tableint, bool>> SelectEdge(int pid, int ql, int qr, int edge_limit, std::vector<std::pair<int, int>> &queryrange, int current_step)
        {
            iRangeGraph::TreeNode *cur_node = nullptr, *nxt_node = tree->root;
            std::vector<std::pair<tableint, bool>> selected_edges;
            do
            {
                cur_node = nxt_node;
                bool contain = false;
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); i++)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, ql, qr) == GetOverLap(nxt_node->lbound, nxt_node->rbound, ql, qr))
                        {
                            cur_node = nxt_node;
                            contain = true;
                        }
                    }
                } while (contain);

                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);

                for (size_t j = 1; j <= size; j++)
                {
                    int neighborId = *(data + j);
                    if (neighborId < ql || neighborId > qr)
                        continue;
                    int prob = 1;
                    int next_step = current_step + 1;
                    bool inrange = CheckInQueryRange(neighborId, queryrange);
                    if (!inrange)
                        prob = ProbFunc(next_step);
                    if (!prob)
                        continue;

                    selected_edges.emplace_back(neighborId, inrange);
                    if (selected_edges.size() == edge_limit)
                        return selected_edges;
                }

            } while (cur_node->lbound < ql || cur_node->rbound > qr);
            return selected_edges;
        }

        std::priority_queue<PFI> TopDown_search(const void *query_data, int ef, int query_k, int QL, int QR, int edge_limit, std::vector<std::pair<int, int>> queryrange, std::vector<iRangeGraph::TreeNode *> &filterednodes)
        {
            visited_tag++;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFII, std::vector<PFII>, std::greater<PFII>> candidate_set;
            std::priority_queue<PFI> top_candidates;

            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                visitedpool[pid] = visited_tag;
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                candidate_set.emplace(std::make_pair(dis, std::make_pair(pid, -1)));
                if (CheckInQueryRange(pid, queryrange))
                    top_candidates.emplace(dis, storage->original_id[pid]);
            }

            float lowerBound = std::numeric_limits<float>::max();

            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                metric_hops++;
                if (current_point_pair.first > lowerBound)
                {
                    break;
                }
                candidate_set.pop();
                int current_pid = current_point_pair.second.first;
                int current_step = current_point_pair.second.second;
                auto selected_edges = SelectEdge(current_pid, QL, QR, edge_limit, queryrange, current_step);

                while (selected_edges.size())
                {
                    auto neighbor_pair = selected_edges.back();
                    selected_edges.pop_back();
                    int neighbor_id = neighbor_pair.first;
                    bool inrange = neighbor_pair.second;

                    if (visitedpool[neighbor_id] == visited_tag)
                        continue;
                    visitedpool[neighbor_id] = visited_tag;
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    metric_distance_computations++;

                    if (top_candidates.size() < ef || dis < lowerBound)
                    {
                        int next_step = current_step + 1;
                        if (inrange)
                        {
                            top_candidates.emplace(dis, storage->original_id[neighbor_id]);
                            next_step = -1;
                        }
                        candidate_set.emplace(std::make_pair(dis, std::make_pair(neighbor_id, next_step)));

                        if (top_candidates.size() > ef)
                            top_candidates.pop();
                        if (top_candidates.size())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }

            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }

        void search(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit = 32)
        {
            for (auto range : storage->query_range)
            {
                std::string domain = range.first;
                std::vector<std::vector<int>> &gt = storage->ground_truth[domain];
                std::string savepath = saveprefix + domain + ".csv";
                CheckPath(savepath);

                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }

                std::vector<int> HOP;
                std::vector<int> DCO;
                std::vector<float> QPS;
                std::vector<float> RECALL;

                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {
                        auto cons = range.second[i];
                        int ql = storage->mapped_queryrange[domain][i].first;
                        int qr = storage->mapped_queryrange[domain][i].second;

                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        auto filterednodes = tree->range_filter(tree->root, ql, qr);
                        auto res = TopDown_search(storage->query_points[i].data(), ef, storage->query_K, ql, qr, edge_limit, cons.attr_constraints, filterednodes);
                        gettimeofday(&t2, NULL);
                        searchtime += GetTime(t1, t2);

                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;

                    HOP.emplace_back(hop);
                    DCO.emplace_back(dco);
                    QPS.emplace_back(qps);
                    RECALL.emplace_back(recall);
                }

                for (int i = 0; i < RECALL.size(); i++)
                {
                    outfile << SearchEF[i] << "," << RECALL[i] << "," << QPS[i] << "," << DCO[i] << "," << HOP[i] << std::endl;
                }
                outfile.close();
            }
        }
    };
}