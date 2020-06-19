#include<iostream>
#include<vector>
#include<string>
#include<random>
#include<time.h>
#include<math.h>
#include<Rcpp.h>

struct counts {
    int c_and_p, c_and_no_p, no_c_and_p;
};

std::vector<std::vector<float> > get_local_scores(std::vector<std::vector<bool> > &DataFrame, int k, float log_epsilon, std::string model);
std::vector<std::vector<int> > get_best_parents(std::vector<std::vector<float> > &best_scores, int p);
std::vector<int> get_best_sinks(int p, std::vector<std::vector<int> > &best_parents, std::vector<std::vector<float> > &local_scores, std::vector<float> &scores);
std::vector<int> get_best_ordering(int p, std::vector<int> &sinks);
std::vector<std::vector<int> > get_best_network(int p, std::vector<int> &ordering, std::vector<std::vector<int> > &best_parents);

float infer_theta(counts count);
counts get_counts(int v, std::vector<int> S, std::vector<std::vector<bool> > &DataFrame, std::string model);
float local_score(int v, int s, std::vector<std::vector<bool> > &DataFrame, float log_epsilon, std::string model);
std::vector<int> int_to_subset(int v, int s, int p);
bool penalty_check(int s, int k, int p);
std::vector<int> get_candidate_sinks(int s, int p);
std::vector<int> get_parents(int v, int s, int p);
int reshifter(int v, int s, int p);


// [[Rcpp::export]]
Rcpp::List dp(Rcpp::List input) {
    Rcpp::Rcout << "In and didnt die \n";
    std::vector<int> dims = input["dims"];
    int N = dims[0];
    int p = dims[1];
    std::vector<bool> M = input["df"];
    std::vector<std::vector<bool> > DataFrame(N, std::vector<bool>(p, 0));
    std::vector<float> options = input["options"];
    int k = (int)options[0];
    std::string model = input["model"];
    bool verbose = input["verbose"];
    float epsilon = options[1];
    Rcpp::List out = Rcpp::List::create();

    //error catching
    if(k < 1) {
        Rcpp::stop("in-degree bound (k) must be positive.");
    }
    if(epsilon < 0.0 || epsilon > 1.0) {
        Rcpp::stop("Penalty must be in (0,1).");
    }

    float logepsilon = log(epsilon);

    if(verbose) {Rcpp::Rcout << "Loading data ... ... \n";}
    int it = 0;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < p; ++j) {
            DataFrame[i][j] = M[it];
            ++it;
        }
    }
    M.erase(M.begin(), M.end());

    clock_t start, end;
    start = clock();

    if(verbose) {Rcpp::Rcout << "Computing local scores ... ... \n";}
    std::vector<std::vector<float> > local_scores = get_local_scores(DataFrame, k, logepsilon, model);
    if(verbose) {Rcpp::Rcout << "Computing best parent sets ... ... \n";}
    std::vector<std::vector<int> > best_parents = get_best_parents(local_scores, p);

    if(verbose) {Rcpp::Rcout << "Finding best sinks ... ...\n";}
    std::vector<float> scores(1 << p, 0.0);
    std::vector<int> sinks = get_best_sinks(p, best_parents, local_scores, scores);
    Rcpp::NumericVector bs = Rcpp::wrap(scores[(1 << p) - 1]);

    if(verbose) {Rcpp::Rcout << "Finding best ordering ... ... \n";}
    std::vector<int> ordering = get_best_ordering(p, sinks);
    if(verbose) {Rcpp::Rcout << "Building optimal network ... ... \n";}
    std::vector<std::vector<int> > parents = get_best_network(p, ordering, best_parents);

    end = clock();
    if(verbose) {Rcpp::Rcout << "Optimal network with score " << bs[0] << " found in " << (float)(end - start)/CLOCKS_PER_SEC << " seconds.\n";}

    Rcpp::IntegerVector pi;
    for(int i = 0; i < p; ++i) {
        if(parents[i].size() == 0) {
            pi.push_back(p + 1);
            pi.push_back(i + 1);
        }
        for(int j = 0; j < parents[i].size(); ++j) {
            pi.push_back(parents[i][j] + 1);
            pi.push_back(i + 1);
        }
    }
    out.push_back(pi);
    out.push_back(bs);
    return out;
}

float infer_theta(counts count) {
    return (float)count.c_and_p/(float)(count.c_and_p + count.no_c_and_p);
}

counts get_counts(int v, std::vector<int> S, std::vector<std::vector<bool> > &DataFrame, std::string model) {
    int c_and_p = 0;
    int c_and_no_p = 0;
    int no_c_and_p = 0;
    bool is_cbn = false;
    if(model == "CBN") {
        is_cbn = true;
    }
    for(int i = 0; i < DataFrame.size(); ++i) {
        if(DataFrame[i][v] == 1) {
            if(S.size() == 0) {
                ++c_and_p;
            }
            else {
                for(int j = 0; j < S.size(); ++j) {
                    if(DataFrame[i][S[j]] == 1 && !is_cbn) {
                        ++c_and_p;
                        break;
                    }
                    else if(DataFrame[i][S[j]] == 0 && is_cbn) {
                        ++c_and_no_p;
                        break;
                    }
                    else if(j == (int)S.size() - 1) {
                        if(is_cbn) {++c_and_p;}
                        else if(!is_cbn) {++c_and_no_p;}
                    }
                }
            }
        }
        else {
            if(S.size() == 0) {
                ++no_c_and_p;
            }
            else {
                for(int j = 0; j < S.size(); ++j) {
                    if(DataFrame[i][S[j]] == 1 && !is_cbn) {
                        ++no_c_and_p;
                        break;
                    }
                    else if(DataFrame[i][S[j]] == 0 && is_cbn) {
                        break;
                    }
                    else if(j == ((int)S.size() - 1) && is_cbn) {
                        ++no_c_and_p;
                    }
                }
            }
        }
    }
    counts count;
    count.c_and_p = c_and_p;
    count.c_and_no_p = c_and_no_p;
    count.no_c_and_p = no_c_and_p;
    return count; 
}

float local_score(int v, int s, std::vector<std::vector<bool> > &DataFrame, float log_epsilon, std::string model) {
    std::vector<int> S = int_to_subset(v, s, DataFrame[0].size());
    counts count = get_counts(v, S, DataFrame, model);

    if(count.c_and_p == 0 || count.no_c_and_p == 0) {
        return -std::numeric_limits<float>::infinity();
    }

    float theta = infer_theta(count);
    return (count.c_and_p)*log(theta) + (count.no_c_and_p)*log(1-theta) + (count.c_and_no_p)*log_epsilon;
}

std::vector<std::vector<float> > get_local_scores(std::vector<std::vector<bool> > &DataFrame, int k, float log_epsilon, std::string model) {
    std::vector<std::vector<float > > local_scores(DataFrame[0].size(), std::vector<float>(1 << (DataFrame[0].size() - 1), 0.0));
    for(int i = 0; i < DataFrame[0].size(); ++i) {
        for(int j = 0; j < local_scores[0].size(); ++j) {
            if(penalty_check(j, k, DataFrame[0].size())) {
                local_scores[i][j] = local_score(i, j, DataFrame, log_epsilon, model);
            }
            else {
                local_scores[i][j] = -std::numeric_limits<float>::infinity();
            }
        }
    }
    return local_scores;
}

std::vector<int> get_best_sinks(int p, std::vector<std::vector<int> > &best_parents, std::vector<std::vector<float> > &local_scores, std::vector<float> &scores) {
    std::vector<int> sinks(1 << p, -1);

    for(int i = 0; i < (1 << p); ++i) {
        std::vector<int> candidate_sinks = get_candidate_sinks(i, p);
        
        for(std::vector<int>::iterator j = candidate_sinks.begin(); j != candidate_sinks.end(); ++j) {
            int upvars = i - (1 << (p-1-*j));
            float score = scores[upvars];
            int bps_i = reshifter(*j, upvars, p);
            score += local_scores[*j][best_parents[*j][bps_i]];
            if(sinks[i] == -1 || score > scores[i]) {
                scores[i] = score;
                sinks[i] = *j;
            }
        }
    }
    return sinks;
}

std::vector<std::vector<int> > get_best_parents(std::vector<std::vector<float> > &best_scores, int p) {
    std::vector<std::vector<int> > best_parents(p, std::vector<int> (1 << (p-1), 0));

    for(int i = 0; i < p; ++i) {
        for(int j = 0; j < best_parents[0].size(); ++j) {
            best_parents[i][j] = j;
            for(int k = 0; k < p-1; ++k) {
                int a = 1 << k;
                if((a & j) != 0) {
                    if(best_scores[i][j-a] >= best_scores[i][j]) {
                        best_scores[i][j] = best_scores[i][j-a];
                        best_parents[i][j] = best_parents[i][j-a];
                    }
                }
            }
        }
    }
    return best_parents;
}

std::vector<int> get_best_ordering(int p, std::vector<int> &sinks) {
    std::vector<int> ordering(p, 0);
    int left = (1 << p) - 1;
    for(int i = p-1; i >= 0; --i) {
        ordering[i] = sinks[left];
        left -= (1 << (p-1-ordering[i]));
    }
    return(ordering);
}

std::vector<std::vector<int> > get_best_network(int p, std::vector<int> &ordering, std::vector<std::vector<int> > &best_parents) {
    int predecs = 0;
    std::vector<std::vector<int> > parents;
    for(int i = 0; i < p; ++i) {
        std::vector<int> empty_vec;
        parents.push_back(empty_vec);
    }

    for(std::vector<int>::iterator i = ordering.begin(); i != ordering.end(); ++i) {
        int predecs_index = reshifter(*i, predecs, p);
        int pi = best_parents[*i][predecs_index];
        predecs += (1 << (p-1-*i));

        std::vector<int> pi_v = get_parents(*i, pi, p);
        for(std::vector<int>::iterator j = pi_v.begin(); j != pi_v.end(); ++j) {
            parents[*i].push_back(*j);
        }
    }
    return parents;
}

std::vector<int> int_to_subset(int v, int s, int p) {
    std::vector<bool> S_bool(p - 1, 0); 
    for(int i = p-2; i >= 0; --i) {
        int a = 1 << i;
        if((a & s) != 0) {
            S_bool[p-2-i] = 1;
        }
    }
    std::vector<int> S; 
    for(int i = 0; i < p-1; ++i) {
        if(S_bool[i] == 1 && i < v) {
            S.push_back(i);
        }
        else if(S_bool[i] == 1) {
            S.push_back(i + 1);
        }
    }
    return S;
}

bool penalty_check(int s, int k, int p) {
    int on = 0;
    for(int i = 0; i < p-1; ++i) {
        int a = 1 << i;
        if((a & s) != 0) {
            ++on;
        }
        if(on == k+1) {
            return 0;
        }
    }
    return 1;
}

std::vector<int> get_candidate_sinks(int s, int p) {
    std::vector<int> candidate_sinks;

    for(int i = 0; i < p; ++i) {
        int a = 1 << i;
        if((a & s) != 0) {
            candidate_sinks.push_back(p-1-i);
        }
    }

    return candidate_sinks;
}

std::vector<int> get_parents(int v, int s, int p) {
    std::vector<int> parents; 
    for(int i = 0; i < p-1; ++i) {
        int a = 1 << i;
        if((a & s) != 0) {
            if(i < p-1-v) {
                parents.push_back(p-1-i);
            }
            else {
                parents.push_back(p-2-i);
            }
        }
    }
    return parents;
}

int reshifter(int v, int s, int p) {
    int val = 0;
    for(int i = 0; i < p-v-1; ++i) {
        int a = 1 << i;
        if((a & s) != 0) {
            val += a;
        }
    }
    for(int i = p-v; i < p; ++i) {
        int a = 1 << i;
        if((a & s) != 0) {
            val += (1 << (i-1));
        }
    }
    return val;
}