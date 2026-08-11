#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <mex.h>
using namespace std;

double eps=1e-10;

struct Vertex
{
    enum class Type
    {
        None,
        R_T,
        C_Q,
        C,
    } type;
    int id;
};
    
std::vector<Vertex> matching; // C -> R_T \cup C_Q

std::vector<int> dual;
vector<vector<double>> P;
    int n, m_Q, m_T;
    
    std::vector<bool> matched_T;  // R_T -> {true, false}
    std::vector<int>  base;       // R_Q -> C
    std::vector<int>  base_inv;   // C   -> R_Q
    
    std::vector<Vertex> queue;
    std::vector<Vertex> prev_C;
    std::vector<Vertex> prev_C_Q;
    std::vector<Vertex> prev_R_T;
    
    Vertex search_shortest_path(vector<vector<double>> &Q, vector<vector<double>> &T)
    {
        queue.resize(2*n + m_T);
        int r = 0;
        // push start vertices in R_T
        prev_R_T.assign(m_T, {Vertex::Type::None, -1});
        for(int i=0; i<m_T; ++i){
            if(!matched_T[i]){
                prev_R_T[i].id = 0;
                queue[r] = {Vertex::Type::R_T, i};
                ++r;
            }
        }
        // push start vertices in C_Q
        prev_C_Q.assign(n, {Vertex::Type::None, -1});
        for(int h=0; h<m_Q; ++h){
            if(base[h] != -1) continue;
            for(int j=0;j<n;j++){
                if(prev_C_Q[j].id == -1 && abs(P[h][j])>eps){
                    prev_C_Q[j].id = 0;
                    queue[r] = {Vertex::Type::C_Q, j};
                    ++r;
                }
            }
        }
        
        prev_C.assign(n, {Vertex::Type::None, -1});
        // BFS
        for(int l=0; l<r; ++l){
            const Vertex v = queue[l];
            
            const auto add_C = [&](const int j, const Vertex& prev) -> bool {
                prev_C[j] = prev;
                if(matching[j].type == Vertex::Type::None){
                    return true;
                }
                queue[r] = {Vertex::Type::C, j};
                ++r;
                return false;
            };
            
            if(v.type == Vertex::Type::R_T){
                for(int j=0;j<n;j++){
                    if(T[v.id][j]!=0 && prev_C[j].id == -1){
                        if(add_C(j, v)) return {Vertex::Type::C, j};
                    }
                }
            }
            else if(v.type == Vertex::Type::C_Q){
                if(prev_C[v.id].id == -1){
                    if(add_C(v.id, v)) return {Vertex::Type::C, v.id};
                }
                
                // hige
                const int h = base_inv[v.id];
                if(h != -1){
                    for(int j=0;j<n;j++){
                        if(abs(P[h][j])>eps && prev_C_Q[j].id == -1){
                            prev_C_Q[j] = v;
                            queue[r] = {Vertex::Type::C_Q, j};
                            ++r;
                        }
                    }
                }
            }
            else if(v.type == Vertex::Type::C){
                const auto type = matching[v.id].type;
                const auto id   = matching[v.id].id;
                if(type != Vertex::Type::None){
                    if(type == Vertex::Type::R_T){
                        if(prev_R_T[id].id == -1){
                            prev_R_T[id] = v;
                            queue[r] = {type, id};
                            ++r;
                        }
                    }
                    else if(type == Vertex::Type::C_Q){
                        if(prev_C_Q[id].id == -1){
                            prev_C_Q[id] = v;
                            queue[r] = {type, id};
                            ++r;
                        }
                    }
                    else assert(false);
                }
            }
            else assert(false);
        }
        
        return {Vertex::Type::None, -1};
    }
    
    void pivot(const int h, const int j, const double& piv)
    {
        assert(abs(piv)>eps);
        for(int k=0; k<m_Q; ++k){
            if(k == h) continue;
            double z=-P[k][j]/piv;
            for(int i=0;i<P[k].size();i++){
                P[k][i]+=P[h][i]*z;
            }
        }
    }
    

    pair<int,vector<int>> compute(vector<vector<double>> Q, vector<vector<double>> T)
    {
        m_Q=Q.size();
        m_T=T.size();
        n=Q[0].size();
        
        matching.assign(n, {Vertex::Type::None, -1});
        matched_T.assign(m_T, false);
        base.assign(m_Q, -1);
        base_inv.assign(n, -1);
        P = Q;
        
        // the main loop
        for(;;){
            Vertex v = search_shortest_path(Q,T);
            if(v.type == Vertex::Type::None) break;
            
            while(v.type != Vertex::Type::None){
                //assert(v.type == Vertex::Type::C);
                const Vertex u = prev_C[v.id];
                matching[v.id] = u;
                
                if(u.type == Vertex::Type::R_T){
                    matched_T[u.id] = true;
                    const Vertex t = prev_R_T[u.id];
                    if(t.type == Vertex::Type::None) break; 
                    v = t;
                }
                else if(u.type == Vertex::Type::C_Q){
                    const int j = u.id;
                    const Vertex t = prev_C_Q[j];
                    
                    if(t.type == Vertex::Type::None){
                        // add the j-th column to the base
                        for(int h=0; h<m_Q; ++h){
                            const auto piv = P[h][j];
                            if(base[h] == -1 && abs(piv)>eps){
                                base[h] = j;
                                base_inv[j] = h;
                                pivot(h, j, piv);
                                break;
                            }
                        }
                        break;
                    }
                    else if(t.type == Vertex::Type::C_Q){
                        // remove the i-th column and add the j-th column
                        const int i = t.id;
                        const int h = base_inv[i];
                        base[h] = j;
                        base_inv[i] = -1;
                        base_inv[j] = h;
                        pivot(h, j,P[h][j]);
                        v = prev_C_Q[i];
                    }
                    else assert(false);
                }
                else assert(false);
            }
        }
        
        // accumulate the result
        int rank = 0;
        for(const Vertex& v : matching){
            if(v.type != Vertex::Type::None) ++rank;
        }
        
        // compute the dual
        dual.clear();
        for(int j=0; j<n; ++j){
            if(prev_C[j].id == -1) {
                dual.push_back(j+1);
            }
        }
        return make_pair(rank,dual);
    }

void matlabArrayToCppVector(const mxArray *matlabArray, std::vector<std::vector<double>> &cppVector) {
    mwSize numRows = mxGetM(matlabArray);
    mwSize numCols = mxGetN(matlabArray);

    cppVector.resize(numRows, std::vector<double>(numCols, 0.0));

    double *matlabData = mxGetPr(matlabArray);

    for (mwIndex row = 0; row < numRows; ++row) {
        for (mwIndex col = 0; col < numCols; ++col) {
            cppVector[row][col] = matlabData[row + numRows * col];
        }
    }
}

void cppVectorToMatlabArray(const std::vector<int> &cppVector, mxArray *matlabArray) {
    mwSize numElements = cppVector.size();
    mxSetM(matlabArray, 1);
    mxSetN(matlabArray, numElements);

    double *data = mxGetPr(matlabArray);
    for (mwIndex i = 0; i < numElements; ++i) {
        data[i] = cppVector[i];
    }
}

void mexFunction(int nlhs,  mxArray *plhs[],int nrhs, const mxArray *prhs[]){
    vector<vector<double>> Q,T;
    matlabArrayToCppVector(prhs[0],Q);
    matlabArrayToCppVector(prhs[1],T);
    pair<int,vector<int>> ret=compute(Q,T);
    plhs[0]=mxCreateDoubleScalar(ret.first);
    plhs[1]=mxCreateDoubleMatrix(1, ret.second.size(), mxREAL);
    cppVectorToMatlabArray(ret.second,plhs[1]);
}