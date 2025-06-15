/*
 * FSPBWt.h
 *
 *  Created on: May 20, 2024
 *      Author: Cui Rongyue
 */

#include <chrono>
#include<iostream>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <sstream>
#include <vector>
#include<cmath>
#include <ctime>
#include <string>
#include <unistd.h>

using namespace std;
struct multiPBWT {
    int M = 0;
    int N = 0;
    int maxSite = 0;
    int t=0;
    double readPaneltime = 0;
    double makePanelTime = 0;
    double inPanelQuerytime = 0;
    double readQuerytime = 0;
    double outPanelQuerytime = 0;
    u_long inPanelMatchNum = 0;
    u_long outPanelMatchNum = 0;
    vector<string> IDs;
    vector<vector<uint8_t> > X; // MN bits
    vector<vector<int> > array; // 32MN/B bits
    vector<vector<int> > divergence; // 32MN/B bits
    int *u;

    int Q = 0;
    vector<vector<uint8_t> > Z;
    vector<string> qIDs;

    int readMacsPanel(string txt_file);

    int makePanel();
    
    int inPanelLongMatchQuery(int L, string inPanelOutput_file);

    int outPanelLongMatchQuery(int L, string outPanelOutput_file, bool even);
};

int multiPBWT::readMacsPanel(string panel_file) {
    clock_t start, end;
    start = clock();
    std::ifstream in(panel_file);
    if (in.fail()) {
        std::cerr << "无法打开文件: " << panel_file << std::endl;
        return 1;
    }

    std::string line;

    // Step 1: 计算单倍型数 (M)
    M = 0;
    bool found_site = false;
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) == 0) {
            found_site = true;
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, '\t')) {
                tokens.push_back(token);
            }
            if (tokens.size() < 5) {
                std::cerr << "SITE行格式错误: 需要至少5个字段，实际为 " << tokens.size() << std::endl;
                return 2;
            }
            M = tokens[4].size();
            break;
        }
    }
    if (!found_site) {
        std::cerr << "未找到SITE行" << std::endl;
        return 2;
    }
    if (M < 1) {
        std::cerr << "无效的M: " << M << std::endl;
        return 3;
    }
    std::cerr << "M = " << M << std::endl;

    // Step 2: 设置IDs
    IDs.resize(M);
    for (int i = 0; i < M; i++) {
        IDs[i] = std::to_string(i);
    }

    // Step 3: 计算位点数 (N)
    N = 0;
    in.clear();
    in.seekg(0);
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) == 0) {
            N++;
        }
    }
    if (N < 1) {
        std::cerr << "无效的N: " << N << std::endl;
        return 4;
    }

    // Step 4: 初始化数据结构
    try {
        X.resize(M,vector<uint8_t>(N));
        array.resize(N + 1, std::vector<int>(M));
        std::iota(array[0].begin(), array[0].end(), 0);
        divergence.resize(N + 1, std::vector<int>(M, 0));

    } catch (const std::bad_alloc& e) {
        std::cerr << "内存分配失败: " << e.what() << std::endl;
        return -1;
    }

    // Step 5: 处理SITE行
    in.clear();
    in.seekg(0);


    int K = 0;
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) != 0) {
            continue;
        }

        if (K >= N) {
            std::cerr << "SITE行数过多: K=" << K << ", 预期N=" << N << std::endl;
            return 10;
        }


        std::stringstream ss(line);
        std::string token;
        std::getline(ss, token, '\t'); // Skip "SITE:"
        std::getline(ss, token, '\t'); // Skip index
        std::getline(ss, token, '\t'); // Skip physLoc
        std::getline(ss, token, '\t'); // Skip other column
        std::getline(ss, token, '\t'); // Get haplotype data

        if (token.size() != M) {
            std::cerr << "单倍型数据长度不匹配: 预期 " << M << ", 实际 " << token.size() << ", K=" << K << std::endl;
            return 6;
        }

        int index = 0;

        for (char c : token) {
            if (index >= M) {
                std::cerr << "索引越界: index=" << index << ", M=" << M << ", K=" << K << std::endl;
                return 9;
            }
            int site = c-'0';
            if (site > maxSite)
            {
                maxSite = site;
            }
            X[index][K]=site;
            index++;
        }
        if (index != M) {
            std::cerr << "处理了 " << index << " 个单倍型，预期 " << M << ", K=" << K << std::endl;
            return 9;
        }

        K++;
    }

    if (K != N) {
        std::cerr << "处理了 " << K << " 个位点，预期 " << N << std::endl;
        return 10;
    }
    t = maxSite+1;
    u = new int[(unsigned long)N * M * t];

    end = clock();
    readPaneltime = ((double)(end - start)) / CLOCKS_PER_SEC;

    return 0;
}


int multiPBWT::makePanel() // fuzzy way :overall
{

    clock_t start, end;
    start = clock();

    int a_count[t] = {0};
    int d_count[t] = {0};
    int a[t][M];
    int d[t][M];
    int p[t];
    for (int k = 0; k < N; k++) //each position
    {
        for (int _ = 0; _ < t; _++) {
            p[_] = k + 1;
        }

        for (int i = 0; i < M; i++) //each hapolotype
        {
            //update u[][][]
            for (int _ = 0; _ < t; _++) {
                u[k * (M * t) + i * t + _]=a_count[_];
                if (divergence[k][i] > p[_]) {
                    p[_] = divergence[k][i];
                }
            }
            int index = array[k][i];

            int site = X[index][k];

            a[site][a_count[site]] = index;
            d[site][d_count[site]] = p[site];
            a_count[site]++;
            d_count[site]++;
            p[site] = 0;
        } //end j in  M
        //update array divergence u[][][]
        int m = 0;
        for (int _ = 0; _ < t; _++) {
            for (int w = 0; w < a_count[_]; w++) {
                array[k + 1][m] = a[_][w];
                divergence[k + 1][m] = d[_][w];
                ++m;
            }
        }
        if (m != M) {
            return 4;
        }
        //put plus into u for w()
        for (int _ = 1; _ < t; _++) {
            int plus = 0;
            for (int j = 0; j < _; j++) {
                plus += a_count[j];
            }
            for (int j = 0; j < M; j++) {
                u[k * (M * t) + j * t + _] += plus;
            }
        }
        for (int _ = 0; _ < t; _++) {
            a_count[_] = 0;
            d_count[_] = 0;
        }
    }
    end = clock();
    makePanelTime = ((double) (end - start)) / CLOCKS_PER_SEC;
    return 0;
}

int multiPBWT::inPanelLongMatchQuery(int L, string inPanelOutput_file) {
    clock_t start, end;
    start = clock();

    ofstream out(inPanelOutput_file);
    if (out.fail())
        return 2;

    int k;
    for (k = 0; k < N - 1; k++) {
        bool m[t];
        for (int _ = 0; _ < t; _++) {
            m[_] = false;
        }
        int top = 0;
        bool report = false;
        for (int i = 0; i < M; i++) {
            if (divergence[k][i] > k - L) {
                for (int w = 0; w < t - 1; w++) {
                    for (int v = w + 1; v < t; v++) {
                        if (m[w] == true && m[v] == true) {
                            report = true;
                            break;
                        }
                    }
                }
                if (report == true) {
                    for (int i_a = top; i_a < i - 1; i_a++) {
                        int maxDivergence = 0;
                        for (int i_b = i_a + 1; i_b < i; i_b++) {
                            if (divergence[k][i_b] > maxDivergence) {
                                maxDivergence = divergence[k][i_b];
                            }
                            uint32_t temp1, temp2, fuzzy1, fuzzy2;
                            int index_a = array[k][i_a], index_b =
                                    array[k][i_b];
                            int site1 = X[index_a][k];
                            int site2 = X[index_b][k];

                            if (site1 != site2) {
                                out << IDs[index_a] << '\t' << IDs[index_b] << '\t' << maxDivergence << '\t'
           << k-1 << '\n';
                                ++this->inPanelMatchNum;
                            }
                        }
                    }
                    report = false;
                } //end if(report==true)

                top = i;
                for (int _ = 0; _ < t; _++) {
                    m[_] = false;
                }
                //report = false;
            } // end if(divergence[k][i]>k-l)
            //change m[]
            
            int site = X[ array[k][i] ][k];

            m[site] = true;
        } //for i from 0 to M-1
        //cheak bottom block
        for (int w = 0; w < t - 1; w++) {
            for (int v = w + 1; v < t; v++) {
                if (m[w] == true && m[v] == true) {
                    report = true;
                }
            }
        }
        if (report == true) {
            for (int i_a = top; i_a < M - 1; i_a++) {
                int maxDivergence = 0;
                for (int i_b = i_a + 1; i_b < M; i_b++) {
                    if (divergence[k][i_b] > maxDivergence) {
                        maxDivergence = divergence[k][i_b];
                    }
                    uint32_t temp1, temp2, fuzzy1, fuzzy2;
                    int index_a = array[k][i_a];
                    int index_b = array[k][i_b];
                    int site1 = X[index_a][k];
                    int site2 = X[index_b][k];
                    
                    if (site1 != site2) {
                        out << IDs[index_a] << '\t' << IDs[index_b] << '\t' << maxDivergence << '\t'
            << k-1 << '\n';
                    }
                }
            }
        }
    }

    //late site     
    int top = 0;
    for (int i = 0; i < M; i++) {
        if (divergence[k][i] > k - L + 1) {
            for (int i_a = top; i_a < i - 1; i_a++) {
                int maxDivergence = 0;
                for (int i_b = i_a + 1; i_b < i; i_b++) {
                    if (divergence[k][i_b] > maxDivergence) {
                        maxDivergence = divergence[k][i_b];
                    }

                    int index_a = array[k][i_a];
                    int index_b = array[k][i_b];

                    int site1 = X[index_a][k];
                    int site2 = X[index_b][k];
                    




                    if (site1 == site2) {
                        out << IDs[index_a] << '\t' << IDs[index_b] << '\t' << maxDivergence << '\t'
           << k<< '\n';
                    } else if (site1 != site2) {
                        if (k - maxDivergence >= L) {
                            out << IDs[index_a] << '\t' << IDs[index_b] << '\t' << maxDivergence << '\t'
            << k << '\n';
                        }
                    }
                }
            }
            top = i;
        }
    }
    //bottom block
    for (int i_a = top; i_a < M - 1; i_a++) {
        int maxDivergence = 0;
        for (int i_b = i_a + 1; i_b < M; i_b++) {
            int index_a = array[k][i_a];
            int index_b = array[k][i_b];
            if (divergence[k][i_b] > maxDivergence) {
                maxDivergence = divergence[k][i_b];
            }
            out << IDs[index_a] << '\t' << IDs[index_b] << '\t' << maxDivergence << '\t'
           << N<< '\n';
        }
    }

    end = clock();

     this->inPanelQuerytime = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    return 0;
    out.close();
    cout << "matches has been put into " << inPanelOutput_file << endl;
}