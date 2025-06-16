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

    int readMacsQuery(string txt_file);

    int makePanel();
    
    int inPanelLongMatchQuery(int L, string inPanelOutput_file);

    int outPanelLongMatchQuery(int L, string outPanelOutput_file);

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

int multiPBWT::readMacsQuery(string txt_file)
{
    clock_t start, end;
    start = clock();

    // 打开查询文件
    std::ifstream in(txt_file);
    if (in.fail()) {
        std::cerr << "无法打开查询文件: " << txt_file << std::endl;
        return 1;
    }

    std::string line;

    // Step 1: 计算查询单倍型数 (Q)
    Q = 0;
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
            Q = tokens[4].size();
            break;
        }
    }
    if (!found_site) {
        std::cerr << "未找到SITE行" << std::endl;
        return 2;
    }
    if (Q < 1) {
        std::cerr << "无效的Q: " << Q << std::endl;
        return 3;
    }
    std::cerr << "Q = " << Q << std::endl;

    // Step 2: 设置查询 IDs
    qIDs.resize(Q);
    for (int i = 0; i < Q; i++) {
        qIDs[i] = std::to_string(i);
    }

    // Step 3: 计算位点数并验证
    int query_N = 0;
    in.clear();
    in.seekg(0);
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) == 0) {
            query_N++;
        }
    }
    if (query_N < 1) {
        std::cerr << "无效的查询位点数: " << query_N << std::endl;
        return 4;
    }
    if (N > 0 && query_N != N) {
        std::cerr << "查询位点数 " << query_N << " 与面板位点数 " << N << " 不匹配" << std::endl;
        return 5;
    }
    if (N == 0) {
        N = query_N; // 若未调用 readMacsPanel，设置 N
    }

    // Step 4: 初始化数据结构
    try {
        Z.resize(Q, std::vector<uint8_t>(N));
    } catch (const std::bad_alloc& e) {
        std::cerr << "内存分配失败: " << e.what() << std::endl;
        return -1;
    }

    // Step 5: 处理 SITE 行
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

        if (token.size() != Q) {
            std::cerr << "查询单倍型数据长度不匹配: 预期 " << Q << ", 实际 " << token.size() << ", K=" << K << std::endl;
            return 6;
        }

        int index = 0;
        for (char c : token) {
            if (index >= Q) {
                std::cerr << "索引越界: index=" << index << ", Q=" << Q << ", K=" << K << std::endl;
                return 9;
            }
            int site = c - '0';
            if (site < 0 || site > 9) {
                std::cerr << "无效的位点值: '" << c << "' 在 K=" << K << ", index=" << index << std::endl;
                return 7;
            }
            Z[index][K] = site;
            index++;
        }
        if (index != Q) {
            std::cerr << "处理了 " << index << " 个查询单倍型，预期 " << Q << ", K=" << K << std::endl;
            return 9;
        }

        K++;
    }

    if (K != N) {
        std::cerr << "处理了 " << K << " 个位点，预期 " << N << std::endl;
        return 10;
    }


    // 记录时间
    end = clock();
    readQuerytime = ((double)(end - start)) / CLOCKS_PER_SEC;

    in.close();
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

int multiPBWT::outPanelLongMatchQuery(int L, string outPanelOutput_file)
{
    ofstream out(outPanelOutput_file);
    if (out.fail())
        return 2;

    vector<int> dZ(M); //match start
    dZ.shrink_to_fit();
    vector<int> fakeLocation(N + 1); //fake location
    fakeLocation.shrink_to_fit();
    vector<int> Zdivergence(N + 2); //divergence of Z	因为要从n到0计算Zdivergence和belowZdivergence时要用到dZ[n+1]=n
    Zdivergence.shrink_to_fit();
    vector<int> belowZdivergence(N + 2); //divergence of Z.below
    belowZdivergence.shrink_to_fit();
    for (int q = 0; q < Q; q++) {
        fill(dZ.begin(), dZ.end(), 0);
        fill(fakeLocation.begin(), fakeLocation.end(), 0);
        fill(Zdivergence.begin(), Zdivergence.end(), 0);
        fill(belowZdivergence.begin(), belowZdivergence.end(), 0);

        string qID = qIDs[q >> 1] + "-" + to_string(q & 1);
        fakeLocation[0] = 0;

        //fake location
        for (int k = 0; k < N; k++) {

            int site = Z[q][k];
            if (fakeLocation[k] != M) {
                fakeLocation[k + 1] =
                    u[k * (M * t) + fakeLocation[k] * t + site];
          } else {
                if (site < t - 1) {
                    fakeLocation[k + 1] =u[k * (M * t) + 0 * t + site+1];
                } else if (site == t - 1) {
                    fakeLocation[k + 1] = M;
                } else {
                    return 3;
                }
            }
        }

        Zdivergence[N + 1] = belowZdivergence[N + 1] = N;
        for (int k = N; k >= 0; --k) {
            Zdivergence[k] = std::min(Zdivergence[k + 1], k);
            belowZdivergence[k] = std::min(belowZdivergence[k + 1], k);
            if (fakeLocation[k] != 0) {
                int index = array[k][fakeLocation[k] - 1];
                int panelSite=X[index][(Zdivergence[k] - 1)];
                int querySite=Z[q][(Zdivergence[k] - 1)];

                //向前更新Zdivergence
                while (Zdivergence[k] > 0 && panelSite == querySite) {
                    --Zdivergence[k];
                   panelSite=X[index][(Zdivergence[k] - 1)];
                    querySite=Z[q][(Zdivergence[k] - 1)];
                }
            } else {
                //t[k]==0
                Zdivergence[k] = k;
            }
            if (fakeLocation[k] < M) {
                int index = array[k][fakeLocation[k]]; //hapolotype below query
                int panelSite = X[index][belowZdivergence[k] - 1];
                int querySite = Z[q][belowZdivergence[k]-1];
                //向前更新belowZdivergence
                while (belowZdivergence[k] > 0 && panelSite == querySite) {
                    belowZdivergence[k]--;
                    panelSite = X[index][belowZdivergence[k] - 1];
                    querySite = Z[q][belowZdivergence[k]-1];
                }
            } else {
                belowZdivergence[k] = k;
            }
        }

        int f, g;
        f = g = fakeLocation[0];
        vector<int> ftemp, gtemp;
        ftemp.resize(t);
        gtemp.resize(t);

        for (int k = 0; k < N; k++) {
            int querySite = Z[q][k];
            if (g == M) {
                if (f == M) {
                    //update ftemp
                    for (int i = 0; i < t; i++) {
                        if (querySite != i) {
                            if (i != t - 1) {
                                ftemp[i] = u[k * (M * t) + 0 * t + i+1];
                                   // u[k][0][i + 1];
                            } else {
                                ftemp[i] = M;
                            }
                        }
                    }

                    if (querySite != t - 1) {
                        f =u[k * (M * t) + 0 * t + querySite+1];
                            //u[k][0][fuzzyQ + 1];
                    } else {
                        f = M;
                    }
                } else //f!=M
                {
                    for (int i = 0; i < t; i++) {
                        if (querySite != i) {
                            ftemp[i] =u[k * (M * t) + f * t + i];
                                //u[k][f][i];
                        }
                    }
                    f = u[k * (M * t) + f * t + querySite];
                        //u[k][f][fuzzyQ];
                }
                //update gtemp and g
                for (int i = 0; i < t; i++) {
                    if (querySite != i) {
                        if (i < t - 1) {
                            gtemp[i] =u[k * (M * t) + 0 * t + i + 1];
                                //u[k][0][i + 1];
                        } else {
                            gtemp[i] = M;
                        }
                    }
                }
                if (querySite < t - 1) {
                    g =u[k * (M * t) + 0 * t + querySite+1];
                        //u[k][0][fuzzyQ + 1];
                } else {
                    g = M;
                }
            } else //g!=M
            {
                for (int i = 0; i < t; i++) {
                    if (i != querySite) {
                        ftemp[i] = u[k * (M * t) + f * t + i];
                            //u[k][f][i];
                        gtemp[i] = u[k * (M * t) + g * t + i];
                            //u[k][g][i];
                    }
                }
                f =u[k * (M * t) + f * t + querySite];
                    //u[k][f][fuzzyQ];
                g =u[k * (M * t) + g * t + querySite];
                    //u[k][g][fuzzyQ];
            }

            //output matches
            for (int i = 0; i < t; i++) {
                if (i != querySite) {
                    while (ftemp[i] != gtemp[i]) {
                        //output Match
                        //int start = 0, end = 0;
                        int index = array[k + 1][ftemp[i]];
                        out << IDs[index] << '\t' << this->qIDs[q] << '\t' << dZ[index] << '\t'
<< k-1 << '\n';
                        ++ftemp[i];
                    }
                }
            }

            if (f == g) {
                if (k + 1 - Zdivergence[k + 1] == L) {
                    --f;
                    dZ[array[k + 1][f]] = k + 1 - L;
                    //store divergence
                }

                //if (k + 1 - belowZdivergence[k + 1] == l)
                if (k + 1 - belowZdivergence[k + 1] == L) {
                    //store divergence
                    dZ[array[k + 1][g]] = k + 1 - L;
                    ++g;
                }
            }
            if (f != g) {
                while (divergence[k + 1][f] <= k + 1 - L) {
                    --f;
                    dZ[array[k + 1][f]] = k + 1 - L;
                }
                while (g < M && divergence[k + 1][g] <= k + 1 - L) {
                    dZ[array[k + 1][g]] = k + 1 - L;

                    ++g;
                }
            }
        }

        //mathces no ending at
        while (f != g) {
            //output Match
            //	int start = 0, end = 0;
            int index = array[N][f];
            out << IDs[index] << '\t' << this->qIDs[q] << '\t' << dZ[index] << '\t'
<< N-1 << '\n';
            ++f;
        }

        std::fill(fakeLocation.begin(), fakeLocation.end(), 0);
        std::fill(Zdivergence.begin(), Zdivergence.end(), 0);
        std::fill(belowZdivergence.begin(), belowZdivergence.end(), 0);
    }

    cout << "matches has been put into " << outPanelOutput_file << endl;
    return 0;
}
