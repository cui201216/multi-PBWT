#include "multiPBWT.h" // 假设 multiPBWT 类定义在此头文件中

// 打印帮助信息
void printHelp(const char* programName) {
    std::cout << "用法: " << programName << " [选项]\n"
              << "选项:\n"
              << "  -i <file>  指定输入 MaCS 格式单倍型面板文件 (默认: sites.txt)\n"
              << "  -q <file>  指定面板外查询的 MaCS 格式查询文件 (可选，面板外查询时必须)\n"
              << "  -o <file>  指定输出文件 (默认: <输入面板文件>.out)\n"
              << "  -l <int>   指定最小匹配长度 (默认: 100)\n"
              << "  -t <type>  指定查询类型: 'in' (面板内查询) 或 'out' (面板外查询) (默认: in)\n"
              << "  -h         显示此帮助信息\n"
              << "示例:\n"
              << "  面板内查询: " << programName << " -i panel.txt -l 100 -o output.txt -t in\n"
              << "  面板外查询: " << programName << " -i panel.txt -q query.txt -l 100 -o output.txt -t out\n";
}

// 验证文件有效性
bool validateFiles(const std::string& panel, const std::string& query, const std::string& output, bool isExternalQuery) {
    // 检查面板文件
    std::ifstream panelFile(panel);
    if (!panelFile.good()) {
        std::cerr << "错误: 无法打开输入面板文件 '" << panel << "'\n";
        return false;
    }
    panelFile.close();

    // 如果是面板外查询，检查查询文件
    if (isExternalQuery && !query.empty()) {
        std::ifstream queryFile(query);
        if (!queryFile.good()) {
            std::cerr << "错误: 无法打开查询文件 '" << query << "'\n";
            return false;
        }
        queryFile.close();
    }

    // 检查输出文件是否可写
    std::ofstream outFile(output);
    if (!outFile.good()) {
        std::cerr << "错误: 无法写入输出文件 '" << output << "'\n";
        return false;
    }
    outFile.close();
    return true;
}

int main(int argc, char* argv[]) {
    std::string panel = "sites.txt";      // 默认 MaCS 输入文件
    std::string query;                    // 查询文件（面板外查询时使用）
    std::string outputFile;               // 输出文件动态生成
    int queryLength = 100;                // 默认最小匹配长度
    std::string queryType = "in";         // 默认查询类型为面板内查询

    // 打印命令行参数（用于调试）
    for (int i = 0; i < argc; ++i) {
        std::cerr << argv[i] << " ";
    }
    std::cerr << "\n";

    // 解析命令行参数
    int opt;
    while ((opt = getopt(argc, argv, "i:I:q:Q:o:O:l:L:t:T:hH")) != -1) {
        try {
            switch (opt) {
                case 'i':
                case 'I':
                    panel = optarg;
                    break;
                case 'q':
                case 'Q':
                    query = optarg;
                    break;
                case 'o':
                case 'O':
                    outputFile = optarg;
                    break;
                case 'l':
                case 'L':
                    queryLength = std::stoi(optarg);
                    break;
                case 't':
                case 'T':
                    queryType = optarg;
                    if (queryType != "in" && queryType != "out") {
                        std::cerr << "错误: 查询类型必须为 'in' 或 'out'\n";
                        return 1;
                    }
                    break;
                case 'h':
                case 'H':
                    printHelp(argv[0]);
                    return 0;
                default:
                    std::cerr << "错误: 未知选项 '" << (char)opt << "'，使用 -h 查看用法\n";
                    return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "错误: 无效的参数值 '" << optarg << "'\n";
            return 1;
        }
    }

    // 自动生成输出文件名
    if (outputFile.empty()) {
        outputFile = panel + ".out";
    }

    // 验证参数
    if (panel.empty() || queryLength <= 0) {
        std::cerr << "错误: 输入面板文件和查询长度必须有效\n";
        return 1;
    }
    if (queryType == "out" && query.empty()) {
        std::cerr << "错误: 面板外查询必须提供查询文件 (-q)\n";
        return 1;
    }
    if (!validateFiles(panel, query, outputFile, queryType == "out")) {
        return 1;
    }

    // 输出参数信息
    std::cout << "参数:\n"
              << "输入面板文件: " << panel << "\n"
              << "查询文件: " << (query.empty() ? "无（面板内查询）" : query) << "\n"
              << "输出文件: " << outputFile << "\n"
              << "查询长度: " << queryLength << "\n"
              << "查询类型: " << (queryType == "in" ? "面板内查询" : "面板外查询") << "\n";

    // 创建 PBWT 处理器
    multiPBWT haplotypeMatcher;
    int a = haplotypeMatcher.readMacsPanel(panel);
    std::cout << "读取面板: " << a << "\n";
    if (a != 0) return a;

    // 读取查询文件（仅面板外查询）
    if (queryType == "out") {
        int d = haplotypeMatcher.readMacsQuery(query);
        std::cout << "读取查询文件: " << d << "\n";
        if (d != 0) return d;
    }

    int b = haplotypeMatcher.makePanel();
    std::cout << "生成面板: " << b << "\n";
    if (b != 0) return b;

    // 根据查询类型执行查询
    int c;
    if (queryType == "in") {
        c = haplotypeMatcher.inPanelLongMatchQuery(queryLength, outputFile);
        std::cout << "面板内查询完成: " << c << "\n";
    } else {
        c = haplotypeMatcher.outPanelLongMatchQuery(queryLength, outputFile);
        std::cout << "面板外查询完成: " << c << "\n";
    }
    return c;
}