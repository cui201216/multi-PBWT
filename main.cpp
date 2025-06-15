
#include <getopt.h>

// 假设 multiPBWT 定义在以下头文件中
#include "multiPBWT.h"

// 显示帮助信息
void printHelp(const char* programName) {
    std::cout << "用法: " << programName << " [选项]\n"
              << "处理 MaCS 格式单倍型数据，查询长匹配区间。\n"
              << "选项（大小写不敏感）:\n"
              << "  -i/-I <file>  输入 MaCS 文件 (默认: sites.txt)\n"
              << "  -o/-O <file>  输出文件，保存匹配区间 (默认: <input>.out)\n"
              << "  -l/-L <int>   最小匹配长度 (默认: 100)\n"
              << "  -h/-H         显示帮助\n"
              << "示例: " << programName << " -i data.txt -l 200 -o result.txt\n"
              << "输入格式: MaCS 文件，每行以 'SITE:' 开头，第五字段为单倍型数据。\n";
}

// 验证文件路径
bool validateFiles(const std::string& inputFile, const std::string& outputFile) {
    std::ifstream in(inputFile);
    if (!in.good()) {
        std::cerr << "错误: 输入文件 '" << inputFile << "' 不存在\n";
        return false;
    }
    in.close();
    std::ofstream out(outputFile, std::ios::app);
    if (!out.good()) {
        std::cerr << "错误: 无法写入输出文件 '" << outputFile << "'\n";
        return false;
    }
    out.close();
    return true;
}

// 处理 MaCS 格式单倍型数据，查询长匹配区间
int main(int argc, char* argv[]) {
    std::string panel = "sites.txt";      // 默认 MaCS 输入文件
    std::string outputFile;               // 输出文件动态生成
    int queryLength = 100;                // 默认最小匹配长度

    for (int i = 0; i < argc; ++i) {
        std::cerr << argv[i] << " ";
    }
    std::cerr << "\n";

    // 解析命令行参数
    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:l:L:hH")) != -1) {
        try {
            switch (opt) {
                case 'i':
                case 'I':
                    panel = optarg;
                    break;
                case 'o':
                case 'O':
                    outputFile = optarg;
                    break;
                case 'l':
                case 'L':
                    queryLength = std::stoi(optarg);
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
        std::cerr << "错误: 输入文件和查询长度必须有效\n";
        return 1;
    }
    if (!validateFiles(panel, outputFile)) {
        return 1;
    }

    // 输出参数信息
    std::cout << "参数:\n"
              << "输入文件: " << panel << "\n"
              << "输出文件: " << outputFile << "\n"
              << "查询长度: " << queryLength << "\n";

    // 创建 PBWT 处理器
    multiPBWT haplotypeMatcher;
    int a = haplotypeMatcher.readMacsPanel(panel);
    std::cout << "读取面板: " << a << "\n";
    if (a != 0) return a;

    int b = haplotypeMatcher.makePanel();
    std::cout << "生成面板: " << b << "\n";
    if (a != 0) return b;

    int c = haplotypeMatcher.inPanelLongMatchQuery(queryLength, outputFile);
    std::cout << "查询完成: " << c << "\n";
    return c;
}