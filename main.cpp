#include <iostream>
#include"multiPBWT.h"
// TIP 要<b>Run</b>代码，请按 <shortcut actionId="Run"/> 或点击装订区域中的 <icon src="AllIcons.Actions.Execute"/> 图标。
int main(int argc, char* argv[]) {
    std::string panel = "sites.txt";  // 默认输入文件名
    int queryLength = 100;            // 默认查询长度
    std::string outputFile = "outputFile";  // 默认输出文件名

    // 解析命令行参数
    int opt;
    while ((opt = getopt(argc, argv, "i:L:o:")) != -1) {
        switch (opt) {
        case 'i':
        case 'I':  // 不区分大小写
            panel = optarg;
            break;
        case 'L':
        case 'l':
            queryLength = std::stoi(optarg);
            break;
        case 'o':
        case 'O':
            outputFile = optarg;
            break;
        default:
            std::cerr << "用法: " << argv[0]
                      << " -i <input_file> -L <query_length> -o <output_file>" << std::endl;
            return 1;
        }
    }

    // 验证必要参数
    if (panel.empty() || queryLength <= 0 || outputFile.empty()) {
        std::cerr << "错误: 请提供有效的输入文件、查询长度和输出文件名。" << std::endl;
        return 1;
    }

    // 执行程序逻辑
    multiPBWT mpbwt;
    int a = mpbwt.readMacsPanel(panel);
    std::cout << "读取面板文件完成: " << a << std::endl;

    int b = mpbwt.makePanel();
    std::cout << "生成面板完成: " << b << std::endl;

    int c = mpbwt.inPanelLongMatchQuery(queryLength, outputFile);
    std::cout << "查询完成: " << c << std::endl;

    return 0;
}