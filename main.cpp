#include <iostream>
#include"multiPBWT.h"
// TIP 要<b>Run</b>代码，请按 <shortcut actionId="Run"/> 或点击装订区域中的 <icon src="AllIcons.Actions.Execute"/> 图标。
int main()
{
    string panel = "sites.txt";
    multiPBWT mpbwt;
    int a= mpbwt.readMacsPanel(panel);
    std::cout << "read panel file done: " << a << endl;
    int b =mpbwt.makePanel();
    std::cout << "make panel file done: " << b << endl;
    int c= mpbwt.inPanelLongMatchQuery(100,"outputFile");
    std::cout << "query done: " << c << endl;
    return 0;
}