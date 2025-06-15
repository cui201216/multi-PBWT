#include <iostream>
#include"multiPBWT.h"
// TIP 要<b>Run</b>代码，请按 <shortcut actionId="Run"/> 或点击装订区域中的 <icon src="AllIcons.Actions.Execute"/> 图标。
int main()
{
    string panel = "out";
    multiPBWT mpbwt;
    mpbwt.readMacsPanel(panel);
    int a=0;
    while (mpbwt.X[a][0]==0){a++;}
    int c=a;
    mpbwt.makePanel();
    return 0;
}