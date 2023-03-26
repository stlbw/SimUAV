#include <iostream>
#include "declaredFun.h"
int main() {
    std::cout << "Hello, World!" << std::endl;
    getDba("dba.ini");
    openFile("battery.ini");
    return 0;
}
