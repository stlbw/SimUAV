#include <iostream>
#include "declaredFun.h"
int main() {
    std::cout << "Hello, World!" << std::endl;
    getDba("test.db"); // function implemented by Mario
    AeroDB db0; // create db object
    db0 = readData("dba_2000.ini"); // Open database and print it to the screen

    return 0;
}
