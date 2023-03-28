#include <iostream>
#include "declaredFun.h"
int main() {

    // getDba("dba.ini"); // function implemented by Mario
    AeroDB db0; // create db object
    db0 = readData("dba_2000.ini"); // Open database and print it to the screen

    return 0;
}
