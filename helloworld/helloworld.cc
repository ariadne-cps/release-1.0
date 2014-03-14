#include <ariadne/ariadne.h>

using namespace Ariadne;

int main() {
    std::cout << "This is Ariadne!" << std::endl;
#ifdef HAVE_BUDDY_H
    std::cout << "You have BDD support" << std::endl;
#else
    std::cout << "You don't have BDD support" << std::endl;
#endif

    Interval i1(1,2);
    Interval i2(3,4);

    std::cout << i1 << " + " << i2 << " = " << i1+i2 << std::endl;

    return 0;
}


