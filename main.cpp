#include <iostream>
#include "strassen.h"

int main() {
    matrix_complex m1;
    std::vector<complex> row1 = {complex(10, 2), complex(5, 3), complex(6, 1), complex(1, 1)};
    std::vector<complex> row2 = {complex(5, 2), complex(2, 3), complex(10, 5), complex(20, 2)};
    std::vector<complex> row3 = {complex(7, 2), complex(4, 3), complex(6, 5), complex(2,3)};
    std::vector<complex> row4 = {complex(4, 3), complex(2, 3), complex(10, 2), complex(2, 5)};
    m1.push_back(row1);
    m1.push_back(row2);
    m1.push_back(row3);
    m1.push_back(row4);

    matrix_complex m2;
    m2.push_back(row1);
    m2.push_back(row2);
    m2.push_back(row3);
    m2.push_back(row4);

    strassen(m1, m2);



    return 0;
}
