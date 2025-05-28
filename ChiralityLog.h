#include <iostream>

#define CHIRALITY_INFO(s)                                  \
    {                                                      \
        std::cout << "\033[32m[INFO]: " + s + "\033[0m\n"; \
    }

#define CHIRALITY_WARN(s)                                  \
    {                                                      \
        std::cout << "\033[33m[WARN]: " + s + "\033[0m\n"; \
    }

#define CHIRALITY_ERROR(s)                                  \
    {                                                       \
        std::cout << "\033[31m[ERROR]: " + s + "\033[0m\n"; \
    }