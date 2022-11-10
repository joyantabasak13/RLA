#include <iostream>
#include <string>
#include <algorithm>

int main() {
    std::string str = "Hello8$World";

    auto last = std::remove_if(str.begin(), str.end(), [](auto ch) {
        return ::isdigit(ch) || ::ispunct(ch);
    });
    str.erase(last, str.end());

    std::cout << str << std::endl;
    return 0;
}