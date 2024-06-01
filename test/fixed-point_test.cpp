#include <bits/stdint-uintn.h>
#include <iostream>

int main() {
    double number1 = -3.14159265358979323846264338327950288;
    double number2 = 1.23456789876543212345678987654321;
    uint64_t fixed_number1 = (uint64_t(number1 * (1ull << 31)));
    uint64_t fixed_number2 = (uint64_t(number2 * (1ull << 31)));
    uint64_t fixed_number = fixed_number1 * fixed_number2;
    fixed_number >>= 31;
    printf("%.20lf\n", (double)fixed_number / (1ull << 31));
    printf("%.20lf\n", number1 * number2);
}