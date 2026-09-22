// Identical CPU preconditioning before every measured process; never timed.
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <unistd.h>

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    using Clock = std::chrono::steady_clock;
    const auto end = Clock::now() + std::chrono::milliseconds(std::strtoul(argv[1], nullptr, 10));
    uint64_t value = 0x9e3779b97f4a7c15ULL;
    do {
        for (unsigned i = 0; i < 4096; ++i) {
            value ^= value << 7;
            value ^= value >> 9;
            asm volatile("" : "+r"(value));
        }
    } while (Clock::now() < end);
    execvp(argv[2], argv + 2);
    std::cerr << "cannot execute benchmark after preconditioning\n";
    return 2;
}
