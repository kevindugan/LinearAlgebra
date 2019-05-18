#ifndef PARTITIONING_H_VS11
#define PARTITIONING_H_VS11

// Index that is at the beginning of the range and and index one past the end
struct IndexRange {
    IndexRange() : begin(0), end(1), skip(1) {}
    IndexRange(unsigned int a, unsigned int b) : begin(a), end(b), skip(1) {}
    IndexRange(unsigned int a, unsigned int b, unsigned int s) : begin(a), end(b), skip(s) {}
    unsigned int begin, end, skip;
};

#endif // PARTITIONING_H_VS11