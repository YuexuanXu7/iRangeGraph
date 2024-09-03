#pragma once

#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <sys/mman.h>


namespace memory
{
    template <uint32_t alignment>
    inline void *align_mm(size_t nbytes)
    {
        size_t size = (nbytes + alignment - 1) / alignment * alignment;
        void *p = std::aligned_alloc(alignment, size);
        std::memset(p, 0, size);
        return p;
    }


    template <typename T>
    struct align_alloc
    {
        T *ptr = nullptr;
        using value_type = T;
        T *allocate(int n)
        {
            if (n <= 1 << 14)
            {
                int sz = (n * sizeof(T) + 63) >> 6 << 6;
                return ptr = (T *)std::aligned_alloc(64, sz);
            }
            int sz = (n * sizeof(T) + (1 << 21) - 1) >> 21 << 21;
            ptr = (T *)std::aligned_alloc(1 << 21, sz);
            madvise(ptr, sz, MADV_HUGEPAGE);
            return ptr;
        }
        void deallocate(T *, int) { free(ptr); }
        template <typename U>
        struct rebind
        {
            typedef align_alloc<U> other;
        };
        bool operator!=(const align_alloc &rhs) { return ptr != rhs.ptr; }
    };


    inline void prefetch_L1(const void *address)
    {
#if defined(__SSE2__)
        _mm_prefetch((const char *)address, _MM_HINT_T0);
#else
        __builtin_prefetch(address, 0, 3);
#endif
    }


    inline void mem_prefetch_L1(char *ptr, const int num_lines)
    {
        switch (num_lines)
        {
        default:
            [[fallthrough]];
        case 28:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 27:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 26:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 25:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 24:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 23:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 22:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 21:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 20:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 19:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 18:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 17:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 16:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 15:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 14:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 13:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 12:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 11:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 10:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 9:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 8:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 7:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 6:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 5:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 4:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 3:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 2:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 1:
            prefetch_L1(ptr);
            ptr += 64;
            [[fallthrough]];
        case 0:
            break;
        }
    }
}