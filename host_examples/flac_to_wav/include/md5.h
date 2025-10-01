// Simple MD5 implementation for FLAC audio verification
// Based on RFC 1321: https://www.ietf.org/rfc/rfc1321.txt

#pragma once

#include <cstdint>
#include <cstring>
#include <array>

class MD5 {
public:
    MD5() { reset(); }

    void reset() {
        count_[0] = 0;
        count_[1] = 0;
        state_[0] = 0x67452301;
        state_[1] = 0xefcdab89;
        state_[2] = 0x98badcfe;
        state_[3] = 0x10325476;
    }

    void update(const uint8_t* input, size_t length) {
        size_t index = (count_[0] >> 3) & 0x3F;

        count_[0] += length << 3;
        if (count_[0] < (length << 3)) {
            count_[1]++;
        }
        count_[1] += length >> 29;

        size_t part_len = 64 - index;
        size_t i = 0;

        if (length >= part_len) {
            std::memcpy(&buffer_[index], input, part_len);
            transform(buffer_);

            for (i = part_len; i + 63 < length; i += 64) {
                transform(&input[i]);
            }
            index = 0;
        }

        std::memcpy(&buffer_[index], &input[i], length - i);
    }

    std::array<uint8_t, 16> finalize() {
        uint8_t bits[8];
        encode(bits, count_, 8);

        size_t index = (count_[0] >> 3) & 0x3f;
        size_t pad_len = (index < 56) ? (56 - index) : (120 - index);

        static const uint8_t padding[64] = { 0x80 };
        update(padding, pad_len);
        update(bits, 8);

        std::array<uint8_t, 16> digest;
        encode(digest.data(), state_, 16);

        reset();
        return digest;
    }

private:
    uint32_t state_[4];
    uint32_t count_[2];
    uint8_t buffer_[64];

    static inline uint32_t F(uint32_t x, uint32_t y, uint32_t z) {
        return (x & y) | (~x & z);
    }

    static inline uint32_t G(uint32_t x, uint32_t y, uint32_t z) {
        return (x & z) | (y & ~z);
    }

    static inline uint32_t H(uint32_t x, uint32_t y, uint32_t z) {
        return x ^ y ^ z;
    }

    static inline uint32_t I(uint32_t x, uint32_t y, uint32_t z) {
        return y ^ (x | ~z);
    }

    static inline uint32_t rotate_left(uint32_t x, uint32_t n) {
        return (x << n) | (x >> (32 - n));
    }

    static inline void FF(uint32_t &a, uint32_t b, uint32_t c, uint32_t d, uint32_t x, uint32_t s, uint32_t ac) {
        a += F(b, c, d) + x + ac;
        a = rotate_left(a, s);
        a += b;
    }

    static inline void GG(uint32_t &a, uint32_t b, uint32_t c, uint32_t d, uint32_t x, uint32_t s, uint32_t ac) {
        a += G(b, c, d) + x + ac;
        a = rotate_left(a, s);
        a += b;
    }

    static inline void HH(uint32_t &a, uint32_t b, uint32_t c, uint32_t d, uint32_t x, uint32_t s, uint32_t ac) {
        a += H(b, c, d) + x + ac;
        a = rotate_left(a, s);
        a += b;
    }

    static inline void II(uint32_t &a, uint32_t b, uint32_t c, uint32_t d, uint32_t x, uint32_t s, uint32_t ac) {
        a += I(b, c, d) + x + ac;
        a = rotate_left(a, s);
        a += b;
    }

    void transform(const uint8_t block[64]) {
        uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        uint32_t x[16];

        decode(x, block, 64);

        // Round 1
        FF(a, b, c, d, x[ 0],  7, 0xd76aa478);
        FF(d, a, b, c, x[ 1], 12, 0xe8c7b756);
        FF(c, d, a, b, x[ 2], 17, 0x242070db);
        FF(b, c, d, a, x[ 3], 22, 0xc1bdceee);
        FF(a, b, c, d, x[ 4],  7, 0xf57c0faf);
        FF(d, a, b, c, x[ 5], 12, 0x4787c62a);
        FF(c, d, a, b, x[ 6], 17, 0xa8304613);
        FF(b, c, d, a, x[ 7], 22, 0xfd469501);
        FF(a, b, c, d, x[ 8],  7, 0x698098d8);
        FF(d, a, b, c, x[ 9], 12, 0x8b44f7af);
        FF(c, d, a, b, x[10], 17, 0xffff5bb1);
        FF(b, c, d, a, x[11], 22, 0x895cd7be);
        FF(a, b, c, d, x[12],  7, 0x6b901122);
        FF(d, a, b, c, x[13], 12, 0xfd987193);
        FF(c, d, a, b, x[14], 17, 0xa679438e);
        FF(b, c, d, a, x[15], 22, 0x49b40821);

        // Round 2
        GG(a, b, c, d, x[ 1],  5, 0xf61e2562);
        GG(d, a, b, c, x[ 6],  9, 0xc040b340);
        GG(c, d, a, b, x[11], 14, 0x265e5a51);
        GG(b, c, d, a, x[ 0], 20, 0xe9b6c7aa);
        GG(a, b, c, d, x[ 5],  5, 0xd62f105d);
        GG(d, a, b, c, x[10],  9, 0x02441453);
        GG(c, d, a, b, x[15], 14, 0xd8a1e681);
        GG(b, c, d, a, x[ 4], 20, 0xe7d3fbc8);
        GG(a, b, c, d, x[ 9],  5, 0x21e1cde6);
        GG(d, a, b, c, x[14],  9, 0xc33707d6);
        GG(c, d, a, b, x[ 3], 14, 0xf4d50d87);
        GG(b, c, d, a, x[ 8], 20, 0x455a14ed);
        GG(a, b, c, d, x[13],  5, 0xa9e3e905);
        GG(d, a, b, c, x[ 2],  9, 0xfcefa3f8);
        GG(c, d, a, b, x[ 7], 14, 0x676f02d9);
        GG(b, c, d, a, x[12], 20, 0x8d2a4c8a);

        // Round 3
        HH(a, b, c, d, x[ 5],  4, 0xfffa3942);
        HH(d, a, b, c, x[ 8], 11, 0x8771f681);
        HH(c, d, a, b, x[11], 16, 0x6d9d6122);
        HH(b, c, d, a, x[14], 23, 0xfde5380c);
        HH(a, b, c, d, x[ 1],  4, 0xa4beea44);
        HH(d, a, b, c, x[ 4], 11, 0x4bdecfa9);
        HH(c, d, a, b, x[ 7], 16, 0xf6bb4b60);
        HH(b, c, d, a, x[10], 23, 0xbebfbc70);
        HH(a, b, c, d, x[13],  4, 0x289b7ec6);
        HH(d, a, b, c, x[ 0], 11, 0xeaa127fa);
        HH(c, d, a, b, x[ 3], 16, 0xd4ef3085);
        HH(b, c, d, a, x[ 6], 23, 0x04881d05);
        HH(a, b, c, d, x[ 9],  4, 0xd9d4d039);
        HH(d, a, b, c, x[12], 11, 0xe6db99e5);
        HH(c, d, a, b, x[15], 16, 0x1fa27cf8);
        HH(b, c, d, a, x[ 2], 23, 0xc4ac5665);

        // Round 4
        II(a, b, c, d, x[ 0],  6, 0xf4292244);
        II(d, a, b, c, x[ 7], 10, 0x432aff97);
        II(c, d, a, b, x[14], 15, 0xab9423a7);
        II(b, c, d, a, x[ 5], 21, 0xfc93a039);
        II(a, b, c, d, x[12],  6, 0x655b59c3);
        II(d, a, b, c, x[ 3], 10, 0x8f0ccc92);
        II(c, d, a, b, x[10], 15, 0xffeff47d);
        II(b, c, d, a, x[ 1], 21, 0x85845dd1);
        II(a, b, c, d, x[ 8],  6, 0x6fa87e4f);
        II(d, a, b, c, x[15], 10, 0xfe2ce6e0);
        II(c, d, a, b, x[ 6], 15, 0xa3014314);
        II(b, c, d, a, x[13], 21, 0x4e0811a1);
        II(a, b, c, d, x[ 4],  6, 0xf7537e82);
        II(d, a, b, c, x[11], 10, 0xbd3af235);
        II(c, d, a, b, x[ 2], 15, 0x2ad7d2bb);
        II(b, c, d, a, x[ 9], 21, 0xeb86d391);

        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
    }

    static void decode(uint32_t* output, const uint8_t* input, size_t len) {
        for (size_t i = 0, j = 0; j < len; i++, j += 4) {
            output[i] = ((uint32_t)input[j]) | (((uint32_t)input[j+1]) << 8) |
                       (((uint32_t)input[j+2]) << 16) | (((uint32_t)input[j+3]) << 24);
        }
    }

    static void encode(uint8_t* output, const uint32_t* input, size_t len) {
        for (size_t i = 0, j = 0; j < len; i++, j += 4) {
            output[j]   = (uint8_t)(input[i] & 0xff);
            output[j+1] = (uint8_t)((input[i] >> 8) & 0xff);
            output[j+2] = (uint8_t)((input[i] >> 16) & 0xff);
            output[j+3] = (uint8_t)((input[i] >> 24) & 0xff);
        }
    }
};
