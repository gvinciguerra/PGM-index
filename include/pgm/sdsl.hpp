/* sdsl - succinct data structures library
    Copyright (C) 2008-2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file int_vector.hpp
    \brief int_vector.hpp contains the sdsl::int_vector class.
    \author Simon Gog
*/

#pragma once

#ifndef INCLUDED_SDSL_INT_VECTOR
#define INCLUDED_SDSL_INT_VECTOR

/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file bits.hpp
    \brief bits.hpp contains the sdsl::bits class.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_BITS
#define INCLUDED_SDSL_BITS

#include <stdint.h> // for uint64_t uint32_t declaration
#include <iostream>// for cerr
#include <cassert>
#ifdef __BMI2__
#include <immintrin.h>
#endif
#ifdef __SSE4_2__
#include <xmmintrin.h>
#endif

#ifdef WIN32
#include <ciso646>
#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A helper class for bitwise tricks on 64 bit words.
/*!
	bits is a helper class for bitwise tricks and
	techniques. For the basic tricks and techiques we refer to Donald E. Knuth's
	"The Art of Computer Programming", Volume 4A, Chapter 7.1.3 and
	the informative website of Sean E. Anderson about the topic:
	http://www-graphics.stanford.edu/~seander/bithacks.html .

	We have added new functions like: cnt11 and sel11.

	All members of this class are static variables or methods.
	This class cannot be instantiated.

	\author Simon Gog
 */
struct bits {
    bits() = delete;
    //! 64bit mask with all bits set to 1.
    constexpr static uint64_t  all_set {-1ULL};

    //! This constant represents a de Bruijn sequence B(k,n) for k=2 and n=6.
    /*! Details for de Bruijn sequences see
       http://en.wikipedia.org/wiki/De_bruijn_sequence
       deBruijn64 is used in combination with the
       array lt_deBruijn_to_idx.
    */
    constexpr static uint64_t deBruijn64 {0x0218A392CD3D5DBFULL};

    //! This table maps a 6-bit subsequence S[idx...idx+5] of constant deBruijn64 to idx.
    /*! \sa deBruijn64
    */
    static const uint32_t lt_deBruijn_to_idx[64];

    //! Array containing Fibonacci numbers less than \f$2^64\f$.
    static const uint64_t lt_fib[92];

    //! Lookup table for byte popcounts.
    static const uint8_t lt_cnt[256];

    //! Lookup table for most significant set bit in a byte.
    static const uint32_t lt_hi[256];

    //! lo_set[i] is a 64-bit word with the i least significant bits set and the high bits not set.
    /*! lo_set[0] = 0ULL, lo_set[1]=1ULL, lo_set[2]=3ULL...
     */
    static const uint64_t lo_set[65];

    //! lo_unset[i] is a 64-bit word with the i least significant bits not set and the high bits set.
    /*! lo_unset[0] = FFFFFFFFFFFFFFFFULL, lo_unset_set[1]=FFFFFFFFFFFFFFFEULL, ...
     */
    static const uint64_t lo_unset[65];

    //! Lookup table for least significant set bit in a byte.
    static const uint8_t lt_lo[256];

    //! Lookup table for select on bytes.
    /*! Entry at idx = 256*j + i equals the position of the
        (j+1)-th set bit in byte i. Positions lie in the range \f$[0..7]\f$.
     */
    static const uint8_t lt_sel[256*8];

    //! Use to help to decide if a prefix sum stored in a byte overflows.
    static const uint64_t ps_overflow[65];

    //! Counts the number of set bits in x.
    /*! \param  x 64-bit word
        \return Number of set bits.
     */
    static uint64_t cnt(uint64_t x);

    //! Position of the most significant set bit the 64-bit word x
    /*! \param x 64-bit word
        \return The position (in 0..63) of the most significant set bit
                in `x` or 0 if x equals 0.
    	\sa sel, lo
    */
    static uint32_t hi(uint64_t x);

    //! Calculates the position of the rightmost 1-bit in the 64bit integer x if it exists
    /*! \param x 64 bit integer.
    	\return The position (in 0..63) of the rightmost 1-bit in the 64bit integer x if
    	        x>0 and 0 if x equals 0.
    	\sa sel, hi
    */
    static uint32_t lo(uint64_t x);

    //! Counts the number of 1-bits in the 32bit integer x.
    /*! This function is a variant of the method cnt. If
    	32bit multiplication is fast, this method beats the cnt.
    	for 32bit integers.
    	\param x 64bit integer to count the bits.
    	\return The number of 1-bits in x.
     */
    static uint32_t cnt32(uint32_t x);

    //! Count the number of consecutive and distinct 11 in the 64bit integer x.
    /*!
      	\param x 64bit integer to count the terminating sequence 11 of a Fibonacci code.
    	\param c Carry equals msb of the previous 64bit integer.
     */
    static uint32_t cnt11(uint64_t x, uint64_t& c);

    //! Count the number of consecutive and distinct 11 in the 64bit integer x.
    /*!
      	\param x 64bit integer to count the terminating sequence 11 of a Fibonacci code.
     */
    static uint32_t cnt11(uint64_t x);

    //! Count 10 bit pairs in the word x.
    /*!
     * \param x 64bit integer to count the 10 bit pairs.
     * \param c Carry equals msb of the previous 64bit integer.
     */
    static uint32_t cnt10(uint64_t x, uint64_t& c);

    //! Count 01 bit pairs in the word x.
    /*!
     * \param x 64bit integer to count the 01 bit pairs.
     * \param c Carry equals msb of the previous 64bit integer.
     */
    static uint32_t cnt01(uint64_t x, uint64_t& c);

    //! Map all 10 bit pairs to 01 or 1 if c=1 and the lsb=0. All other pairs are mapped to 00.
    static uint64_t map10(uint64_t x, uint64_t c=0);

    //! Map all 01 bit pairs to 01 or 1 if c=1 and the lsb=0. All other pairs are mapped to 00.
    static uint64_t map01(uint64_t x, uint64_t c=1);

    //! Calculate the position of the i-th rightmost 1 bit in the 64bit integer x
    /*!
      	\param x 64bit integer.
    	\param i Argument i must be in the range \f$[1..cnt(x)]\f$.
    	\pre Argument i must be in the range \f$[1..cnt(x)]\f$.
      	\sa hi, lo
     */
    static uint32_t sel(uint64_t x, uint32_t i);
    static uint32_t _sel(uint64_t x, uint32_t i);

    //! Calculates the position of the i-th rightmost 11-bit-pattern which terminates a Fibonacci coded integer in x.
    /*!	\param x 64 bit integer.
        \param i Index of 11-bit-pattern. \f$i \in [1..cnt11(x)]\f$
    	\param c Carry bit from word before
     	\return The position (in 1..63) of the i-th 11-bit-pattern which terminates a Fibonacci coded integer in x if
    	        x contains at least i 11-bit-patterns and a undefined value otherwise.
        \sa cnt11, hi11, sel

     */
    static uint32_t sel11(uint64_t x, uint32_t i, uint32_t c=0);

    //! Calculates the position of the leftmost 11-bit-pattern which terminates a Fibonacci coded integer in x.
    /*! \param x 64 bit integer.
        \return The position (in 1..63) of the leftmost 1 of the leftmost 11-bit-pattern which
    	        terminates a Fibonacci coded integer in x if x contains a 11-bit-pattern
    			and 0 otherwise.
    	\sa cnt11, sel11
    */
    static uint32_t hi11(uint64_t x);

    //! Writes value x to an bit position in an array.
    static void write_int(uint64_t* word, uint64_t x, const uint8_t offset=0, const uint8_t len=64);

    //! Writes value x to an bit position in an array and moves the bit-pointer.
    static void write_int_and_move(uint64_t*& word, uint64_t x, uint8_t& offset, const uint8_t len);

    //! Reads a value from a bit position in an array.
    static uint64_t read_int(const uint64_t* word, uint8_t offset=0, const uint8_t len=64);

    //! Reads a value from a bit position in an array and moved the bit-pointer.
    static uint64_t read_int_and_move(const uint64_t*& word, uint8_t& offset, const uint8_t len=64);

    //! Reads an unary decoded value from a bit position in an array.
    static uint64_t read_unary(const uint64_t* word, uint8_t offset=0);

    //! Reads an unary decoded value from a bit position in an array and moves the bit-pointer.
    static uint64_t read_unary_and_move(const uint64_t*& word, uint8_t& offset);

    //! Move the bit-pointer (=uint64_t word and offset) `len` to the right.
    /*!\param word   64-bit word part of the bit pointer
     * \param offset Offset part of the bit pointer
     * \param len    Move distance. \f$ len \in [0..64] \f$
     * \sa move_left
     */
    static void move_right(const uint64_t*& word, uint8_t& offset, const uint8_t len);

    //! Move the bit-pointer (=uint64_t word and offset) `len` to the left.
    /*!\param word   64-bit word part of the bit pointer
     * \param offset Offset part of the bit pointer
     * \param len    Move distance. \f$ len \in [0..64] \f$
     * \sa move_right
     */
    static void move_left(const uint64_t*& word, uint8_t& offset, const uint8_t len);

    //! Get the first one bit in the interval \f$[idx..\infty )\f$
    static uint64_t next(const uint64_t* word, uint64_t idx);

    //! Get the one bit with the greatest position in the interval \f$[0..idx]\f$
    static uint64_t prev(const uint64_t* word, uint64_t idx);

    //! reverses a given 64 bit word
    static uint64_t rev(uint64_t x);
};

// ============= inline - implementations ================

// see page 11, Knuth TAOCP Vol 4 F1A
inline uint64_t bits::cnt(uint64_t x)
{
#ifdef __SSE4_2__
    return __builtin_popcountll(x);
#else
    #ifdef POPCOUNT_TL
    return lt_cnt[x&0xFFULL] + lt_cnt[(x>>8)&0xFFULL] +
           lt_cnt[(x>>16)&0xFFULL] + lt_cnt[(x>>24)&0xFFULL] +
           lt_cnt[(x>>32)&0xFFULL] + lt_cnt[(x>>40)&0xFFULL] +
           lt_cnt[(x>>48)&0xFFULL] + lt_cnt[(x>>56)&0xFFULL];
#else
    x = x-((x>>1) & 0x5555555555555555ull);
    x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
    return (0x0101010101010101ull*x >> 56);
#endif
#endif
}

inline uint32_t bits::cnt32(uint32_t x)
{
    x = x-((x>>1) & 0x55555555);
    x = (x & 0x33333333) + ((x>>2) & 0x33333333);
    return (0x10101010*x >>28)+(0x01010101*x >>28);
}

inline uint32_t bits::cnt11(uint64_t x, uint64_t& c)
{
    // extract "11" 2bit blocks
    uint64_t ex11 = (x&(x>>1))&0x5555555555555555ULL, t;
    // extract "10" 2bit blocks
    uint64_t ex10or01 = (ex11|(ex11<<1))^x;

    x = ex11 | ((t=(ex11|(ex11<<1))+(((ex10or01<<1)&0x5555555555555555ULL)|c))&(ex10or01&0x5555555555555555ULL));
    c = (ex10or01>>63) or(t < (ex11|(ex11<<1)));

    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (0x0101010101010101ULL*x >> 56);
}

inline uint32_t bits::cnt11(uint64_t x)
{
    // extract "11" 2bit blocks
    uint64_t ex11 = (x&(x>>1))&0x5555555555555555ULL;
    // extract "10" 2bit blocks
    uint64_t ex10or01 = (ex11|(ex11<<1))^x;

    x = ex11 | (((ex11|(ex11<<1))+((ex10or01<<1)&0x5555555555555555ULL))&(ex10or01&0x5555555555555555ULL));

    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (0x0101010101010101ULL*x >> 56);
}

inline uint32_t bits::cnt10(uint64_t x, uint64_t& c)
{
    uint32_t res = cnt((x ^((x<<1) | c)) & (~x));
    c = (x >> 63);
    return res;
}

inline uint64_t bits::map10(uint64_t x, uint64_t c)
{
    return ((x ^((x << 1) | c)) & (~x));
}

inline uint32_t bits::cnt01(uint64_t x, uint64_t& c)
{
    uint32_t res = cnt((x ^((x<<1) | c)) & x);
    c = (x >> 63);
    return res;
}
inline uint64_t bits::map01(uint64_t x, uint64_t c)
{
    return ((x ^((x << 1) | c)) &  x);
}

inline uint32_t bits::sel(uint64_t x, uint32_t i)
{
#ifdef __BMI2__
    // index i is 1-based here, (i-1) changes it to 0-based
    return __builtin_ctzll(_pdep_u64(1ull << (i-1), x));
#elif defined(__SSE4_2__)
    uint64_t s = x, b;
    s = s-((s>>1) & 0x5555555555555555ULL);
    s = (s & 0x3333333333333333ULL) + ((s >> 2) & 0x3333333333333333ULL);
    s = (s + (s >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    s = 0x0101010101010101ULL*s;
// now s contains 8 bytes s[7],...,s[0]; s[j] contains the cumulative sum
// of (j+1)*8 least significant bits of s
    b = (s+ps_overflow[i]) & 0x8080808080808080ULL;
// ps_overflow contains a bit mask x consisting of 8 bytes
// x[7],...,x[0] and x[j] is set to 128-j
// => a byte b[j] in b is >= 128 if cum sum >= j

// __builtin_ctzll returns the number of trailing zeros, if b!=0
    int  byte_nr = __builtin_ctzll(b) >> 3;   // byte nr in [0..7]
    s <<= 8;
    i -= (s >> (byte_nr<<3)) & 0xFFULL;
    return (byte_nr << 3) + lt_sel[((i-1) << 8) + ((x>>(byte_nr<<3))&0xFFULL) ];
#else
    return _sel(x, i);
#endif
}

inline uint32_t bits::_sel(uint64_t x, uint32_t i)
{
    uint64_t s = x, b;  // s = sum
    s = s-((s>>1) & 0x5555555555555555ULL);
    s = (s & 0x3333333333333333ULL) + ((s >> 2) & 0x3333333333333333ULL);
    s = (s + (s >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    s = 0x0101010101010101ULL*s;
    b = (s+ps_overflow[i]);//&0x8080808080808080ULL;// add something to the partial sums to cause overflow
    i = (i-1)<<8;
    if (b&0x0000000080000000ULL) // byte <=3
        if (b&0x0000000000008000ULL) //byte <= 1
            if (b&0x0000000000000080ULL)
                return    lt_sel[(x&0xFFULL) + i];
            else
                return 8 +lt_sel[(((x>>8)&0xFFULL)  + i - ((s&0xFFULL)<<8))&0x7FFULL];//byte 1;
        else//byte >1
        if (b&0x0000000000800000ULL) //byte <=2
            return 16+lt_sel[(((x>>16)&0xFFULL) + i - (s&0xFF00ULL))&0x7FFULL];//byte 2;
        else
            return 24+lt_sel[(((x>>24)&0xFFULL) + i - ((s>>8)&0xFF00ULL))&0x7FFULL];//byte 3;
    else//  byte > 3
    if (b&0x0000800000000000ULL) // byte <=5
        if (b&0x0000008000000000ULL) //byte <=4
            return 32+lt_sel[(((x>>32)&0xFFULL) + i - ((s>>16)&0xFF00ULL))&0x7FFULL];//byte 4;
        else
            return 40+lt_sel[(((x>>40)&0xFFULL) + i - ((s>>24)&0xFF00ULL))&0x7FFULL];//byte 5;
    else// byte >5
    if (b&0x0080000000000000ULL) //byte<=6
        return 48+lt_sel[(((x>>48)&0xFFULL) + i - ((s>>32)&0xFF00ULL))&0x7FFULL];//byte 6;
    else
        return 56+lt_sel[(((x>>56)&0xFFULL) + i - ((s>>40)&0xFF00ULL))&0x7FFULL];//byte 7;
    return 0;
}

// using built-in method or
// 64-bit version of 32-bit proposal of
// http://www-graphics.stanford.edu/~seander/bithacks.html
inline uint32_t bits::hi(uint64_t x)
{
#ifdef __SSE4_2__
    if (x == 0)
        return 0;
    return 63 - __builtin_clzll(x);
#else
    uint64_t t,tt; // temporaries
    if ((tt = x >> 32)) { // hi >= 32
        if ((t = tt >> 16)) { // hi >= 48
            return (tt = t >> 8) ? 56 + lt_hi[tt] : 48 + lt_hi[t];
        } else { // hi < 48
            return (t = tt >> 8) ? 40 + lt_hi[t] : 32 + lt_hi[tt];
        }
    } else { // hi < 32
        if ((t = x >> 16)) { // hi >= 16
            return (tt = t >> 8) ? 24 + lt_hi[tt] : 16 + lt_hi[t];
        } else { // hi < 16
            return (tt = x >> 8) ?  8 + lt_hi[tt] : lt_hi[x];
        }
    }
#endif
}

// details see: http://citeseer.ist.psu.edu/leiserson98using.html
// or page 10, Knuth TAOCP Vol 4 F1A
inline uint32_t bits::lo(uint64_t x)
{
#ifdef __SSE4_2__
    if (x==0)
        return 0;
    return __builtin_ctzll(x);
#else
    if (x&1) return 0;
    if (x&3) return 1;
    if (x&7) return 2;
    if (x&0x7FF) { // in average every second random number x can be answered this way
        return lt_lo[(x&0x7FF)>>3]+3;
    }
    // x&-x equals x with only the lsb set
    return lt_deBruijn_to_idx[((x&-x)*deBruijn64)>>58];
#endif
}

inline uint32_t bits::hi11(uint64_t x)
{
    // extract "11" 2bit blocks
    uint64_t ex11 = (x&(x>>1))&0x5555555555555555ULL;
    // extract "10" 2bit blocks
    uint64_t ex10or01 = (ex11|(ex11<<1))^x;
    // extract "10" 2bit blocks
    ex11 += (((ex11|(ex11<<1))+((ex10or01<<1)&0x5555555555555555ULL)) & ((ex10or01&0x5555555555555555ULL)|ex11));
    return hi(ex11);
}

inline uint32_t bits::sel11(uint64_t x, uint32_t i, uint32_t c)
{
    uint64_t ex11 = (x&(x>>1))&0x5555555555555555ULL;
    uint64_t ex10or01 = (ex11|(ex11<<1))^x;
    ex11 += (((ex11|(ex11<<1))+(((ex10or01<<1)&0x5555555555555555ULL)|c)) & ((ex10or01&0x5555555555555555ULL)|ex11));
    return sel(ex11,i);
}

inline void bits::write_int(uint64_t* word, uint64_t x, uint8_t offset, const uint8_t len)
{
    x &= bits::lo_set[len];
    if (offset + len < 64) {
        *word &=
            ((bits::all_set << (offset+len)) | bits::lo_set[offset]); // mask 1..10..01..1
        *word |= (x << offset);
//		*word ^= ((*word ^ x) & (bits::lo_set[len] << offset) );
//      surprisingly the above line is slower than the lines above
    } else {
        *word &=
            ((bits::lo_set[offset]));  // mask 0....01..1
        *word |= (x << offset);
        if ((offset = (offset+len)&0x3F)) { // offset+len > 64
            *(word+1) &= (~bits::lo_set[offset]); // mask 1...10..0
//			*(word+1) &= bits::lo_unset[offset]; // mask 1...10..0
//          surprisingly the above line is slower than the line above
            *(word+1) |= (x >> (len-offset));
        }
    }
}

inline void bits::write_int_and_move(uint64_t*& word, uint64_t x, uint8_t& offset, const uint8_t len)
{
    x &= bits::lo_set[len];
    if (offset + len < 64) {
        *word &=
            ((bits::all_set << (offset+len)) | bits::lo_set[offset]); // mask 1..10..01..1
        *word |= (x << offset);
        offset += len;
    } else {
        *word &=
            ((bits::lo_set[offset]));  // mask 0....01..1
        *word |= (x << offset);
        if ((offset= (offset+len))>64) {// offset+len >= 64
            offset &= 0x3F;
            *(++word) &= (~bits::lo_set[offset]); // mask 1...10..0
            *word |= (x >> (len-offset));
        } else {
            offset = 0;
            ++word;
        }
    }
}

inline uint64_t bits::read_int(const uint64_t* word, uint8_t offset, const uint8_t len)
{
    uint64_t w1 = (*word)>>offset;
    if ((offset+len) > 64) { // if offset+len > 64
        return w1 |  // w1 or w2 adepted:
            ((*(word+1) & bits::lo_set[(offset+len)&0x3F])   // set higher bits zero
                << (64-offset));  // move bits to the left
    } else {
        return w1 & bits::lo_set[len];
    }
}

inline uint64_t bits::read_int_and_move(const uint64_t*& word, uint8_t& offset, const uint8_t len)
{
    uint64_t w1 = (*word)>>offset;
    if ((offset = (offset+len))>=64) {  // if offset+len > 64
        if (offset==64) {
            offset &= 0x3F;
            ++word;
            return w1;
        } else {
            offset &= 0x3F;
            return w1 |
                (((*(++word)) & bits::lo_set[offset]) << (len-offset));
        }
    } else {
        return w1 & bits::lo_set[len];
    }
}

inline uint64_t bits::read_unary(const uint64_t* word, uint8_t offset)
{
    uint64_t w = *word >> offset;
    if (w) {
        return bits::lo(w);
    } else {
        if (0!=(w=*(++word)))
            return bits::lo(w)+64-offset;
        uint64_t cnt=2;
        while (0==(w=*(++word)))
            ++cnt;
        return bits::lo(w)+(cnt<<6)-offset;
    }
    return 0;
}

inline uint64_t bits::read_unary_and_move(const uint64_t*& word, uint8_t& offset)
{
    uint64_t w = (*word) >> offset; // temporary variable is good for the performance
    if (w) {
        uint8_t r = bits::lo(w);
        offset = (offset + r+1)&0x3F;
        // we know that offset + r +1 <= 64, so if the new offset equals 0 increase word
        word += (offset==0);
        return r;
    } else {
        uint8_t rr=0;
        if (0!=(w=*(++word))) {
            rr = bits::lo(w)+64-offset;
            offset = (offset+rr+1)&0x3F;
            word += (offset==0);
            return rr;
        } else {
            uint64_t cnt_1=1;
            while (0==(w=*(++word)))
                ++cnt_1;
            rr = bits::lo(w)+64-offset;
            offset = (offset+rr+1)&0x3F;
            word += (offset==0);
            return ((cnt_1)<<6) + rr;
        }
    }
    return 0;
}

inline void bits::move_right(const uint64_t*& word, uint8_t& offset, const uint8_t len)
{
    if ((offset+=len)&0xC0) { // if offset >= 65
        offset&=0x3F;
        ++word;
    }
}

inline void bits::move_left(const uint64_t*& word, uint8_t& offset, const uint8_t len)
{
    if ((offset-=len)&0xC0) {  // if offset-len<0
        offset&=0x3F;
        --word;
    }
}

inline uint64_t bits::next(const uint64_t* word, uint64_t idx)
{
    word += (idx>>6);
    if (*word & ~lo_set[idx&0x3F]) {
        return (idx & ~((size_t)0x3F)) + lo(*word & ~lo_set[idx&0x3F]);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    ++word;
    while (*word==0) {
        idx += 64;
        ++word;
    }
    return idx + lo(*word);
}

inline uint64_t bits::prev(const uint64_t* word, uint64_t idx)
{
    word += (idx>>6);
    if (*word & lo_set[(idx&0x3F)+1]) {
        return (idx & ~((size_t)0x3F)) + hi(*word & lo_set[(idx&0x3F)+1]);
    }
    idx = (idx & ~((size_t)0x3F)) - 64;
    --word;
    while (*word==0) {
        idx -= 64;
        --word;
    }
    return idx + hi(*word);
}

inline uint64_t bits::rev(uint64_t x)
{
    x = ((x & 0x5555555555555555ULL) << 1) | ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1);
    x = ((x & 0x3333333333333333ULL) << 2) | ((x & 0xCCCCCCCCCCCCCCCCULL) >> 2);
    x = ((x & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((x & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
    x = ((x & 0x00FF00FF00FF00FFULL) << 8) | ((x & 0xFF00FF00FF00FF00ULL) >> 8);
    x = ((x & 0x0000FFFF0000FFFFULL) <<16) | ((x & 0xFFFF0000FFFF0000ULL) >>16);
    x = ((x & 0x00000000FFFFFFFFULL) <<32) | ((x & 0xFFFFFFFF00000000ULL) >>32);
    return x;
}

} // end namespace sdsl

#endif

/*!\file structure_tree.hpp
   \brief structure_tree.hpp contains a helper class which can represent the memory structure of a class.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_STRUCTURE_TREE
#define INCLUDED_SDSL_STRUCTURE_TREE

#ifndef INCLUDED_SDSL_UINTX_T
#define INCLUDED_SDSL_UINTX_T

#include <cstdint>

using std::int8_t;
using std::int16_t;
using std::int32_t;
using std::int64_t;

using std::uint8_t;
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;

#endif

#include <unordered_map>
#include <string>
#include <iostream>
#include <sstream>
#include <memory>
#ifndef SDSL_CONFIG
#define SDSL_CONFIG

#include <map>
#include <string>

namespace sdsl
{
namespace conf  // namespace for library constant
{
// size of the buffer for reading and writing data in elements (not in bytes)
const uint64_t SDSL_BLOCK_SIZE = (uint64_t)1<<22;

const char KEY_BWT[] 		= "bwt";
const char KEY_BWT_INT[]	= "bwt_int";
const char KEY_SA[] 		= "sa";
const char KEY_CSA[] 		= "csa";
const char KEY_CST[] 		= "cst";
const char KEY_ISA[] 		= "isa";
const char KEY_TEXT[] 		= "text";
const char KEY_TEXT_INT[] 	= "text_int";
const char KEY_PSI[] 		= "psi";
const char KEY_LCP[] 		= "lcp";
const char KEY_SAMPLE_CHAR[]= "sample_char";
}
typedef uint64_t int_vector_size_type;

typedef std::map<std::string, std::string> tMSS;

enum format_type {JSON_FORMAT, R_FORMAT, HTML_FORMAT};

enum byte_sa_algo_type {LIBDIVSUFSORT, SE_SAIS};

//! Helper class for construction process
struct cache_config {
    bool 		delete_files;   // Flag which indicates if all files which were created
    // during construction should be deleted.
    std::string dir;    		// Directory for temporary files.
    std::string id;     		// Identifier is part of temporary file names. If
    // id is the empty string, then it will be replace
    // a concatenation of PID and a unique ID inside the
    // current process.
    tMSS 		file_map;		// Files stored during the construction process.
    cache_config(bool f_delete_files=true, std::string f_dir="./", std::string f_id="", tMSS f_file_map=tMSS());
};

//! Helper classes to transform width=0 and width=8 to corresponding text key
template<uint8_t width>
struct key_text_trait {
    static const char* KEY_TEXT;
};

//! Helper classes to transform width=0 and width=8 to corresponding bwt key
template<uint8_t width>
struct key_bwt_trait {
    static const char* KEY_BWT;
};
}

#endif

//! Namespace for the succinct data structure library
namespace sdsl
{

class structure_tree_node
{
private:
    using map_type = std::unordered_map<std::string,std::unique_ptr<structure_tree_node>>;
    map_type            m_children;
public:
    const map_type& children = m_children;
    size_t              size = 0;
    std::string         name;
    std::string         type;
public:
    structure_tree_node(const std::string& n, const std::string& t) : name(n) , type(t) {}
    structure_tree_node* add_child(const std::string& n, const std::string& t) {
        auto hash = n+t;
        auto child_itr = m_children.find(hash);
        if (child_itr == m_children.end()) {
            // add new child as we don't have one of this type yet
            structure_tree_node* new_node = new structure_tree_node(n,t);
            m_children[hash] = std::unique_ptr<structure_tree_node>(new_node);
            return new_node;
        } else {
            // child of same type and name exists
            return (*child_itr).second.get();
        }
    }
    void add_size(size_t s) { size += s; }
};

class structure_tree
{
public:
    static structure_tree_node* add_child(structure_tree_node* v, const std::string& name, const std::string& type) {
        if (v) return v->add_child(name,type);
        return nullptr;
    };
    static void add_size(structure_tree_node* v, uint64_t value) {
        if (v) v->add_size(value);
    };
};

template<format_type F>
void write_structure_tree(const structure_tree_node* v, std::ostream& out, size_t level = 0);

}
#endif

/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file util.hpp
    \brief util.hpp contains some helper methods for int_vector and other stuff like demangle class names.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_UTIL
#define INCLUDED_SDSL_UTIL

/*!\file sfstream.hpp
   \brief sfstream.hpp contains a two stream class which can be used to read/write from/to files or strings.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_SFSTREAM
#define INCLUDED_SDSL_SFSTREAM

#include <fstream>
#include <sstream>
#include <string>
/*! \file ram_fs.hpp
 * \brief ram_fs.hpp
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_RAM_FS
#define INCLUDED_SDSL_RAM_FS

#include <string>
#include <map>
#include <vector>
#include <mutex>

namespace sdsl
{

class ram_fs_initializer
{
public:
    ram_fs_initializer();
    ~ram_fs_initializer();
};

} // end namespace sdsl

static sdsl::ram_fs_initializer init_ram_fs;

namespace sdsl
{

//! ram_fs is a simple store for RAM-files.
/*!
 * Simple key-value store which maps file names
 * (strings) to file content (content_type).
 */
class ram_fs
{
public:
    typedef std::vector<char> content_type;

private:
    friend class ram_fs_initializer;
    typedef std::map<std::string, content_type> mss_type;
    static mss_type m_map;
    static std::recursive_mutex m_rlock;

public:
    //! Default construct
    ram_fs();
    static void store(const std::string& name, content_type data);
    //! Check if the file exists
    static bool exists(const std::string& name);
    //! Get the file size
    static size_t file_size(const std::string& name);
    //! Get the content
    static content_type& content(const std::string& name);
    //! Remove the file with key `name`
    static int remove(const std::string& name);
    //! Rename the file. Change key `old_filename` into `new_filename`.
    static int rename(const std::string old_filename, const std::string new_filename);
};

//! Determines if the given file is a RAM-file.
bool is_ram_file(const std::string& file);

//! Returns the corresponding RAM-file name for file.
std::string ram_file_name(const std::string& file);

//! Returns for a RAM-file the corresponding disk file name
std::string disk_file_name(const std::string& file);

//! Remove a file.
int remove(const std::string& file);

//! Rename a file
int rename(const std::string& old_filename, const std::string& new_filename);

} // end namespace sdsl
#endif

#ifndef INCLUDED_SDSL_RAM_FSTREAMBUF
#define INCLUDED_SDSL_RAM_FSTREAMBUF

#include <fstream>
#include <vector>

namespace sdsl
{

class ram_filebuf : public std::streambuf
{
private:
    ram_fs::content_type* m_ram_file = nullptr;  // file handle
    void pbump64(std::ptrdiff_t);

public:
    virtual ~ram_filebuf();

    ram_filebuf();
    ram_filebuf(std::vector<char>& ram_file);

    std::streambuf*
    open(const std::string s, std::ios_base::openmode mode);

    bool is_open();

    ram_filebuf*
    close();

    pos_type
    seekpos(pos_type sp,
            std::ios_base::openmode which = std::ios_base::in | std::ios_base::out) override;

    pos_type
    pubseekoff(off_type off, std::ios_base::seekdir way,
               std::ios_base::openmode which = std::ios_base::in | std::ios_base::out);

    pos_type
    pubseekpos(pos_type sp,
               std::ios_base::openmode which = std::ios_base::in | std::ios_base::out);

//    std::streamsize
//    xsputn(const char_type* s, std::streamsize n) override;

    int
    sync() override;

    int_type
    overflow(int_type c = traits_type::eof()) override;
};

}

#endif

namespace sdsl
{

class osfstream : public std::ostream
{
public:
    typedef std::streambuf* buf_ptr_type;
private:
    buf_ptr_type m_streambuf = nullptr;
    std::string  m_file      = "";
public:
    typedef void* voidptr;
    //! Standard constructor.
    osfstream();
    //! Constructor taking a file name and open mode.
    osfstream(const std::string& file, std::ios_base::openmode mode = std::ios_base::out);
    //! Open the stream.
    buf_ptr_type
    open(const std::string& file, std::ios_base::openmode mode = std::ios_base::out);
    //! Is the stream close?
    bool is_open();
    //! Close the stream.
    void close();
    //! Standard destructor
    ~osfstream();
    //! Cast to void*
    operator  voidptr() const;

    osfstream& seekp(pos_type pos);
    osfstream& seekp(off_type off, ios_base::seekdir way);
    std::streampos tellp();
};

class isfstream : public std::istream
{
    typedef std::streambuf* buf_ptr_type;
private:
    buf_ptr_type m_streambuf = nullptr;
    std::string  m_file      = "";
public:
    typedef void* voidptr;
    //! Standard constructor.
    isfstream();
    //! Constructor taking a file name and open mode.
    isfstream(const std::string& file, std::ios_base::openmode mode = std::ios_base::in);
    //! Open the stream.
    buf_ptr_type
    open(const std::string& file, std::ios_base::openmode mode = std::ios_base::in);
    //! Is the stream close?
    bool is_open();
    //! Close the stream.
    void close();
    //! Standard destructor
    ~isfstream();
    //! Cast to void*
    operator  voidptr() const;

    isfstream& seekg(pos_type pos);
    isfstream& seekg(off_type off, ios_base::seekdir way);
    std::streampos tellg();
};

} // end namespace

#endif

#include <iosfwd>      // forward declaration of ostream
#include <stdint.h>    // for uint64_t uint32_t declaration
#include <cassert>
#include <ctime>       // for rand initialization
#include <string>
#include <functional>  // for class_to_hash
#include <string.h>    // for strlen and strdup
#include <cstdlib>
#include <sstream>     // for to_string method
#include <stdexcept>   // for std::logic_error
#include <typeinfo>    // for typeid
#include <iomanip>
#include <numeric>
#include <random>
#include <chrono>
#include <atomic>
#include <mutex>
#include <algorithm>

// macros to transform a defined name to a string
#define SDSL_STR(x) #x
#define SDSL_XSTR(s) SDSL_STR(s)

#ifndef MSVC_COMPILER
#define SDSL_UNUSED __attribute__ ((unused))
#include <sys/time.h>  // for struct timeval
#include <sys/resource.h> // for struct rusage
#include <libgen.h>    // for basename
#include <unistd.h>    // for getpid, file_size, clock_gettime
#else
#include <process.h>
#include <ciso646>
#define SDSL_UNUSED
#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t>
class int_vector;     // forward declaration

//! A namespace for helper functions
namespace util
{

//============= Debug information =========================

SDSL_UNUSED static bool verbose = false;

void set_verbose();

//============ Manipulating int_vectors ===================

//! Sets all bits of the int_vector to pseudo-random bits.
/*! \param v The int_vector whose bits should be set to random bits
 *  \param seed If seed = 0, the time is used to initialize the
 *              pseudo random number generator, otherwise the seed
 *              parameter is used.
 */
template<class t_int_vec>
void set_random_bits(t_int_vec& v, int seed=0);
//! Sets all bits of the int_vector to 0-bits.
template<class t_int_vec>
void _set_zero_bits(t_int_vec& v);
//! Sets all bits of the int_vector to 1-bits.
template<class t_int_vec>
void _set_one_bits(t_int_vec& v);

//! Bit compress the int_vector
/*! Determine the biggest value X and then set the
 *  int_width to the smallest possible so that we
 *  still can represent X
 */
template<class t_int_vec>
void bit_compress(t_int_vec& v);

//! Expands the integer width to new_width >= v.width()
template<class t_int_vec>
void expand_width(t_int_vec& v, uint8_t new_width);

//! All elements of v modulo m
template<class t_int_vec>
void mod(t_int_vec& v, typename t_int_vec::size_type m);

//! Set all entries of int_vector to value k
/*! \param  v The int_vector which should be set
 *  \param  k The value which should be inserted into v.
 *  \par Details
 *   This method pre-calculates the content of at most 64
 *   words and then repeatedly inserts these words into v.
 */
template<class t_int_vec>
void set_to_value(t_int_vec& v, uint64_t k);

//! Sets each entry of the numerical vector v at position \$fi\f$ to value \$fi\$f
template<class t_int_vec>
void set_to_id(t_int_vec& v);

//! Number of set bits in v.
/*! \param v  int_vector object.
      \return The number of 1-bits in v.
 */
template<class t_int_vec>
typename t_int_vec::size_type cnt_one_bits(const t_int_vec& v);

//! Number of occurrences of bit pattern `10` in v.
/*! \sa getOneBits, getOneZeroBits
 */
template<class t_int_vec>
typename t_int_vec::size_type cnt_onezero_bits(const t_int_vec& v);

//! Number of occurrences of bit pattern `01` in v.
/*! \sa getOneBits, getZeroOneBits
 */
template <class t_int_vec>
typename t_int_vec::size_type cnt_zeroone_bits(const t_int_vec& v);

//! Get the smallest position \f$i\geq idx\f$ where a bit is set
/*! \param v The int_vector in which the bit is searched
 *  \param idx The start position for the search \f$ 0\leq idx < v.bit_size()\f$
 *  \return The smallest position greater or equal to idx, where corresponding bit is 1 or v.bit_size() if no such position exists
 *  \par Time complexity
 *      \f$ \Order{n} \f$
 */
template <class t_int_vec>
typename t_int_vec::size_type next_bit(const t_int_vec& v, uint64_t idx);

//! Get the greatest position \f$i\leq idx\f$ where a bit is set
/*! \param v The int_vector in which the bit is searched
 *  \param idx The start position for the search \f$ 0\leq idx < v.bit_size()\f$
 *  \return The greatest position smaller or equal to idx, where corresponding bit is 1 or v.bit_size() if no such position exists
 *  \par Time complexity
 *     \f$ \Order{n} \f$
*/
template <class t_int_vec>
typename t_int_vec::size_type prev_bit(const t_int_vec& v, uint64_t idx);

//============= Handling files =============================

//! Get the size of a file in bytes
/*! \param file  Path to a file.
 *  \returns     Size of the specified file in bytes.
 */
size_t file_size(const std::string& file);

//! Returns the basename of a file
/*! \param file  Path to a file.
 *  \returns     Basename of the specified file.
 */
std::string basename(std::string file);

//! Returns the directory of a file. A trailing `/` will be removed.
/*! \param file  Path to a file.
 *  \returns     Directory name part of the specified path.
 */
std::string dirname(std::string file);

//! Demangle the class name of typeid(...).name()
/*!
 * \param name A pointer to the result of typeid(...).name()
 */
std::string demangle(const std::string& name);

//! Demangle the class name of typeid(...).name() and remove the "sdsl::"-prefix, "unsigned int",...
std::string demangle2(const std::string& name);

//! Convert type to string
template<typename T>
std::string to_string(const T& t, int w=1);

//! Transforms the demangled class name of an object to a hash value.
template<class T>
uint64_t hashvalue_of_classname(const T&)
{
    std::hash<std::string> str_hash;
    return str_hash(sdsl::util::demangle2(typeid(T).name()));
}

//! Transforms the demangled class name of an object to a hash value.
template<class T>
std::string class_to_hash(const T& t)
{
    return to_string(hashvalue_of_classname(t));
}

template<class T>
std::string class_name(const T& t)
{
    std::string result = demangle2(typeid(t).name());
    size_t template_pos = result.find("<");
    if (template_pos != std::string::npos) {
        result = result.erase(template_pos);
    }
    return result;
}

//! Get the process id of the current process
uint64_t pid();

// convert an errno number to a readable msg
char* str_from_errno();

class _id_helper
{
private:
    static uint64_t id;
public:
    static uint64_t getId()
    {
        return id++;
    }
};

//! Get a unique id inside the process
uint64_t id();

template<typename T>
std::string to_latex_string(const T& t);

std::string to_latex_string(unsigned char c);

//! Delete all files of the file_map.
void delete_all_files(tMSS& file_map);

// thanks to Stefan Arnold for the assign functions
//! Assigns the value x of type T to the value of y of type U.
/*!
 * \param x    The assigned variable.
 * \param y    The variable which provides the value that is assigned to x.
 */
template<class T, class U>
void assign(T& x, const U& y)
{
    x = T(y);
}

//! Swaps variables x and y.
/*!
 * \param x Reference to the first variable.
 * \param y Reference to the second variable.
 */
template<class T>
void assign(T& x, T& y)
{
    x.swap(y);
}

//! clear the space used by x
/*!
 * \param x Reference to the data structure.
 */
template<class T>
void clear(T& x)
{
    T y;
    x.swap(y);
}

//! Swap support data structure and assign to new vector
/*! \param s1 First support structure.
 *  \param s2 Second support structure.
 *  \param p1 First supported structure.
 *  \param p2 Second supported structure.
 *  s1 is swapped with s2 and after the execution s1 supports p1 and s2 supports
 *  p2. I.e. if p1 and p2 are members of a complex data structure, we have to
 *  swap p1 and p2 before we use this method.
 */
template<class S, class P>
void swap_support(S& s1, S& s2, const P* p1, const P* p2)
{
    s1.swap(s2);
    s1.set_vector(p1);
    s2.set_vector(p2);
}

//! Initialise support data structure with
/*! \param s Support structure which should be initialized
 *  \param x Pointer to the data structure which should be supported.
 */
template<class S, class X>
void init_support(S& s, const X* x)
{
    S temp(x);       // generate a temporary support object
    s.swap(temp);    // swap its content with the target object
    s.set_vector(x); // set the support object's  pointer to x
}

class spin_lock
{
private:
    std::atomic_flag m_slock;
public:
    spin_lock()
    {
        m_slock.clear();
    }
    void lock()
    {
        while (m_slock.test_and_set(std::memory_order_acquire)) {
            /* spin */
        }
    };
    void unlock()
    {
        m_slock.clear(std::memory_order_release);
    };
};

//! Create 2^{log_s} random integers mod m with seed x
/*
 */
template<class t_int_vec>
t_int_vec rnd_positions(uint8_t log_s, uint64_t& mask, uint64_t mod=0, uint64_t seed=17)
{
    mask = (1<<log_s)-1;
    t_int_vec rands(1<<log_s ,0);
    set_random_bits(rands, seed);
    if (mod > 0) {
        util::mod(rands, mod);
    }
    return rands;
}

//! Checks at compile time whether type is regular or not
/*  static_assert(is_regular<YOUR_TYPE>::value);
 *  Code is from a talk of Aerix Consulting
 */
template<typename T>
struct is_regular : std::integral_constant< bool,
                                            std::is_default_constructible<T>::value&&
                                                std::is_copy_constructible<T>::value&&
                                                std::is_move_constructible<T>::value&&
                                                std::is_copy_assignable<T>::value&&
                                                std::is_move_assignable<T>::value > {};

} // end namespace util

//==================== Template functions ====================

template<class t_int_vec>
void util::set_random_bits(t_int_vec& v, int seed)
{
    std::mt19937_64 rng;
    if (0 == seed) {
        rng.seed(std::chrono::system_clock::now().time_since_epoch().count() + util::id());
    } else
        rng.seed(seed);

    uint64_t* data = v.data();
    if (v.empty())
        return;
    *data = rng();
    for (typename t_int_vec::size_type i=1; i < (v.capacity()>>6); ++i) {
        *(++data) = rng();
    }
}

// all elements of vector v modulo m
template<class t_int_vec>
void util::mod(t_int_vec& v, typename t_int_vec::size_type m)
{
    for (typename t_int_vec::size_type i=0; i < v.size(); ++i) {
        v[i] = v[i] % m;
    }
}

template<class t_int_vec>
void util::bit_compress(t_int_vec& v)
{
    auto max_elem = std::max_element(v.begin(),v.end());
    uint64_t max = 0;
    if (max_elem != v.end()) {
        max = *max_elem;
    }
    uint8_t min_width = bits::hi(max)+1;
    uint8_t old_width = v.width();
    if (old_width > min_width) {
        const uint64_t* read_data = v.data();
        uint64_t* write_data = v.data();
        uint8_t read_offset = 0;
        uint8_t write_offset = 0;
        for (typename t_int_vec::size_type i=0; i < v.size(); ++i) {
            uint64_t x = bits::read_int_and_move(read_data, read_offset, old_width);
            bits::write_int_and_move(write_data,  x, write_offset, min_width);
        }
        v.bit_resize(v.size()*min_width);
        v.width(min_width);
    }
}

template<class t_int_vec>
void util::expand_width(t_int_vec& v, uint8_t new_width)
{
    uint8_t old_width = v.width();
    typename t_int_vec::size_type n = v.size();
    if (new_width > old_width) {
        if (n > 0) {
            typename t_int_vec::size_type i, old_pos, new_pos;
            new_pos = (n-1)*new_width;
            old_pos = (n-1)*old_width;
            v.bit_resize(v.size()*new_width);
            for (i=0; i < n; ++i, new_pos-=new_width, old_pos-=old_width) {
                v.set_int(new_pos, v.get_int(old_pos, old_width), new_width);
            }
        }
        v.width(new_width);
    }
}

template<class t_int_vec>
void util::_set_zero_bits(t_int_vec& v)
{
    uint64_t* data = v.data();
    if (v.empty())
        return;
    // TODO: replace by memset() but take care of size_t in the argument!
    *data = 0ULL;
    for (typename t_int_vec::size_type i=1; i < (v.capacity()>>6); ++i) {
        *(++data) = 0ULL;
    }
}

template<class t_int_vec>
void util::_set_one_bits(t_int_vec& v)
{
    uint64_t* data = v.data();
    if (v.empty())
        return;
    *data = 0xFFFFFFFFFFFFFFFFULL;
    for (typename t_int_vec::size_type i=1; i < (v.capacity()>>6); ++i) {
        *(++data) = 0xFFFFFFFFFFFFFFFFULL;
    }
}

template<class t_int_vec>
void util::set_to_value(t_int_vec& v, uint64_t k)
{
    uint64_t* data = v.data();
    if (v.empty())
        return;
    uint8_t int_width = v.width();
    if (int_width == 0) {
        throw std::logic_error("util::set_to_value can not be performed with int_width=0!");
    }
    if (0 == k) {
        _set_zero_bits(v);
        return;
    }
    if (bits::lo_set[int_width] == k) {
        _set_one_bits(v);
        return;
    }
    k = k & (0xFFFFFFFFFFFFFFFFULL >> (64-int_width));
    uint64_t vec[67] = {0}; // allocate memory for the mask and initialize with zeros
    vec[0] = 0;
    uint8_t offset = 0;
    uint64_t n=0, vals=0;
    do { // loop terminates after at most 64 iterations
        vec[n] = vec[n] | (k << offset);
        offset += int_width;
        vals++;
        if (offset >= 64) {
            vec[n+1] = 0;
            vec[++n] = k >> (int_width-(offset-64));
            offset -= 64;
        }
    } while (offset != 0);

    typename t_int_vec::size_type n64 = v.capacity()/64;
    for (typename t_int_vec::size_type i=0; i < n64;) {
        for (uint64_t ii=0; ii < n and i < n64; ++ii,++i) {
            *(data++) = vec[ii];
        }
    }
}

//! Set v[i] = i for i=[0..v.size()-1]
template<class t_int_vec>
void util::set_to_id(t_int_vec& v)
{
    std::iota(v.begin(), v.end(), 0ULL);
}

template<class t_int_vec>
typename t_int_vec::size_type util::cnt_one_bits(const t_int_vec& v)
{
    const uint64_t* data = v.data();
    if (v.empty())
        return 0;
    typename t_int_vec::size_type result = bits::cnt(*data);
    for (typename t_int_vec::size_type i=1; i < (v.capacity()>>6); ++i) {
        result += bits::cnt(*(++data));
    }
    if (v.bit_size()&0x3F) {
        result -= bits::cnt((*data) & (~bits::lo_set[v.bit_size()&0x3F]));
    }
    return result;
}

template<class t_int_vec>
typename t_int_vec::size_type util::cnt_onezero_bits(const t_int_vec& v)
{
    const uint64_t* data = v.data();
    if (v.empty())
        return 0;
    uint64_t carry = 0, oldcarry=0;
    typename t_int_vec::size_type result = bits::cnt10(*data, carry);
    for (typename t_int_vec::size_type i=1; i < (v.capacity()>>6); ++i) {
        oldcarry = carry;
        result += bits::cnt10(*(++data), carry);
    }
    if (v.bit_size()&0x3F) {// if bit_size is not a multiple of 64, subtract the counts of the additional bits
        result -= bits::cnt(bits::map10(*data, oldcarry) & bits::lo_unset[v.bit_size()&0x3F]);
    }
    return result;
}

template<class t_int_vec>
typename t_int_vec::size_type util::cnt_zeroone_bits(const t_int_vec& v)
{
    const uint64_t* data = v.data();
    if (v.empty())
        return 0;
    uint64_t carry = 1, oldcarry = 1;
    typename t_int_vec::size_type result = bits::cnt01(*data, carry);
    for (typename t_int_vec::size_type i=1; i < (v.capacity()>>6); ++i) {
        oldcarry = carry;
        result += bits::cnt01(*(++data), carry);
    }
    if (v.bit_size()&0x3F) {// if bit_size is not a multiple of 64, subtract the counts of the additional bits
        result -= bits::cnt(bits::map01(*data, oldcarry) & bits::lo_unset[v.bit_size()&0x3F]);
    }
    return result;
}

template <class t_int_vec>
typename t_int_vec::size_type util::next_bit(const t_int_vec& v, uint64_t idx)
{
    uint64_t pos = idx>>6;
    uint64_t node = v.data()[pos];
    node >>= (idx&0x3F);
    if (node) {
        return idx+bits::lo(node);
    } else {
        ++pos;
        while ((pos<<6) < v.bit_size()) {
            if (v.data()[pos]) {
                return (pos<<6)|bits::lo(v.data()[pos]);
            }
            ++pos;
        }
        return v.bit_size();
    }
}

template <class t_int_vec>
typename t_int_vec::size_type util::prev_bit(const t_int_vec& v, uint64_t idx)
{
    uint64_t pos = idx>>6;
    uint64_t node = v.data()[pos];
    node <<= 63-(idx&0x3F);
    if (node) {
        return bits::hi(node)+(pos<<6)-(63-(idx&0x3F));
    } else {
        --pos;
        while ((pos<<6) < v.bit_size()) {
            if (v.data()[pos]) {
                return (pos<<6)|bits::hi(v.data()[pos]);
            }
            --pos;
        }
        return v.bit_size();
    }
}

template<typename T>
std::string util::to_string(const T& t, int w)
{
    std::stringstream ss;
    ss<<std::setw(w)<<t;
    return ss.str();
}

template<typename T>
std::string util::to_latex_string(const T& t)
{
    return to_string(t);
}

}// end namespace sdsl
#endif

/* sdsl - succinct data structures library
    Copyright (C) 2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file io.hpp
    \brief io.hpp contains some methods for reading/writing sdsl structures.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_IO
#define INCLUDED_SDSL_IO

/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file sdsl_concepts.hpp
    \brief Contains declarations and definitions of data structure concepts.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CONCEPTS
#define INCLUDED_SDSL_CONCEPTS

namespace sdsl
{

struct bv_tag {}; // bitvector tag
struct iv_tag {}; // int_vector tag

struct csa_tag {}; // compressed suffix array (CSAs) tag
struct cst_tag {}; // compressed suffix tree (CST) tag
struct wt_tag {};  // wavelet tree tag

struct psi_tag {}; // tag for CSAs based on the psi function
struct lf_tag {}; // tag for CSAs based on the LF function

struct csa_member_tag {}; // tag for text, bwt, LF, \Psi members of CSA

struct lcp_tag {};
struct lcp_plain_tag {};
struct lcp_permuted_tag {};
struct lcp_tree_compressed_tag {};
struct lcp_tree_and_lf_compressed_tag {};

struct alphabet_tag {};
struct byte_alphabet_tag { static const uint8_t WIDTH=8; };
struct int_alphabet_tag { static const uint8_t WIDTH=0; };

struct sa_sampling_tag {};
struct isa_sampling_tag {};

template<class t_T, class t_r = void>
struct enable_if_type {
    typedef t_r type;
};

template<class t_idx, class t_enable = void>
struct index_tag {
    typedef t_enable type;
};

template<class t_idx>
struct index_tag<t_idx, typename enable_if_type<typename t_idx::index_category>::type> {
    using type = typename t_idx::index_category;
};

template<class t_sampling, class t_enable = void>
struct sampling_tag {
    typedef t_enable type;
};

template<class t_sampling>
struct sampling_tag<t_sampling, typename enable_if_type<typename t_sampling::sampling_category>::type> {
    using type = typename t_sampling::sampling_category;
};

template<class t_enc_vec, class t_enable = void>
struct is_enc_vec {
    static const bool value = false;
};

template<class t_enc_vec>
struct is_enc_vec<t_enc_vec, typename enable_if_type<typename t_enc_vec::enc_vec_type>::type> {
    static const bool value = true;
};

template<class t_alphabet, class t_enable = void>
struct is_alphabet {
    static const bool value = false;
};

template<class t_alphabet>
struct is_alphabet<t_alphabet, typename enable_if_type<typename t_alphabet::alphabet_category>::type> {
    static const bool value = true;
};

} // end namespace sdsl

#endif

#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <cctype>

namespace sdsl
{

template<typename T>
void load_vector(std::vector<T>&, std::istream&);

template<class T>
uint64_t
serialize_vector(const std::vector<T>&, std::ostream&,
                 sdsl::structure_tree_node* v=nullptr, std::string="");

// has_serialize<X>::value is true if class X has
// implement method serialize
// Adapted solution from jrok's proposal:
// http://stackoverflow.com/questions/87372/check-if-a-class-has-a-member-function-of-a-given-signature
template<typename X>
struct has_serialize {
    template<typename T>
    static constexpr auto check(T*)
    -> typename
    std::is_same<
        decltype(std::declval<T>().serialize(
            std::declval<std::ostream&>(),
            std::declval<structure_tree_node*>(),
            std::declval<std::string>()
        )),
        typename T::size_type>::type {return std::true_type();}
    template<typename>
    static constexpr std::false_type check(...) {return std::false_type();}
    typedef decltype(check<X>(nullptr)) type;
    static constexpr bool value = type::value;
};

// has_load<X>::value is true if class X has
// implement method load
template<typename X>
struct has_load {
    template<typename T>
    static constexpr auto check(T*)
    -> typename
    std::is_same<
        decltype(std::declval<T>().load(
            std::declval<std::istream&>()
        )),
        void>::type {return std::true_type();}
    template<typename>
    static constexpr std::false_type check(...) {return std::false_type();}
    typedef decltype(check<X>(nullptr)) type;
    static constexpr bool value = type::value;
};

// Writes primitive-typed variable t to stream out
template<class T>
size_t write_member(const T& t, std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")
{
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, util::class_name(t));
    out.write((char*)&t, sizeof(t));
    size_t written_bytes = sizeof(t);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

// Specialization for std::string
template<>
size_t write_member<std::string>(const std::string& t, std::ostream& out, sdsl::structure_tree_node* v, std::string name);

// Writes primitive-typed variable t to stream out
template<class T>
void read_member(T& t, std::istream& in)
{
    in.read((char*)&t, sizeof(t));
}

// Specialization for std::string
template<>
void read_member<std::string>(std::string& t, std::istream& in);

template<typename X>
typename std::enable_if<has_serialize<X>::value,typename X::size_type>::type
serialize(const X& x,
          std::ostream& out, structure_tree_node* v=nullptr,
          std::string name="")
{
    return x.serialize(out, v, name);
}

template<typename X>
typename std::enable_if<std::is_pod<X>::value,uint64_t>::type
serialize(const X& x,
          std::ostream& out, structure_tree_node* v=nullptr,
          std::string name="")
{
    return write_member(x, out, v, name);
}

template<typename X>
uint64_t
serialize(const std::vector<X>& x,
          std::ostream& out, structure_tree_node* v=nullptr,
          std::string name="")
{

    return serialize(x.size(), out, v, name)
        + serialize_vector(x, out, v, name);
}

template<typename X>
typename std::enable_if<has_load<X>::value,void>::type
load(X& x, std::istream& in)
{
    x.load(in);
}

template<typename X>
typename std::enable_if<std::is_pod<X>::value,void>::type
load(X& x, std::istream& in)
{
    read_member(x, in);
}

template<typename X>
void load(std::vector<X>& x, std::istream& in)
{
    typename std::vector<X>::size_type size;
    load(size, in);
    x.resize(size);
    load_vector(x, in);
}

//! Load sdsl-object v from a file.
/*!
 * \param v sdsl-
 * \param file Name of the serialized file.
 */
template<class T>
bool load_from_file(T& v, const std::string& file);

//! Load an int_vector from a plain array of `num_bytes`-byte integers with X in \{0, 1,2,4,8\} from disk.
// TODO: Remove ENDIAN dependency.
template<class t_int_vec>
bool load_vector_from_file(t_int_vec& v, const std::string& file, uint8_t num_bytes=1, uint8_t max_int_width=64)
{
    if ((uint8_t)0 == num_bytes) {  // if byte size is variable read int_vector<0> from file
        return load_from_file(v, file);
    } else if (num_bytes == 'd') {
        uint64_t x = 0, max_x = 0;
        isfstream in(file, std::ios::in | std::ios::binary);
        if (!in) {
            return false;
        } else {
            std::vector<uint64_t> tmp;
            while (in >> x) {
                tmp.push_back(x);
                max_x = std::max(x, max_x);
            }
            v.width(bits::hi(max_x)+1); v.resize(tmp.size());
            for (size_t i=0; i < tmp.size(); ++i) {
                v[i] = tmp[i];
            }
            return true;
        }
    } else {
        off_t file_size = util::file_size(file);
        if (file_size == 0) {
            v.resize(0);
            return true;
        }
        if (file_size % num_bytes != 0) {
            throw std::logic_error("file size "+util::to_string(file_size)+" of \""+ file
                                       +"\" is not a multiple of "+util::to_string(num_bytes));
            return false;
        }
        isfstream in(file, std::ios::in | std::ios::binary);
        if (in) {
            v.width(std::min((int)8*num_bytes, (int)max_int_width));
            v.resize(file_size / num_bytes);
            if (8 == t_int_vec::fixed_int_width and 1 == num_bytes) {  // if int_vector<8> is created from byte alphabet file
                in.read((char*)v.data(), file_size);
            } else {
                size_t idx=0;
                const size_t block_size = conf::SDSL_BLOCK_SIZE*num_bytes;
                std::vector<uint8_t> buf(block_size);
                // TODO: check for larger alphabets with num_bytes*8 = v::fixed_int_width

                uint64_t x = 0; // value
                uint8_t  cur_byte = 0;
                do {
                    in.read((char*)buf.data(), block_size);
                    size_t read = in.gcount();
                    uint8_t* begin = buf.data();
                    uint8_t* end   = begin+read;
                    while (begin < end) {
                        x |= ((uint64_t)(*begin)) << (cur_byte*8);
                        ++cur_byte;
                        if (cur_byte == num_bytes) {
                            v[idx++] = x;
                            cur_byte = 0;
                            x = 0ULL;
                        }
                        ++begin;
                    }
                } while (idx < v.size());
                in.close();
            }
            return true;
        } else {
            return false;
        }
    }
}

//! Store a data structure to a file.
/*! The data structure has to provide a serialize function.
 *  \param v Data structure to store.
 *  \param file Name of the file where to store the data structure.
 *  \param Return if the data structure was stored successfully
 */
template<class T>
bool store_to_file(const T& v, const std::string& file);

//! Specialization of store_to_file for a char array
bool store_to_file(const char* v, const std::string& file);

//! Specialization of store_to_file for int_vector
template<uint8_t t_width>
bool store_to_file(const int_vector<t_width>& v, const std::string& file, bool write_fixed_as_variable=false);

//! Store an int_vector as plain int_type array to disk
template<class int_type, class t_int_vec>
bool store_to_plain_array(t_int_vec& v, const std::string& file)
{
    osfstream out(file, std::ios::out | std::ios::binary);
    if (out) {
        for (typename t_int_vec::size_type i=0; i<v.size(); ++i) {
            int_type x = v[i];
            out.write((char*)&x, sizeof(int_type));
        }
        return true;
    } else {
        return false;
    }
}

template<class T>
size_t serialize_empty_object(std::ostream&, structure_tree_node* v=nullptr, std::string name="", const T* t=nullptr)
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*t));
    size_t written_bytes = 0;
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

//! Get the size of a data structure in bytes.
/*!
 *  \param v A reference to the data structure for which the size in bytes should be calculated.
 */
template<class T>
typename T::size_type size_in_bytes(const T& t);

//! Get the size of a data structure in mega bytes (MiB).
/*!
 *  \param t A reference to the data structure for which the size in bytes should be calculated.
 */
template<class T>
double size_in_mega_bytes(const T& t);

struct nullstream : std::ostream {
    struct nullbuf: std::streambuf {
        int overflow(int c)
        {
            return traits_type::not_eof(c);
        }
        int xputc(int) { return 0; }
        std::streamsize xsputn(char const*, std::streamsize n) { return n; }
        int sync() { return 0; }
    } m_sbuf;
    nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf), m_sbuf() {}
};

//! Serialize each element of an std::vector
/*!
 * \param vec The vector which should be serialized.
 * \param out Output stream to which should be written.
 * \param v   Structure tree node. Note: If all elements have the same
 *            structure, then it is tried to combine all elements (i.e.
 *            make one node w with size set to the cumulative sum of all
 *           sizes of the children)
 */
template<class T>
uint64_t
serialize_vector(const std::vector<T>& vec, std::ostream& out, sdsl::structure_tree_node* v, std::string name)
{
    if (vec.size() > 0) {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, "std::vector<"+util::class_name(vec[0])+">");
        size_t written_bytes = 0;
        for (const auto& x : vec) {
            written_bytes += serialize(x, out, child, "[]");
        }
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    } else {
        return 0;
    }
}

//! Load all elements of a vector from a input stream
/*! \param vec  Vector whose elements should be loaded.
 *  \param in   Input stream.
 *  \par Note
 *   The vector has to be resized prior the loading
 *   of its elements.
 */
template<class T>
void load_vector(std::vector<T>& vec, std::istream& in)
{
    for (typename std::vector<T>::size_type i = 0; i < vec.size(); ++i) {
        load(vec[i], in);
    }
}

template<format_type F, typename X>
void write_structure(const X& x, std::ostream& out)
{
    std::unique_ptr<structure_tree_node> st_node(new structure_tree_node("name","type"));
    nullstream ns;
    serialize(x, ns, st_node.get(), "");
    if (st_node.get()->children.size() > 0) {
        for (const auto& child: st_node.get()->children) {
            sdsl::write_structure_tree<F>(child.second.get(), out);
        }
    }
}

template<format_type F, typename X>
void write_structure(const X& x, std::string file)
{
    std::ofstream out(file);
    write_structure<F>(x, out);
}

template<format_type F, typename... Xs>
void write_structure(std::ostream& out, Xs... xs)
{
    typedef std::unique_ptr<structure_tree_node> up_stn_type;
    up_stn_type st_node(new structure_tree_node("name","type"));
    _write_structure(st_node, xs...);
    sdsl::write_structure_tree<F>(st_node.get(), out);
}

template<typename X, typename... Xs>
void _write_structure(std::unique_ptr<structure_tree_node>& st_node, X x, Xs... xs)
{
    nullstream ns;
    serialize(x, ns, st_node.get(), "");
    _write_structure(st_node, xs...);
}

inline
void _write_structure(std::unique_ptr<structure_tree_node>&) {}

//! Internal function used by csXprintf
uint64_t _parse_number(std::string::const_iterator& c, const std::string::const_iterator& end);

//! Internal function used by csXprintf
template<class t_csa>
const t_csa& _idx_csa(const t_csa& t, csa_tag)
{
    return t;
}

//! Internal function used by csXprintf
template<class t_cst>
const typename t_cst::csa_type& _idx_csa(const t_cst& t, cst_tag)
{
    return t.csa;
}

//! Internal function used by csXprintf
template<class t_csa>
std::string _idx_lcp_val(const t_csa&, uint64_t, uint64_t, csa_tag)
{
    return "";
}

//! Internal function used by csXprintf
template<class t_cst>
std::string _idx_lcp_val(const t_cst& t, uint64_t i, uint64_t w, cst_tag)
{
    return util::to_string(t.lcp[i], w);
}

template<class t_csx, class t_alph=typename t_csx::alphabet_category>
struct default_sentinel {
    static const char value = '$';
};

template<class t_csx>
struct default_sentinel<t_csx, byte_alphabet_tag> {
    static const char value = '$';
};

template<class t_csx>
struct default_sentinel<t_csx, int_alphabet_tag> {
    static const char value = '0';
};

//! Prints members of CSAs and CSTs
/*! This is a printf like method to write members of CSAs and CSTs into an outstream.
 * \tparam t_idx   Type of the index. Class should be of concept csa_tag or cst_tag.
 * \param out      Output stream.
 * \param format   Format string. See explanation below.
 * \param idx      CSA or CST object.
 * \param sentinel Character which should replace the 0-symbol in BWT/ TEXT.
 *
 * \par Format string
 *   Each line of the output will be formatted according to the format string.
 *   All content, except tokens which start with `%` will be copied. Tokens
 *   which start with `%` will be replaced as follows (let w be a positive
 *    number. setw(w) is used to format single numbers):
 *
 *      Token      |  Replacement | Comment
 *      -----------------------------------------------------------------------
 *       %[w]I     | Row index i.                           |
 *       %[w]S     | SA[i]                                  |
 *       %[w]s     | ISA[i]                                 |
 *       %[w]P     | PSI[i]                                 |
 *       %[w]p     | LF[i]                                  |
 *       %[w]L     | LCP[i]                                 | only for CSTs
 *       %[w]B     | BWT[i]                                 |
 *       %[w[:W]]T | Print min(idx.size(),w) chars of each  |
 *                 | suffix, each char formatted by setw(W).|
 *       %%        | %                                      |
 */
template<class t_idx>
void
csXprintf(std::ostream& out, const std::string& format,
          const t_idx& idx, char sentinel=default_sentinel<t_idx>::value)
{
    typename t_idx::index_category cat;
    const typename t_idx::csa_type& csa = _idx_csa(idx, cat);
    std::vector<std::string> res(csa.size());
    for (std::string::const_iterator c = format.begin(), s=c; c != format.end(); s=c) {
        while (c != format.end() and* c != '%') ++c;   // string before the next `%`
        if (c > s) {  // copy format string part
            std::vector<std::string> to_copy(csa.size(), std::string(s, c));
            transform(res.begin(), res.end(), to_copy.begin(), res.begin(), std::plus<std::string>());
        }
        if (c == format.end()) break;
        ++c; // skip `%`
        uint64_t w = _parse_number(c, format.end());  // element width
        if (c == format.end()) break;
        uint64_t W = 0; // character width
        if (':' == *c) {
            ++c;
            W = _parse_number(c, format.end());
        }
        if (c == format.end()) break;
        for (uint64_t i=0; i<csa.size(); ++i) {
            switch (*c) {
                case 'I': res[i] += util::to_string(i,w); break;
                case 'S': res[i] += util::to_string(csa[i],w); break;
                case 's': res[i] += util::to_string(csa.isa[i],w); break;
                case 'P': res[i] += util::to_string(csa.psi[i],w); break;
                case 'p': res[i] += util::to_string(csa.lf[i],w); break;
                case 'L': res[i] += _idx_lcp_val(idx,i,w, cat); break;
                case 'B': if (0 == csa.bwt[i]) {
                        res[i] += util::to_string(sentinel,w);
                    } else {
                        res[i] += util::to_string(csa.bwt[i],w);
                    }
                    break;
                case 'T': for (uint64_t k=0; (w>0 and k < w) or(0==w and k < csa.size()); ++k) {
                        if (0 == csa.text[(csa[i]+k)%csa.size()]) {
                            res[i] += util::to_string(sentinel, W);
                        } else {
                            res[i] += util::to_string(csa.text[(csa[i]+k)%csa.size()], W);
                        }
                    }
                    break;
                case '%': res[i] += "%"; break;
            }
        }
        ++c;
    }
    for (size_t i=0; i<res.size(); ++i) out << res[i] << std::endl;
}

//! Returns the file name of the resource.
/*!
 * \param  key        Resource key.
 * \param  config    Cache configuration.
 * \return The file name of the resource.
 */
std::string cache_file_name(const std::string& key,const cache_config& config);

//! Returns the file name of the resource.
/*!
 * \param  key        Resource key.
 * \param  config    Cache configuration.
 * \return The file name of the resource.
 */
template<typename T>
std::string cache_file_name(const std::string& key,const cache_config& config)
{
    return cache_file_name(key+"_"+util::class_to_hash(T()), config);
}

//! Register the existing resource specified by the key to the cache
/*!
 *  \param key        Resource key.
 *  \param config    Cache configuration.
 *
 *  Note: If the resource does not exist under the given key,
 *  it will be not added to the cache configuration.
 */
void register_cache_file(const std::string& key, cache_config& config);

//! Checks if the resource specified by the key exists in the cache.
/*!
  \param key    Resource key.
  \param config Cache configuration.
  \return True, if the file exists, false otherwise.
*/
bool cache_file_exists(const std::string& key, const cache_config& config);

//! Checks if the resource specified by the key and type exists in the cache.
/*!
  \tparam T     Type.
  \param key    Resource key.
  \param config Cache configuration.
  \return True, if the file exists, false otherwise.
*/
template<typename T>
bool cache_file_exists(const std::string& key, const cache_config& config)
{
    return cache_file_exists(key+"_"+util::class_to_hash(T()), config);
}

//! Returns a name for a temporary file. I.e. the name was not used before.
std::string tmp_file(const cache_config& config, std::string name_part="");

//! Returns a name for a temporary file. I.e. the name was not used before.
std::string tmp_file(const std::string& filename, std::string name_part="");

template<class T>
bool load_from_cache(T& v, const std::string& key, const cache_config& config, bool add_type_hash=false)
{
    std::string file;
    if (add_type_hash) {
        file = cache_file_name<T>(key, config);
    } else {
        file = cache_file_name(key, config);
    }
    if (load_from_file(v, file)) {
        if (util::verbose) {
            std::cerr << "Load `" << file << std::endl;
        }
        return true;
    } else {
        std::cerr << "WARNING: Could not load file '";
        std::cerr << file << "'" << std::endl;
        return false;
    }
}

//! Stores the object v as a resource in the cache.
/*!
 *  \param
 */
template<class T>
bool store_to_cache(const T& v, const std::string& key, cache_config& config, bool add_type_hash=false)
{
    std::string file;
    if (add_type_hash) {
        file = cache_file_name<T>(key, config);
    } else {
        file = cache_file_name(key, config);
    }
    if (store_to_file(v, file)) {
        config.file_map[std::string(key)] = file;
        return true;
    } else {
        std::cerr<<"WARNING: store_to_cache: could not store file `"<< file <<"`" << std::endl;
        return false;
    }
}

//==================== Template functions ====================

template<class T>
typename T::size_type size_in_bytes(const T& t)
{
    nullstream ns;
    return serialize(t, ns);
}

template<class T>
double size_in_mega_bytes(const T& t)
{
    return size_in_bytes(t)/(1024.0*1024.0);
}

template<class T>
void add_hash(const T& t, std::ostream& out)
{
    uint64_t hash_value = util::hashvalue_of_classname(t);
    write_member(hash_value, out);
}

template<class T>
bool store_to_file(const T& t, const std::string& file)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        if (util::verbose) {
            std::cerr<<"ERROR: store_to_file not successful for: `"<<file<<"`"<<std::endl;
        }
        return false;
    }
    serialize(t,out);
    out.close();
    if (util::verbose) {
        std::cerr<<"INFO: store_to_file: `"<<file<<"`"<<std::endl;
    }
    return true;
}

template<class T>
bool store_to_checked_file(const T& t, const std::string& file)
{
    std::string checkfile = file+"_check";
    osfstream out(checkfile, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        if (util::verbose) {
            std::cerr<<"ERROR: store_to_checked_file not successful for: `"<<checkfile<<"`"<<std::endl;
        }
        return false;
    }
    add_hash(t, out);
    out.close();
    return store_to_file(t, file);
}

bool store_to_file(const char* v, const std::string& file);

bool store_to_file(const std::string& v, const std::string& file);

template<uint8_t t_width>
bool store_to_file(const int_vector<t_width>& v, const std::string& file, bool write_fixed_as_variable)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        std::cerr<<"ERROR: util::store_to_file:: Could not open file `"<<file<<"`"<<std::endl;
        return false;
    } else {
        if (util::verbose) {
            std::cerr<<"INFO: store_to_file: `"<<file<<"`"<<std::endl;
        }
    }
    v.serialize(out, nullptr, "", write_fixed_as_variable);
    out.close();
    return true;
}

template<uint8_t t_width>
bool store_to_checked_file(const int_vector<t_width>& v, const std::string& file, bool write_fixed_as_variable)
{
    std::string checkfile = file+"_check";
    osfstream out(checkfile, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        std::cerr<<"ERROR: util::store_to_checked_file: Could not open check file `"<<checkfile<<"`"<<std::endl;
        return false;
    } else {
        if (util::verbose) {
            std::cerr<<"INFO: store_to_checked_file: `"<<checkfile<<"`"<<std::endl;
        }
    }
    add_hash(v, out);
    out.close();
    return store_to_file(v, file, write_fixed_as_variable);
}

template<class T>
bool load_from_file(T& v, const std::string& file)
{
    isfstream in(file, std::ios::binary | std::ios::in);
    if (!in) {
        if (util::verbose) {
            std::cerr << "Could not load file `" << file << "`" << std::endl;
        }
        return false;
    }
    load(v, in);
    in.close();
    if (util::verbose) {
        std::cerr << "Load file `" << file << "`" << std::endl;
    }
    return true;
}

template<class T>
bool load_from_checked_file(T& v, const std::string& file)
{
    isfstream in(file+"_check", std::ios::binary | std::ios::in);
    if (!in) {
        if (util::verbose) {
            std::cerr << "Could not load check file `" << file << "_check`" << std::endl;
        }
        return false;
    }
    uint64_t hash_value;
    read_member(hash_value, in);
    if (hash_value != util::hashvalue_of_classname(v)) {
        if (util::verbose) {
            std::cerr << "File `" << file << "` is not an instance of the class `" << sdsl::util::demangle2(typeid(T).name()) << "`" << std::endl;
        }
        return false;
    }
    return load_from_file(v, file);
}

template<class t_iv>
inline typename std::enable_if<
    std::is_same<typename t_iv::index_category, iv_tag>::value or
        std::is_same<typename t_iv::index_category, csa_tag>::value or
        std::is_same<typename t_iv::index_category, lcp_tag>::value
    , std::ostream&>::type
operator<<(std::ostream& os, const t_iv& v)
{
    for (auto it=v.begin(), end = v.end(); it != end; ++it) {
        os << *it;
        if (it+1 != end) os << " ";
    }
    return os;
}

template<class t_iv>
inline typename std::enable_if<std::is_same<typename t_iv::index_category ,wt_tag>::value, std::ostream&>::type
operator<<(std::ostream& os, const t_iv& v)
{
    for (auto it=v.begin(), end = v.end(); it != end; ++it) {
        os << *it;
        if (it+1 != end and std::is_same<typename t_iv::alphabet_category,int_alphabet_tag>::value) os << " ";
    }
    return os;
}

template<class t_int>
inline typename std::enable_if<std::is_integral<t_int>::value, std::ostream&>::type
operator<<(std::ostream& os, const std::vector<t_int>& v)
{
    for (auto it=v.begin(), end = v.end(); it != end; ++it) {
        os << *it;
        if (it+1 != end) os << " ";
    }
    return os;
}

template<class t_iv>
inline typename std::enable_if<std::is_same<typename t_iv::category ,csa_member_tag>::value, std::ostream&>::type
operator<<(std::ostream& os, const t_iv& v)
{
    for (auto it=v.begin(), end = v.end(); it != end; ++it) {
        os << *it;
        if (it+1 != end and std::is_same<typename t_iv::alphabet_category,int_alphabet_tag>::value) os << " ";
    }
    return os;
}

}
#endif

/*!\file memory_management.hpp
\brief memory_management.hpp contains two function for allocating and deallocating memory
\author Simon Gog
*/
#ifndef INCLUDED_SDSL_MEMORY_MANAGEMENT
#define INCLUDED_SDSL_MEMORY_MANAGEMENT

#include <map>
#include <iostream>
#include <cstdlib>
#include <mutex>
#include <chrono>
#include <cstring>
#include <set>
#include <cstddef>
#include <stack>
#include <vector>

#include <fcntl.h>

#ifdef MSVC_COMPILER
// windows.h has min/max macro which causes problems when using std::min/max
#define NOMINMAX
#include <windows.h>
#include <io.h>
#else
#include <sys/mman.h>
#endif

namespace sdsl
{

class memory_monitor;

template<format_type F>
void write_mem_log(std::ostream& out, const memory_monitor& m);

class memory_monitor
{
public:
    using timer = std::chrono::high_resolution_clock;
    struct mm_alloc {
        timer::time_point timestamp;
        int64_t usage;
        mm_alloc(timer::time_point t, int64_t u) : timestamp(t), usage(u) {};
    };
    struct mm_event {
        std::string name;
        std::vector<mm_alloc> allocations;
        mm_event(std::string n, int64_t usage) : name(n)
        {
            allocations.emplace_back(timer::now(), usage);
        };
        bool operator< (const mm_event& a) const
        {
            if (a.allocations.size() && this->allocations.size()) {
                if (this->allocations[0].timestamp == a.allocations[0].timestamp) {
                    return this->allocations.back().timestamp < a.allocations.back().timestamp;
                } else {
                    return this->allocations[0].timestamp < a.allocations[0].timestamp;
                }
            }
            return true;
        }
    };
    struct mm_event_proxy {
        bool add;
        timer::time_point created;
        mm_event_proxy(const std::string& name, int64_t usage, bool a) : add(a)
        {
            if (add) {
                auto& m = the_monitor();
                std::lock_guard<util::spin_lock> lock(m.spinlock);
                m.event_stack.emplace(name, usage);
            }
        }
        ~mm_event_proxy()
        {
            if (add) {
                auto& m = the_monitor();
                std::lock_guard<util::spin_lock> lock(m.spinlock);
                auto& cur = m.event_stack.top();
                auto cur_time = timer::now();
                cur.allocations.emplace_back(cur_time, m.current_usage);
                m.completed_events.emplace_back(std::move(cur));
                m.event_stack.pop();
                // add a point to the new "top" with the same memory
                // as before but just ahead in time
                if (!m.event_stack.empty()) {
                    if (m.event_stack.top().allocations.size()) {
                        auto last_usage = m.event_stack.top().allocations.back().usage;
                        m.event_stack.top().allocations.emplace_back(cur_time, last_usage);
                    }
                }
            }
        }
    };
    std::chrono::milliseconds log_granularity = std::chrono::milliseconds(20ULL);
    int64_t current_usage = 0;
    bool track_usage = false;
    std::vector<mm_event> completed_events;
    std::stack<mm_event> event_stack;
    timer::time_point start_log;
    timer::time_point last_event;
    util::spin_lock spinlock;
private:
    // disable construction of the object
    memory_monitor() {};
    ~memory_monitor()
    {
        if (track_usage) {
            stop();
        }
    }
    memory_monitor(const memory_monitor&) = delete;
    memory_monitor& operator=(const memory_monitor&) = delete;
private:
    static memory_monitor& the_monitor()
    {
        static memory_monitor m;
        return m;
    }
public:
    static void granularity(std::chrono::milliseconds ms)
    {
        auto& m = the_monitor();
        m.log_granularity = ms;
    }
    static int64_t peak()
    {
        auto& m = the_monitor();
        int64_t max = 0;
        for (auto events : m.completed_events) {
            for (auto alloc : events.allocations) {
                if (max < alloc.usage) {
                    max = alloc.usage;
                }
            }
        }
        return max;
    }

    static void start()
    {
        auto& m = the_monitor();
        m.track_usage = true;
        // clear if there is something there
        if (m.completed_events.size()) {
            m.completed_events.clear();
        }
        while (m.event_stack.size()) {
            m.event_stack.pop();
        }
        m.start_log = timer::now();
        m.current_usage = 0;
        m.last_event = m.start_log;
        m.event_stack.emplace("unknown", 0);
    }
    static void stop()
    {
        auto& m = the_monitor();
        while (!m.event_stack.empty()) {
            m.completed_events.emplace_back(std::move(m.event_stack.top()));
            m.event_stack.pop();
        }
        m.track_usage = false;
    }
    static void record(int64_t delta)
    {
        auto& m = the_monitor();
        if (m.track_usage) {
            std::lock_guard<util::spin_lock> lock(m.spinlock);
            auto cur = timer::now();
            if (m.last_event + m.log_granularity < cur) {
                m.event_stack.top().allocations.emplace_back(cur, m.current_usage);
                m.current_usage = m.current_usage + delta;
                m.event_stack.top().allocations.emplace_back(cur, m.current_usage);
                m.last_event = cur;
            } else {
                if (m.event_stack.top().allocations.size()) {
                    m.current_usage = m.current_usage + delta;
                    m.event_stack.top().allocations.back().usage = m.current_usage;
                    m.event_stack.top().allocations.back().timestamp = cur;
                }
            }
        }
    }
    static mm_event_proxy event(const std::string& name)
    {
        auto& m = the_monitor();
        if (m.track_usage) {
            return mm_event_proxy(name, m.current_usage, true);
        }
        return mm_event_proxy(name, m.current_usage, false);
    }
    template<format_type F>
    static void write_memory_log(std::ostream& out)
    {
        write_mem_log<F>(out, the_monitor());
    }
};

#pragma pack(push, 1)
typedef struct mm_block {
    size_t size;
    struct mm_block* next;
    struct mm_block* prev;
} mm_block_t;

typedef struct bfoot {
    size_t size;
} mm_block_foot_t;
#pragma pack(pop)

#ifndef MSVC_COMPILER
class hugepage_allocator
{
private:
    uint8_t* m_base = nullptr;
    mm_block_t* m_first_block = nullptr;
    uint8_t* m_top = nullptr;
    size_t m_total_size = 0;
    std::multimap<size_t, mm_block_t*> m_free_large;
private:
    size_t determine_available_hugepage_memory();
    void coalesce_block(mm_block_t* block);
    void split_block(mm_block_t* bptr, size_t size);
    uint8_t* hsbrk(size_t size);
    mm_block_t* new_block(size_t size);
    void remove_from_free_set(mm_block_t* block);
    void insert_into_free_set(mm_block_t* block);
    mm_block_t* find_free_block(size_t size_in_bytes);
    mm_block_t* last_block();
    void print_heap();
public:
    void init(SDSL_UNUSED size_t size_in_bytes = 0)
    {
#ifdef MAP_HUGETLB
        if (size_in_bytes == 0) {
                size_in_bytes = determine_available_hugepage_memory();
            }

            m_total_size = size_in_bytes;
            m_base = (uint8_t*)mmap(nullptr, m_total_size,
                                    (PROT_READ | PROT_WRITE),
                                    (MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE), 0, 0);
            if (m_base == MAP_FAILED) {
                throw std::system_error(ENOMEM, std::system_category(),
                                        "hugepage_allocator could not allocate hugepages");
            } else {
                // init the allocator
                m_top = m_base;
                m_first_block = (mm_block_t*)m_base;
            }
#else
        throw std::system_error(ENOMEM, std::system_category(),
                                "hugepage_allocator: MAP_HUGETLB / hugepage support not available");
#endif
    }
    void* mm_realloc(void* ptr, size_t size);
    void* mm_alloc(size_t size_in_bytes);
    void mm_free(void* ptr);
    bool in_address_space(void* ptr)
    {
        // check if ptr is in the hugepage address space
        if (ptr == nullptr) {
            return true;
        }
        if (ptr >= m_base && ptr < m_top) {
            return true;
        }
        return false;
    }
    static hugepage_allocator& the_allocator()
    {
        static hugepage_allocator a;
        return a;
    }
};
#endif

class memory_manager
{
private:
    bool hugepages = false;
private:
    static memory_manager& the_manager()
    {
        static memory_manager m;
        return m;
    }
public:
    static uint64_t* alloc_mem(size_t size_in_bytes)
    {
#ifndef MSVC_COMPILER
        auto& m = the_manager();
        if (m.hugepages) {
            return (uint64_t*)hugepage_allocator::the_allocator().mm_alloc(size_in_bytes);
        }
#endif
        return (uint64_t*)calloc(size_in_bytes, 1);
    }
    static void free_mem(uint64_t* ptr)
    {
#ifndef MSVC_COMPILER
        auto& m = the_manager();
        if (m.hugepages and hugepage_allocator::the_allocator().in_address_space(ptr)) {
            hugepage_allocator::the_allocator().mm_free(ptr);
            return;
        }
#endif
        std::free(ptr);
    }
    static uint64_t* realloc_mem(uint64_t* ptr, size_t size)
    {
#ifndef MSVC_COMPILER
        auto& m = the_manager();
        if (m.hugepages and hugepage_allocator::the_allocator().in_address_space(ptr)) {
            return (uint64_t*)hugepage_allocator::the_allocator().mm_realloc(ptr, size);
        }
#endif
        uint64_t* temp = (uint64_t*)realloc(ptr, size);
        if (temp == NULL) {
            throw std::bad_alloc();
        }
        return temp;
    }
public:
    static void use_hugepages(size_t bytes = 0)
    {
#ifndef MSVC_COMPILER
        auto& m = the_manager();
        hugepage_allocator::the_allocator().init(bytes);
        m.hugepages = true;
#else
        throw std::runtime_error("hugepages not support on MSVC_COMPILER");
#endif
    }
    template<class t_vec>
    static void resize(t_vec& v, const typename t_vec::size_type size)
    {
        uint64_t old_size_in_bytes = ((v.m_size + 63) >> 6) << 3;
        uint64_t new_size_in_bytes = ((size + 63) >> 6) << 3;
        bool do_realloc = old_size_in_bytes != new_size_in_bytes;
        v.m_size = size;
        if (do_realloc || v.m_data == nullptr) {
            // Note that we allocate 8 additional bytes if m_size % 64 == 0.
            // We need this padding since rank data structures do a memory
            // access to this padding to answer rank(size()) if size()%64 ==0.
            // Note that this padding is not counted in the serialize method!
            size_t allocated_bytes = (size_t)(((size + 64) >> 6) << 3);
            v.m_data = memory_manager::realloc_mem(v.m_data, allocated_bytes);
            if (allocated_bytes != 0 && v.m_data == nullptr) {
                throw std::bad_alloc();
            }
            // update and fill with 0s
            if (v.bit_size() < v.capacity()) {
                uint8_t len = (uint8_t)(v.capacity() - v.bit_size());
                uint8_t in_word_offset = (uint8_t)(v.bit_size() & 0x3F);
                bits::write_int(v.m_data + (v.bit_size() >> 6), 0, in_word_offset, len);
            }
            if (((v.m_size) % 64) == 0) {  // initialize unreachable bits with 0
                v.m_data[v.m_size / 64] = 0;
            }

            // update stats
            if (do_realloc) {
                memory_monitor::record((int64_t)new_size_in_bytes - (int64_t)old_size_in_bytes);
            }
        }
    }
    template<class t_vec>
    static void clear(t_vec& v)
    {
        int64_t size_in_bytes = ((v.m_size + 63) >> 6) << 3;
        // remove mem
        memory_manager::free_mem(v.m_data);
        v.m_data = nullptr;

        // update stats
        if (size_in_bytes) {
            memory_monitor::record(size_in_bytes*-1);
        }
    }

    static int open_file_for_mmap(std::string& filename, std::ios_base::openmode mode) {
#ifdef MSVC_COMPILER
        int fd = -1;
            if (!(mode&std::ios_base::out)) _sopen_s(&fd,filename.c_str(), _O_BINARY| _O_RDONLY, _SH_DENYNO, _S_IREAD);
            else _sopen_s(&fd, filename.c_str(), _O_BINARY | _O_RDWR, _SH_DENYNO, _S_IREAD | _S_IWRITE);
            return fd;
#else
        if (!(mode&std::ios_base::out)) return open(filename.c_str(), O_RDONLY);
        else return open(filename.c_str(), O_RDWR);
#endif
        return -1;
    }

    static void* mmap_file(int fd,uint64_t file_size, std::ios_base::openmode mode) {
#ifdef MSVC_COMPILER
        HANDLE fh = (HANDLE)_get_osfhandle(fd);
            if (fh == INVALID_HANDLE_VALUE) {
                return nullptr;
            }
            HANDLE fm;
            if (!(mode&std::ios_base::out)) { // read only?
                fm = CreateFileMapping(fh, NULL, PAGE_READONLY, 0, 0, NULL);
            } else fm = CreateFileMapping(fh, NULL, PAGE_READWRITE, 0, 0, NULL);
            if (fm == NULL) {
                return nullptr;
            }
            void* map = nullptr;
            if (!(mode&std::ios_base::out)) { // read only?
                map = MapViewOfFile(fm, FILE_MAP_READ, 0, 0, file_size);
            } else map = MapViewOfFile(fm, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, file_size);
            // we can close the file handle before we unmap the view: (see UnmapViewOfFile Doc)
            // Although an application may close the file handle used to create a file mapping object,
            // the system holds the corresponding file open until the last view of the file is unmapped.
            // Files for which the last view has not yet been unmapped are held open with no sharing restrictions.
            CloseHandle(fm);
            return map;
#else
        void* map = nullptr;
        if (!(mode&std::ios_base::out)) map = mmap(NULL,file_size,PROT_READ,MAP_SHARED,fd, 0);
        else map = mmap(NULL,file_size,PROT_READ | PROT_WRITE,MAP_SHARED,fd, 0);
        if(map == MAP_FAILED) map = nullptr; // unify windows and unix error behaviour
        return map;
#endif
        return nullptr;
    }

    static int mem_unmap(void* addr,const uint64_t size) {
#ifdef MSVC_COMPILER
        if (UnmapViewOfFile(addr)) return 0;
            return -1;
#else
        return munmap(addr, size);
#endif
        return -1;
    }

    static int close_file_for_mmap(int fd) {
#ifdef MSVC_COMPILER
        return _close(fd);
#else
        return close(fd);
#endif
        return -1;
    }

    static int truncate_file_mmap(int fd,const uint64_t new_size) {
#ifdef MSVC_COMPILER
        auto ret = _chsize_s(fd,new_size);
            if(ret != 0) ret = -1;
            return ret;
#else
        return ftruncate(fd,new_size);
#endif
        return -1;
    }

};

} // end namespace

#endif

#include <iosfwd>    // forward declaration of ostream
#include <stdexcept> // for exceptions
#include <iostream>  // for cerr
#include <typeinfo>
#include <cassert>
#include <iterator>
#include <cstdlib>
#include <cstddef>
#include <ctime>    // for rand initialization
#include <cstring>  // for memcpy
#include <ostream>
#include <istream>
#include <string>
#include <initializer_list>
#include <type_traits>
#include <vector>
#include <ios>

//! Namespace for the succinct data structure library.
namespace sdsl
{

typedef uint64_t std_size_type_for_int_vector;

template<uint8_t t_width=0>
class int_vector;

template<class int_vector_type>
class mm_item;

namespace algorithm
{
template<uint8_t t_width>
static void calculate_sa(const unsigned char* c,
                         typename int_vector<t_width>::size_type len,
                         int_vector<t_width>& sa);
}

//! bit_vector is a specialization of the int_vector.
typedef int_vector<1> bit_vector;

template<class t_int_vector>
class int_vector_reference;

template<class t_int_vector>
class int_vector_iterator_base;

template<class t_int_vector>
class int_vector_iterator;

template<class t_int_vector>
class int_vector_const_iterator;

template<uint8_t t_width,std::ios_base::openmode t_mode>
class int_vector_mapper;

template<uint8_t b, uint8_t t_patter_len>  // forward declaration
class rank_support_v;

class rank_support;

class select_support;

template<uint8_t t_bit_pattern, uint8_t t_pattern_len>
class select_support_mcl;

namespace coder
{
class fibonacci;
class elias_delta;
class elias_gamma;
template<uint8_t t_width> class comma;
}

template<uint8_t t_width>
struct int_vec_category_trait {
    typedef iv_tag type;
};

template<>
struct int_vec_category_trait<1> {
    typedef bv_tag type;
};

template<uint8_t t_width>
struct int_vector_trait {
    typedef uint64_t                                    value_type;
    typedef int_vector<t_width>                         int_vector_type;
    typedef int_vector_reference<int_vector_type>       reference;
    typedef const uint64_t                              const_reference;
    typedef uint8_t                                     int_width_type;
    typedef int_vector_iterator<int_vector_type>        iterator;
    typedef int_vector_const_iterator<int_vector_type>  const_iterator;

    static iterator begin(int_vector_type* v, uint64_t*)
    {
        return iterator(v, 0);
    }
    static iterator end(int_vector_type* v, uint64_t*, int_vector_size_type)
    {
        return iterator(v, v->size()*v->width()) ;
    }
    static const_iterator begin(const int_vector_type* v, const uint64_t*)
    {
        return const_iterator(v, 0);
    }
    static const_iterator end(const int_vector_type* v, const uint64_t*, int_vector_size_type)
    {
        return const_iterator(v, v->size()*v->width());
    }

    static void set_width(uint8_t new_width, int_width_type& width)
    {
        if (t_width == 0) {
            if (0 < new_width and new_width <= 64)
                width = new_width;
            else
                width = 64;
        }
    }
};

template<>
struct int_vector_trait<64> {
    typedef uint64_t        value_type;
    typedef int_vector<64>  int_vector_type;
    typedef uint64_t&       reference;
    typedef const uint64_t  const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint64_t*       iterator;
    typedef const uint64_t* const_iterator;

    static iterator begin(int_vector_type*, uint64_t* begin)
    {
        return begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size)
    {
        return begin+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin)
    {
        return begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size)
    {
        return begin+size;
    }

    static void set_width(uint8_t, int_width_type) {}
};

template<>
struct int_vector_trait<32> {
    typedef uint32_t        value_type;
    typedef int_vector<32>  int_vector_type;
    typedef uint32_t&       reference;
    typedef const uint32_t  const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint32_t*       iterator;
    typedef const uint32_t* const_iterator;

    static iterator begin(int_vector_type*, uint64_t* begin)
    {
        return (uint32_t*)begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size)
    {
        return ((uint32_t*)begin)+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin)
    {
        return (uint32_t*)begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size)
    {
        return ((uint32_t*)begin)+size;
    }

    static void set_width(uint8_t, int_width_type) {}
};

template<>
struct int_vector_trait<16> {
    typedef uint16_t        value_type;
    typedef int_vector<16>  int_vector_type;
    typedef uint16_t&       reference;
    typedef const uint16_t  const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint16_t*       iterator;
    typedef const uint16_t* const_iterator;

    static iterator begin(int_vector_type*, uint64_t* begin)
    {
        return (uint16_t*)begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size)
    {
        return ((uint16_t*)begin)+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin)
    {
        return (uint16_t*)begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size)
    {
        return ((uint16_t*)begin)+size;
    }

    static void set_width(uint8_t, int_width_type) {}
};

template<>
struct int_vector_trait<8> {
    typedef uint8_t         value_type;
    typedef int_vector<8>   int_vector_type;
    typedef uint8_t&        reference;
    typedef const uint8_t   const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint8_t*        iterator;
    typedef const uint8_t*  const_iterator;

    static iterator begin(int_vector_type*, uint64_t* begin)
    {
        return (uint8_t*)begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size)
    {
        return ((uint8_t*)begin)+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin)
    {
        return (uint8_t*)begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size)
    {
        return ((uint8_t*)begin)+size;
    }

    static void set_width(uint8_t, int_width_type) {}
};

//! A generic vector class for integers of width \f$w\in [1..64]\f$.
/*! \author Simon Gog
 *
 *    This generic vector class could be used to generate a vector
 *    that contains integers of fixed width \f$w\in [1..64]\f$.
 *
 *  \tparam t_width Width of the integer. If set to `0` it is variable
 *          during runtime, otherwise fixed at compile time.
 *  @ingroup int_vector
 */
template<uint8_t t_width>
class int_vector
{
private:
    static_assert(t_width <= 64 , "int_vector: width of must be at most 64bits.");
public:
    typedef typename int_vector_trait<t_width>::value_type      value_type;
    typedef typename int_vector_trait<t_width>::iterator        iterator;
    typedef typename int_vector_trait<t_width>::const_iterator  const_iterator;
    typedef typename int_vector_trait<t_width>::reference       reference;
    typedef typename int_vector_trait<t_width>::const_reference const_reference;
    typedef int_vector_reference<int_vector>*                   pointer;
    typedef const value_type*                                   const_pointer;
    typedef ptrdiff_t                                           difference_type;
    typedef int_vector_size_type                                size_type;
    typedef typename int_vector_trait<t_width>::int_width_type  int_width_type;
    typedef rank_support_v<1,1>                                 rank_1_type;
    typedef rank_support_v<0,1>                                 rank_0_type;
    typedef select_support_mcl<1,1>                             select_1_type;
    typedef select_support_mcl<0,1>                             select_0_type;
    typedef typename int_vec_category_trait<t_width>::type      index_category;

    friend struct int_vector_trait<t_width>;
    friend class  int_vector_iterator_base<int_vector>;
    friend class  int_vector_iterator<int_vector>;
    friend class  int_vector_const_iterator<int_vector>;
    template<uint8_t,std::ios_base::openmode> friend class int_vector_mapper;
    friend class  coder::elias_delta;
    friend class  coder::elias_gamma;
    friend class  coder::fibonacci;
    template<uint8_t> friend class coder::comma;
    friend class  memory_manager;

    enum { fixed_int_width = t_width }; // make template parameter accessible

private:

    size_type      m_size;  //!< Number of bits needed to store int_vector.
    uint64_t*      m_data;  //!< Pointer to the memory for the bits.
    int_width_type m_width; //!< Width of the integers.

public:

    //! Constructor for int_vector.
    /*! \param size          Number of elements. Default value is 0.
        \param default_value Initialize all value to `default value`.
        \param int_width     The width of each integer.
        \sa resize, width
     */

    int_vector(size_type size, value_type default_value,
               uint8_t int_width = t_width);

    //! Constructor to fix possible comparison with integeres issue.
    explicit int_vector(size_type size = 0) : int_vector(size, static_cast<value_type>(0), t_width) {

    }
    //! Constructor for initializer_list.
    template<class t_T>
    int_vector(std::initializer_list<t_T> il) : int_vector(0,0)
    {
        resize(il.size());
        size_type idx = 0;
        for (auto x : il) {
            (*this)[idx++] = x;
        }
    }

    //! Move constructor.
    int_vector(int_vector&& v);

    //! Copy constructor.
    int_vector(const int_vector& v);

    //! Destructor.
    ~int_vector();

    //! Equivalent to size() == 0.
    bool empty() const
    {
        return 0==m_size;
    }

    //! Swap method for int_vector.
    void swap(int_vector& v);

    //! Resize the int_vector in terms of elements.
    /*! \param size The size to resize the int_vector in terms of elements.
     */
    void resize(const size_type size)
    {
        bit_resize(size * width());
    }

    //! Resize the int_vector in terms of bits.
    /*! \param size The size to resize the int_vector in terms of bits.
     */
    void bit_resize(const size_type size);

    //! The number of elements in the int_vector.
    /*! \sa max_size, bit_size, capacity
     */
    size_type size() const
    {
        return m_size/m_width;
    }

    //! Maximum size of the int_vector.
    /*! \sa size, bit_size, capacity
    */
    static size_type max_size()
    {
        return ((size_type)1)<<(sizeof(size_type)*8-6);
    }

    //! The number of bits in the int_vector.
    /*!  \sa size, max_size, bit_size, capacity
     */
    size_type bit_size() const
    {
        return m_size;
    }

    //! Returns the size of the occupied bits of the int_vector.
    /*! The capacity of a int_vector is greater or equal to the
        bit_size of the vector: capacity() >= bit_size().
        \sa size, bit_size, max_size, capacity
     */
    size_type capacity() const
    {
        return ((m_size+63)>>6)<<6;
    }

    //! Pointer to the raw data of the int_vector
    /*! \returns Const pointer to the raw data of the int_vector
     */
    const uint64_t* data() const
    {
        return m_data;
    }

    //! Pointer to the raw data of the int_vector
    /*! \returns pointer to the raw data of the int_vector
     */
    uint64_t* data()
    {
        return m_data;
    }

    //! Get the integer value of the binary string of length len starting at position idx in the int_vector.
    /*! \param idx Starting index of the binary representation of the integer.
        \param len Length of the binary representation of the integer. Default value is 64.
        \returns The integer value of the binary string of length len starting at position idx.
        \sa setInt, getBit, setBit
    */
    value_type get_int(size_type idx, const uint8_t len=64) const;

    //! Set the bits from position idx to idx+len-1 to the binary representation of integer x.
    /*! The bit at position idx represents the least significant bit(lsb), and the bit at
        position idx+len-1 the most significant bit (msb) of x.
        \param idx Starting index of the binary representation of x.
        \param x   The integer to store in the int_vector.
        \param len The length used to store x in the int_vector. Default value is 64.
        \sa getInt, getBit, setBit
    */
    void set_int(size_type idx, value_type x, const uint8_t len=64);

    //! Returns the width of the integers which are accessed via the [] operator.
    /*! \returns The width of the integers which are accessed via the [] operator.
        \sa width
    */
    uint8_t width() const
    {
        return m_width;
    }

    //! Sets the width of the integers which are accessed via the [] operator, if t_width equals 0.
    /*! \param intWidth New width of the integers accessed via the [] operator.
        \note This method has no effect if t_width is in the range [1..64].
          \sa width
    */
    void width(uint8_t new_width)
    {
        int_vector_trait<t_width>::set_width(new_width, m_width);
    }

    // Write data (without header) to a stream.
    size_type write_data(std::ostream& out) const;

    //! Serializes the int_vector to a stream.
    /*! \return The number of bytes written to out.
     *  \sa load
     */
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                        std::string name = "", bool write_fixed_as_variable=false) const;

    //! Load the int_vector for a stream.
    void load(std::istream& in);

    //! non const version of [] operator
    /*! \param i Index the i-th integer of length width().
     *  \return A reference to the i-th integer of length width().
     */
    inline reference operator[](const size_type& i);

    //! const version of [] operator
    /*! \param i Index the i-th integer of length width().
     *  \return The value of the i-th integer of length width().
     */
    inline const_reference operator[](const size_type& i) const;

    //! Assignment operator.
    /*! \param v The vector v which should be assigned
     */
    int_vector& operator=(const int_vector& v);

    //! Move assignment operator.
    int_vector& operator=(int_vector&& v);

    //! Equality operator for two int_vectors.
    /*! Two int_vectors are equal if
     *    - capacities and sizes are equal and
     *    - width are equal and
     *    - the bits in the range [0..bit_size()-1] are equal.
     */
    bool operator==(const int_vector& v) const;

    //! Inequality operator for two int_vectors.
    /*! Two int_vectors are not equal if
     *    - capacities and sizes are not equal or
     *    - int widths are not equal or
     *    - the bits in the range [0..bit_size()-1] are not equal.
     */
    bool operator!=(const int_vector& v) const;

    //! Less operator for two int_vectors
    /*! int_vector w is less than v if
     *    - w[i]==v[i] for i<j and w[j]<v[j] with j in [0, min(w.size(), v.size()) )
     *    - or w[i]==v[i] for all i < min(w.size(), v.size()) and w.size()<v.size().
     *  \sa operator>
    */
    bool operator<(const int_vector& v) const;

    //! Greater operator for two int_vectors
    /*! int_vector w is greater than v if
     *    - w[i]==v[i] for i<j and w[j]>v[j] with j in [0, min(w.size(), v.size()) )
     *    - or w[i]==v[i] for all i < min(w.size(), v.size()) and w.size()>v.size().
    */
    bool operator>(const int_vector& v) const;

    //! Less or equal operator
    bool operator<=(const int_vector& v) const;

    //! Greater of equal operator
    bool operator>=(const int_vector& v) const;

    //! bitwise-and-update operator
    int_vector& operator&=(const int_vector& v);

    //! bitwise-or-update equal operator
    int_vector& operator|=(const int_vector& v);

    //! bitwise-xor-update operator
    int_vector& operator^=(const int_vector& v);

    //! Iterator that points to the first element of the int_vector.
    /*!  Time complexity guaranty is O(1).
     */
    const iterator begin()
    {
        return int_vector_trait<t_width>::begin(this, m_data);
    }

    //! Iterator that points to the element after the last element of int_vector.
    /*! Time complexity guaranty is O(1).
     */
    const iterator end()
    {
        return int_vector_trait<t_width>::end(this, m_data, (m_size/m_width));
    }

    //! Const iterator that points to the first element of the int_vector.
    const const_iterator begin() const
    {
        return int_vector_trait<t_width>::begin(this, m_data);
    }

    //! Const iterator that points to the element after the last element of int_vector.
    const const_iterator end() const
    {
        return int_vector_trait<t_width>::end(this, m_data, (m_size/m_width));
    }

    //! Flip all bits of bit_vector
    void flip()
    {
        static_assert(1 == t_width, "int_vector: flip() is available only for bit_vector.");
        if (!empty()) {
            for (uint64_t i=0; i<(capacity()>>6); ++i) {
                m_data[i] = ~m_data[i];
            }
        }
    }

    //! Read the size and int_width of a int_vector
    static void read_header(int_vector_size_type& size, int_width_type& int_width, std::istream& in)
    {
        read_member(size, in);
        if (0 == t_width) {
            read_member(int_width, in);
        }
    }

    //! Write the size and int_width of a int_vector
    static uint64_t write_header(uint64_t size, uint8_t int_width, std::ostream& out)
    {
        uint64_t written_bytes = write_member(size, out);
        if (0 == t_width) {
            written_bytes += write_member(int_width, out);
        }
        return written_bytes;
    }

    struct raw_wrapper {
        const int_vector& vec;
        raw_wrapper() = delete;
        raw_wrapper(const int_vector& _vec) : vec(_vec) {}

        size_type
        serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            auto written_bytes = vec.write_data(out);
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
    };

    const raw_wrapper raw = raw_wrapper(*this);
};

//! A proxy class that acts as a reference to an integer of length \p len bits in a int_vector.
/*! \tparam t_int_vector The specific int_vector class.
 */
template<class t_int_vector>
class int_vector_reference
{
public:
    typedef typename t_int_vector::value_type value_type;
private:
    typename t_int_vector::value_type* const m_word;
    const uint8_t m_offset;
    const uint8_t m_len; //!< Length of the integer referred to in bits.
public:
    //! Constructor for the reference class
    /*! \param word Pointer to the corresponding 64bit word in the int_vector.
        \param offset Offset to the starting bit (offset in [0..63])
        \param len length of the integer, should be v->width()!!!
    */
    int_vector_reference(value_type* word, uint8_t offset, uint8_t len):
        m_word(word),m_offset(offset),m_len(len) {};

    int_vector_reference() = delete;
    int_vector_reference(const int_vector_reference &) noexcept = default;
    int_vector_reference(int_vector_reference &&) noexcept = default;

    //! Assignment operator for the proxy class
    /*!
        The integer x is assign to the referenced
        position in the t_int_vector with the specified width
        of the int_vector
        \param x 64bit integer to assign
        \return A const_reference to the assigned reference
     */
    int_vector_reference& operator=(value_type x)
    {
        bits::write_int(m_word, x, m_offset, m_len);
        return *this;
    };

    int_vector_reference& operator=(const int_vector_reference& x)
    {
        return *this = value_type(x);
    };

    //! Cast the reference to a int_vector<>::value_type
    operator value_type()const
    {
        return bits::read_int(m_word, m_offset, m_len);
    }

    //! Prefix increment of the proxy object
    int_vector_reference& operator++()
    {
        value_type x = bits::read_int(m_word, m_offset, m_len);
        bits::write_int(m_word, x+1, m_offset, m_len);
        return *this;
    }

    //! Postfix increment of the proxy object
    value_type operator++(int)
    {
        value_type val = (typename t_int_vector::value_type)*this;
        ++(*this);
        return val;
    }

    //! Prefix decrement of the proxy object
    int_vector_reference& operator--()
    {
        value_type x = bits::read_int(m_word, m_offset, m_len);
        bits::write_int(m_word, x-1, m_offset, m_len);
        return *this;
    }

    //! Postfix decrement of the proxy object
    value_type operator--(int)
    {
        value_type val = (value_type)*this;
        --(*this);
        return val;
    }

    //! Add assign from the proxy object
    int_vector_reference& operator+=(const value_type x)
    {
        value_type w = bits::read_int(m_word, m_offset, m_len);
        bits::write_int(m_word, w+x, m_offset, m_len);
        return *this;
    }

    //! Subtract assign from the proxy object
    int_vector_reference& operator-=(const value_type x)
    {
        value_type w = bits::read_int(m_word, m_offset, m_len);
        bits::write_int(m_word, w-x, m_offset, m_len);
        return *this;
    }

    bool operator==(const int_vector_reference& x)const
    {
        return value_type(*this) == value_type(x);
    }

    bool operator<(const int_vector_reference& x)const
    {
        return value_type(*this) < value_type(x);
    }
};

// For C++11
template<class t_int_vector>
inline void swap(int_vector_reference<t_int_vector> x,
                 int_vector_reference<t_int_vector> y)
{
    // TODO: more efficient solution?
    typename int_vector_reference<t_int_vector>::value_type tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<class t_int_vector>
inline void swap(typename int_vector_reference<t_int_vector>::value_type& x,
                 int_vector_reference<t_int_vector> y)
{
    // TODO: more efficient solution?
    typename int_vector_reference<t_int_vector>::value_type tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<class t_int_vector>
inline void swap(int_vector_reference<t_int_vector> x,
                 typename int_vector_reference<t_int_vector>::value_type& y)
{
    // TODO: more efficient solution?
    typename int_vector_reference<t_int_vector>::value_type tmp = x;
    x = y;
    y = tmp;
}

// specialization for int_vector_reference for int_vector == bit_vector
// special thanks to Timo Beller, who pointed out that the specialization is missing
// Same implementation as in stl_bvector.h.
template<>
class int_vector_reference<bit_vector>
{
public:
    typedef bool value_type;
private:
    uint64_t* const m_word;
    uint64_t m_mask;
public:
    //! Constructor for the reference class
    /*! \param word Pointer to the corresponding 64bit word in the int_vector.
        \param offset Offset to the starting bit (offset in [0..63])
    */
    int_vector_reference(uint64_t* word, uint8_t offset, uint8_t):
        m_word(word),m_mask(1ULL<<offset) {};

    int_vector_reference() = delete;
    int_vector_reference(const int_vector_reference &) noexcept = default;
    int_vector_reference(int_vector_reference &&) noexcept = default;

    //! Assignment operator for the proxy class
    int_vector_reference& operator=(bool x)
    {
        if (x)
            *m_word |= m_mask;
        else
            *m_word &= ~m_mask;
        return *this;
    };

    int_vector_reference& operator=(const int_vector_reference& x)
    {
        return *this = bool(x);
    };

    //! Cast the reference to a bool
    operator bool()const
    {
        return !!(*m_word & m_mask);
    }

    bool operator==(const int_vector_reference& x)const
    {
        return bool(*this) == bool(x);
    }

    bool operator<(const int_vector_reference& x)const
    {
        return !bool(*this) && bool(x);
    }
};

// For C++11
template<>
inline void swap(int_vector_reference<bit_vector> x,
                 int_vector_reference<bit_vector> y)
{
    // TODO: more efficient solution?
    bool tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<>
inline void swap(bool& x,
                 int_vector_reference<bit_vector> y)
{
    // TODO: more efficient solution?
    bool tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<>
inline void swap(int_vector_reference<bit_vector> x,
                 bool& y)
{
    // TODO: more efficient solution?
    bool tmp = x;
    x = y;
    y = tmp;
}

template<class t_int_vector>
class int_vector_iterator_base: public std::iterator<std::random_access_iterator_tag, typename t_int_vector::value_type, typename t_int_vector::difference_type>
{
public:
    typedef uint64_t  size_type;
protected:
    uint8_t           m_offset;
    uint8_t           m_len;

public:
    int_vector_iterator_base(uint8_t offset, uint8_t len):
        m_offset(offset),m_len(len) {}

    int_vector_iterator_base(const t_int_vector* v=nullptr, size_type idx=0):
        m_offset(idx&0x3F), m_len(v==nullptr ? 0 : v->m_width) {}
};

template<class t_int_vector>
class int_vector_iterator : public int_vector_iterator_base<t_int_vector>
{
public:

    typedef int_vector_reference<t_int_vector>     reference;
    typedef uint64_t                               value_type;
    typedef int_vector_iterator                    iterator;
    typedef reference*                             pointer;
    typedef typename t_int_vector::size_type       size_type;
    typedef typename t_int_vector::difference_type difference_type;

    friend class int_vector_const_iterator<t_int_vector>;
private:

    using int_vector_iterator_base<t_int_vector>::m_offset; // make m_offset easy usable
    using int_vector_iterator_base<t_int_vector>::m_len;    // make m_len easy usable

    typename t_int_vector::value_type* m_word;

public:

    int_vector_iterator(t_int_vector* v=nullptr, size_type idx=0):
        int_vector_iterator_base<t_int_vector>(v, idx),
        m_word((v != nullptr) ? v->m_data + (idx>>6) : nullptr) {}

    int_vector_iterator(const int_vector_iterator<t_int_vector>& it) :
        int_vector_iterator_base<t_int_vector>(it), m_word(it.m_word)
    {
        m_offset = it.m_offset;
        m_len = it.m_len;
    }

    reference operator*() const
    {
        return reference(m_word, m_offset, m_len);
    }

    //! Prefix increment of the Iterator
    iterator& operator++()
    {
        m_offset+=m_len;
        if (m_offset >= 64) {
            m_offset &= 0x3F;
            ++m_word;
        }
        return *this;
    }

    //! Postfix increment of the Iterator
    iterator operator++(int)
    {
        int_vector_iterator it = *this;
        ++(*this);
        return it;
    }

    //! Prefix decrement of the Iterator
    iterator& operator--()
    {
        m_offset-=m_len;
        if (m_offset >= 64) {
            m_offset &= 0x3F;
            --m_word;
        }
        return *this;
    }

    //! Postfix decrement of the Iterator
    iterator operator--(int)
    {
        int_vector_iterator it = *this;
        --(*this);
        return it;
    }

    iterator& operator+=(difference_type i)
    {
        if (i<0)
            return *this -= (-i);
        difference_type t = i*m_len;
        m_word += (t>>6);
        if ((m_offset+=(t&0x3F))&~0x3F) {  // if new offset is >= 64
            ++m_word;       // add one to the word
            m_offset&=0x3F; // offset = offset mod 64
        }
        return *this;
    }

    iterator& operator-=(difference_type i)
    {
        if (i<0)
            return *this += (-i);
        difference_type t = i*m_len;
        m_word -= (t>>6);
        if ((m_offset-=(t&0x3F))&~0x3F) {  // if new offset is < 0
            --m_word;
            m_offset&=0x3F;
        }
        return *this;
    }

    iterator& operator=(const int_vector_iterator<t_int_vector>& it)
    {
        if (this != &it) {
            m_word   = it.m_word;
            m_offset = it.m_offset;
            m_len    = it.m_len;
        }
        return *this;
    }

    iterator operator+(difference_type i) const
    {
        iterator it = *this;
        return it += i;
    }

    iterator operator-(difference_type i) const
    {
        iterator it = *this;
        return it -= i;
    }

    reference operator[](difference_type i) const
    {
        return *(*this + i);
    }

    bool operator==(const int_vector_iterator& it)const
    {
        return it.m_word == m_word && it.m_offset == m_offset;
    }

    bool operator!=(const int_vector_iterator& it)const
    {
        return !(*this==it);
    }

    bool operator<(const int_vector_iterator& it)const
    {
        if (m_word == it.m_word)
            return m_offset < it.m_offset;
        return m_word < it.m_word;
    }

    bool operator>(const int_vector_iterator& it)const
    {
        if (m_word == it.m_word)
            return m_offset > it.m_offset;
        return m_word > it.m_word;
    }

    bool operator>=(const int_vector_iterator& it)const
    {
        return !(*this < it);
    }

    bool operator<=(const int_vector_iterator& it)const
    {
        return !(*this > it);
    }
    inline difference_type operator-(const int_vector_iterator& it)
    {
        return (((m_word - it.m_word)<<6) + m_offset - it.m_offset) / m_len;
    }
};

//template<class t_int_vector>
//void swap(const int_vector_iterator<t_int_vector> &x, const int_vector_iterator<t_int_vector> &y){
//  x->swap(*y);
//}

template<class t_int_vector>
inline int_vector_iterator<t_int_vector> operator+(typename int_vector_iterator<t_int_vector>::difference_type n, const int_vector_iterator<t_int_vector>& it)
{
    return it+n;
}

template<class t_int_vector>
class int_vector_const_iterator : public int_vector_iterator_base<t_int_vector>
{
public:

    typedef typename t_int_vector::value_type        const_reference;
    typedef const typename t_int_vector::value_type* pointer;
    typedef int_vector_const_iterator                const_iterator;
    typedef typename t_int_vector::size_type         size_type;
    typedef typename t_int_vector::difference_type   difference_type;

    template<class X>
    friend typename int_vector_const_iterator<X>::difference_type
    operator-(const int_vector_const_iterator<X>& x, const int_vector_const_iterator<X>& y);
    friend class int_vector_iterator<t_int_vector>;
    friend class int_vector_iterator_base<t_int_vector>;

private:

    using int_vector_iterator_base<t_int_vector>::m_offset; // make m_offset easy usable
    using int_vector_iterator_base<t_int_vector>::m_len;    // make m_len easy usable

    const typename t_int_vector::value_type* m_word;

public:

    int_vector_const_iterator(const t_int_vector* v=nullptr, size_type idx=0):
        int_vector_iterator_base<t_int_vector>(v, idx),
        m_word((v != nullptr) ? v->m_data + (idx>>6) : nullptr) {}

    int_vector_const_iterator(const int_vector_const_iterator& it):
        int_vector_iterator_base<t_int_vector>(it), m_word(it.m_word)
    {
        m_offset = it.m_offset;
        m_len = it.m_len;
    }

    int_vector_const_iterator(const int_vector_iterator<t_int_vector>& it):
        m_word(it.m_word)
    {
        m_offset = it.m_offset;
        m_len = it.m_len;
    }

    const_reference operator*() const
    {
        if (m_offset+m_len <= 64) {
            return ((*m_word)>>m_offset)&bits::lo_set[m_len];
        }
        return ((*m_word)>>m_offset) |
            ((*(m_word+1) & bits::lo_set[(m_offset+m_len)&0x3F])<<(64-m_offset));
    }

    //! Prefix increment of the Iterator
    const_iterator& operator++()
    {
        m_offset+=m_len;
        if (m_offset >= 64) {
            m_offset &= 0x3F;
            ++m_word;
        }
        return *this;
    }

    //! Postfix increment of the Iterator
    const_iterator operator++(int)
    {
        int_vector_const_iterator it = *this;
        ++(*this);
        return it;
    }

    //! Prefix decrement of the Iterator
    const_iterator& operator--()
    {
        m_offset-=m_len;
        if (m_offset >= 64) {
            m_offset &= 0x3F;
            --m_word;
        }
        return *this;
    }

    //! Postfix decrement of the Iterator
    const_iterator operator--(int)
    {
        int_vector_const_iterator it = *this;
        --(*this);
        return it;
    }

    const_iterator& operator+=(difference_type i)
    {
        if (i<0)
            return *this -= (-i);
        difference_type t = i*m_len;
        m_word += (t>>6);
        if ((m_offset+=(t&0x3F))&~0x3F) {// if new offset >= 64
            ++m_word;       // add one to the word
            m_offset&=0x3F; // offset = offset mod 64
        }
        return *this;
    }

    const_iterator& operator-=(difference_type i)
    {
        if (i<0)
            return *this += (-i);
        difference_type t = i*m_len;
        m_word -= (t>>6);
        if ((m_offset-=(t&0x3F))&~0x3F) {// if new offset is < 0
            --m_word;
            m_offset&=0x3F;
        }
        return *this;
    }

    const_iterator operator+(difference_type i) const
    {
        const_iterator it = *this;
        return it += i;
    }

    const_iterator operator-(difference_type i) const
    {
        const_iterator it = *this;
        return it -= i;
    }

    const_reference operator[](difference_type i) const
    {
        return *(*this + i);
    }

    bool operator==(const int_vector_const_iterator& it)const
    {
        return it.m_word == m_word && it.m_offset == m_offset;
    }

    bool operator!=(const int_vector_const_iterator& it)const
    {
        return !(*this==it);
    }

    bool operator<(const int_vector_const_iterator& it)const
    {
        if (m_word == it.m_word)
            return m_offset < it.m_offset;
        return m_word < it.m_word;
    }

    bool operator>(const int_vector_const_iterator& it)const
    {
        if (m_word == it.m_word)
            return m_offset > it.m_offset;
        return m_word > it.m_word;
    }

    bool operator>=(const int_vector_const_iterator& it)const
    {
        return !(*this < it);
    }

    bool operator<=(const int_vector_const_iterator& it)const
    {
        return !(*this > it);
    }

};

template<class t_int_vector>
inline typename int_vector_const_iterator<t_int_vector>::difference_type
operator-(const int_vector_const_iterator<t_int_vector>& x,
          const int_vector_const_iterator<t_int_vector>& y)
{
    return (((x.m_word - y.m_word)<<6) + x.m_offset - y.m_offset) / x.m_len;
}

template<class t_int_vector>
inline int_vector_const_iterator<t_int_vector>
operator+(typename int_vector_const_iterator<t_int_vector>::difference_type n,
          const int_vector_const_iterator<t_int_vector>& it)
{
    return it + n;
}

template<class t_bv>
inline typename std::enable_if<std::is_same<typename t_bv::index_category ,bv_tag>::value, std::ostream&>::type
operator<<(std::ostream& os, const t_bv& bv)
{
    for (auto b : bv) {
        os << b;
    }
    return os;
}

// ==== int_vector implementation  ====

template<uint8_t t_width>
inline int_vector<t_width>::int_vector(size_type size, value_type default_value, uint8_t intWidth):
    m_size(0), m_data(nullptr), m_width(t_width)
{
    width(intWidth);
    resize(size);
    util::set_to_value(*this, default_value); // new initialization
}

template<uint8_t t_width>
inline int_vector<t_width>::int_vector(int_vector&& v) :
    m_size(v.m_size), m_data(v.m_data), m_width(v.m_width)
{
    v.m_data = nullptr; // ownership of v.m_data now transferred
    v.m_size = 0;
}

template<uint8_t t_width>
inline int_vector<t_width>::int_vector(const int_vector& v):
    m_size(0), m_data(nullptr), m_width(v.m_width)
{
    bit_resize(v.bit_size());
    if (v.capacity() > 0) {
        if (memcpy(m_data, v.data() ,v.capacity()/8)==nullptr) {
            throw std::bad_alloc(); // LCOV_EXCL_LINE
        }
    }
    width(v.m_width);
}

template<uint8_t t_width>
int_vector<t_width>& int_vector<t_width>::operator=(const int_vector& v)
{
    if (this != &v) {// if v is not the same object
        bit_resize(v.bit_size());
        if (v.bit_size()>0) {
            if (memcpy(m_data, v.data() ,v.capacity()/8)==nullptr) {
                throw std::bad_alloc(); // LCOV_EXCL_LINE
            }
        }
        width(v.width());
    }
    return *this;
}

template<uint8_t t_width>
int_vector<t_width>& int_vector<t_width>::operator=(int_vector&& v)
{
    swap(v);
    return *this;
}

// Destructor
template<uint8_t t_width>
int_vector<t_width>::~int_vector()
{
    memory_manager::clear(*this);
}

template<uint8_t t_width>
void int_vector<t_width>::swap(int_vector& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        size_type size     = m_size;
        uint64_t* data     = m_data;
        uint8_t  int_width = m_width;
        m_size   = v.m_size;
        m_data   = v.m_data;
        width(v.m_width);
        v.m_size = size;
        v.m_data = data;
        v.width(int_width);
    }
}

template<uint8_t t_width>
void int_vector<t_width>::bit_resize(const size_type size)
{
    memory_manager::resize(*this, size);
}

template<uint8_t t_width>
auto int_vector<t_width>::get_int(size_type idx, const uint8_t len)const -> value_type
{
#ifdef SDSL_DEBUG
    if (idx+len > m_size) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::get_int(size_type, uint8_t); idx+len > size()!");
    }
    if (len > 64) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::get_int(size_type, uint8_t); len>64!");
    }
#endif
    return bits::read_int(m_data+(idx>>6), idx&0x3F, len);
}

template<uint8_t t_width>
inline void int_vector<t_width>::set_int(size_type idx, value_type x, const uint8_t len)
{
#ifdef SDSL_DEBUG
    if (idx+len > m_size) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::set_int(size_type, uint8_t); idx+len > size()!");
    }
    if (len > 64) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::set_int(size_type, uint8_t); len>64!");
    }
#endif
    bits::write_int(m_data+(idx>>6), x, idx&0x3F, len);
}

template<uint8_t t_width>
inline auto int_vector<t_width>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    size_type j = idx * m_width;
    return reference(this->m_data + (j>>6), j&0x3F, m_width);
}

// specialized [] operator for 64 bit access.
template<>
inline auto int_vector<64>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(this->m_data+idx);
}

// specialized [] operator for 32 bit access.
template<>
inline auto int_vector<32>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(((uint32_t*)(this->m_data))+idx);
}

// specialized [] operator for 16 bit access.
template<>
inline auto int_vector<16>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(((uint16_t*)(this->m_data))+idx);
}

// specialized [] operator for 8 bit access.
template<>
inline auto int_vector<8>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(((uint8_t*)(this->m_data))+idx);
}

template<uint8_t t_width>
inline auto
int_vector<t_width>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return get_int(idx * t_width, t_width);
}

template<>
inline auto
int_vector<0>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return get_int(idx * m_width, m_width);
}

template<>
inline auto
int_vector<64>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(this->m_data+idx);
}

template<>
inline auto
int_vector<32>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(((uint32_t*)this->m_data)+idx);
}

template<>
inline auto
int_vector<16>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(((uint16_t*)this->m_data)+idx);
}

template<>
inline auto
int_vector<8>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(((uint8_t*)this->m_data)+idx);
}

template<>
inline auto
int_vector<1>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return ((*(m_data+(idx>>6)))>>(idx&0x3F))&1;
}
template<uint8_t t_width>
bool int_vector<t_width>::operator==(const int_vector& v)const
{
    if (capacity() != v.capacity())
        return false;
    if (bit_size() != v.bit_size())
        return false;
    if (empty())
        return true;
    const uint64_t* data1 = v.data();
    const uint64_t* data2 = data();
    for (size_type i=0; i < (capacity()>>6)-1; ++i) {
        if (*(data1++) != *(data2++))
            return false;
    }
    int8_t l = 64-(capacity()-bit_size());
    return ((*data1)&bits::lo_set[l])==((*data2)&bits::lo_set[l]);
}

template<uint8_t t_width>
bool int_vector<t_width>::operator<(const int_vector& v)const
{
    size_type min_size = size();
    if (min_size > v.size())
        min_size = v.size();
    for (auto it = begin(), end = begin()+min_size, it_v = v.begin(); it!=end; ++it, ++it_v) {
        if (*it == *it_v)
            continue;
        else
            return *it < *it_v;
    }
    return  size() < v.size();
}

template<uint8_t t_width>
bool int_vector<t_width>::operator>(const int_vector& v)const
{
    size_type min_size = size();
    if (min_size > v.size())
        min_size = v.size();
    for (auto it = begin(), end = begin()+min_size, it_v = v.begin(); it!=end; ++it, ++it_v) {
        if (*it == *it_v)
            continue;
        else
            return *it > *it_v;
    }
    return  size() > v.size();
}

template<uint8_t t_width>
bool int_vector<t_width>::operator<=(const int_vector& v)const
{
    return *this==v or *this<v;
}

template<uint8_t t_width>
bool int_vector<t_width>::operator>=(const int_vector& v)const
{
    return *this==v or *this>v;
}

template<uint8_t t_width>
bool int_vector<t_width>::operator!=(const int_vector& v)const
{
    return !(*this==v);
}

template<uint8_t t_width>
int_vector<t_width>& int_vector<t_width>::operator&=(const int_vector& v)
{
    assert(bit_size() == v.bit_size());
    assert(v.capacity() <= capacity());
    for (uint64_t i=0; i<(v.capacity()>>6); ++i)
        m_data[i] &= v.m_data[i];
    return *this;
}

template<uint8_t t_width>
int_vector<t_width>& int_vector<t_width>::operator|=(const int_vector& v)
{
    assert(bit_size() == v.bit_size());
    assert(v.capacity() <= capacity());
    for (uint64_t i=0; i<(v.capacity()>>6); ++i)
        m_data[i] |= v.m_data[i];
    return *this;
}

template<uint8_t t_width>
int_vector<t_width>& int_vector<t_width>::operator^=(const int_vector& v)
{
    assert(bit_size() == v.bit_size());
    assert(v.capacity() <= capacity());
    for (uint64_t i=0; i<(v.capacity()>>6); ++i)
        m_data[i] ^= v.m_data[i];
    return *this;
}

template<uint8_t t_width>
typename int_vector<t_width>::size_type int_vector<t_width>::write_data(std::ostream& out) const
{
    size_type written_bytes = 0;
    uint64_t* p = m_data;
    size_type idx = 0;
    while (idx+conf::SDSL_BLOCK_SIZE < (capacity()>>6)) {
        out.write((char*) p, conf::SDSL_BLOCK_SIZE*sizeof(uint64_t));
        written_bytes += conf::SDSL_BLOCK_SIZE*sizeof(uint64_t);
        p     += conf::SDSL_BLOCK_SIZE;
        idx    += conf::SDSL_BLOCK_SIZE;
    }
    out.write((char*) p, ((capacity()>>6)-idx)*sizeof(uint64_t));
    written_bytes += ((capacity()>>6)-idx)*sizeof(uint64_t);
    return written_bytes;
}

template<uint8_t t_width>
typename int_vector<t_width>::size_type int_vector<t_width>::serialize(std::ostream& out,
                                                                       structure_tree_node* v,
                                                                       std::string name,
                                                                       bool write_fixed_as_variable) const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    if (t_width > 0 and write_fixed_as_variable) {
        written_bytes += int_vector<0>::write_header(m_size, t_width, out);
    } else {
        written_bytes += int_vector<t_width>::write_header(m_size, m_width, out);
    }
    written_bytes += write_data(out);
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t t_width>
void int_vector<t_width>::load(std::istream& in)
{
    size_type size;
    int_vector<t_width>::read_header(size, m_width, in);

    bit_resize(size);
    uint64_t* p = m_data;
    size_type idx = 0;
    while (idx+conf::SDSL_BLOCK_SIZE < (capacity()>>6)) {
        in.read((char*) p, conf::SDSL_BLOCK_SIZE*sizeof(uint64_t));
        p     += conf::SDSL_BLOCK_SIZE;
        idx += conf::SDSL_BLOCK_SIZE;
    }
    in.read((char*) p, ((capacity()>>6)-idx)*sizeof(uint64_t));
}

}// end namespace sdsl

/* sdsl - succinct data structures library
    Copyright (C) 2008-2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file int_vector_buffer.hpp
    \brief int_vector_buffer.hpp contains the sdsl::int_vector_buffer class.
    \author Maike Zwerger, Timo Beller and Simon Gog
*/
#ifndef INCLUDED_INT_VECTOR_BUFFER
#define INCLUDED_INT_VECTOR_BUFFER

/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file iterators.hpp
    \brief iterators.hpp contains an generic iterator for random access containers.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_ITERATORS
#define INCLUDED_SDSL_ITERATORS

#include <iterator>

namespace sdsl
{

//! Generic iterator for a random access container
/*! \tparam t_rac Type of random access container.
 */
template<class t_rac>
class random_access_const_iterator: public std::iterator<std::random_access_iterator_tag, typename t_rac::value_type, typename t_rac::difference_type>
{
public:
    typedef const typename t_rac::value_type  const_reference;
    typedef typename t_rac::size_type size_type;
    typedef random_access_const_iterator<t_rac> iterator;
    typedef typename t_rac::difference_type difference_type;

private:
    const t_rac* m_rac;// pointer to the random access container
    typename t_rac::size_type m_idx;

    template<class t_RAC>
    friend typename random_access_const_iterator<t_RAC>::difference_type operator-(const random_access_const_iterator<t_RAC>& x,
                                                                                   const random_access_const_iterator<t_RAC>& y);

public:

    //! Default Constructor
    random_access_const_iterator() : m_rac(nullptr), m_idx(0) { }

    //! Constructor
    random_access_const_iterator(const t_rac* rac, size_type idx = 0) : m_rac(rac), m_idx(idx) { }

    //! Dereference operator for the Iterator.
    const_reference operator*()const
    {
        return (*m_rac)[m_idx];
    }

    //! Prefix increment of the Iterator.
    iterator& operator++()
    {
        ++m_idx;
        return *this;
    }

    //! Postfix increment of the Iterator.
    iterator operator++(int)
    {
        random_access_const_iterator it = *this;
        ++(*this);
        return it;
    }

    //! Prefix decrement of the Iterator.
    iterator& operator--()
    {
        --m_idx;
        return *this;
    }

    //! Postfix decrement of the Iterator.
    iterator operator--(int)
    {
        random_access_const_iterator it = *this;
        --(*this);
        return it;
    }

    iterator& operator+=(difference_type i)
    {
        if (i<0)
            return *this -= (-i);
        m_idx += i;
        return *this;
    }

    iterator& operator-=(difference_type i)
    {
        if (i<0)
            return *this += (-i);
        m_idx -= i;
        return *this;
    }

    iterator operator+(difference_type i) const
    {
        iterator it = *this;
        return it += i;
    }

    iterator operator-(difference_type i) const
    {
        iterator it = *this;
        return it -= i;
    }

    const_reference operator[](difference_type i) const
    {
        return *(*this + i);
    }

    bool operator==(const iterator& it)const
    {
        return it.m_rac == m_rac && it.m_idx == m_idx;
    }

    bool operator!=(const iterator& it)const
    {
        return !(*this==it);
    }

    bool operator<(const iterator& it)const
    {
        return m_idx < it.m_idx;
    }

    bool operator>(const iterator& it)const
    {
        return m_idx > it.m_idx;
    }

    bool operator>=(const iterator& it)const
    {
        return !(*this < it);
    }

    bool operator<=(const iterator& it)const
    {
        return !(*this > it);
    }

};

template<class t_rac>
inline typename random_access_const_iterator<t_rac>::difference_type operator-(const random_access_const_iterator<t_rac>& x, const random_access_const_iterator<t_rac>& y)
{
    return (typename random_access_const_iterator<t_rac>::difference_type)x.m_idx
        - (typename random_access_const_iterator<t_rac>::difference_type)y.m_idx;
}

template<class t_rac>
inline random_access_const_iterator<t_rac> operator+(typename random_access_const_iterator<t_rac>::difference_type n, const random_access_const_iterator<t_rac>& it)
{
    return it+n;
}

template<typename t_F>
struct random_access_container {
    typedef int_vector<>::size_type                               size_type;
    typedef int_vector<>::difference_type                         difference_type;
    typedef typename std::result_of<t_F(size_type)>::type         value_type;
    typedef random_access_const_iterator<random_access_container> iterator_type;

    t_F f;
    size_type m_size;

    random_access_container() {};
    random_access_container(t_F ff, size_type size) : f(ff), m_size(size) { }

    value_type operator[](size_type i) const { return f(i); }

    size_type size() const { return m_size; }

    iterator_type begin() const
    {
        return iterator_type(this, 0);
    }

    iterator_type end() const
    {
        return iterator_type(this, size());
    }
};

} // end namespace sdsl
#endif

#include <cassert>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

namespace sdsl
{

template<uint8_t t_width=0>
class int_vector_buffer
{
public:
    class iterator;
    typedef typename int_vector<t_width>::difference_type difference_type;
    typedef typename int_vector<t_width>::value_type      value_type;

private:
    static_assert(t_width <= 64 , "int_vector_buffer: width must be at most 64 bits.");
    sdsl::isfstream     m_ifile;
    sdsl::osfstream     m_ofile;
    std::string         m_filename;
    int_vector<t_width> m_buffer;
    bool                m_need_to_write = false;
    // length of int_vector header in bytes: 0 for plain, 8 for int_vector<t_width> (0 < t_width), 9 for int_vector<0>
    uint64_t            m_offset     = 0;
    uint64_t            m_buffersize = 8;    // in elements! m_buffersize*width() must be a multiple of 8!
    uint64_t            m_size       = 0;    // size of int_vector_buffer
    uint64_t            m_begin      = 0;    // number in elements

    //! Read block containing element at index idx.
    void read_block(const uint64_t idx)
    {
        m_begin = (idx/m_buffersize)*m_buffersize;
        if (m_begin >= m_size) {
            util::set_to_value(m_buffer, 0);
        } else {
            m_ifile.seekg(m_offset+(m_begin*width())/8);
            assert(m_ifile.good());
            m_ifile.read((char*) m_buffer.data(), (m_buffersize*width())/8);
            if ((uint64_t)m_ifile.gcount() < (m_buffersize*width())/8) {
                m_ifile.clear();
            }
            assert(m_ifile.good());
            for (uint64_t i=m_size-m_begin; i<m_buffersize; ++i) {
                m_buffer[i] = 0;
            }
        }
    }

    //! Write current block to file.
    void write_block()
    {
        if (m_need_to_write) {
            m_ofile.seekp(m_offset+(m_begin*width())/8);
            assert(m_ofile.good());
            if (m_begin+m_buffersize >= m_size) {
                //last block in file
                uint64_t wb = ((m_size-m_begin)*width()+7)/8;
                m_ofile.write((char*) m_buffer.data(), wb);
            } else {
                m_ofile.write((char*) m_buffer.data(), (m_buffersize*width())/8);
            }
            m_ofile.flush();
            assert(m_ofile.good());
            m_need_to_write = false;
        }
    }

    //! Read value from idx.
    uint64_t read(const uint64_t idx)
    {
        assert(is_open());
        assert(idx < m_size);
        if (idx < m_begin or m_begin+m_buffersize <= idx) {
            write_block();
            read_block(idx);
        }
        return m_buffer[idx-m_begin];
    }

    //! Write value to idx.
    void write(const uint64_t idx, const uint64_t value)
    {
        assert(is_open());
        // If idx is not in current block, write current block and load needed block
        if (idx < m_begin or m_begin+m_buffersize <= idx) {
            write_block();
            read_block(idx);
        }
        if (m_size <= idx) {
            m_size = idx+1;
        }
        m_need_to_write = true;
        m_buffer[idx-m_begin] = value;
    }

public:

    //! Constructor.
    int_vector_buffer()
    {
        m_buffer = int_vector<t_width>();
    }

    //! Constructor for int_vector_buffer.
    /*! \param filename   File that contains the data read from / written to.
     *  \param mode       Openmode:
     *                    std::ios::in opens an existing file (that must exist already),
     *                    std::ios::out creates a new file (that may exist already).
     *  \param buffersize Buffersize in bytes. This has to be a multiple of 8, if not the next multiple of 8 will be taken
     *  \param int_width  The width of each integer.
     *  \param is_plain   If false (default) the file will be interpreted as int_vector.
     *                    If true the file will be interpreted as plain array with t_width bits per integer.
     *                    In second case (is_plain==true), t_width must be 8, 16, 32 or 64.
     */
    int_vector_buffer(const std::string filename, std::ios::openmode mode=std::ios::in, const uint64_t buffer_size=1024*1024, const uint8_t int_width=t_width, const bool is_plain=false)
    {
        m_filename = filename;
        assert(!(mode&std::ios::app));
        mode &= ~std::ios::app;
        m_buffer.width(int_width);
        if (is_plain) {
            // is_plain is only allowed with width() in {8, 16, 32, 64}
            assert(8==width() or 16==width() or 32==width() or 64==width());
        } else {
            m_offset = t_width ? 8 : 9;
        }

        // Open file for IO
        m_ofile.open(m_filename, mode|std::ios::out|std::ios::binary);
        assert(m_ofile.good());
        m_ifile.open(m_filename, std::ios::in|std::ios::binary);
        assert(m_ifile.good());
        if (mode & std::ios::in) {
            uint64_t size  = 0;
            if (is_plain) {
                m_ifile.seekg(0, std::ios_base::end);
                size = m_ifile.tellg()*8;
            } else {
                uint8_t width = 0;
                int_vector<t_width>::read_header(size, width, m_ifile);
                m_buffer.width(width);
            }
            assert(m_ifile.good());
            m_size = size/width();
        }
        buffersize(buffer_size);
    }

    //! Move constructor.
    int_vector_buffer(int_vector_buffer&& ivb) :
        m_filename(std::move(ivb.m_filename)),
        m_buffer(std::move(ivb.m_buffer)),
        m_need_to_write(ivb.m_need_to_write),
        m_offset(ivb.m_offset),
        m_buffersize(ivb.m_buffersize),
        m_size(ivb.m_size),
        m_begin(ivb.m_begin)
    {
        ivb.m_ifile.close();
        ivb.m_ofile.close();
        m_ifile.open(m_filename, std::ios::in|std::ios::binary);
        m_ofile.open(m_filename, std::ios::in|std::ios::out|std::ios::binary);
        assert(m_ifile.good());
        assert(m_ofile.good());
        // set ivb to default-constructor state
        ivb.m_filename = "";
        ivb.m_buffer = int_vector<t_width>();
        ivb.m_need_to_write = false;
        ivb.m_offset = 0;
        ivb.m_buffersize = 8;
        ivb.m_size = 0;
        ivb.m_begin = 0;
    }

    //! Destructor.
    ~int_vector_buffer()
    {
        close();
    }

    //! Move assignment operator.
    int_vector_buffer<t_width>& operator=(int_vector_buffer&& ivb)
    {
        close();
        ivb.m_ifile.close();
        ivb.m_ofile.close();
        m_filename = ivb.m_filename;
        m_ifile.open(m_filename, std::ios::in|std::ios::binary);
        m_ofile.open(m_filename, std::ios::in|std::ios::out|std::ios::binary);
        assert(m_ifile.good());
        assert(m_ofile.good());
        // assign the values of ivb to this
        m_buffer = (int_vector<t_width>&&)ivb.m_buffer;
        m_need_to_write = ivb.m_need_to_write;
        m_offset = ivb.m_offset;
        m_buffersize = ivb.m_buffersize;
        m_size = ivb.m_size;
        m_begin = ivb.m_begin;
        // set ivb to default-constructor state
        ivb.m_filename = "";
        ivb.m_buffer = int_vector<t_width>();
        ivb.m_need_to_write = false;
        ivb.m_offset = 0;
        ivb.m_buffersize = 8;
        ivb.m_size = 0;
        ivb.m_begin = 0;
        return *this;
    }

    //! Returns the width of the integers which are accessed via the [] operator.
    uint8_t width() const
    {
        return m_buffer.width();
    }

    //! Returns the number of elements currently stored.
    uint64_t size() const
    {
        return m_size;
    }

    //! Returns the filename.
    std::string filename() const
    {
        return m_filename;
    }

    //! Returns the buffersize in bytes
    uint64_t buffersize() const
    {
        assert(m_buffersize*width()%8==0);
        return (m_buffersize*width())/8;
    }

    //! Set the buffersize in bytes
    void buffersize(uint64_t buffersize)
    {
        if (0ULL == buffersize)
            buffersize = 8;
        write_block();
        if (0==(buffersize*8)%width()) {
            m_buffersize = buffersize*8/width(); // m_buffersize might not be multiple of 8, but m_buffersize*width() is.
        } else {
            uint64_t element_buffersize = (buffersize*8)/width()+1; // one more element than fits into given buffersize in byte
            m_buffersize = element_buffersize+7 - (element_buffersize+7)%8; // take next multiple of 8
        }
        m_buffer = int_vector<t_width>(m_buffersize, 0, width());
        if (0!=m_buffersize) read_block(0);
    }

    //! Returns whether state of underlying streams are good
    bool good()
    {
        return m_ifile.good() and m_ofile.good();
    }

    //! Returns whether underlying streams are currently associated to a file
    bool is_open()
    {
        return m_ifile.is_open() and m_ofile.is_open();;
    }

    //! Delete all content and set size to 0
    void reset()
    {
        // reset file
        assert(m_ifile.good());
        assert(m_ofile.good());
        m_ifile.close();
        m_ofile.close();
        m_ofile.open(m_filename, std::ios::out|std::ios::binary);
        assert(m_ofile.good());
        m_ifile.open(m_filename, std::ios::in|std::ios::binary);
        assert(m_ifile.good());
        assert(m_ofile.good());
        // reset member variables
        m_need_to_write = false;
        m_size = 0;
        // reset buffer
        read_block(0);
    }

    // Forward declaration
    class reference;

    //! [] operator
    /*! \param i Index the i-th integer of length width().
     *  \return A reference to the i-th integer of length width().
     */
    reference operator[](uint64_t idx)
    {
        return reference(this, idx);
    }

    //! Appends the given element value to the end of the int_vector_buffer
    void push_back(const uint64_t value)
    {
        write(m_size, value);
    }

    //! Close the int_vector_buffer.
    /*! It is not possible to read from / write into the int_vector_buffer after calling this method
     *  \param remove_file If true, the underlying file will be removed on closing.
     */
    void close(bool remove_file=false)
    {
        if (is_open()) {
            if (!remove_file) {
                write_block();
                if (0 < m_offset) { // in case of int_vector, write header and trailing zeros
                    uint64_t size = m_size*width();
                    m_ofile.seekp(0, std::ios::beg);
                    int_vector<t_width>::write_header(size, width(), m_ofile);
                    assert(m_ofile.good());
                    uint64_t wb = (size+7)/8;
                    if (wb%8) {
                        m_ofile.seekp(m_offset+wb);
                        assert(m_ofile.good());
                        m_ofile.write("\0\0\0\0\0\0\0\0", 8-wb%8);
                        assert(m_ofile.good());
                    }
                }
            }
            m_ifile.close();
            assert(m_ifile.good());
            m_ofile.close();
            assert(m_ofile.good());
            if (remove_file) {
                sdsl::remove(m_filename);
            }
        }
    }

    iterator begin()
    {
        return iterator(*this, 0);
    }

    iterator end()
    {
        return iterator(*this, size());
    }

    //! Swap method for int_vector_buffer.
    void swap(int_vector_buffer<t_width>& ivb)
    {
        if (this != &ivb) {
            m_ifile.close();
            ivb.m_ifile.close();
            m_ofile.close();
            ivb.m_ofile.close();
            std::swap(m_filename, ivb.m_filename);
            m_ifile.open(m_filename, std::ios::in|std::ios::binary);
            assert(m_ifile.good());
            m_ofile.open(m_filename, std::ios::in|std::ios::out|std::ios::binary);
            assert(m_ofile.good());
            ivb.m_ifile.open(ivb.m_filename, std::ios::in|std::ios::binary);
            assert(ivb.m_ifile.good());
            ivb.m_ofile.open(ivb.m_filename, std::ios::in|std::ios::out|std::ios::binary);
            assert(ivb.m_ofile.good());
            std::swap(m_buffer, ivb.m_buffer);
            std::swap(m_need_to_write, ivb.m_need_to_write);
            std::swap(m_offset, ivb.m_offset);
            std::swap(m_buffersize, ivb.m_buffersize);
            std::swap(m_size, ivb.m_size);
            std::swap(m_begin, ivb.m_begin);
        }
    }

    class reference
    {
        friend class int_vector_buffer<t_width>;
    private:
        int_vector_buffer<t_width>* const m_int_vector_buffer = nullptr;
        uint64_t m_idx = 0;

        reference() {}

        reference(int_vector_buffer<t_width>* _int_vector_buffer, uint64_t _idx) :
            m_int_vector_buffer(_int_vector_buffer), m_idx(_idx) {}

    public:

        //! Conversion to int for read operations
        operator uint64_t ()const
        {
            return m_int_vector_buffer->read(m_idx);
        }

        //! Assignment operator for write operations
        reference& operator=(const uint64_t& val)
        {
            m_int_vector_buffer->write(m_idx, val);
            return *this;
        }

        //! Assignment operator
        reference& operator=(reference& x)
        {
            return *this = (uint64_t)(x);
        };

        //! Prefix increment of the proxy object
        reference& operator++()
        {
            uint64_t x = m_int_vector_buffer->read(m_idx);
            m_int_vector_buffer->write(m_idx, x+1);
            return *this;
        }

        //! Postfix increment of the proxy object
        uint64_t operator++(int)
        {
            uint64_t val = (uint64_t)*this;
            ++(*this);
            return val;
        }

        //! Prefix decrement of the proxy object
        reference& operator--()
        {
            uint64_t x = m_int_vector_buffer->read(m_idx);
            m_int_vector_buffer->write(m_idx, x-1);
            return *this;
        }

        //! Postfix decrement of the proxy object
        uint64_t operator--(int)
        {
            uint64_t val = (uint64_t)*this;
            --(*this);
            return val;
        }

        //! Add assign from the proxy object
        reference& operator+=(const uint64_t x)
        {
            uint64_t w = m_int_vector_buffer->read(m_idx);
            m_int_vector_buffer->write(m_idx, w+x);
            return *this;
        }

        //! Subtract assign from the proxy object
        reference& operator-=(const uint64_t x)
        {
            uint64_t w = m_int_vector_buffer->read(m_idx);
            m_int_vector_buffer->write(m_idx, w-x);
            return *this;
        }

        bool operator==(const reference& x)const
        {
            return (uint64_t)*this == (uint64_t)x;
        }

        bool operator<(const reference& x)const
        {
            return (uint64_t)*this < (uint64_t)x;
        }
    };

    class iterator: public std::iterator<std::random_access_iterator_tag, value_type, difference_type, value_type*, reference>
    {
    private:
        int_vector_buffer<t_width>& m_ivb;
        uint64_t m_idx = 0;
    public:

        iterator() = delete;
        iterator(int_vector_buffer<t_width>& ivb, uint64_t idx=0) : m_ivb(ivb), m_idx(idx) {}

        iterator& operator++()
        {
            ++m_idx;
            return *this;
        }

        iterator operator++(int)
        {
            iterator it = *this;
            ++(*this);
            return it;
        }

        iterator& operator--()
        {
            --m_idx;
            return *this;
        }

        iterator operator--(int)
        {
            iterator it = *this;
            --(*this);
            return it;
        }

        reference operator*()const
        {
            return m_ivb[m_idx];
        }

        iterator& operator+=(difference_type i)
        {
            if (i<0)
                return *this -= (-i);
            m_idx += i;
            return *this;
        }

        iterator& operator-=(difference_type i)
        {
            if (i<0)
                return *this += (-i);
            m_idx -= i;
            return *this;
        }

        iterator operator+(difference_type i) const
        {
            iterator it = *this;
            return it += i;
        }

        iterator& operator-(difference_type i) const
        {
            iterator it = *this;
            return it -= i;
        }

        bool operator==(const iterator& it) const
        {
            return &m_ivb == &(it.m_ivb) and m_idx == it.m_idx;
        }

        bool operator!=(const iterator& it) const
        {
            return !(*this == it);
        }
        inline difference_type operator-(const iterator& it)
        {
            return (m_idx - it.m_idx);
        }
    };
};

} // end of namespace

#endif // include guard

#endif

/* sdsl - succinct data structures library
    Copyright (C) 2012-2014 Simon Gog
    Copyright (C) 2015 Genome Research Ltd.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!\file sd_vector.hpp
   \brief sd_vector.hpp contains the sdsl::sd_vector class, and
          classes which support rank and select for sd_vector.
   \author Simon Gog, Jouni Siren
*/
#ifndef INCLUDED_SDSL_SD_VECTOR
#define INCLUDED_SDSL_SD_VECTOR

/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file select_support_mcl.hpp
    \brief select_support_mcl.hpp contains classes that support a sdsl::bit_vector with constant time select information.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_MCL
#define INCLUDED_SDSL_SELECT_SUPPORT_MCL

/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file select_support.hpp
    \brief select_support.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT
#define INCLUDED_SDSL_SELECT_SUPPORT

/** \defgroup select_support_group Select Support (SCS)
 * This group contains data structures which support an sdsl::bit_vector with the select method.
 */

/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file rank_support.hpp
    \brief rank_support.hpp contains classes that support a sdsl::bit_vector with constant time rank information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT
#define INCLUDED_SDSL_RANK_SUPPORT

/** \defgroup rank_support_group Rank Support (RS)
 * This group contains data structures which support an sdsl::bit_vector with the rank method.
 */

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! The base class of classes supporting rank_queries for a sdsl::bit_vector in constant time.
/*!
*/
class rank_support
{
protected:
    const bit_vector* m_v; //!< Pointer to the rank supported bit_vector
public:
    typedef bit_vector::size_type size_type;

    //! Constructor
    /*! \param v The supported bit_vector.
     */
    rank_support(const bit_vector* v = nullptr);
    //! Copy constructor
    rank_support(const rank_support&) = default;
    rank_support(rank_support&&) = default;
    rank_support& operator=(const rank_support&) = default;
    rank_support& operator=(rank_support&&) = default;
    //! Destructor
    virtual ~rank_support() {}

    //! Answers rank queries for the supported bit_vector.
    /*!	\param i Argument for the length of the prefix v[0..i-1].
        \returns Number of 1-bits in the prefix [0..i-1] of the supported bit_vector.
        \note Method init has to be called before the first call of rank.
        \sa init
     */
    virtual size_type rank(size_type i) const = 0;
    //! Alias for rank(i)
    virtual size_type operator()(size_type idx) const = 0;
    //! Serializes rank_support.
    /*! \param out Out-Stream to serialize the data to.
    */
    virtual size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const = 0;
    //! Loads the rank_support.
    /*! \param in In-Stream to load the rank_support data from.
        \param v The supported bit_vector.
     */
    virtual void load(std::istream& in, const bit_vector* v=nullptr) = 0;
    //! Sets the supported bit_vector to the given pointer.
    /*! \param v The new bit_vector to support.
     *  \note Method init has to be called before the next call of rank.
     *  \sa init, rank
     */
    virtual void set_vector(const bit_vector* v=nullptr) = 0;
};

inline rank_support::rank_support(const bit_vector* v)
{
    m_v = v;
}

//----------------------------------------------------------------------

template<uint8_t bit_pattern, uint8_t pattern_len>
struct rank_support_trait {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t, uint64_t&)
    {
        return 0;
    }

    static uint32_t word_rank(const uint64_t*, size_type)
    {
        return 0;
    }

    static uint32_t full_word_rank(const uint64_t*, size_type)
    {
        return 0;
    }

    static uint64_t init_carry()
    {
        return 0;
    }
};

template<>
struct rank_support_trait<0,1> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t&)
    {
        return bits::cnt(~w);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx)
    {
        return	bits::cnt((~*(data+(idx>>6))) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx)
    {
        return	bits::cnt((~*(data+(idx>>6))));
    }

    static uint64_t init_carry()
    {
        return 0;
    }
};

template<>
struct rank_support_trait<1,1> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t&)
    {
        return bits::cnt(w);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx)
    {
        return	bits::cnt(*(data+(idx>>6)) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx)
    {
        return	bits::cnt(*(data+(idx>>6)));
    }

    static uint64_t init_carry()
    {
        return 0;
    }
};

template<>
struct rank_support_trait<10,2> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        return bits::cnt10(w, carry);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return	bits::cnt(bits::map10(*data, carry) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return	bits::cnt(bits::map10(*data, carry));
    }

    static uint64_t init_carry()
    {
        return 0;
    }
};

template<>
struct rank_support_trait<01,2> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        return bits::cnt01(w, carry);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 1;
        return	bits::cnt(bits::map01(*data, carry) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 1;
        return	bits::cnt(bits::map01(*data, carry));
    }

    static uint64_t init_carry()
    {
        return 1;
    }
};

template<>
struct rank_support_trait<00,2> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        size_type res = bits::cnt(~(w | (w<<1 | carry)));
        carry = (w >> 63);
        return res;
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 1;
        return bits::cnt((~(*data | ((*data)<<1 | carry))) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 1;
        return bits::cnt(~(*data | ((*data)<<1 | carry)));
    }

    static uint64_t init_carry()
    {
        return 1;
    }
};

template<>
struct rank_support_trait<11,2> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        size_type res = bits::cnt(w & (w<<1 | carry));
        carry = (w >> 63);
        return res;
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return bits::cnt((*data & ((*data)<<1 | carry)) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx)
    {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return bits::cnt(*data & ((*data)<<1 | carry));
    }

    static uint64_t init_carry()
    {
        return 0;
    }
};

}// end namespace sdsl

/* sdsl - succinct data structures library
    Copyright (C) 2008-2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file rank_support_v.hpp
    \brief rank_support_v.hpp contains rank_support_v.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_V
#define INCLUDED_SDSL_RANK_SUPPORT_V

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A rank structure proposed by Sebastiano Vigna
/*! \par Space complexity
 *  \f$ 0.25n\f$ for a bit vector of length n bits.
 *
 * The superblock size is 512. Each superblock is subdivided into 512/64 = 8
 * blocks. So absolute counts for the superblock add 64/512 bits on top of each
 * supported bit. Since the first of the 8 relative count values is 0, we can
 * fit the remaining 7 (each of width log(512)=9) in a 64bit word. The relative
 * counts add another 64/512 bits on top of each supported bit.
 * In total this results in 128/512=25% overhead.
 *
 * \tparam t_b       Bit pattern `0`,`1`,`10`,`01` which should be ranked.
 * \tparam t_pat_len Length of the bit pattern.
 *
 * \par Reference
 *    Sebastiano Vigna:
 *    Broadword Implementation of Rank/Select Queries.
 *    WEA 2008: 154-168
 *
 * @ingroup rank_support_group
 */
template<uint8_t t_b=1, uint8_t t_pat_len=1>
class rank_support_v : public rank_support
{
private:
    static_assert(t_b == 1u or t_b == 0u or t_b == 10u or t_b == 11, "rank_support_v: bit pattern must be `0`,`1`,`10` or `01`");
    static_assert(t_pat_len == 1u or t_pat_len == 2u , "rank_support_v: bit pattern length must be 1 or 2");
public:
    typedef bit_vector                          bit_vector_type;
    typedef rank_support_trait<t_b, t_pat_len>  trait_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = t_pat_len };
private:
    // basic block for interleaved storage of superblockrank and blockrank
    int_vector<64> m_basic_block;
public:
    explicit rank_support_v(const bit_vector* v = nullptr) {
        set_vector(v);
        if (v == nullptr) {
            return;
        } else if (v->empty()) {
            m_basic_block = int_vector<64>(2,0);   // resize structure for basic_blocks
            return;
        }
        size_type basic_block_size = ((v->capacity() >> 9)+1)<<1;
        m_basic_block.resize(basic_block_size);   // resize structure for basic_blocks
        if (m_basic_block.empty())
            return;
        const uint64_t* data = m_v->data();
        size_type i, j=0;
        m_basic_block[0] = m_basic_block[1] = 0;

        uint64_t carry = trait_type::init_carry();
        uint64_t sum   = trait_type::args_in_the_word(*data, carry);
        uint64_t second_level_cnt = 0;
        for (i = 1; i < (m_v->capacity()>>6) ; ++i) {
            if (!(i&0x7)) {// if i%8==0
                j += 2;
                m_basic_block[j-1] = second_level_cnt;
                m_basic_block[j] 	= m_basic_block[j-2] + sum;
                second_level_cnt = sum = 0;
            } else {
                second_level_cnt |= sum<<(63-9*(i&0x7));//  54, 45, 36, 27, 18, 9, 0
            }
            sum += trait_type::args_in_the_word(*(++data), carry);
        }
        if (i&0x7) { // if i%8 != 0
            second_level_cnt |= sum << (63-9*(i&0x7));
            m_basic_block[j+1] = second_level_cnt;
        } else { // if i%8 == 0
            j += 2;
            m_basic_block[j-1] = second_level_cnt;
            m_basic_block[j]   = m_basic_block[j-2] + sum;
            m_basic_block[j+1] = 0;
        }
    }

    rank_support_v(const rank_support_v&)  = default;
    rank_support_v(rank_support_v&&) = default;
    rank_support_v& operator=(const rank_support_v&) = default;
    rank_support_v& operator=(rank_support_v&&) = default;

    size_type rank(size_type idx) const {
        assert(m_v != nullptr);
        assert(idx <= m_v->size());
        const uint64_t* p = m_basic_block.data()
            + ((idx>>8)&0xFFFFFFFFFFFFFFFEULL); // (idx/512)*2
        if (idx&0x3F)  // if (idx%64)!=0
            return  *p + ((*(p+1)>>(63 - 9*((idx&0x1FF)>>6)))&0x1FF) +
                trait_type::word_rank(m_v->data(), idx);
        else
            return  *p + ((*(p+1)>>(63 - 9*((idx&0x1FF)>>6)))&0x1FF);
    }

    inline size_type operator()(size_type idx)const {
        return rank(idx);
    }

    size_type size()const {
        return m_v->size();
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                        std::string name="")const {
        size_type written_bytes = 0;
        structure_tree_node* child = structure_tree::add_child(v, name,
                                                               util::class_name(*this));
        written_bytes += m_basic_block.serialize(out, child,
                                                 "cumulative_counts");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in, const int_vector<1>* v=nullptr) {
        set_vector(v);
        m_basic_block.load(in);
    }

    void set_vector(const bit_vector* v=nullptr) {
        m_v = v;
    }

    void swap(rank_support_v& rs) {
        if (this != &rs) { // if rs and _this_ are not the same object
            m_basic_block.swap(rs.m_basic_block);
        }
    }
};

}// end namespace sds

#endif // end file

/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file rank_support_v5.hpp
    \brief rank_support_v5.hpp contains rank_support_v5.5
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_VFIVE
#define INCLUDED_SDSL_RANK_SUPPORT_VFIVE

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t, uint8_t>
struct rank_support_trait;

//! A class supporting rank queries in constant time.
/*! \par Space complexity
 *  \f$ 0.0625n\f$ bits for a bit vector of length n bits.
 *
 * The superblock size is 2048. Each superblock is subdivided into
 * 2048/(6*64) = 5 blocks (with some bit remaining). So absolute counts for
 * the superblock add 64/2048 bits on top of each supported bit. Since the
 * first of the 6 relative count values is 0, we can fit the remaining 5
 * (each of width log(2048)=11) in a 64 bit word. The relative counts add
 * another 64/2048 bits bits on top of each supported bit. In total this
 * results in 128/2048= 6.25% overhead.
 *
 * \tparam t_b       Bit pattern `0`,`1`,`10`,`01` which should be ranked.
 * \tparam t_pat_len Length of the bit pattern.
 *
 * @ingroup rank_support_group
 */
template<uint8_t t_b=1, uint8_t t_pat_len=1>
class rank_support_v5 : public rank_support
{
private:
    static_assert(t_b == 1u or t_b == 0u or t_b == 10u or t_b == 11u, "rank_support_v5: bit pattern must be `0`,`1`,`10` or `01` or `11`");
    static_assert(t_pat_len == 1u or t_pat_len == 2u , "rank_support_v5: bit pattern length must be 1 or 2");
public:
    typedef bit_vector bit_vector_type;
    typedef rank_support_trait<t_b, t_pat_len>  trait_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = t_pat_len };
private:
//      basic block for interleaved storage of superblockrank and blockrank
    int_vector<64> m_basic_block;
public:
    explicit rank_support_v5(const bit_vector* v = nullptr) {
        set_vector(v);
        if (v == nullptr) {
            return;
        } else if (v->empty()) {
            m_basic_block = int_vector<64>(2,0);
            return;
        }
        size_type basic_block_size = ((v->capacity() >> 11)+1)<<1;
        m_basic_block.resize(basic_block_size);   // resize structure for basic_blocks
        if (m_basic_block.empty())
            return;
        const uint64_t* data = m_v->data();
        size_type i, j=0;
        m_basic_block[0] = m_basic_block[1] = 0;

        uint64_t carry = trait_type::init_carry();
        uint64_t sum   = trait_type::args_in_the_word(*data, carry);
        uint64_t second_level_cnt = 0;
        uint64_t cnt_words=1;
        for (i = 1; i < (m_v->capacity()>>6) ; ++i, ++cnt_words) {
            if (cnt_words == 32) {
                j += 2;
                m_basic_block[j-1] = second_level_cnt;
                m_basic_block[j]     = m_basic_block[j-2] + sum;
                second_level_cnt = sum = cnt_words = 0;
            } else if ((cnt_words%6)==0) {
                // pack the prefix sum for each 6x64bit block into the second_level_cnt
                second_level_cnt |= sum<<(60-12*(cnt_words/6));//  48, 36, 24, 12, 0
            }
            sum += trait_type::args_in_the_word(*(++data), carry);
        }

        if ((cnt_words%6)==0) {
            second_level_cnt |= sum<<(60-12*(cnt_words/6));
        }
        if (cnt_words == 32) {
            j += 2;
            m_basic_block[j-1] = second_level_cnt;
            m_basic_block[j]   = m_basic_block[j-2] + sum;
            m_basic_block[j+1] = 0;
        } else {
            m_basic_block[j+1] = second_level_cnt;
        }
    }

    rank_support_v5(const rank_support_v5&) = default;
    rank_support_v5(rank_support_v5&&) = default;
    rank_support_v5& operator=(const rank_support_v5&) = default;
    rank_support_v5& operator=(rank_support_v5&&) = default;

    size_type rank(size_type idx) const {
        assert(m_v != nullptr);
        assert(idx <= m_v->size());
        const uint64_t* p = m_basic_block.data()
            + ((idx>>10)&0xFFFFFFFFFFFFFFFEULL);// (idx/2048)*2
//                     ( prefix sum of the 6x64bit blocks | (idx%2048)/(64*6) )
        size_type result = *p
            + ((*(p+1)>>(60-12*((idx&0x7FF)/(64*6))))&0x7FFULL)
            + trait_type::word_rank(m_v->data(), idx);
        idx -= (idx&0x3F);
        uint8_t to_do = ((idx>>6)&0x1FULL)%6;
        --idx;
        while (to_do) {
            result += trait_type::full_word_rank(m_v->data(), idx);
            --to_do;
            idx-=64;
        }
        return result;
    }

    inline size_type operator()(size_type idx)const {
        return rank(idx);
    }
    size_type size()const {
        return m_v->size();
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        size_type written_bytes = 0;
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        written_bytes += m_basic_block.serialize(out, child, "cumulative_counts");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in, const bit_vector* v=nullptr) {
        set_vector(v);
        assert(m_v != nullptr); // supported bit vector should be known
        m_basic_block.load(in);
    }

    void set_vector(const bit_vector* v=nullptr) {
        m_v = v;
    }
    //! swap Operator
    void swap(rank_support_v5& rs) {
        if (this != &rs) { // if rs and _this_ are not the same object
            m_basic_block.swap(rs.m_basic_block);
        }
    }
};

}// end namespace sds

#endif // end file

/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file rank_support_scan.hpp
    \brief rank_support_scan.hpp contains rank_support_scan that support a sdsl::bit_vector with linear time rank information.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_SCAN
#define INCLUDED_SDSL_RANK_SUPPORT_SCAN

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A class supporting rank queries in linear time.
/*! \par Space complexity
 *       Constant.
 *  \par Time complexity
 *       Linear in the size of the supported vector.
 *
 *  \tparam t_b       Bit pattern which should be supported. Either `0`,`1`,`10`,`01`.
 *  \tparam t_pat_len Length of the bit pattern.
 * @ingroup rank_support_group
 */
template<uint8_t t_b=1, uint8_t t_pat_len=1>
class rank_support_scan : public rank_support
{
private:
    static_assert(t_b == 1u or t_b == 0u or t_b == 10u or t_b == 11u, "rank_support_scan: bit pattern must be `0`,`1`,`10` or `01`");
    static_assert(t_pat_len == 1u or t_pat_len == 2u , "rank_support_scan: bit pattern length must be 1 or 2");
public:
    typedef bit_vector bit_vector_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = t_pat_len };
public:
    explicit rank_support_scan(const bit_vector* v = nullptr)
    {
        set_vector(v);
    }
    rank_support_scan(const rank_support_scan& rs)
    {
        set_vector(rs.m_v);
    }
    size_type rank(size_type idx) const;
    size_type operator()(size_type idx)const
    {
        return rank(idx);
    };
    size_type size()const
    {
        return m_v->size();
    };
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        return serialize_empty_object(out, v, name, this);
    }
    void load(std::istream&, const int_vector<1>* v=nullptr)
    {
        set_vector(v);
    }
    void set_vector(const bit_vector* v=nullptr)
    {
        m_v=v;
    }

    //! Assign Operator
    rank_support_scan& operator=(const rank_support_scan& rs)
    {
        set_vector(rs.m_v);
        return *this;
    }

    //! swap Operator
    void swap(rank_support_scan&) {}
};

template<uint8_t t_b, uint8_t t_pat_len>
inline typename rank_support_scan<t_b, t_pat_len>::size_type rank_support_scan<t_b, t_pat_len>::rank(size_type idx)const
{
    assert(m_v != nullptr);
    assert(idx <= m_v->size());
    const uint64_t* p   = m_v->data();
    size_type       i   = 0;
    size_type   result  = 0;
    while (i+64 <= idx) {
        result += rank_support_trait<t_b, t_pat_len>::full_word_rank(p, i);
        i += 64;
    }
    return  result+rank_support_trait<t_b, t_pat_len>::word_rank(p, idx);
}

}// end namespace sds

#endif // end file

#endif // end file

//! Namespace for the succinct data structure library.
namespace sdsl
{
//! The base class of classes supporting select queries for a sdsl::bit_vector in constant time.
/*! Abstract base class for classes supporting select queries.
 */
class select_support
{
protected:
    const int_vector<1>* m_v; //!< Pointer to the select supported sdsl::bit_vector.
public:
    typedef int_vector<1>::size_type size_type;
    const bit_vector* vv;

    //! Constructor of select_support.
    /*! \param v The bit_vector to support rank queries.
     */
    select_support(const int_vector<1>* f_v=nullptr):vv(f_v)
    {
        m_v = f_v;
    }
    //! Copy constructor
    /*! Copy the whole select_support including the  pointer
     *  to the supported bit_vector.
     */
    select_support(const select_support& f_v);
    //! Destructor of select_support.
    virtual ~select_support() {};

    //! Select returns the index of the i-th 1-bit in the supported bit_vector.
    /*!	\param i Argument to calculate the index of the i-th 1-bit in the supported bit_vector.
        \return The index \f$\in [0..v.size()-1]\f$ of the i-th 1-bit in the supported bit_vector.
        Call init or load to initialize the data structure before the first call of this method.
         \sa init, load.
     */
    virtual size_type select(size_type i) const = 0;

    //! Alias for select
    virtual size_type operator()(size_type i) const = 0;
    //! Serialize the select_support to an out file stream.
    virtual size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const = 0;
    //! Load the select_support from an in file stream.
    /*!	Load an previously serialized select_support from a std::istream.
        \param in The std::istream to load the select_support.
        \param v The bit_vector to be supported.
        \sa init, select.
     */
    virtual void load(std::istream& in, const int_vector<1>* v=nullptr) = 0;

    //! This method sets the supported bit_vector
    virtual void set_vector(const int_vector<1>* v=nullptr) = 0;
};

template<uint8_t bit_pattern, uint8_t pattern_len>
struct select_support_trait {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector&)
    {
        return 0;
    }

    static size_type args_in_the_first_word(uint64_t, uint8_t, uint64_t)
    {
        return 0;
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t, size_type, uint8_t, uint64_t)
    {
        return 0;
    }

    static size_type args_in_the_word(uint64_t, uint64_t&)
    {
        return 0;
    }

    static size_type ith_arg_pos_in_the_word(uint64_t, size_type, uint64_t)
    {
        return 0;
    }

    static bool found_arg(size_type, const bit_vector&)
    {
        return 0;
    }

    static uint64_t init_carry(const uint64_t*, size_type)
    {
        return 0;
    }

    static uint64_t get_carry(uint64_t)
    {
        return 0;
    }
};

template<>
struct select_support_trait<0,1> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v)
    {
        return v.bit_size()-util::cnt_one_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t)
    {
        return bits::cnt((~w)& bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t)
    {
        return bits::sel(~w & bits::lo_unset[offset], (uint32_t) i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t&)
    {
        return bits::cnt(~w);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t)
    {
        return bits::sel(~w, (uint32_t) i);
    }
    static bool found_arg(size_type i, const bit_vector& v)
    {
        return !v[i];
    }
    static uint64_t init_carry(const uint64_t*, size_type)
    {
        return 0;
    }
    static uint64_t get_carry(uint64_t)
    {
        return 0;
    }
};

template<>
struct select_support_trait<1,1> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v)
    {
        return util::cnt_one_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t)
    {
        return bits::cnt(w & bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t)
    {
        return bits::sel(w & bits::lo_unset[offset], (uint32_t) i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t&)
    {
        return bits::cnt(w);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t)
    {
        return bits::sel(w, (uint32_t) i);
    }
    static bool found_arg(size_type i, const bit_vector& v)
    {
        return v[i] == 1;
    }
    static uint64_t init_carry(const uint64_t*, size_type)
    {
        return 0;
    }
    static uint64_t get_carry(uint64_t)
    {
        return 0;
    }
};

template<>
struct select_support_trait<10,2> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v)
    {
        return util::cnt_onezero_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry)
    {
        return bits::cnt(bits::map10(w, carry) & bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry)
    {
        return bits::sel(bits::map10(w, carry) & bits::lo_unset[offset], (uint32_t) i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        return bits::cnt10(w, carry);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry)
    {
        return bits::sel(bits::map10(w, carry), (uint32_t) i);
    }
    static bool found_arg(size_type i, const bit_vector& v)
    {
        if (i > 0 and v[i-1] and !v[i])
            return true;
        return false;
    }
    static uint64_t init_carry(const uint64_t* data, size_type word_pos)
    {
        return word_pos ? (*(data-1)>>63) : 0;
    }
    static uint64_t get_carry(uint64_t w)
    {
        return w>>63;
    }
};

template<>
struct select_support_trait<01,2> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v)
    {
        return util::cnt_zeroone_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry)
    {
        return bits::cnt(bits::map01(w, carry) & bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry)
    {
        return bits::sel(bits::map01(w, carry) & bits::lo_unset[offset], (uint32_t) i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        return bits::cnt01(w, carry);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry)
    {
        return bits::sel(bits::map01(w, carry), (uint32_t) i);
    }
    static bool found_arg(size_type i, const bit_vector& v)
    {
        if (i > 0 and !v[i-1] and v[i])
            return true;
        return false;
    }
    static uint64_t init_carry(const uint64_t* data, size_type word_pos)
    {
        return word_pos ? (*(data-1)>>63) : 1;
    }
    static uint64_t get_carry(uint64_t w)
    {
        return w>>63;
    }
};

template<>
struct select_support_trait<00,2> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v)
    {
        const uint64_t* data = v.data();
        if (v.empty())
            return 0;
        uint64_t carry = rank_support_trait<00,2>::init_carry();
        size_type result = 0;
        for (auto end = v.data() + (v.size() >> 6); data < end; ++data) {
            result += rank_support_trait<00,2>::args_in_the_word(*data, carry);
        }
        if (v.bit_size()&0x3F) {   // if bit_size is not a multiple of 64, add count of the last word
            result += rank_support_trait<00,2>::args_in_the_word((*data)|bits::lo_unset[v.bit_size()&0x3F], carry);
        }
        return result;
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry)
    {
        size_type res = 0;
        if (offset == 0)
            res = rank_support_trait<00,2>::args_in_the_word(w, carry);
        else {
            res = bits::cnt((~(w | (w<<1))) & bits::lo_unset[offset]);
        }
        return res;
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry)
    {
        return bits::sel((~(((w << 1) | carry) | w)) & bits::lo_unset[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        return rank_support_trait<00,2>::args_in_the_word(w, carry);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry)
    {
        return bits::sel(~(((w << 1) | carry) | w), i);
    }
    static bool found_arg(size_type i, const bit_vector& v)
    {
        return i > 0 and !v[i-1] and !v[i];
    }
    static uint64_t init_carry(const uint64_t* data, size_type word_pos)
    {
        return word_pos ? (*(data-1)>>63) : 1;
    }
    static uint64_t get_carry(uint64_t w)
    {
        return w>>63;
    }
};

template<>
struct select_support_trait<11,2> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v)
    {
        const uint64_t* data = v.data();
        if (v.empty())
            return 0;
        uint64_t carry = rank_support_trait<11,2>::init_carry();
        size_type result = 0;
        for (auto end = v.data() + (v.size() >> 6); data < end; ++data) {
            result += rank_support_trait<11,2>::args_in_the_word(*data, carry);
        }
        if (v.bit_size()&0x3F) {   // if bit_size is not a multiple of 64, add count of the last word
            result += rank_support_trait<11,2>::args_in_the_word((*data)&bits::lo_set[v.bit_size()&0x3F], carry);
        }
        return result;
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry)
    {
        size_type res = 0;
        if (offset == 0)
            res = rank_support_trait<11,2>::args_in_the_word(w, carry);
        else {
            res = bits::cnt((w>>(offset-1)) & (w>>offset));
        }
        return res;
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry)
    {
        return bits::sel((((w << 1) | carry) & w) & bits::lo_unset[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t& carry)
    {
        return rank_support_trait<11,2>::args_in_the_word(w, carry);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry)
    {
        return bits::sel(((w << 1) | carry) & w, i);
    }
    static bool found_arg(size_type i, const bit_vector& v)
    {
        if (i > 0 and v[i-1] and v[i])
            return true;
        return false;
    }
    static uint64_t init_carry(const uint64_t* data, size_type word_pos)
    {
        return word_pos ? (*(data-1)>>63) : 0;
    }
    static uint64_t get_carry(uint64_t w)
    {
        return w>>63;
    }
};

} // end namespace sdsl

/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file select_support_scan.hpp
    \brief select_support_scan.hpp contains classes that support a sdsl::bit_vector with linear time select.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_SCAN
#define INCLUDED_SDSL_SELECT_SUPPORT_SCAN

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A class supporting linear time select queries.
/*! \par Space complexity
 *       Constant.
 *  \par Time complexity
 *       Linear in the size of the supported vector.
 *
 *  \tparam t_b       Bit pattern which should be supported. Either `0`,`1`,`10`,`01`.
 *  \tparam t_pat_len Length of the bit pattern.
 * @ingroup select_support_group
 */
template<uint8_t t_b=1, uint8_t t_pat_len=1>
class select_support_scan : public select_support
{
private:
    static_assert(t_b == 1u or t_b == 0u or t_b == 10u , "select_support_scan: bit pattern must be `0`,`1`,`10` or `01`");
    static_assert(t_pat_len == 1u or t_pat_len == 2u , "select_support_scan: bit pattern length must be 1 or 2");
public:
    typedef bit_vector bit_vector_type;
    enum { bit_pat = t_b };
public:
    explicit select_support_scan(const bit_vector* v=nullptr) : select_support(v) {}
    select_support_scan(const select_support_scan<t_b,t_pat_len>& ss) : select_support(ss.m_v) {}

    inline size_type select(size_type i) const;
    inline size_type operator()(size_type i)const
    {
        return select(i);
    }
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        return serialize_empty_object(out, v, name, this);
    }
    void load(std::istream&, SDSL_UNUSED const bit_vector* v=nullptr)
    {
        set_vector(v);
    }

    void set_vector(const bit_vector* v=nullptr)
    {
        m_v = v;
    }
    select_support_scan<t_b, t_pat_len>& operator=(const select_support_scan& ss)
    {
        set_vector(ss.m_v);
        return *this;
    }
    void swap(select_support_scan<t_b, t_pat_len>&) {}
};

template<uint8_t t_b, uint8_t t_pat_len>
inline typename select_support_scan<t_b,t_pat_len>::size_type select_support_scan<t_b,t_pat_len>::select(size_type i)const
{
    const uint64_t* data = m_v->data();
    size_type word_pos = 0;
    size_type word_off = 0;
    uint64_t carry = select_support_trait<t_b,t_pat_len>::init_carry(data, word_pos);
    size_type args = select_support_trait<t_b,t_pat_len>::args_in_the_first_word(*data, word_off, carry);
    if (args >= i) {
        return (word_pos<<6)+select_support_trait<t_b,t_pat_len>::ith_arg_pos_in_the_first_word(*data, i, word_off, carry);
    }
    word_pos+=1;
    size_type sum_args = args;
    carry = select_support_trait<t_b,t_pat_len>::get_carry(*data);
    uint64_t old_carry = carry;
    args = select_support_trait<t_b,t_pat_len>::args_in_the_word(*(++data), carry);
    while (sum_args + args < i) {
        sum_args += args;
        assert(data+1 < m_v->data() + (m_v->capacity()>>6));
        old_carry = carry;
        args = select_support_trait<t_b,t_pat_len>::args_in_the_word(*(++data), carry);
        word_pos+=1;
    }
    return (word_pos<<6) + select_support_trait<t_b,t_pat_len>::ith_arg_pos_in_the_word(*data, i-sum_args, old_carry);
}

} // end namespace
#endif

#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A class supporting constant time select queries.
/*!
 * \par Space usage
 *      The space usage of the data structure depends on the number of \f$ m \f$ of ones in the
 *      original bitvector $b$. We store the position of every $4096$th set bit
 *      (called L1-sampled bits) of $b$.
 *      This takes in the worst case \f$\frac{m}{4096} \log{n} \leq \frac{n}{64}\f$ bits.
 *      Next,
 *      (1) if the distance of two adjacent L1-sampled bits $b[i]$ and $b[j]$
 *      is greater or equal than $\log^4 n$, then
 *      we store each of the 4096 positions of the set $b$ in [i..j-1] with
 *      $\log{n}$ bits. This results in at most
 *      \$ \frac{4096\cdot \log n}{\log^4 n}=\frac{4096}{\log^3 n}\$ bits per bit.
 *      For a bitvector of 4GB, i.e. \f$ \log n = 35 \f$ we get about 0.01 bits per bit.
 *      If the $j-i+1 < \log^4 n$ then
 *      (2) we store the relative position of every $64$th set bit (called L2-sampled bits)
 *      in b[i..j-1] in at most $4\log\log n$ bits per L2-sampled bits.
 *      An pessimistic upper bound for the space would be
 *      \f$ \frac{4\log\log n}{64} \leq \frac{24}{64} = 0.375\f$ bit per
 *      bit (since $\log\log n\leq 6$. It is very pessimistic, since we store
 *      the relative position in $\log\log(j-i+1)\leq \log\log n$ bits.
 *
 * \tparam t_b       Bit pattern `0`,`1`,`10`,`01` which should be ranked.
 * \tparam t_pat_len Length of the bit pattern.
 *
 * The implementation is a practical variant of the following reference:
 *
 * \par Reference
 *      David Clark:
 *      PhD Thesis: Compact Pat Trees
 *      University of Waterloo, 1996 (Section 2.2.2).
 *      http://www.nlc-bnc.ca/obj/s4/f2/dsk3/ftp04/nq21335.pdf
 *
 * @ingroup select_support_group
 */
template<uint8_t t_b=1, uint8_t t_pat_len=1>
class select_support_mcl : public select_support
{
private:
    static_assert(t_b == 1u or t_b == 0u or t_b == 10u or t_b == 11u, "select_support_mcl: bit pattern must be `0`,`1`,`10`, `01`, or `11`");
    static_assert(t_pat_len == 1u or t_pat_len == 2u , "select_support_mcl: bit pattern length must be 1 or 2");
public:
    typedef bit_vector bit_vector_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = t_pat_len };
private:
    uint32_t m_logn                 = 0,     // \f$ log(size) \f$
    m_logn2                = 0,     // \f$ log^2(size) \f$
    m_logn4                = 0;     // \f$ log^4(size) \f$
    // entry i of m_superblock equals the answer to select_1(B,i*4096)
    int_vector<0> m_superblock;
    int_vector<0>* m_longsuperblock = nullptr;
    int_vector<0>* m_miniblock      = nullptr;
    size_type m_arg_cnt             = 0;
    void copy(const select_support_mcl<t_b, t_pat_len>& ss);
    void initData();
    void init_fast(const bit_vector* v=nullptr);
public:
    explicit select_support_mcl(const bit_vector* v=nullptr);
    select_support_mcl(const select_support_mcl<t_b,t_pat_len>& ss);
    select_support_mcl(select_support_mcl<t_b,t_pat_len>&& ss);
    ~select_support_mcl();
    void init_slow(const bit_vector* v=nullptr);
    //! Select function
    inline size_type select(size_type i) const;
    //! Alias for select(i).
    inline size_type operator()(size_type i)const;
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;
    void load(std::istream& in, const bit_vector* v=nullptr);
    void set_vector(const bit_vector* v=nullptr);
    select_support_mcl<t_b, t_pat_len>& operator=(const select_support_mcl& ss);
    select_support_mcl<t_b, t_pat_len>& operator=(select_support_mcl&&);
    void swap(select_support_mcl<t_b, t_pat_len>& ss);
};

template<uint8_t t_b, uint8_t t_pat_len>
select_support_mcl<t_b,t_pat_len>::select_support_mcl(const bit_vector* f_v):select_support(f_v)
{
    if (t_pat_len>1 or(vv!=nullptr and  vv->size() < 100000))
        init_slow(vv);
    else
        init_fast(vv);
    return;
}

template<uint8_t t_b, uint8_t t_pat_len>
select_support_mcl<t_b,t_pat_len>::select_support_mcl(const select_support_mcl& ss):select_support(ss.m_v)
{
    copy(ss);
}

template<uint8_t t_b, uint8_t t_pat_len>
select_support_mcl<t_b,t_pat_len>::select_support_mcl(select_support_mcl&& ss) : select_support(ss.m_v)
{
    *this = std::move(ss);
}

template<uint8_t t_b, uint8_t t_pat_len>
select_support_mcl<t_b, t_pat_len>& select_support_mcl<t_b,t_pat_len>::operator=(const select_support_mcl& ss)
{
    if (this != &ss) {
        copy(ss);
    }
    return *this;
}

template<uint8_t t_b, uint8_t t_pat_len>
select_support_mcl<t_b, t_pat_len>& select_support_mcl<t_b,t_pat_len>::operator=(select_support_mcl&& ss)
{
    if (this != &ss) {
        m_logn       = ss.m_logn;      // copy log n
        m_logn2      = ss.m_logn2;      // copy (logn)^2
        m_logn4      = ss.m_logn4;      // copy (logn)^4
        m_superblock = std::move(ss.m_superblock); // move long superblock
        m_arg_cnt    = ss.m_arg_cnt;    // copy count of 1-bits
        m_v          = ss.m_v;          // copy pointer to the supported bit vector

        delete [] m_longsuperblock;
        m_longsuperblock = ss.m_longsuperblock;
        ss.m_longsuperblock = nullptr;

        delete [] m_miniblock;
        m_miniblock = ss.m_miniblock;
        ss.m_miniblock = nullptr;
    }
    return *this;
}

template<uint8_t t_b, uint8_t t_pat_len>
void select_support_mcl<t_b,t_pat_len>::swap(select_support_mcl& ss)
{
    std::swap(m_logn, ss.m_logn);
    std::swap(m_logn2, ss.m_logn2);
    std::swap(m_logn4, ss.m_logn4);
    m_superblock.swap(ss.m_superblock);
    std::swap(m_longsuperblock, ss.m_longsuperblock);
    std::swap(m_miniblock, ss.m_miniblock);
    std::swap(m_arg_cnt, ss.m_arg_cnt);
}

template<uint8_t t_b, uint8_t t_pat_len>
void select_support_mcl<t_b,t_pat_len>::copy(const select_support_mcl<t_b, t_pat_len>& ss)
{
    m_logn        = ss.m_logn;      // copy log n
    m_logn2      = ss.m_logn2;      // copy (logn)^2
    m_logn4      = ss.m_logn4;      // copy (logn)^4
    m_superblock = ss.m_superblock; // copy long superblock
    m_arg_cnt    = ss.m_arg_cnt;    // copy count of 1-bits
    m_v          = ss.m_v;          // copy pointer to the supported bit vector
    size_type sb = (m_arg_cnt+4095)>>12;
    delete [] m_longsuperblock;
    m_longsuperblock = nullptr;
    if (ss.m_longsuperblock!=nullptr) {
        m_longsuperblock = new int_vector<0>[sb]; //copy longsuperblocks
        for (size_type i=0; i<sb; ++i) {
            m_longsuperblock[i] = ss.m_longsuperblock[i];
        }
    }
    delete [] m_miniblock;
    m_miniblock = nullptr;
    if (ss.m_miniblock!=nullptr) {
        m_miniblock = new int_vector<0>[sb]; // copy miniblocks
        for (size_type i=0; i<sb; ++i) {
            m_miniblock[i] = ss.m_miniblock[i];
        }
    }
}

template<uint8_t t_b, uint8_t t_pat_len>
select_support_mcl<t_b,t_pat_len>::~select_support_mcl()
{
    delete[] m_longsuperblock;
    delete[] m_miniblock;
}

template<uint8_t t_b, uint8_t t_pat_len>
void select_support_mcl<t_b,t_pat_len>::init_slow(const bit_vector* v)
{
    set_vector(v);
    initData();
    if (m_v==nullptr)
        return;
    // Count the number of arguments in the bit vector
    m_arg_cnt = select_support_trait<t_b,t_pat_len>::arg_cnt(*v);

    const size_type SUPER_BLOCK_SIZE = 4096;

    if (m_arg_cnt==0) // if there are no arguments in the vector we are done...
        return;

    size_type sb = (m_arg_cnt+SUPER_BLOCK_SIZE-1)/SUPER_BLOCK_SIZE; // number of superblocks
    delete [] m_miniblock;
    m_miniblock = new int_vector<0>[sb];

    m_superblock = int_vector<0>(sb, 0, m_logn);

    size_type arg_position[SUPER_BLOCK_SIZE], arg_cnt=0;
    size_type sb_cnt=0;
    for (size_type i=0; i < v->size(); ++i) {
        if (select_support_trait<t_b,t_pat_len>::found_arg(i, *v)) {
            arg_position[ arg_cnt%SUPER_BLOCK_SIZE ] = i;
            assert(arg_position[arg_cnt%SUPER_BLOCK_SIZE] == i);
            ++arg_cnt;
            if (arg_cnt % SUPER_BLOCK_SIZE == 0 or arg_cnt == m_arg_cnt) { //
                assert(sb_cnt < sb);
                m_superblock[sb_cnt] = arg_position[0];

                size_type pos_diff = arg_position[(arg_cnt-1)%SUPER_BLOCK_SIZE]-arg_position[0];
                if (pos_diff > m_logn4) { // longblock
                    if (m_longsuperblock == nullptr) m_longsuperblock = new int_vector<0>[sb]; // create longsuperblock
                    m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0, bits::hi(arg_position[(arg_cnt-1)%SUPER_BLOCK_SIZE]) + 1);

                    for (size_type j=0; j <= (arg_cnt-1)%SUPER_BLOCK_SIZE ; ++j) m_longsuperblock[sb_cnt][j] = arg_position[j]; // copy argument positions to longsuperblock
                } else { // short block
                    m_miniblock[sb_cnt] = int_vector<0>(64, 0, bits::hi(pos_diff)+1);
                    for (size_type j=0; j <= (arg_cnt-1)%SUPER_BLOCK_SIZE; j+=64) {
                        m_miniblock[sb_cnt][j/64] = arg_position[j]-arg_position[0];
                    }
                }
                ++sb_cnt;
            }
        }
    }
}

// TODO: find bug, detected by valgrind
template<uint8_t t_b, uint8_t t_pat_len>
void select_support_mcl<t_b,t_pat_len>::init_fast(const bit_vector* v)
{
    set_vector(v);
    initData();
    if (m_v==nullptr)
        return;
    // Count the number of arguments in the bit vector
    m_arg_cnt = select_support_trait<t_b,t_pat_len>::arg_cnt(*v);

    const size_type SUPER_BLOCK_SIZE = 64*64;

    if (m_arg_cnt==0) // if there are no arguments in the vector we are done...
        return;

//    size_type sb = (m_arg_cnt+63+SUPER_BLOCK_SIZE-1)/SUPER_BLOCK_SIZE; // number of superblocks, add 63 as the last block could contain 63 uninitialized bits
    size_type sb = (m_arg_cnt+SUPER_BLOCK_SIZE-1)/SUPER_BLOCK_SIZE; // number of superblocks
    delete [] m_miniblock;
    m_miniblock = new int_vector<0>[sb];

    m_superblock = int_vector<0>(sb, 0, m_logn);// TODO: hier koennte man logn noch optimieren...s

    bit_vector::size_type arg_position[SUPER_BLOCK_SIZE];
    const uint64_t* data = v->data();
    uint64_t carry_new=0;
    size_type last_k64 = 1, sb_cnt=0;
    for (size_type i=0, cnt_old=0, cnt_new=0, last_k64_sum=1; i < v->capacity(); i+=64, ++data) {
        cnt_new += select_support_trait<t_b, t_pat_len>::args_in_the_word(*data, carry_new);
        if (cnt_new >= last_k64_sum) {
            arg_position[last_k64-1] = i + select_support_trait<t_b, t_pat_len>::ith_arg_pos_in_the_word(*data, last_k64_sum  - cnt_old, carry_new);
            last_k64 += 64;
            last_k64_sum += 64;

            if (last_k64 == SUPER_BLOCK_SIZE+1) {
                m_superblock[sb_cnt] = arg_position[0];
                size_type pos_of_last_arg_in_the_block = arg_position[last_k64-65];

                for (size_type ii=arg_position[last_k64-65]+1, j=last_k64-65; ii < v->size() and j < SUPER_BLOCK_SIZE; ++ii)
                    if (select_support_trait<t_b,t_pat_len>::found_arg(ii, *v)) {
                        pos_of_last_arg_in_the_block = ii;
                        ++j;
                    }
                size_type pos_diff = pos_of_last_arg_in_the_block - arg_position[0];
                if (pos_diff > m_logn4) { // long block
                    if (m_longsuperblock == nullptr) m_longsuperblock = new int_vector<0>[sb+1]; // create longsuperblock
                    // GEANDERT am 2010-07-17 +1 nach pos_of_last_arg..
                    m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0, bits::hi(pos_of_last_arg_in_the_block) + 1);
                    for (size_type j=arg_position[0], k=0; k < SUPER_BLOCK_SIZE and j <= pos_of_last_arg_in_the_block; ++j)
                        if (select_support_trait<t_b, t_pat_len>::found_arg(j, *v)) {
                            if (k>=SUPER_BLOCK_SIZE) {
                                for (size_type ii=0; ii < SUPER_BLOCK_SIZE; ++ii) {
                                    std::cout<<"("<<ii<<","<<m_longsuperblock[sb_cnt][ii]<<") ";
                                }
                                std::cout << std::endl;
                                std::cout<<"k="<<k<<" SUPER_BLOCK_SIZE="<<SUPER_BLOCK_SIZE<<std::endl;
                                std::cout<<"pos_of_last_arg_in_the_block"<< pos_of_last_arg_in_the_block<<std::endl;
                                std::cout.flush();
                            }
                            m_longsuperblock[sb_cnt][k++] = j;
                        }
                } else {
                    m_miniblock[sb_cnt] = int_vector<0>(64, 0, bits::hi(pos_diff)+1);
                    for (size_type j=0; j < SUPER_BLOCK_SIZE; j+=64) {
                        m_miniblock[sb_cnt][j/64] = arg_position[j]-arg_position[0];
                    }
                }
                ++sb_cnt;
                last_k64 = 1;
            }
        }
        cnt_old = cnt_new;
    }
    // handle last block: append long superblock
    if (last_k64 > 1) {
        if (m_longsuperblock == nullptr) m_longsuperblock = new int_vector<0>[sb+1]; // create longsuperblock
        m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0, bits::hi(v->size()-1) + 1);
        for (size_type i=arg_position[0],k=0; i < v->size(); ++i) {
            if (select_support_trait<t_b, t_pat_len>::found_arg(i, *v)) {
                m_longsuperblock[sb_cnt][k++] = i;
            }
        }
        ++sb_cnt;
    }
}

template<uint8_t t_b, uint8_t t_pat_len>
inline auto select_support_mcl<t_b,t_pat_len>::select(size_type i)const -> size_type
{
    assert(i > 0 and i <= m_arg_cnt);

    i = i-1;
    size_type sb_idx = i>>12;   // i/4096
    size_type offset = i&0xFFF; // i%4096
    if (m_longsuperblock!=nullptr and !m_longsuperblock[sb_idx].empty()) {
        return m_longsuperblock[sb_idx][offset];
    } else {
        if ((offset&0x3F)==0) {
            assert(sb_idx < m_superblock.size());
            assert((offset>>6) < m_miniblock[sb_idx].size());
            return m_superblock[sb_idx] + m_miniblock[sb_idx][offset>>6/*/64*/];
        } else {
            i = i-(sb_idx<<12)-((offset>>6)<<6);
            // now i > 0 and i <= 64
            assert(i > 0);
            size_type pos = m_superblock[sb_idx] + m_miniblock[sb_idx][offset>>6] + 1;

            // now pos is the position from where we search for the ith argument
            size_type word_pos = pos>>6;
            size_type word_off = pos&0x3F;
            const uint64_t* data = m_v->data() + word_pos;
            uint64_t carry = select_support_trait<t_b,t_pat_len>::init_carry(data, word_pos);
            size_type args = select_support_trait<t_b,t_pat_len>::args_in_the_first_word(*data, word_off, carry);

            if (args >= i) {
                return (word_pos<<6)+select_support_trait<t_b,t_pat_len>::ith_arg_pos_in_the_first_word(*data, i, word_off, carry);
            }
            word_pos+=1;
            size_type sum_args = args;
            carry = select_support_trait<t_b,t_pat_len>::get_carry(*data);
            uint64_t old_carry = carry;
            args = select_support_trait<t_b,t_pat_len>::args_in_the_word(*(++data), carry);
            while (sum_args + args < i) {
                sum_args += args;
                assert(data+1 < m_v->data() + (m_v->capacity()>>6));
                old_carry = carry;
                args = select_support_trait<t_b,t_pat_len>::args_in_the_word(*(++data), carry);
                word_pos+=1;
            }
            return (word_pos<<6) +
                select_support_trait<t_b,t_pat_len>::ith_arg_pos_in_the_word(*data, i-sum_args, old_carry);
        }
    }
}

template<uint8_t t_b, uint8_t t_pat_len>
inline auto select_support_mcl<t_b,t_pat_len>::operator()(size_type i)const -> size_type
{
    return select(i);
}

template<uint8_t t_b, uint8_t t_pat_len>
void select_support_mcl<t_b,t_pat_len>::initData()
{
    m_arg_cnt = 0;
    if (nullptr == m_v) {
        m_logn = m_logn2 = m_logn4 = 0;
    } else {
        m_logn = bits::hi(m_v->capacity())+1; // TODO maybe it's better here to take a max(...,12)
        m_logn2 = m_logn*m_logn;
        m_logn4 = m_logn2*m_logn2;
    }
    delete[] m_longsuperblock;
    m_longsuperblock = nullptr;
    delete[] m_miniblock;
    m_miniblock = nullptr;
}

template<uint8_t t_b, uint8_t t_pat_len>
void select_support_mcl<t_b,t_pat_len>::set_vector(const bit_vector* v)
{
    m_v = v;
}

template<uint8_t t_b, uint8_t t_pat_len>
auto select_support_mcl<t_b,t_pat_len>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const -> size_type
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    // write the number of 1-bits in the supported bit_vector
    out.write((char*) &m_arg_cnt, sizeof(size_type)/sizeof(char));
    written_bytes = sizeof(size_type)/sizeof(char);
    // number of superblocks in the data structure
    size_type sb = (m_arg_cnt+4095)>>12;

    if (m_arg_cnt) { // if there exists 1-bits to be supported
        written_bytes += m_superblock.serialize(out, child, "superblock"); // serialize superblocks
        bit_vector mini_or_long;// Helper vector: mini or long block?
        if (m_longsuperblock!=nullptr) {
            mini_or_long.resize(sb); // resize indicator bit_vector to the number of superblocks
            for (size_type i=0; i< sb; ++i)
                mini_or_long[i] = !m_miniblock[i].empty();
        }
        written_bytes += mini_or_long.serialize(out, child, "mini_or_long");
        size_type written_bytes_long = 0;
        size_type written_bytes_mini = 0;
        for (size_type i=0; i < sb; ++i)
            if (!mini_or_long.empty() and !mini_or_long[i]) {
                written_bytes_long += m_longsuperblock[i].serialize(out);
            } else {
                written_bytes_mini += m_miniblock[i].serialize(out);
            }
        written_bytes += written_bytes_long;
        written_bytes += written_bytes_mini;
        structure_tree_node* child_long = structure_tree::add_child(child, "longsuperblock", util::class_name(m_longsuperblock));
        structure_tree::add_size(child_long, written_bytes_long);
        structure_tree_node* child_mini = structure_tree::add_child(child, "minisuperblock", util::class_name(m_miniblock));
        structure_tree::add_size(child_mini, written_bytes_mini);
    }
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t t_b, uint8_t t_pat_len>
void select_support_mcl<t_b,t_pat_len>::load(std::istream& in, const bit_vector* v)
{
    set_vector(v);
    initData();
    // read the number of 1-bits in the supported bit_vector
    in.read((char*) &m_arg_cnt, sizeof(size_type)/sizeof(char));
    size_type sb = (m_arg_cnt+4095)>>12;

    if (m_arg_cnt) { // if there exists 1-bits to be supported
        m_superblock.load(in); // load superblocks

        delete[] m_miniblock;
        m_miniblock = nullptr;
        delete[] m_longsuperblock;
        m_longsuperblock = nullptr;

        bit_vector mini_or_long;// Helper vector: mini or long block?
        mini_or_long.load(in); // Load the helper vector
        m_miniblock = new int_vector<0>[sb]; // Create miniblock int_vector<0>
        if (!mini_or_long.empty())
            m_longsuperblock = new int_vector<0>[sb]; // Create longsuperblock int_vector<0>

        for (size_type i=0; i < sb; ++i)
            if (!mini_or_long.empty() and not mini_or_long[i]) {
                m_longsuperblock[i].load(in);
            } else {
                m_miniblock[i].load(in);
            }
    }
}

}

#endif

//! Namespace for the succinct data structure library
namespace sdsl
{

// forward declaration needed for friend declaration
template<uint8_t t_b          = 1,
    class t_hi_bit_vector= bit_vector,
    class t_select_1     = typename t_hi_bit_vector::select_1_type,
    class t_select_0     = typename t_hi_bit_vector::select_0_type>
class rank_support_sd;  // in sd_vector

// forward declaration needed for friend declaration
template<uint8_t t_b          = 1,
    class t_hi_bit_vector= bit_vector,
    class t_select_1     = typename t_hi_bit_vector::select_1_type,
    class t_select_0     = typename t_hi_bit_vector::select_0_type>
class select_support_sd;  // in sd_vector

// forward declaration needed for friend declaration
template<typename, typename, typename>
class sd_vector;  // in sd_vector

//! Class for in-place construction of sd_vector from a strictly increasing sequence
/*! \par Building an sd_vector will clear the builder.
 */
class sd_vector_builder
{
    template<typename, typename, typename>
    friend class sd_vector;

public:
    typedef bit_vector::size_type size_type;

private:
    size_type m_size, m_capacity;
    size_type m_wl;
    size_type m_tail, m_items;
    size_type m_last_high, m_highpos;

    int_vector<> m_low;
    bit_vector   m_high;

public:
    sd_vector_builder();

    //! Constructor
    /*! \param n Vector size.
     *  \param m The number of 1-bits.
     */
    sd_vector_builder(size_type n, size_type m);

    inline size_type size() const { return m_size; }
    inline size_type capacity() const { return m_capacity; }
    inline size_type tail() const { return m_tail; }
    inline size_type items() const { return m_items; }

    //! Set a bit to 1.
    /*! \param i The position of the bit.
     *  \par The position must be strictly greater than for the previous call.
     */
    inline void set(size_type i)
    {
        assert(i >= m_tail && i < m_size);
        assert(m_items < m_capacity);

        size_type cur_high = i >> m_wl;
        m_highpos += (cur_high - m_last_high);
        m_last_high = cur_high;
        m_low[m_items++] = i; // int_vector truncates the most significant logm bits
        m_high[m_highpos++] = 1;  // write 1 for the entry
        m_tail = i + 1;
    }

    //! Swap method
    void swap(sd_vector_builder& sdb);
};

//! A bit vector which compresses very sparse populated bit vectors by
// representing the positions of 1 by the Elias-Fano representation for non-decreasing sequences
/*!
 * \par Other implementations of this data structure:
 *  - the sdarray of Okanohara and Sadakane
 *  - Sebastiano Vigna implemented a elias_fano class in this sux library.
 *
 * \par References
 *  - P. Elias: ,,Efficient storage and retrieval by content and address of static files'',
 *              Journal of the ACM, 1974
 *  - R. Fano: ,,On the number of bits required to implement an associative memory''.
 *             Memorandum 61. Computer Structures Group, Project MAC, MIT, 1971
 *  - D. Okanohara, K. Sadakane: ,,Practical Entropy-Compressed Rank/Select Dictionary'',
 *             Proceedings of ALENEX 2007.
 *
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<class t_hi_bit_vector = bit_vector,
    class t_select_1     = typename t_hi_bit_vector::select_1_type,
    class t_select_0     = typename t_hi_bit_vector::select_0_type>
class sd_vector
{
public:
    typedef bit_vector::size_type                   size_type;
    typedef size_type                               value_type;
    typedef bit_vector::difference_type             difference_type;
    typedef random_access_const_iterator<sd_vector> iterator;
    typedef iterator                                const_iterator;
    typedef bv_tag                                  index_category;
    typedef t_select_0                              select_0_support_type;
    typedef t_select_1                              select_1_support_type;

    typedef rank_support_sd<0, t_hi_bit_vector, select_1_support_type, select_0_support_type> rank_0_type;
    typedef rank_support_sd<1, t_hi_bit_vector, select_1_support_type, select_0_support_type> rank_1_type;
    typedef select_support_sd<0, t_hi_bit_vector, select_1_support_type, select_0_support_type> select_0_type;
    typedef select_support_sd<1, t_hi_bit_vector, select_1_support_type, select_0_support_type> select_1_type;

    typedef t_hi_bit_vector hi_bit_vector_type;
private:
    // we need this variables to represent the m ones of the original bit vector of size n
    size_type m_size = 0;  // length of the original bit vector
    uint8_t   m_wl   = 0;  // log n - log m, where n is the length of the original bit vector
    // and m is the number of ones in the bit vector, wl is the abbreviation
    // for ,,width (of) low (part)''

    int_vector<>          m_low;           // vector for the least significant bits of the positions of the m ones
    hi_bit_vector_type    m_high;          // bit vector that represents the most significant bit in permuted order
    select_1_support_type m_high_1_select; // select support for the ones in m_high
    select_0_support_type m_high_0_select; // select support for the zeros in m_high

    void copy(const sd_vector& v)
    {
        m_size = v.m_size;
        m_wl   = v.m_wl;
        m_low  = v.m_low;
        m_high = v.m_high;
        m_high_1_select = v.m_high_1_select;
        m_high_1_select.set_vector(&m_high);
        m_high_0_select = v.m_high_0_select;
        m_high_0_select.set_vector(&m_high);
    }

public:
    const uint8_t&               wl            = m_wl;
    const hi_bit_vector_type&    high          = m_high;
    const int_vector<>&          low           = m_low;
    const select_1_support_type& high_1_select = m_high_1_select;
    const select_0_support_type& high_0_select = m_high_0_select;

    sd_vector() { }

    sd_vector(const sd_vector& sd)
    {
        copy(sd);
    }

    sd_vector(sd_vector&& sd)
    {
        *this = std::move(sd);
    }

    sd_vector(const bit_vector& bv)
    {
        m_size = bv.size();
        size_type m = util::cnt_one_bits(bv);
        uint8_t logm = bits::hi(m)+1;
        uint8_t logn = bits::hi(m_size)+1;
        if (logm == logn) {
            --logm;    // to ensure logn-logm > 0
        }
        m_wl    = logn - logm;
        m_low = int_vector<>(m, 0, m_wl);
        bit_vector high = bit_vector(m + (1ULL<<logm), 0); //
        const uint64_t* bvp = bv.data();
        for (size_type i=0, mm=0,last_high=0,highpos=0; i < (bv.size()+63)/64; ++i, ++bvp) {
            size_type position = 64*i;
            uint64_t  w = *bvp;
            while (w) {  // process bit_vector word by word
                uint8_t offset = bits::lo(w);
                w >>= offset;   // note:  w >>= (offset+1) can not be applied for offset=63!
                position += offset;
                if (position >= bv.size()) // check that we have not reached the end of the bitvector
                    break;
                // (1) handle high part
                size_type cur_high = position >> m_wl;
                highpos += (cur_high - last_high);   // write cur_high-last_high 0s
                last_high = cur_high;
                // (2) handle low part
                m_low[mm++] = position; // int_vector truncates the most significant logm bits
                high[highpos++] = 1;     // write 1 for the entry
                position += 1;
                w >>= 1;
            }
        }
        util::assign(m_high, high);
        util::init_support(m_high_1_select, &m_high);
        util::init_support(m_high_0_select, &m_high);
    }

    template<class t_itr>
    sd_vector(const t_itr begin,const t_itr end)
    {
        if (begin == end) {
            return;
        }
        if (! std::is_sorted(begin,end)) {
            throw std::runtime_error("sd_vector: source list is not sorted.");
        }
        size_type m = std::distance(begin,end);
        m_size = *(end-1)+1;
        uint8_t logm = bits::hi(m)+1;
        uint8_t logn = bits::hi(m_size)+1;
        if (logm == logn) {
            --logm;    // to ensure logn-logm > 0
        }
        m_wl    = logn - logm;
        m_low = int_vector<>(m, 0, m_wl);
        bit_vector high = bit_vector(m + (1ULL<<logm), 0);
        auto itr = begin;
        size_type mm=0,last_high=0,highpos=0;
        while (itr != end) {
            auto position = *itr;
            // (1) handle high part
            size_type cur_high = position >> m_wl;
            highpos += (cur_high - last_high);   // write cur_high-last_high 0s
            last_high = cur_high;
            // (2) handle low part
            m_low[mm++] = position; // int_vector truncates the most significant logm bits
            high[highpos++] = 1;     // write 1 for the entry
            ++itr;
        }

        util::assign(m_high, high);
        util::init_support(m_high_1_select, &m_high);
        util::init_support(m_high_0_select, &m_high);
    }

    sd_vector(sd_vector_builder& builder)
    {
        if (builder.items() != builder.capacity()) {
            throw std::runtime_error("sd_vector: builder is not at full capacity.");
        }

        m_size = builder.m_size;
        m_wl = builder.m_wl;
        m_low.swap(builder.m_low);
        util::assign(m_high, builder.m_high);
        util::init_support(m_high_1_select, &(this->m_high));
        util::init_support(m_high_0_select, &(this->m_high));

        builder = sd_vector_builder();
    }

    //! Accessing the i-th element of the original bit_vector
    /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
    *   \return The i-th bit of the original bit_vector
    *   \par Time complexity
    *           \f$ \Order{t_{select0} + n/m} \f$, where m equals the number of zeros
    *    \par Remark
     *         The time complexity can be easily improved to
    *            \f$\Order{t_{select0}+\log(n/m)}\f$
    *        by using binary search in the second step.
    */
    value_type operator[](size_type i)const
    {
        size_type high_val = (i >> (m_wl));
        size_type sel_high = m_high_0_select(high_val + 1);
        size_type rank_low = sel_high - high_val;
        if (0 == rank_low)
            return 0;
        size_type val_low = i & bits::lo_set[ m_wl ]; // extract the low m_wl = log n -log m bits
        --sel_high; --rank_low;
        while (m_high[sel_high] and m_low[rank_low] > val_low) {
            if (sel_high > 0) {
                --sel_high; --rank_low;
            } else
                return 0;
        }
        return m_high[sel_high] and m_low[rank_low] == val_low;
    }

    //! Get the integer value of the binary string of length len starting at position idx.
    /*! \param idx Starting index of the binary representation of the integer.
     *  \param len Length of the binary representation of the integer. Default value is 64.
     *  \returns The integer value of the binary string of length len starting at position idx.
     *
     *  \pre idx+len-1 in [0..size()-1]
     *  \pre len in [1..64]
     */
    uint64_t get_int(size_type idx, const uint8_t len=64) const
    {
        uint64_t i = idx+len-1;
        uint64_t high_val = (i >> (m_wl));
        uint64_t sel_high = m_high_0_select(high_val + 1);
        uint64_t rank_low = sel_high - high_val;
        if (0 == rank_low)
            return 0;
        size_type val_low = i & bits::lo_set[ m_wl ]; // extract the low m_wl = log n -log m bits
        --sel_high; --rank_low;
        while (m_high[sel_high] and m_low[rank_low] > val_low) {
            if (sel_high > 0) {
                --sel_high; --rank_low;
            } else
                return 0;
        }
        uint64_t res = 0;
        while (true) {
            while (!m_high[sel_high]) {
                if (sel_high > 0 and(high_val << m_wl) >=idx) {
                    --sel_high; --high_val;
                } else {
                    return res;
                }
            }
            while (m_high[sel_high]) {
                uint64_t val = (high_val << m_wl) + m_low[rank_low];
                if (val >= idx) {
                    res |= 1ULL<<(val-idx);
                } else {
                    return res;
                }
                if (sel_high > 0) {
                    --sel_high; --rank_low;
                } else {
                    return res;
                }
            }
        }
    }

    //! Swap method
    void swap(sd_vector& v)
    {
        if (this != &v) {
            std::swap(m_size, v.m_size);
            std::swap(m_wl, v.m_wl);
            m_low.swap(v.m_low);
            m_high.swap(v.m_high);
            util::swap_support(m_high_1_select, v.m_high_1_select, &m_high, &v.m_high);
            util::swap_support(m_high_0_select, v.m_high_0_select, &m_high, &v.m_high);
        }
    }

    //! Returns the size of the original bit vector.
    size_type size()const
    {
        return m_size;
    }

    sd_vector& operator=(const sd_vector& v)
    {
        if (this != &v) {
            copy(v);
        }
        return *this;
    }

    sd_vector& operator=(sd_vector&& v)
    {
        if (this != &v) {
            m_size = v.m_size;
            m_wl   = v.m_wl;
            m_low  = std::move(v.m_low);
            m_high = std::move(v.m_high);
            m_high_1_select = std::move(v.m_high_1_select);
            m_high_1_select.set_vector(&m_high);
            m_high_0_select = std::move(v.m_high_0_select);
            m_high_0_select.set_vector(&m_high);
        }
        return *this;
    }

    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_size, out, child, "size");
        written_bytes += write_member(m_wl, out, child, "wl");
        written_bytes += m_low.serialize(out, child, "low");
        written_bytes += m_high.serialize(out, child, "high");
        written_bytes += m_high_1_select.serialize(out, child, "high_1_select");
        written_bytes += m_high_0_select.serialize(out, child, "high_0_select");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream& in)
    {
        read_member(m_size, in);
        read_member(m_wl, in);
        m_low.load(in);
        m_high.load(in);
        m_high_1_select.load(in, &m_high);
        m_high_0_select.load(in, &m_high);
    }

    iterator begin() const
    {
        return iterator(this, 0);
    }

    iterator end() const
    {
        return iterator(this, size());
    }
};

//! Specialized constructor that is a bit more space-efficient than the default.
template<> sd_vector<>::sd_vector(sd_vector_builder& builder);

template<uint8_t t_b>
struct rank_support_sd_trait {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r,size_type)
    {
        return r;
    }
};

template<>
struct rank_support_sd_trait<0> {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type n)
    {
        return n - r;
    }
};

//! Rank data structure for sd_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
class rank_support_sd
{
    static_assert(t_b == 1u or t_b == 0u , "rank_support_sd: bit pattern must be `0` or `1`");
public:
    typedef bit_vector::size_type size_type;
    typedef sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = (uint8_t)1 };
private:
    const bit_vector_type* m_v;

public:

    explicit rank_support_sd(const bit_vector_type* v=nullptr)
    {
        set_vector(v);
    }

    size_type rank(size_type i)const
    {
        assert(m_v != nullptr);
        assert(i <= m_v->size());
        // split problem in two parts:
        // (1) find  >=
        size_type high_val = (i >> (m_v->wl));
        size_type sel_high = m_v->high_0_select(high_val + 1);
        size_type rank_low = sel_high - high_val; //
        if (0 == rank_low)
            return rank_support_sd_trait<t_b>::adjust_rank(0, i);
        size_type val_low = i & bits::lo_set[ m_v->wl ];
        // now since rank_low > 0 => sel_high > 0
        do {
            if (!sel_high)
                return rank_support_sd_trait<t_b>::adjust_rank(0, i);
            --sel_high; --rank_low;
        } while (m_v->high[sel_high] and m_v->low[rank_low] >= val_low);
        return rank_support_sd_trait<t_b>::adjust_rank(rank_low+1, i);
    }

    size_type operator()(size_type i)const
    {
        return rank(i);
    }

    size_type size()const
    {
        return m_v->size();
    }

    void set_vector(const bit_vector_type* v=nullptr)
    {
        m_v = v;
    }

    rank_support_sd& operator=(const rank_support_sd& rs)
    {
        if (this != &rs) {
            set_vector(rs.m_v);
        }
        return *this;
    }

    void swap(rank_support_sd&) { }

    void load(std::istream&, const bit_vector_type* v=nullptr)
    {
        set_vector(v);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        return serialize_empty_object(out, v, name, this);
    }
};

template<uint8_t t_b, class t_sd_vec>
struct select_support_sd_trait {
    typedef bit_vector::size_type size_type;
    static size_type select(size_type i, const t_sd_vec* v)
    {
        return v->low[i-1] +  // lower part of the number
            ((v->high_1_select(i) + 1 - i)  << (v->wl));  // upper part
        //^-number of 0 before the i-th 1-^    ^-shift by wl
    }
};

template<class t_sd_vec>
struct select_support_sd_trait<0, t_sd_vec> {
    typedef bit_vector::size_type size_type;
    static size_type select(size_type i, const t_sd_vec* v)
    {
        auto ones  = v->low.size();
        assert(0 < i and i <= v->size() - ones);
        size_type lb = 1, rb = ones+1;
        size_type r0 = 0;
        size_type pos = (size_type)-1;
        // rb exclusive
        // invariant: rank0(select_1(rb)) >= i
        while (lb < rb) {
            auto mid = lb + (rb-lb)/2;
            auto x = select_support_sd_trait<1, t_sd_vec>::select(mid, v);
            auto rank0 = x + 1 - mid;
            if (rank0 >= i) {
                rb = mid;
            } else {
                r0 = rank0;
                pos = x;
                lb = mid + 1;
            }
        }
        return pos + i - r0;
    }
};

//! Select data structure for sd_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
class select_support_sd
{
public:
    typedef bit_vector::size_type size_type;
    typedef sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = (uint8_t)1 };
private:
    const bit_vector_type* m_v;
public:

    explicit select_support_sd(const bit_vector_type* v=nullptr)
    {
        set_vector(v);
    }

    //! Returns the position of the i-th occurrence in the bit vector.
    size_type select(size_type i)const
    {
        return select_support_sd_trait<t_b, bit_vector_type>::select(i, m_v);
    }

    size_type operator()(size_type i)const
    {
        return select(i);
    }

    size_type size()const
    {
        return m_v->size();
    }

    void set_vector(const bit_vector_type* v=nullptr)
    {
        m_v = v;
    }

    void swap(select_support_sd&) { }

    void load(std::istream&, const bit_vector_type* v=nullptr)
    {
        set_vector(v);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        return serialize_empty_object(out, v, name, this);
    }
};

//! Select_0 data structure for sd_vector
/*! \tparam t_sd_vector sd_vector type
 *  \tparam t_rank_1    Rank support for high part of sd_vector
 */
template<typename t_sd_vector=sd_vector<>>
class select_0_support_sd
{
public:
    typedef bit_vector::size_type size_type;
    typedef t_sd_vector           bit_vector_type;
    using rank_1 = typename t_sd_vector::rank_1_type;
    using sel0_type = typename t_sd_vector::select_0_type;
    typedef bit_vector           y_high_type;
    enum { bit_pat = 0 };
    enum { bit_pat_len = (uint8_t)1 };
private:
    const bit_vector_type* m_v;
    int_vector<>           m_pointer;
    int_vector<>           m_rank1;
public:

    explicit select_0_support_sd(const bit_vector_type* v=nullptr)
    {
        set_vector(v);
        if (nullptr != m_v) {
            size_type rank_0 = 0; // rank0 in H
            const size_type bs = 1ULL << (m_v->wl);
            size_type z = 0;
            size_type rank1 = 0;// rank1 in H
            size_type zeros = m_v->size() - rank_1(m_v)(m_v->size()); // zeros in B
            m_pointer = int_vector<>(zeros/(64*bs)+1, 0, bits::hi(m_v->high.size()/64)+1);
            m_rank1   = int_vector<>(m_pointer.size(), 0, bits::hi(m_v->high.size())+1);
            uint64_t w=0;
            for (size_type i=0, sel0=1; i < m_v->high.size(); i+=64) {
                size_type old_rank1 = rank1;
                w = m_v->high.get_int(i, 64);
                rank1 += bits::cnt(w);
                rank_0 = (i+64)-rank1;
                if (rank1 > 0 and (w>>63)&1) {
                    uint64_t pos = rank_0*bs + m_v->low[rank1-1]; // pos of last one (of previous block in B
                    z = pos + 1 - rank1;
                } else {
                    z = rank_0*bs  - rank1;
                }
                while (sel0 <= z and sel0 <= zeros) {
                    m_pointer[(sel0-1)/(64*bs)] = i/64;
                    m_rank1[(sel0-1)/(64*bs)]   = old_rank1;
                    sel0 += 64*bs;
                }
            }
        }
    }

    //! Returns the position of the i-th occurrence in the bit vector.
    size_type select(size_type i)const
    {
        const size_type bs = 1ULL << (m_v->wl);
        size_type j = m_pointer[(i-1)/(64*bs)]*64;// index into m_high
        size_type rank1 = m_rank1[(i-1)/(64*bs)]; // rank_1(j*bs*64) in B
        size_type pos = 0;
        size_type rank0 = 0;

        if (rank1 > 0 and (m_v->high[j-1])&1) {
            pos  = (j-rank1)*bs + m_v->low[rank1-1]; // starting position of current block
            rank0 = pos+1-rank1;
        } else {
            pos  = (j-rank1)*bs;// starting position of current block
            rank0 = pos-rank1;
        }
        uint64_t w = m_v->high.get_int(j, 64);
        do {
            uint64_t _rank1 = rank1 + bits::cnt(w);
            uint64_t _rank0 = 0;
            if (_rank1 > 0 and (w>>63)&1) {
                pos = (j+64-_rank1)*bs + m_v->low[_rank1-1];
                _rank0 = pos+1-_rank1;
            } else {
                pos = (j+64-_rank1)*bs;
                _rank0 = pos-_rank1;
            }
            if (_rank0 < i) {
                j+=64;
                w = m_v->high.get_int(j, 64);
                rank1 = _rank1;
            } else {
                break;
            }
        } while (true);
        // invariant i >zeros
        do {
            uint64_t _rank1 = rank1 + bits::lt_cnt[w&0xFFULL];
            uint64_t _rank0 = 0;
            if (_rank1 > 0 and (w>>7)&1) {
                pos = (j+8-_rank1)*bs + m_v->low[_rank1-1];
                _rank0 = pos+1-_rank1;
            } else {
                pos = (j+8-_rank1)*bs;
                _rank0 = pos-_rank1;
            }
            if (_rank0 < i) {
                j+=8;
                w >>= 8;
                rank1 = _rank1;
            } else {
                break;
            }
        } while (true);

        do {
            bool b = w&1ULL;
            w >>= 1; // zeros are shifted in
            ++j;
            if (0 == b) {
                pos = (j-rank1)*bs;
                size_type zeros = pos-rank1;
                if (zeros >= i) {
                    pos = pos - (zeros-i) - 1;
                    break;
                }
            } else {
                pos = (j-1-rank1)*bs;
                size_type one_pos = pos + m_v->low[rank1];
                ++rank1;
                size_type zeros = one_pos + 1 - rank1;
                if (zeros >= i) {
                    pos = one_pos - (zeros-i) - 1;
                    break;
                }
            }
            if (j%64==0) {
                w = m_v->high.get_int(j,64);
            }
        } while (true);
        return pos;
    }

    size_type operator()(size_type i)const
    {
        return select(i);
    }

    size_type size()const
    {
        return m_v->size();
    }

    void set_vector(const bit_vector_type* v=nullptr)
    {
        m_v = v;
    }

    select_0_support_sd& operator=(const select_0_support_sd& ss)
    {
        if (this != &ss) {
            m_pointer = ss.m_pointer;
            m_rank1   = ss.m_rank1;
            set_vector(ss.m_v);
        }
        return *this;
    }

    void swap(select_0_support_sd& ss)
    {
        m_pointer.swap(ss.m_pointer);
        m_rank1.swap(ss.m_rank1);
    }

    void load(std::istream& in, const bit_vector_type* v=nullptr)
    {
        m_pointer.load(in);
        m_rank1.load(in);
        set_vector(v);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_pointer.serialize(out, child, "pointer");
        written_bytes += m_rank1.serialize(out, child, "rank1");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

} // end namespace
#endif

#include <chrono>
#include <algorithm>

using namespace std::chrono;

namespace sdsl
{

void output_event_json(std::ostream& out,const memory_monitor::mm_event& ev,const memory_monitor& m)
{
    out << "\t\t" << "\"name\" : " << "\"" << ev.name << "\",\n";
    out << "\t\t" << "\"usage\" : [" << "\n";
    for (size_t j=0; j<ev.allocations.size(); j++)  {
        out << "\t\t\t[" << duration_cast<milliseconds>(ev.allocations[j].timestamp-m.start_log).count()
            << "," << ev.allocations[j].usage << "]";
        if (j+1<ev.allocations.size()) {
            out << ",\n";
        } else {
            out << "\n";
        }
    }
    out << "\t\t" << "]\n";
}

template<>
void write_mem_log<JSON_FORMAT>(std::ostream& out,const memory_monitor& m)
{
    auto events = m.completed_events;
    std::sort(events.begin(),events.end());

    // output
    out << "[\n";
    for (size_t i=0; i<events.size(); i++) {
        out << "\t{\n";
        output_event_json(out,events[i],m);
        if (i<events.size()-1) {
            out << "\t},\n";
        } else {
            out << "\t}\n";
        }
    }
    out << "]\n";
}

std::string create_mem_html_header()
{
    std::stringstream jsonheader;
    jsonheader
        << "<html>\n"
        << "<head>\n"
        << "<meta charset=\"utf-8\">\n"
        << "<style>\n"
        << "    body { font: 11px sans-serif; }\n"
        << "    .rule { height: 90%; position: absolute; border-right: 1px dotted #000; text-align: right; }\n"
        << "</style>\n"
        << "<title>sdsl memory usage visualization</title>\n"
        << "<script src=\"http://d3js.org/d3.v3.js\"></script>\n"
        << "</head>\n"
        << "<body marginwidth=\"0\" marginheight=\"0\">\n"
        << "<button><a id=\"download\">Save as SVG</a></button>\n"
        << "<div class=\"chart\"><div id=\"visualization\"></div></div><script>\n";
    return jsonheader.str();
}

std::string create_mem_js_body(const std::string& jsonObject)
{
    std::stringstream jsonbody;
    jsonbody
        << "var events = " << jsonObject << ";\n"
        << "var w = window,d = document,e = d.documentElement,g = d.getElementsByTagName('body')[0],\n"
        << "  xw = w.innerWidth || e.clientWidth || g.clientWidth,\n"
        << "  yh = w.innerHeight || e.clientHeight || g.clientHeight;\n\n"
        << "var margin = {top: 20,right: 80,bottom: 120,left: 120},\n"
        << "  width = xw - margin.left - margin.right,height = yh - margin.top - margin.bottom;\n"
        << "var x = d3.scale.linear().range([0, width]);\n"
        << "var y = d3.scale.linear().range([height, 0]);\n"
        << "var xAxis = d3.svg.axis().scale(x).orient(\"bottom\");\n"
        << "var yAxis = d3.svg.axis().scale(y).orient(\"left\").ticks(5);\n"
        << "var color = d3.scale.category10();\n"
        << "var x_max = d3.max(events, function (d) { return d3.max(d.usage, function (u) { return u[0] / 1000;})})\n"
        << "var y_max = d3.max(events, function (d) { return d3.max(d.usage, function (u) { return 1.1 * u[1] / (1024 * 1024);})})\n"
        << "var peak = d3.max(events, function (d) { return d3.max(d.usage, function (u) { return u[1]; })})\n"
        << "var data = []\nevents.forEach(function (d) { data = data.concat(d.usage); });\n"
        << "var peakelem = data.filter(function (a) { return a[1] == peak; });\n"
        << "var peakelem = peakelem.splice(0,1);\n"
        << "x.domain([0, x_max]);\n y.domain([0, y_max]);\n"
        << "var svg = d3.select(\"#visualization\").append(\"svg\")\n"
        << "  .attr(\"width\", width + margin.left + margin.right)\n"
        << "  .attr(\"height\", height + margin.top + margin.bottom)\n"
        << "  .attr(\"xmlns\", \"http://www.w3.org/2000/svg\")\n"
        << "  .append(\"g\").attr(\"transform\",\"translate(\" + margin.left + \",\" + margin.top + \")\");\n\n"
        << "  svg.append(\"g\").attr(\"class\", \"xaxis\").attr(\"transform\", \"translate(0,\" + height + \")\")\n"
        << "  .call(xAxis).append(\"text\").attr(\"text-anchor\", \"end\")\n"
        << "  .attr(\"shape-rendering\", \"crispEdges\").attr(\"x\", width / 2 + 50).attr(\"y\", 70).attr(\"shape-rendering\", \"crispEdges\")\n"
        << "  .attr(\"font-family\", \"sans-serif\").attr(\"font-size\", \"20px\").text(\"Time (seconds)\");\n\n"
        << "svg.append(\"g\").attr(\"class\", \"yaxis\").call(yAxis).append(\"text\").attr(\"transform\", \"rotate(-90)\").attr(\"x\", -height / 2 + 50)\n"
        << "  .attr(\"y\", -80).attr(\"shape-rendering\", \"crispEdges\").attr(\"font-family\", \"sans-serif\").attr(\"font-size\", \"20px\").style(\"text-anchor\", \"end\")\n"
        << "  .text(\"Memory Usage (MiB)\");\n\n"
        << "svg.selectAll(\".tick text\").style(\"font-size\", \"20px\");\n"
        << "svg.selectAll(\".xaxis .tick text\").attr(\"dy\", 23);\nsvg.selectAll(\".yaxis .tick text\").attr(\"dx\", -10);\n"
        << "svg.selectAll(\"line\").attr(\"fill\", \"none\").attr(\"stroke\", \"black\")\nsvg.selectAll(\"path\").attr(\"fill\", \"none\").attr(\"stroke\", \"black\")\n\n"
        << "svg.selectAll(\"line.horizontalGrid\").data(y.ticks(5)).enter().append(\"line\")\n"
        << "  .attr({\"class\": \"horizontalGrid\",\"x1\": 0,\"x2\": width,\"y1\": function (d) { return y(d);},\n"
        << "     \"y2\": function (d) { return y(d); }, \"fill\": \"none\", \"shape-rendering\": \"crispEdges\",\n"
        << "     \"stroke\": \"lightgrey\",\"stroke-dasharray\": \"10,10\",\"stroke-width\": \"1.5px\"});\n\n"
        << "var area = d3.svg.area().x(function (d) { return x(d[0] / 1000);}).y0(height).y1(function (d) { return y(d[1] / (1024 * 1024))});\n\n"
        << "var ev = svg.selectAll(\".event\").data(events).enter().append(\"svg:path\").attr(\"class\", \"area\")\n"
        << "  .attr(\"fill\", function (d) { return d3.rgb(color(d.name)); })\n"
        << "  .attr(\"d\", function (d) { return area(d.usage) })\n"
        << "  .style(\"stroke\", function (d) { return d3.rgb(color(d.name)).darker(2);}).style(\"stroke-width\", \"2px\")\n\n"
        << "svg.selectAll(\".dot\").data(peakelem).enter().append(\"circle\").attr(\"r\", 3).attr(\"fill\", \"red\")\n"
        << "  .attr(\"cx\", function (d) {return x(d[0] / 1000)})\n"
        << "  .attr(\"cy\", function (d) {return y(d[1] / (1024 * 1024))})\n"
        << "  .attr(\"fill\", \"red\").attr(\"stroke-width\", 2).attr(\"stroke\", \"#cc0000\")\n\n"
        << "svg.selectAll(\".dot\").data(peakelem).enter().append(\"svg:text\")\n"
        << "  .attr(\"x\", function (d) {return x(d[0] / 1000)}).attr(\"y\", function (d) {return y(d[1] / (1024 * 1024) * 1.025)})\n"
        << "  .text(function (d) {return \"Peak Usage: \" + Math.round(d[1] / (1024 * 1024)) + \" MB\"})\n"
        << "  .attr(\"font-size\", 12).attr(\"fill\", \"red\");\n\n"
        << "svg.selectAll(\".dot\").data(peakelem).enter().append(\"circle\")\n"
        << "  .attr(\"r\", 5).attr(\"fill\", \"red\")\n"
        << "  .attr(\"cx\", function (d) {return x(d[0] / 1000)})\n"
        << "  .attr(\"cy\", function (d) {return y(d[1] / (1024 * 1024))})\n"
        << "  .attr(\"fill\", \"none\").attr(\"stroke-width\", 2).attr(\"stroke\", \"#cc0000\").each(pulsepeak());\n\n"
        << "function pulsepeak() { return function (d, i, j) {\n"
        << "  d3.select(this).attr(\"r\", 5).style(\"stroke-opacity\", 1.0).transition()\n"
        << "    .ease(\"linear\").duration(1000).attr(\"r\", 10).style(\"stroke-opacity\", 0.0).each(\"end\", pulsepeak());};}\n\n"
        << "var vertical = d3.select(\".chart\").append(\"div\").attr(\"class\", \"remove\")\n"
        << "  .style(\"position\", \"absolute\").style(\"z-index\", \"19\").style(\"width\", \"1px\")\n"
        << "  .style(\"height\", height - margin).style(\"top\", \"30px\").style(\"bottom\", \"50px\")\n"
        << "  .style(\"left\", \"0px\").style(\"opacity\", \"0.4\").style(\"background\", \"black\");\n\n"
        << "var tooltip = d3.select(\".chart\").append(\"div\").attr(\"class\", \"remove\")\n"
        << "  .style(\"position\", \"absolute\").style(\"z-index\", \"20\").style(\"visibility\", \"hidden\").style(\"top\", \"10px\");\n\n"
        << "var circle = svg.append(\"circle\").attr(\"cx\", 100).attr(\"cy\", 350).attr(\"r\", 3).attr(\"fill\", \"black\").style(\"opacity\", \"0\")\n\n"
        << "d3.select(\"svg\").on(\"mousemove\", function () {\n"
        << "  mousex = d3.mouse(this);\n"
        << "  if (mousex[0] < margin.left + 3 || mousex[0] >= xw - margin.right) {\n"
        << "    vertical.style(\"opacity\", \"0\"); tooltip.style(\"opacity\", \"0\"); circle.style(\"opacity\", \"0\")\n"
        << "  } else {\n"
        << "    var xvalue = x.invert(mousex[0] - margin.left); var pos = findPosition(xvalue)\n"
        << "    vertical.style(\"opacity\", \"0.4\"); tooltip.style(\"opacity\", \"1\"); circle.style(\"opacity\", \"1\")\n"
        << "    circle.attr(\"cx\", pos.x).attr(\"cy\", pos.y); vertical.style(\"left\", mousex[0] + \"px\");tooltip.style(\"left\", mousex[0] + 15 + \"px\")\n"
        << "    tooltip.html(\"<p>\" + xvalue.toFixed(2) + \" Seconds <br>\" + Math.round(pos.mem) + \" MiB <br> \" + pos.name + "
        << "  \"<br> Phase Time: \" + pos.ptime + \" Seconds </p>\").style(\"visibility\", \"visible\");\n"
        << "  }\n})"
        << ".on(\"mouseover\", function () {\n"
        << "  mousex = d3.mouse(this);\n  if (mousex[0] < margin.left + 3 || mousex[0] > xw - margin.right) {\n"
        << "    vertical.style(\"opacity\", \"0\")\n  } else {\n    vertical.style(\"opacity\", \"0.4\");vertical.style(\"left\", mousex[0] + 7 + \"px\")\n}})\n"
        << "d3.select(\"#download\").on(\"click\", function () {\n"
        << "d3.select(this).attr(\"href\", 'data:application/octet-stream;base64,' + btoa(d3.select(\"#visualization\").html())).attr(\"download\", \"viz.svg\")})\n\n"
        << "function findPosition(e){correctArea=d3.selectAll(\".area\").filter(function(t){if(t.usage[0][0]<=e*1e3&&t.usage[t.usage.length-1][0]>=e*1e3){return true}"
        << "return false});if(correctArea.empty()){return 0}var t=new Array;correctArea[0].forEach(function(n){t.push(findYValueinArea(n,e))});"
        << "max_elem=d3.max(t,function(e){return e.mem});var n=t.filter(function(e){return e.mem==max_elem});return n[0]}"
        << "function findYValueinArea(e,t){len=e.getTotalLength();var n=0;var r=len;for(var i=0;i<=len;i+=50){var s=e.getPointAtLength(i);"
        << "var o=x.invert(s.x);var u=y.invert(s.y);if(u>0&&o>t){n=Math.max(0,i-50);r=i;break}}var a=e.getPointAtLength(0);"
        << "var f=1;while(n<r){var l=(r+n)/2;a=e.getPointAtLength(l);target_x=x.invert(a.x);if((l==n||l==r)&&Math.abs(target_x-t)>.01){break}if(target_x>t)r=l;"
        << "else if(target_x<t)n=l;else{break}if(f>50){break}f++}var c=new function(){this.mem=y.invert(a.y);this.name=e.__data__.name;"
        << "this.min=d3.min(e.__data__.usage,function(e){return e[0]/1e3});this.max=d3.max(e.__data__.usage,function(e){return e[0]/1e3});"
        << "this.ptime=Math.round(this.max-this.min);this.x=a.x;this.y=a.y};return c}\n</script></body></html>";
    return jsonbody.str();
}

template<>
void write_mem_log<HTML_FORMAT>(std::ostream& out,const memory_monitor& m)
{
    std::stringstream json_data;
    write_mem_log<JSON_FORMAT>(json_data,m);

    out << create_mem_html_header();
    out << create_mem_js_body(json_data.str());
}

#define ALIGNMENT             sizeof(uint64_t)
#define ALIGNSPLIT(size)      (((size)) & ~0x7)
#define ALIGN(size)           (((size) + (ALIGNMENT-1)) & ~0x7)
#define MM_BLOCK_OVERHEAD     (sizeof(size_t)+sizeof(size_t))
#define MIN_BLOCKSIZE         (ALIGN(sizeof(mm_block_t)+sizeof(mm_block_foot_t)))
#define UNMASK_SIZE(size)     ((size)&~1)
#define ISFREE(size)          ((size)&1)
#define SETFREE(size)         ((size)|1)
#define SPLIT_THRESHOLD       (MIN_BLOCKSIZE)

/* from a memory location get the corresponding block header */
using namespace sdsl;

mm_block_t*
block_cur(void* ptr)
{
    mm_block_t* bptr = (mm_block_t*)((uint8_t*)ptr - sizeof(size_t));
    return bptr;
}

/* given a block retrieve the previous block if any. nullptr otherwise */
mm_block_t*
block_prev(mm_block_t* cur_bptr,mm_block_t* first)
{
    /* start of the heap? */
    if (cur_bptr == first) return nullptr;
    mm_block_foot_t* prev_foot = (mm_block_foot_t*)((uint8_t*)cur_bptr - sizeof(mm_block_foot_t));
    mm_block_t* prev_bptr = (mm_block_t*)((uint8_t*)cur_bptr - UNMASK_SIZE(prev_foot->size));
    return prev_bptr;
}

/* given a block retrieve the next block if any. nullptr otherwise */
mm_block_t*
block_next(mm_block_t* cur_bptr,uint8_t* top)
{
    /* end of the heap? */
    if ((uint8_t*)((uint8_t*)cur_bptr+UNMASK_SIZE(cur_bptr->size)) >= top) return nullptr;

    mm_block_t* next_bptr = (mm_block_t*)((uint8_t*)cur_bptr + UNMASK_SIZE(cur_bptr->size));
    return next_bptr;
}

/* calculate the size of a memory block */
size_t
block_size(void* ptr)
{
    mm_block_t* bptr = block_cur(ptr);
    return UNMASK_SIZE(bptr->size);
}

bool
block_isfree(mm_block_t* ptr)
{
    ;
    return ((ptr->size)&1ULL);
}

/* is the next block free */
bool
block_nextfree(mm_block_t* ptr,uint8_t* top)
{
    mm_block_t* next = block_next(ptr,top);
    if (next && block_isfree(next)) return true;
    return false;
}

/* is the prev block free */
bool
block_prevfree(mm_block_t* ptr,mm_block_t* begin)
{
    mm_block_t* prev = block_prev(ptr,begin);
    if (prev && block_isfree(prev)) return 1;
    return 0;
}

/* update the footer with a new size */
void
foot_update(mm_block_t* ptr,size_t size)
{
    mm_block_foot_t* fptr = (mm_block_foot_t*)((uint8_t*)ptr+
        UNMASK_SIZE(size)-sizeof(mm_block_foot_t));
    fptr->size = size;
}

/* update the block with a new size */
void
block_update(mm_block_t* ptr,size_t size)
{
    ptr->size = size;
    foot_update(ptr,size);
}

/* return the pointer to the "data" */
void*
block_data(mm_block_t* ptr)
{
    return (void*)((uint8_t*)ptr+sizeof(size_t));
}

/* return size of the data that can be stored in the block */
size_t
block_getdatasize(mm_block_t* ptr)
{
    size_t blocksize = UNMASK_SIZE(ptr->size);
    return blocksize - sizeof(size_t) - sizeof(mm_block_foot_t);
}

/* mark the block as free */
void
block_markfree(mm_block_t* ptr)
{
    block_update(ptr,SETFREE(ptr->size));
}

/* mark the block as used */
void
block_markused(mm_block_t* ptr)
{
    block_update(ptr,UNMASK_SIZE(ptr->size));
}

#ifndef MSVC_COMPILER
void
hugepage_allocator::coalesce_block(mm_block_t* block)
{
    //std::cout << "coalesce_block()" << std::endl;
    mm_block_t* newblock = block;
    if (block_nextfree(block,m_top)) {
        mm_block_t* next = block_next(block,m_top);
        /* remove the "next" block from the free list */
        remove_from_free_set(next);
        /* add the size of our block */
        block_update(block,UNMASK_SIZE(block->size)+UNMASK_SIZE(next->size));
    }
    if (block_prevfree(block,m_first_block)) {
        mm_block_t* prev = block_prev(block,m_first_block);
        /* we remove the old prev block and read it to the correct
           size list if necessary */
        remove_from_free_set(prev);
        newblock = prev;
        block_update(prev,UNMASK_SIZE(prev->size)+UNMASK_SIZE(block->size));
    }
    if (newblock) {
        block_markfree(newblock);
        insert_into_free_set(newblock);
    }
}

void
hugepage_allocator::split_block(mm_block_t* bptr,size_t size)
{
    //std::cout << "split_block("<< (void*)bptr << ")" << std::endl;
    size_t blocksize = UNMASK_SIZE(bptr->size);
    //std::cout << "cur_block_size = " << blocksize << std::endl;
    /* only split if we get at least a small block
       out of it */
    int64_t newblocksize = ALIGNSPLIT(blocksize - ALIGN(size+MM_BLOCK_OVERHEAD));
    //std::cout << "new_block_size = " << newblocksize << std::endl;
    if (newblocksize >= (int64_t)SPLIT_THRESHOLD) {
        /* update blocksize of old block */
        //std::cout << "block_update = " << blocksize-newblocksize << std::endl;
        block_update(bptr,blocksize-newblocksize);
        mm_block_t* newblock = (mm_block_t*)((char*)bptr+(blocksize-newblocksize));
        //std::cout << "new block ptr = " << (void*)newblock << std::endl;
        block_update(newblock,newblocksize);
        coalesce_block(newblock);
    }
}

uint8_t*
hugepage_allocator::hsbrk(size_t size)
{
    ptrdiff_t left = (ptrdiff_t) m_total_size - (m_top - m_base);
    if (left < (ptrdiff_t) size) {  // enough space left?
        throw std::system_error(ENOMEM,std::system_category(),
                                "hugepage_allocator: not enough hugepage memory available");
    }
    uint8_t* new_mem = m_top;
    m_top += size;
    return new_mem;
}

mm_block_t*
hugepage_allocator::new_block(size_t size)
{
    //std::cout << "new_block(" << size << ")" << std::endl;
    size = ALIGN(size+MM_BLOCK_OVERHEAD);
    if (size < MIN_BLOCKSIZE) size = MIN_BLOCKSIZE;
    mm_block_t* ptr = (mm_block_t*) hsbrk(size);
    block_update(ptr,size);
    return ptr;
}

mm_block_t*
hugepage_allocator::last_block()
{
    mm_block_t* last = nullptr;
    //std::cout << "m_top = " << (void*)m_top << std::endl;
    //std::cout << "m_base = " << (void*)m_base << std::endl;
    if (m_top != m_base) {
        mm_block_foot_t* fptr = (mm_block_foot_t*)(m_top - sizeof(size_t));
        //std::cout << "foot of last = " << (void*)fptr << std::endl;
        //std::cout << "size of last = " << UNMASK_SIZE(fptr->size) << std::endl;
        last = (mm_block_t*)(((uint8_t*)fptr) - UNMASK_SIZE(fptr->size) + sizeof(size_t));
        //std::cout << "last = " << (void*)last << std::endl;
    }
    return last;
}

void
block_print(int id,mm_block_t* bptr)
{
    fprintf(stdout, "%d addr=%p size=%lu (%lu) free=%d\n",id,((void*)bptr),
            UNMASK_SIZE(bptr->size),bptr->size,block_isfree(bptr));
    fflush(stdout);
}

void
hugepage_allocator::print_heap()
{
    mm_block_t* bptr = m_first_block;
    size_t id = 0;
    while (bptr) {
        block_print(id,bptr);
        id++;
        bptr = block_next(bptr,m_top);
    }
}

void
hugepage_allocator::remove_from_free_set(mm_block_t* block)
{
    //std::cout << "remove_from_free_set()" << std::endl;
    auto eq_range = m_free_large.equal_range(block->size);
    // find the block amoung the blocks with equal size
    auto itr = eq_range.first;
    auto last = eq_range.second;
    auto found = m_free_large.end();
    while (itr != last) {
        if (itr->second == block) {
            found = itr;
        }
        ++itr;
    }
    if (found == m_free_large.end()) {
        found = last;
    }
    m_free_large.erase(found);
}

void
hugepage_allocator::insert_into_free_set(mm_block_t* block)
{
    //std::cout << "insert_into_free_set("<< (void*)block << "," << UNMASK_SIZE(block->size) << ")" << std::endl;
    //std::cout << "insert_into_free_set("<< (void*)block << "," << block->size << ")" << std::endl;
    m_free_large.insert({block->size,block});
}

mm_block_t*
hugepage_allocator::find_free_block(size_t size_in_bytes)
{
    //std::cout << "find_free_block(" << size_in_bytes << ")" << std::endl;

    mm_block_t* bptr = nullptr;
    auto free_block = m_free_large.lower_bound(size_in_bytes);
    if (free_block != m_free_large.end()) {
        bptr = free_block->second;
        m_free_large.erase(free_block);
    }
    return bptr;
}

void*
hugepage_allocator::mm_alloc(size_t size_in_bytes)
{
    //std::cout << "ALLOC(" << size_in_bytes << ")" << std::endl;
    mm_block_t* bptr = nullptr;
    if ((bptr=find_free_block(size_in_bytes + MM_BLOCK_OVERHEAD)) != nullptr) {
        //std::cout << "found free block = " << (void*)bptr << std::endl;
        block_markused(bptr);
        /* split if we have a block too large for us? */
        split_block(bptr,size_in_bytes);
    } else {
        //std::cout << "no free block found that is big enough!" << std::endl;
        // check if last block is free
        //std::cout << "check last block" << std::endl;
        bptr = last_block();
        if (bptr && block_isfree(bptr)) {
            //std::cout << "last block is free. -> extend!" << std::endl;
            // extent last block as it is free
            size_t blockdatasize = block_getdatasize(bptr);
            size_t needed = ALIGN(size_in_bytes - blockdatasize);
            hsbrk(needed);
            remove_from_free_set(bptr);
            block_update(bptr,blockdatasize+needed+sizeof(size_t)+sizeof(mm_block_foot_t));
            //insert_into_free_set(bptr);
            block_markused(bptr);
        } else {
            bptr = new_block(size_in_bytes);
        }
    }
    //print_heap();
    //void* ptr = block_data(bptr);
    //std::cout << "return ptr = " << ptr << std::endl;
    return block_data(bptr);
}

void
hugepage_allocator::mm_free(void* ptr)
{
    //print_heap();
    //std::cout << "FREE(" << ptr << ")" << std::endl;
    if (ptr) {
        mm_block_t* bptr = block_cur(ptr);
        block_markfree(bptr);
        /* coalesce if needed. otherwise just add */
        coalesce_block(bptr);
    }
    //print_heap();
}

void*
hugepage_allocator::mm_realloc(void* ptr, size_t size)
{
    //print_heap();
    //std::cout << "REALLOC(" << ptr << "," << size << ")" << std::endl;
    /* handle special cases first */
    if (nullptr==ptr) return mm_alloc(size);
    if (size==0) {
        mm_free(ptr);
        return nullptr;
    }
    mm_block_t* bptr = block_cur(ptr);

    bool need_malloc = 0;
    size_t blockdatasize = block_getdatasize(bptr);
    /* we do nothing if the size is equal to the block */
    if (size == blockdatasize) {
        //std::cout << "return ptr = " << ptr << std::endl;
        return ptr; /* do nothing if size fits already */
    }
    if (size < blockdatasize) {
        /* we shrink */
        /* do we shrink enough to perform a split? */
        //std::cout << "shrink!" << std::endl;
        split_block(bptr,size);
    } else {
        //std::cout << "expand!" << std::endl;
        /* we expand */
        /* if the next block is free we could use it! */
        mm_block_t* next = block_next(bptr,m_top);
        if (!next) {
            //std::cout << "no next! -> expand!" << std::endl;
            // we are the last block so we just expand
            blockdatasize = block_getdatasize(bptr);
            size_t needed = ALIGN(size - blockdatasize);
            hsbrk(needed);
            block_update(bptr,UNMASK_SIZE(bptr->size)+needed);
            return block_data(bptr);
        } else {
            // we are not the last block
            //std::cout << "try combine next" << std::endl;
            if (next && block_isfree(next)) {
                /* do we have enough space if we use the next block */
                if (blockdatasize + UNMASK_SIZE(next->size) >= size) {
                    /* the next block is enough! */
                    /* remove the "next" block from the free list */
                    remove_from_free_set(next);
                    /* add the size of our block */
                    block_update(bptr,UNMASK_SIZE(bptr->size)+UNMASK_SIZE(next->size));
                } else {
                    /* the next block is not enough. we allocate a new one instead */
                    need_malloc = true;
                }
            } else {
                /* try combing the previous block if free */
                //std::cout << "try combine prev" << std::endl;
                mm_block_t* prev = block_prev(bptr,m_first_block);
                if (prev && block_isfree(prev)) {
                    if (blockdatasize + UNMASK_SIZE(prev->size) >= size) {
                        remove_from_free_set(prev);
                        size_t newsize = UNMASK_SIZE(prev->size)+UNMASK_SIZE(bptr->size);
                        block_update(prev,newsize);
                        block_markused(prev);
                        /* move the data into the previous block */
                        ptr = memmove(block_data(prev),ptr,blockdatasize);
                    } else {
                        /* not enough in the prev block */
                        need_malloc = true;
                    }
                } else {
                    /* prev block not free. get more memory */
                    need_malloc = true;
                }
            }
        }
    }
    if (need_malloc) {
        //std::cout << "need_alloc in REALLOC!" << std::endl;
        void* newptr = mm_alloc(size);
        memcpy(newptr,ptr,size);
        mm_free(ptr);
        ptr = newptr;
    }
    //print_heap();
    //std::cout << "return ptr = " << ptr << std::endl;
    return ptr;
}

uint64_t extract_number(std::string& line)
{
    std::string num_str;
    for (size_t i=line.size()-1; i+1>=1; i--) {
        if (isdigit(line[i])) {
            num_str.insert(num_str.begin(),line[i]);
        } else {
            if (num_str.size() > 0) {
                break;
            }
        }
    }
    return std::strtoull(num_str.c_str(),nullptr,10);
}

uint64_t extract_multiplier(std::string& line)
{
    uint64_t num = 1;
    if (line[line.size()-2] == 'k' || line[line.size()-2] == 'K') {
        num = 1024;
    }
    if (line[line.size()-2] == 'm' || line[line.size()-2] == 'M') {
        num = 1024*1024;
    }
    if (line[line.size()-2] == 'g' || line[line.size()-2] == 'G') {
        num = 1024*1024*1024;
    }
    return num;
}

size_t
hugepage_allocator::determine_available_hugepage_memory()
{
    size_t size_in_bytes = 0;
    size_t page_size_in_bytes = 0;
    size_t num_free_pages = 0;
    const std::string meminfo_file = "/proc/meminfo";
    const std::string ps_str = "Hugepagesize:";
    const std::string pf_str = "HugePages_Free:";
    std::ifstream mifs(meminfo_file);
    if (mifs.is_open()) {
        // find size of one page
        std::string line;
        while (std::getline(mifs, line)) {
            auto ps = std::mismatch(ps_str.begin(),ps_str.end(), line.begin());
            if (ps.first == ps_str.end()) {
                page_size_in_bytes = extract_number(line) * extract_multiplier(line);
            }
            auto pf = std::mismatch(pf_str.begin(),pf_str.end(), line.begin());
            if (pf.first == pf_str.end()) {
                num_free_pages = extract_number(line);
            }
        }
        size_in_bytes = page_size_in_bytes*num_free_pages;
    } else {
        throw std::system_error(ENOMEM,std::system_category(),
                                "hugepage_allocator could not automatically determine available hugepages");
    }
    return size_in_bytes;
}
#endif

}

/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/

namespace sdsl
{

const uint8_t bits::lt_cnt[] = {
    0, 1, 1, 2, 1, 2, 2, 3,
    1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7,
    5, 6, 6, 7, 6, 7, 7, 8
};

const uint32_t bits::lt_deBruijn_to_idx[] = {
    0, 1, 2, 7, 3,13, 8,19,
    4,25,14,28, 9,34,20,40,
    5,17,26,38,15,46,29,48,
    10,31,35,54,21,50,41,57,
    63, 6,12,18,24,27,33,39,
    16,37,45,47,30,53,49,56,
    62,11,23,32,36,44,52,55,
    61,22,43,51,60,42,59,58
};

const uint32_t bits::lt_hi[] = {
    0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};

const uint64_t bits::lo_set[] = {
    0x0000000000000000ULL,
    0x0000000000000001ULL,
    0x0000000000000003ULL,
    0x0000000000000007ULL,
    0x000000000000000FULL,
    0x000000000000001FULL,
    0x000000000000003FULL,
    0x000000000000007FULL,
    0x00000000000000FFULL,
    0x00000000000001FFULL,
    0x00000000000003FFULL,
    0x00000000000007FFULL,
    0x0000000000000FFFULL,
    0x0000000000001FFFULL,
    0x0000000000003FFFULL,
    0x0000000000007FFFULL,
    0x000000000000FFFFULL,
    0x000000000001FFFFULL,
    0x000000000003FFFFULL,
    0x000000000007FFFFULL,
    0x00000000000FFFFFULL,
    0x00000000001FFFFFULL,
    0x00000000003FFFFFULL,
    0x00000000007FFFFFULL,
    0x0000000000FFFFFFULL,
    0x0000000001FFFFFFULL,
    0x0000000003FFFFFFULL,
    0x0000000007FFFFFFULL,
    0x000000000FFFFFFFULL,
    0x000000001FFFFFFFULL,
    0x000000003FFFFFFFULL,
    0x000000007FFFFFFFULL,
    0x00000000FFFFFFFFULL,
    0x00000001FFFFFFFFULL,
    0x00000003FFFFFFFFULL,
    0x00000007FFFFFFFFULL,
    0x0000000FFFFFFFFFULL,
    0x0000001FFFFFFFFFULL,
    0x0000003FFFFFFFFFULL,
    0x0000007FFFFFFFFFULL,
    0x000000FFFFFFFFFFULL,
    0x000001FFFFFFFFFFULL,
    0x000003FFFFFFFFFFULL,
    0x000007FFFFFFFFFFULL,
    0x00000FFFFFFFFFFFULL,
    0x00001FFFFFFFFFFFULL,
    0x00003FFFFFFFFFFFULL,
    0x00007FFFFFFFFFFFULL,
    0x0000FFFFFFFFFFFFULL,
    0x0001FFFFFFFFFFFFULL,
    0x0003FFFFFFFFFFFFULL,
    0x0007FFFFFFFFFFFFULL,
    0x000FFFFFFFFFFFFFULL,
    0x001FFFFFFFFFFFFFULL,
    0x003FFFFFFFFFFFFFULL,
    0x007FFFFFFFFFFFFFULL,
    0x00FFFFFFFFFFFFFFULL,
    0x01FFFFFFFFFFFFFFULL,
    0x03FFFFFFFFFFFFFFULL,
    0x07FFFFFFFFFFFFFFULL,
    0x0FFFFFFFFFFFFFFFULL,
    0x1FFFFFFFFFFFFFFFULL,
    0x3FFFFFFFFFFFFFFFULL,
    0x7FFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFFULL
};

const uint64_t bits::lo_unset[] = {
    0xFFFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFEULL,
    0xFFFFFFFFFFFFFFFCULL,
    0xFFFFFFFFFFFFFFF8ULL,
    0xFFFFFFFFFFFFFFF0ULL,
    0xFFFFFFFFFFFFFFE0ULL,
    0xFFFFFFFFFFFFFFC0ULL,
    0xFFFFFFFFFFFFFF80ULL,
    0xFFFFFFFFFFFFFF00ULL,
    0xFFFFFFFFFFFFFE00ULL,
    0xFFFFFFFFFFFFFC00ULL,
    0xFFFFFFFFFFFFF800ULL,
    0xFFFFFFFFFFFFF000ULL,
    0xFFFFFFFFFFFFE000ULL,
    0xFFFFFFFFFFFFC000ULL,
    0xFFFFFFFFFFFF8000ULL,
    0xFFFFFFFFFFFF0000ULL,
    0xFFFFFFFFFFFE0000ULL,
    0xFFFFFFFFFFFC0000ULL,
    0xFFFFFFFFFFF80000ULL,
    0xFFFFFFFFFFF00000ULL,
    0xFFFFFFFFFFE00000ULL,
    0xFFFFFFFFFFC00000ULL,
    0xFFFFFFFFFF800000ULL,
    0xFFFFFFFFFF000000ULL,
    0xFFFFFFFFFE000000ULL,
    0xFFFFFFFFFC000000ULL,
    0xFFFFFFFFF8000000ULL,
    0xFFFFFFFFF0000000ULL,
    0xFFFFFFFFE0000000ULL,
    0xFFFFFFFFC0000000ULL,
    0xFFFFFFFF80000000ULL,
    0xFFFFFFFF00000000ULL,
    0xFFFFFFFE00000000ULL,
    0xFFFFFFFC00000000ULL,
    0xFFFFFFF800000000ULL,
    0xFFFFFFF000000000ULL,
    0xFFFFFFE000000000ULL,
    0xFFFFFFC000000000ULL,
    0xFFFFFF8000000000ULL,
    0xFFFFFF0000000000ULL,
    0xFFFFFE0000000000ULL,
    0xFFFFFC0000000000ULL,
    0xFFFFF80000000000ULL,
    0xFFFFF00000000000ULL,
    0xFFFFE00000000000ULL,
    0xFFFFC00000000000ULL,
    0xFFFF800000000000ULL,
    0xFFFF000000000000ULL,
    0xFFFE000000000000ULL,
    0xFFFC000000000000ULL,
    0xFFF8000000000000ULL,
    0xFFF0000000000000ULL,
    0xFFE0000000000000ULL,
    0xFFC0000000000000ULL,
    0xFF80000000000000ULL,
    0xFF00000000000000ULL,
    0xFE00000000000000ULL,
    0xFC00000000000000ULL,
    0xF800000000000000ULL,
    0xF000000000000000ULL,
    0xE000000000000000ULL,
    0xC000000000000000ULL,
    0x8000000000000000ULL,
    0x0000000000000000ULL
};

const uint64_t bits::ps_overflow[] = {
    0x8080808080808080ULL,
    0x7f7f7f7f7f7f7f7fULL,
    0x7e7e7e7e7e7e7e7eULL,
    0x7d7d7d7d7d7d7d7dULL,
    0x7c7c7c7c7c7c7c7cULL,
    0x7b7b7b7b7b7b7b7bULL,
    0x7a7a7a7a7a7a7a7aULL,
    0x7979797979797979ULL,
    0x7878787878787878ULL,
    0x7777777777777777ULL,
    0x7676767676767676ULL,
    0x7575757575757575ULL,
    0x7474747474747474ULL,
    0x7373737373737373ULL,
    0x7272727272727272ULL,
    0x7171717171717171ULL,
    0x7070707070707070ULL,
    0x6f6f6f6f6f6f6f6fULL,
    0x6e6e6e6e6e6e6e6eULL,
    0x6d6d6d6d6d6d6d6dULL,
    0x6c6c6c6c6c6c6c6cULL,
    0x6b6b6b6b6b6b6b6bULL,
    0x6a6a6a6a6a6a6a6aULL,
    0x6969696969696969ULL,
    0x6868686868686868ULL,
    0x6767676767676767ULL,
    0x6666666666666666ULL,
    0x6565656565656565ULL,
    0x6464646464646464ULL,
    0x6363636363636363ULL,
    0x6262626262626262ULL,
    0x6161616161616161ULL,
    0x6060606060606060ULL,
    0x5f5f5f5f5f5f5f5fULL,
    0x5e5e5e5e5e5e5e5eULL,
    0x5d5d5d5d5d5d5d5dULL,
    0x5c5c5c5c5c5c5c5cULL,
    0x5b5b5b5b5b5b5b5bULL,
    0x5a5a5a5a5a5a5a5aULL,
    0x5959595959595959ULL,
    0x5858585858585858ULL,
    0x5757575757575757ULL,
    0x5656565656565656ULL,
    0x5555555555555555ULL,
    0x5454545454545454ULL,
    0x5353535353535353ULL,
    0x5252525252525252ULL,
    0x5151515151515151ULL,
    0x5050505050505050ULL,
    0x4f4f4f4f4f4f4f4fULL,
    0x4e4e4e4e4e4e4e4eULL,
    0x4d4d4d4d4d4d4d4dULL,
    0x4c4c4c4c4c4c4c4cULL,
    0x4b4b4b4b4b4b4b4bULL,
    0x4a4a4a4a4a4a4a4aULL,
    0x4949494949494949ULL,
    0x4848484848484848ULL,
    0x4747474747474747ULL,
    0x4646464646464646ULL,
    0x4545454545454545ULL,
    0x4444444444444444ULL,
    0x4343434343434343ULL,
    0x4242424242424242ULL,
    0x4141414141414141ULL,
    0x4040404040404040ULL
};

const uint8_t bits::lt_sel[] = {
    0,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,

    0,0,0,1,0,2,2,1,0,3,3,1,3,2,2,1,
    0,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    0,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    0,6,6,1,6,2,2,1,6,3,3,1,3,2,2,1,
    6,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    6,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    0,7,7,1,7,2,2,1,7,3,3,1,3,2,2,1,
    7,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    7,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    7,6,6,1,6,2,2,1,6,3,3,1,3,2,2,1,
    6,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    6,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,

    0,0,0,0,0,0,0,2,0,0,0,3,0,3,3,2,
    0,0,0,4,0,4,4,2,0,4,4,3,4,3,3,2,
    0,0,0,5,0,5,5,2,0,5,5,3,5,3,3,2,
    0,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
    0,0,0,6,0,6,6,2,0,6,6,3,6,3,3,2,
    0,6,6,4,6,4,4,2,6,4,4,3,4,3,3,2,
    0,6,6,5,6,5,5,2,6,5,5,3,5,3,3,2,
    6,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
    0,0,0,7,0,7,7,2,0,7,7,3,7,3,3,2,
    0,7,7,4,7,4,4,2,7,4,4,3,4,3,3,2,
    0,7,7,5,7,5,5,2,7,5,5,3,5,3,3,2,
    7,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
    0,7,7,6,7,6,6,2,7,6,6,3,6,3,3,2,
    7,6,6,4,6,4,4,2,6,4,4,3,4,3,3,2,
    7,6,6,5,6,5,5,2,6,5,5,3,5,3,3,2,
    6,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
    0,0,0,0,0,0,0,4,0,0,0,4,0,4,4,3,
    0,0,0,0,0,0,0,5,0,0,0,5,0,5,5,3,
    0,0,0,5,0,5,5,4,0,5,5,4,5,4,4,3,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,3,
    0,0,0,6,0,6,6,4,0,6,6,4,6,4,4,3,
    0,0,0,6,0,6,6,5,0,6,6,5,6,5,5,3,
    0,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,3,
    0,0,0,7,0,7,7,4,0,7,7,4,7,4,4,3,
    0,0,0,7,0,7,7,5,0,7,7,5,7,5,5,3,
    0,7,7,5,7,5,5,4,7,5,5,4,5,4,4,3,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,3,
    0,7,7,6,7,6,6,4,7,6,6,4,6,4,4,3,
    0,7,7,6,7,6,6,5,7,6,6,5,6,5,5,3,
    7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
    0,0,0,0,0,0,0,5,0,0,0,5,0,5,5,4,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,4,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,5,
    0,0,0,6,0,6,6,5,0,6,6,5,6,5,5,4,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,4,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,5,
    0,0,0,7,0,7,7,5,0,7,7,5,7,5,5,4,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,4,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,5,
    0,7,7,6,7,6,6,5,7,6,6,5,6,5,5,4,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,5,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,5,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,5,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
};

const uint64_t bits::lt_fib[] = {
    1,
    2,
    3,
    5,
    8,
    13,
    21,
    34,
    55,
    89,
    144,
    233,
    377,
    610,
    987,
    1597,
    2584,
    4181,
    6765,
    10946,
    17711,
    28657,
    46368,
    75025,
    121393,
    196418,
    317811,
    514229,
    832040,
    1346269,
    2178309,
    3524578,
    5702887,
    9227465,
    14930352,
    24157817,
    39088169,
    63245986,
    102334155,
    165580141,
    267914296,
    433494437,
    701408733,
    1134903170,
    1836311903,
    2971215073ULL,
    0x11e8d0a40ULL,
    0x1cfa62f21ULL,
    0x2ee333961ULL,
    0x4bdd96882ULL,
    0x7ac0ca1e3ULL,
    0xc69e60a65ULL,
    0x1415f2ac48ULL,
    0x207fd8b6adULL,
    0x3495cb62f5ULL,
    0x5515a419a2ULL,
    0x89ab6f7c97ULL,
    0xdec1139639ULL,
    0x1686c8312d0ULL,
    0x2472d96a909ULL,
    0x3af9a19bbd9ULL,
    0x5f6c7b064e2ULL,
    0x9a661ca20bbULL,
    0xf9d297a859dULL,
    0x19438b44a658ULL,
    0x28e0b4bf2bf5ULL,
    0x42244003d24dULL,
    0x6b04f4c2fe42ULL,
    0xad2934c6d08fULL,
    0x1182e2989ced1ULL,
    0x1c5575e509f60ULL,
    0x2dd8587da6e31ULL,
    0x4a2dce62b0d91ULL,
    0x780626e057bc2ULL,
    0xc233f54308953ULL,
    0x13a3a1c2360515ULL,
    0x1fc6e116668e68ULL,
    0x336a82d89c937dULL,
    0x533163ef0321e5ULL,
    0x869be6c79fb562ULL,
    0xd9cd4ab6a2d747ULL,
    0x16069317e428ca9ULL,
    0x23a367c34e563f0ULL,
    0x39a9fadb327f099ULL,
    0x5d4d629e80d5489ULL,
    0x96f75d79b354522ULL,
    0xf444c01834299abULL,
    0x18b3c1d91e77decdULL,
    0x27f80ddaa1ba7878ULL,
    0x40abcfb3c0325745ULL,
    0x68a3dd8e61eccfbdULL,
    0xa94fad42221f2702ULL
};

const uint8_t bits::lt_lo[]= {
    0x00,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x05,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x06,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x05,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x07,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x05,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x06,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x05,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x04,0x00,0x01,0x00,0x02,0x00,0x01,0x00,
    0x03,0x00,0x01,0x00,0x02,0x00,0x01,0x00
};

} // end namespace sdsl

/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/

#include <sys/types.h> // for file_size
#include <sys/stat.h>  // for file_size
#include <iomanip>
#include <vector>
#include <string>

#include <type_traits>
#include <typeinfo>
#ifndef MSVC_COMPILER
#include <cxxabi.h>
#endif

namespace sdsl
{

namespace util
{

uint64_t _id_helper::id = 0;

std::string basename(std::string file)
{
    file = disk_file_name(file); // remove RAM-prefix
#ifdef MSVC_COMPILER
    char* c = _strdup((const char*)file.c_str());
    char file_name[_MAX_FNAME] = { 0 };
    ::_splitpath_s(c, NULL, 0, NULL, NULL, file_name, _MAX_FNAME, NULL, 0);
    std::string res(file_name);
#else
    char* c = strdup((const char*)file.c_str());
    std::string res = std::string(::basename(c));
#endif
    free(c);
    return res;
}

std::string dirname(std::string file)
{
    bool ram_file = is_ram_file(file);
    file = disk_file_name(file); // remove RAM-prefix
#ifdef MSVC_COMPILER
    char* c = _strdup((const char*)file.c_str());
    char dir_name[_MAX_DIR] = { 0 };
    char drive[_MAX_DRIVE] = {0};
    ::_splitpath_s(c, drive, _MAX_DRIVE, dir_name, _MAX_DIR, NULL,0, NULL,0);
    std::string res = std::string(drive) + std::string(dir_name);
#else
    char* c = strdup((const char*)file.c_str());
    std::string res = std::string(::dirname(c));
#endif
    free(c);
    if (ram_file) {
        if ("." == res) {
            res = ram_file_name("");
        } else if ("/" ==res) {
            res = ram_file_name(res);
        }
    }
    return res;
}

uint64_t pid()
{
#ifdef MSVC_COMPILER
    return _getpid();
#else
    return getpid();
#endif
}

char* str_from_errno()
{
#ifdef MSVC_COMPILER
    #pragma warning(disable:4996)
    return strerror(errno);
#pragma warning(default:4996)
#else
    return strerror(errno);
#endif
}

uint64_t id()
{
    return _id_helper::getId();
}

std::string demangle(const std::string& name)
{
#ifdef HAVE_CXA_DEMANGLE
    char buf[4096];
    size_t size = 4096;
    int status = 0;
    abi::__cxa_demangle(name.c_str(), buf, &size, &status);
    if (status==0)
        return std::string(buf);
    return name;
#else
    return name;
#endif
}

std::string demangle2(const std::string& name)
{
    std::string result = demangle(name);
    std::vector<std::string> words_to_delete;
    words_to_delete.push_back("sdsl::");
    words_to_delete.push_back("(unsigned char)");
    words_to_delete.push_back(", unsigned long");

    for (size_t k=0; k<words_to_delete.size(); ++k) {
        std::string w = words_to_delete[k];
        for (size_t i = result.find(w); i != std::string::npos; i = result.find(w, i)) {
            result.erase(i, w.length());
            ++i;
        }
    }
    size_t index = 0;
    std::string to_replace = "int_vector<1>";
    while ((index = result.find(to_replace, index)) != std::string::npos) {
        result.replace(index, to_replace.size(), "bit_vector");
    }
    return result;
}

void delete_all_files(tMSS& file_map)
{
    for (auto file_pair : file_map) {
        sdsl::remove(file_pair.second);
    }
    file_map.clear();
}

std::string to_latex_string(unsigned char c)
{
    if (c == '_')
        return "\\_";
    else if (c == '\0')
        return "\\$";
    else
        return to_string(c);
}

void set_verbose()
{
    verbose = true;
}

size_t file_size(const std::string& file)
{
    if (is_ram_file(file)) {
        return ram_fs::file_size(file);
    } else {
        struct stat fs;
        stat(file.c_str(), &fs);
        return fs.st_size;
    }
}

}// end namespace util

}// end namespace sdsl

#include <cstdio>
#include <iostream>
#include <algorithm>

static int nifty_counter = 0;

sdsl::ram_fs::mss_type sdsl::ram_fs::m_map;
std::recursive_mutex sdsl::ram_fs::m_rlock;

sdsl::ram_fs_initializer::ram_fs_initializer()
{
    if (0 == nifty_counter++) {
        if (!ram_fs::m_map.empty()) {
            throw std::logic_error("Static preinitialized object is not empty.");
        }
    }
}

sdsl::ram_fs_initializer::~ram_fs_initializer()
{
    if (0 == --nifty_counter) {
        // clean up
    }
}

namespace sdsl
{

ram_fs::ram_fs() {}

void
ram_fs::store(const std::string& name, content_type data)
{
    std::lock_guard<std::recursive_mutex> lock(m_rlock);
    if (!exists(name)) {
        std::string cname = name;
        m_map.insert(std::make_pair(std::move(cname), std::move(data)));
    } else {
        m_map[name] = std::move(data);
    }
}

bool
ram_fs::exists(const std::string& name)
{
    std::lock_guard<std::recursive_mutex> lock(m_rlock);
    return m_map.find(name) != m_map.end();
}

ram_fs::content_type&
ram_fs::content(const std::string& name)
{
    std::lock_guard<std::recursive_mutex> lock(m_rlock);
    return m_map[name];
}

size_t
ram_fs::file_size(const std::string& name)
{
    std::lock_guard<std::recursive_mutex> lock(m_rlock);
    if (exists(name)) {
        return m_map[name].size();
    } else {
        return 0;
    }
}

int
ram_fs::remove(const std::string& name)
{
    std::lock_guard<std::recursive_mutex> lock(m_rlock);
    m_map.erase(name);
    return 0;
}

int
ram_fs::rename(const std::string old_filename, const std::string new_filename)
{
    std::lock_guard<std::recursive_mutex> lock(m_rlock);
    m_map[new_filename] = std::move(m_map[old_filename]);
    remove(old_filename);
    return 0;
}

bool is_ram_file(const std::string& file)
{
    if (file.size() > 0) {
        if (file[0]=='@') {
            return true;
        }
    }
    return false;
}

std::string ram_file_name(const std::string& file)
{
    if (is_ram_file(file)) {
        return file;
    } else {
        return "@" + file;
    }
}

std::string disk_file_name(const std::string& file)
{
    if (!is_ram_file(file)) {
        return file;
    } else {
        return file.substr(1);
    }
}

int remove(const std::string& file)
{
    if (is_ram_file(file)) {
        return ram_fs::remove(file);
    } else {
        return std::remove(file.c_str());
    }
}

int rename(const std::string& old_filename, const std::string& new_filename)
{
    if (is_ram_file(old_filename)) {
        if (!is_ram_file(new_filename)) {  // error, if new file is not also RAM-file
            return -1;
        }
        return ram_fs::rename(old_filename, new_filename);
    } else {
        return std::rename(old_filename.c_str(), new_filename.c_str());
    }
}

} // end namespace sdsl
