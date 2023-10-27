/*
 * random.h
 *
 *  Created on: 3 oct. 2022
 *      Author: Romain Simon
 */

#ifndef RANDOM_MT_H_
#define RANDOM_MT_H_

#include <chrono>
#include <random>

// This header-only Random namespace implements a self-seeding Mersenne Twister
// It can be included into as many code files as needed (The inline keyword avoids ODR violations)
// Freely redistributable, courtesy of learncpp.com
namespace Random
{
    // Returns a seeded Mersenne Twister
    // Note: we'd prefer to return a std::seed_seq (to initialize a std::mt19937), but std::seed can't be copied, so it can't be returned by value.
    // Instead, we'll create a std::mt19937, seed it, and then return the std::mt19937 (which can be copied).
    inline std::mt19937 generate()
    {
        std::random_device rd{};

        // Create seed_seq with clock and 7 random numbers from std::random_device
    //std::seed_seq ss{
                    //static_cast<std::seed_seq::result_type>(std::chrono::steady_clock::now().time_since_epoch().count()),
      //          rd(), rd(), rd(), rd(), rd(), rd(), rd() };
        std::cout << "The PNRG seed is: ";
        std::cout << std::to_string(rd()) << "\n";

        //return std::mt19937{ rd() } ;
        return std::mt19937{ rd() } ;
    }

    // Here's our global std::mt19937 object.
    // The inline keyword means we only have one global instance for our whole program.
    inline std::mt19937 mt{ generate() }; // generates a seeded std::mt19937 and copies it into our global object

    // Generate a random int between [min, max] (inclusive)
    inline int intGenerator(int min, int max)
    {
        return std::uniform_int_distribution{min, max}(mt);
    }

    inline double doubleGenerator(double min, double max)
    {
        return std::uniform_real_distribution<double> {min, max}(mt);
    }

    inline std::vector<double> vectorDoubleGenerator(int vectorSize, double min, double max)
    {
        std::vector<double> randomVector;
        randomVector.reserve(vectorSize);
        for (int i=0; i<vectorSize; ++i)
        {
            randomVector.push_back(std::uniform_real_distribution<double> {min, max}(mt));
        }
        return randomVector;
    }

    // The following function templates can be used to generate random numbers
    // when min and/or max are not type int
    // See https://www.learncpp.com/cpp-tutorial/function-template-instantiation/
    // You can ignore these if you don't understand them

    // Generate a random value between [min, max] (inclusive)
    // * min and max have same type
    // * Return value has same type as min and max
    // * Supported types:
    // *    short, int, long, long long
    // *    unsigned short, unsigned int, unsigned long, or unsigned long long
    // Sample call: Random::get(1L, 6L);             // returns long
    // Sample call: Random::get(1u, 6u);             // returns unsigned int

}
#endif /* RANDOM_MT_H_ */
