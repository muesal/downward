#ifndef TEST_POT_H
#define TEST_POT_H

TEST(testPotSetUp);

#define P(N, P) \
    TEST(testPotFDR_##N);
#define NCE(N, P) \
    TEST(testPotFDR_##N##_noce_strips);
#include "pot_prob.h"

TEST_SUITE(TSPot) {
    TEST_ADD(testPotSetUp),
#define P(N, P) \
    TEST_ADD(testPotFDR_##N),
#define NCE(N, P) \
    TEST_ADD(testPotFDR_##N##_noce_strips),
#include "pot_prob.h"
    TEST_SUITE_CLOSURE
};
#endif
