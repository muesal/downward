#ifndef TEST_MG_STRIPS_H
#define TEST_MG_STRIPS_H

TEST(testMGStripsSetUp);

#define P(N, P) \
    TEST(testMGStrips_##N);
#include "mg_strips_prob.h"

TEST_SUITE(TSMGStrips) {
    TEST_ADD(testMGStripsSetUp),
#define P(N, P) \
    TEST_ADD(testMGStrips_##N),
#include "mg_strips_prob.h"
    TEST_SUITE_CLOSURE
};
#endif
