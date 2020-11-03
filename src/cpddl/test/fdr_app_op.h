#ifndef TEST_FDR_APP_OP_H
#define TEST_FDR_APP_OP_H

TEST(testFDRAppOpSetUp);

#define P(N, P) \
    TEST(testFDRAppOp_##N);
#include "fdr_app_op_prob.h"

TEST_SUITE(TSFDRAppOp) {
    TEST_ADD(testFDRAppOpSetUp),
#define P(N, P) \
    TEST_ADD(testFDRAppOp_##N),
#include "fdr_app_op_prob.h"
    TEST_SUITE_CLOSURE
};
#endif
