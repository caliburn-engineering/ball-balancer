// tests/test_model_library.cpp
#include "test_helpers.h"
#include "analysis/model_library.h"
#include <cstdio>

void test_parse_2x2() {
    auto m = caliburn::parseMatrix("0 1; -2 -3");
    ASSERT_TRUE(m.has_value());
    ASSERT_EQ(m->rows(), 2);
    ASSERT_EQ(m->cols(), 2);
    ASSERT_NEAR((*m)(0, 0), 0.0, 1e-12);
    ASSERT_NEAR((*m)(0, 1), 1.0, 1e-12);
    ASSERT_NEAR((*m)(1, 0), -2.0, 1e-12);
    ASSERT_NEAR((*m)(1, 1), -3.0, 1e-12);
}

void test_parse_1x1() {
    auto m = caliburn::parseMatrix("5");
    ASSERT_TRUE(m.has_value());
    ASSERT_EQ(m->rows(), 1);
    ASSERT_EQ(m->cols(), 1);
    ASSERT_NEAR((*m)(0, 0), 5.0, 1e-12);
}

void test_parse_row_vector() {
    auto m = caliburn::parseMatrix("1 2 3");
    ASSERT_TRUE(m.has_value());
    ASSERT_EQ(m->rows(), 1);
    ASSERT_EQ(m->cols(), 3);
}

void test_parse_with_decimals() {
    auto m = caliburn::parseMatrix("0.5 -1.2; 3.14 0");
    ASSERT_TRUE(m.has_value());
    ASSERT_NEAR((*m)(0, 0), 0.5, 1e-12);
    ASSERT_NEAR((*m)(0, 1), -1.2, 1e-12);
    ASSERT_NEAR((*m)(1, 0), 3.14, 1e-12);
}

void test_parse_invalid_text() {
    auto m = caliburn::parseMatrix("abc");
    ASSERT_TRUE(!m.has_value());
}

void test_parse_ragged_rows() {
    auto m = caliburn::parseMatrix("1 2; 3");
    ASSERT_TRUE(!m.has_value());
}

void test_parse_empty() {
    auto m = caliburn::parseMatrix("");
    ASSERT_TRUE(!m.has_value());
}

void test_builtin_models() {
    auto models = caliburn::getBuiltinModels();
    ASSERT_TRUE(models.size() >= 7);

    // First-Order: 1 state, 1 input, 1 output
    ASSERT_EQ(models[0].system.states(), 1);
    ASSERT_EQ(models[0].system.inputs(), 1);
    ASSERT_EQ(models[0].system.outputs(), 1);

    // Second-Order: 2 states, 1 input, 1 output
    ASSERT_EQ(models[1].system.states(), 2);
    ASSERT_EQ(models[1].system.inputs(), 1);
    ASSERT_EQ(models[1].system.outputs(), 1);

    // Ball-Balancer: 4 states, 2 inputs, 2 outputs
    ASSERT_EQ(models[2].system.states(), 4);
    ASSERT_EQ(models[2].system.inputs(), 2);
    ASSERT_EQ(models[2].system.outputs(), 2);

    // Ball-Balancer Cascade: 7 states, 3 leg inputs, 5 outputs
    ASSERT_EQ(models[3].system.states(), 7);
    ASSERT_EQ(models[3].system.inputs(), 3);
    ASSERT_EQ(models[3].system.outputs(), 5);

    // Inverted Pendulum: 4 states, 1 input, 2 outputs
    ASSERT_EQ(models[4].system.states(), 4);
    ASSERT_EQ(models[4].system.inputs(), 1);
    ASSERT_EQ(models[4].system.outputs(), 2);

    // Quarter-Car: 4 states, 1 input, 2 outputs
    ASSERT_EQ(models[5].system.states(), 4);
    ASSERT_EQ(models[5].system.inputs(), 1);
    ASSERT_EQ(models[5].system.outputs(), 2);

    // Double Mass-Spring: 4 states, 1 input, 2 outputs
    ASSERT_EQ(models[6].system.states(), 4);
    ASSERT_EQ(models[6].system.inputs(), 1);
    ASSERT_EQ(models[6].system.outputs(), 2);
}

// The startup preset is resolved by name, so inserting a model cannot
// silently change which plant the application opens on.
void test_default_model_is_the_cascade() {
    auto models = caliburn::getBuiltinModels();
    const int idx = caliburn::defaultModelIndex(models);
    ASSERT_TRUE(idx >= 0 && idx < (int)models.size());
    ASSERT_TRUE(models[idx].name == "Ball-Balancer Cascade");
}

int main() {
    test_parse_2x2();
    test_parse_1x1();
    test_parse_row_vector();
    test_parse_with_decimals();
    test_parse_invalid_text();
    test_parse_ragged_rows();
    test_parse_empty();
    test_builtin_models();
    test_default_model_is_the_cascade();
    std::printf("All model_library tests passed.\n");
    return 0;
}
