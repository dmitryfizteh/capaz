#include "defines.h"

// Тестовая функция
double Factorial (double x)
{
	return x;
}
//
#ifdef MY_TEST
namespace {

	// Проверить факториал от 0.
	TEST(FactorialTest, HandlesZeroInput) {
		EXPECT_EQ(0, Factorial(0));
	}

	// Проверить факториал некоторых положительных значений.
	TEST(FactorialTest, HandlesPositiveInput) {
		EXPECT_EQ(1, Factorial(1));
		EXPECT_EQ(2, Factorial(2));
		EXPECT_EQ(3, Factorial(3));
		EXPECT_EQ(8, Factorial(8));
	}

}  // namespace
#endif