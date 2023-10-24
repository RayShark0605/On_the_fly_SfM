#include <gtest/gtest.h>

#include "../../src/Base/LRUCache.h"
#include "../../src/Base/Math.h"
#include "../../src/Base/Base.h"

using namespace std;

TEST(CLRUCache, Empty)
{
	CLRUCache<int, int> cache(5, [](const int key) { return key; });
	EXPECT_EQ(cache.GetNumElements(), 0);
	EXPECT_EQ(cache.GetMaxNumElements(), 5);
}

TEST(CLRUCache, Get) 
{
    CLRUCache<int, int> cache(5, [](const int key) { return key; });
    EXPECT_EQ(cache.GetNumElements(), 0);
    for (int i = 0; i < 5; i++) 
    {
        EXPECT_EQ(cache.Get(i), i);
        EXPECT_EQ(cache.GetNumElements(), i + 1);
        EXPECT_TRUE(cache.IsExists(i));
    }

    EXPECT_EQ(cache.Get(5), 5);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_TRUE(cache.IsExists(5));

    EXPECT_EQ(cache.Get(5), 5);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_TRUE(cache.IsExists(5));

    EXPECT_EQ(cache.Get(6), 6);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_FALSE(cache.IsExists(1));
    EXPECT_TRUE(cache.IsExists(6));
}

TEST(CLRUCache, GetRefer1) 
{
    CLRUCache<int, int> cache(5, [](const int key) { return key; });
    EXPECT_EQ(cache.GetNumElements(), 0);
    for (int i = 0; i < 5; i++) 
    {
        EXPECT_EQ(cache.GetRefer(i), i);
        EXPECT_EQ(cache.GetNumElements(), i + 1);
        EXPECT_TRUE(cache.IsExists(i));
    }

    EXPECT_EQ(cache.GetRefer(5), 5);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_TRUE(cache.IsExists(5));

    EXPECT_EQ(cache.GetRefer(5), 5);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_TRUE(cache.IsExists(5));

    EXPECT_EQ(cache.GetRefer(6), 6);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_FALSE(cache.IsExists(1));
    EXPECT_TRUE(cache.IsExists(6));

    cache.GetRefer(6) = 66;
    EXPECT_EQ(cache.GetRefer(6), 66);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_FALSE(cache.IsExists(1));
    EXPECT_TRUE(cache.IsExists(6));
}

TEST(CLRUCache, Set) 
{
    CLRUCache<int, int> cache(5, [](const int key) { return -1; });
    EXPECT_EQ(cache.GetNumElements(), 0);
    for (int i = 0; i < 5; i++) 
    {
        cache.Set(i, i);
        EXPECT_EQ(cache.GetNumElements(), i + 1);
        EXPECT_TRUE(cache.IsExists(i));
    }

    EXPECT_EQ(cache.Get(5), -1);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_TRUE(cache.IsExists(5));

    EXPECT_EQ(cache.Get(6), -1);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_FALSE(cache.IsExists(1));
    EXPECT_TRUE(cache.IsExists(6));
}

TEST(CLRUCache, Pop) 
{
    CLRUCache<int, int> cache(5, [](const int key) { return key; });
    EXPECT_EQ(cache.GetNumElements(), 0);
    for (int i = 0; i < 5; i++) 
    {
        EXPECT_EQ(cache.Get(i), i);
        EXPECT_EQ(cache.GetNumElements(), i + 1);
        EXPECT_TRUE(cache.IsExists(i));
    }

    EXPECT_EQ(cache.Get(5), 5);
    EXPECT_EQ(cache.GetNumElements(), 5);
    EXPECT_FALSE(cache.IsExists(0));
    EXPECT_TRUE(cache.IsExists(5));

    cache.Pop();
    EXPECT_EQ(cache.GetNumElements(), 4);
    cache.Pop();
    EXPECT_EQ(cache.GetNumElements(), 3);
    cache.Pop();
    EXPECT_EQ(cache.GetNumElements(), 2);
    cache.Pop();
    EXPECT_EQ(cache.GetNumElements(), 1);
    cache.Pop();
    EXPECT_EQ(cache.GetNumElements(), 0);
    cache.Pop();
    EXPECT_EQ(cache.GetNumElements(), 0);
}

TEST(CLRUCache, Clear) 
{
    CLRUCache<int, int> cache(5, [](const int key) { return key; });
    EXPECT_EQ(cache.GetNumElements(), 0);
    for (int i = 0; i < 5; i++) 
    {
        EXPECT_EQ(cache.Get(i), i);
        EXPECT_EQ(cache.GetNumElements(), i + 1);
        EXPECT_TRUE(cache.IsExists(i));
    }

    cache.Clear();
    EXPECT_EQ(cache.GetNumElements(), 0);

    EXPECT_EQ(cache.Get(0), 0);
    EXPECT_EQ(cache.GetNumElements(), 1);
    EXPECT_TRUE(cache.IsExists(0));
}

TEST(EnsureTrailingSlash, Nominal) 
{
    EXPECT_EQ(EnsureTrailingSlash(""), "/");
    EXPECT_EQ(EnsureTrailingSlash("/"), "/");
    EXPECT_EQ(EnsureTrailingSlash("////"), "////");
    EXPECT_EQ(EnsureTrailingSlash("test"), "test/");
    EXPECT_EQ(EnsureTrailingSlash("/test"), "/test/");
}

TEST(HasFileExtension, Nominal) 
{
    EXPECT_FALSE(HasFileExtension("", ".jpg"));
    EXPECT_FALSE(HasFileExtension("testjpg", ".jpg"));
    EXPECT_TRUE(HasFileExtension("test.jpg", ".jpg"));
    EXPECT_TRUE(HasFileExtension("test.jpg", ".Jpg"));
    EXPECT_TRUE(HasFileExtension("test.jpg", ".JPG"));
    EXPECT_TRUE(HasFileExtension("test.", "."));
}

TEST(SplitFileExtension, Nominal) 
{
    string root;
    string ext;
    SplitFileExtension("", root, ext);
    EXPECT_EQ(root, "");
    EXPECT_EQ(ext, "");
    SplitFileExtension(".", root, ext);
    EXPECT_EQ(root, "");
    EXPECT_EQ(ext, "");
    SplitFileExtension("file", root, ext);
    EXPECT_EQ(root, "file");
    EXPECT_EQ(ext, "");
    SplitFileExtension("file.", root, ext);
    EXPECT_EQ(root, "file");
    EXPECT_EQ(ext, "");
    SplitFileExtension("file.jpg", root, ext);
    EXPECT_EQ(root, "file");
    EXPECT_EQ(ext, ".jpg");
    SplitFileExtension("dir/file.jpg", root, ext);
    EXPECT_EQ(root, "dir/file");
    EXPECT_EQ(ext, ".jpg");
    SplitFileExtension("/dir/file.jpg", root, ext);
    EXPECT_EQ(root, "/dir/file");
    EXPECT_EQ(ext, ".jpg");
    SplitFileExtension("dir/file.suffix.jpg", root, ext);
    EXPECT_EQ(root, "dir/file.suffix");
    EXPECT_EQ(ext, ".jpg");
    SplitFileExtension("dir.suffix/file.suffix.jpg", root, ext);
    EXPECT_EQ(root, "dir.suffix/file.suffix");
    EXPECT_EQ(ext, ".jpg");
    SplitFileExtension("dir.suffix/file.", root, ext);
    EXPECT_EQ(root, "dir.suffix/file");
    EXPECT_EQ(ext, "");
    SplitFileExtension("./dir.suffix/file.", root, ext);
    EXPECT_EQ(root, "./dir.suffix/file");
    EXPECT_EQ(ext, "");
}

TEST(GetPathBaseName, Nominal) 
{
    EXPECT_EQ(GetPathBaseName(""), "");
    EXPECT_EQ(GetPathBaseName("test"), "test");
    EXPECT_EQ(GetPathBaseName("/test"), "test");
    EXPECT_EQ(GetPathBaseName("test/"), "test");
    EXPECT_EQ(GetPathBaseName("/test/"), "test");
    EXPECT_EQ(GetPathBaseName("test1/test2"), "test2");
    EXPECT_EQ(GetPathBaseName("/test1/test2"), "test2");
    EXPECT_EQ(GetPathBaseName("/test1/test2/"), "test2");
    EXPECT_EQ(GetPathBaseName("/test1/test2/"), "test2");
    EXPECT_EQ(GetPathBaseName("\\test1/test2/"), "test2");
    EXPECT_EQ(GetPathBaseName("\\test1\\test2\\"), "test2");
    EXPECT_EQ(GetPathBaseName("/test1/test2/test3.ext"), "test3.ext");
}

TEST(GetParentDir, Nominal) 
{
    EXPECT_EQ(GetParentDir(""), "");
    EXPECT_EQ(GetParentDir("test"), "");
    EXPECT_EQ(GetParentDir("/test"), "/");
    EXPECT_EQ(GetParentDir("/"), "");
    EXPECT_EQ(GetParentDir("test/test"), "test");
}

TEST(JoinPaths, Nominal) 
{
    EXPECT_EQ(JoinPaths(""), "");
    EXPECT_EQ(JoinPaths("test"), "test");
    EXPECT_EQ(JoinPaths("/test"), "/test");
    EXPECT_EQ(JoinPaths("test/"), "test/");
    EXPECT_EQ(JoinPaths("/test/"), "/test/");
    EXPECT_EQ(JoinPaths("test1/test2"), "test1/test2");
    EXPECT_EQ(JoinPaths("/test1/test2"), "/test1/test2");
    EXPECT_EQ(JoinPaths("/test1/test2/"), "/test1/test2/");
    EXPECT_EQ(JoinPaths("/test1/test2/"), "/test1/test2/");
    EXPECT_EQ(JoinPaths("\\test1/test2/"), "\\test1/test2/");
    EXPECT_EQ(JoinPaths("\\test1\\test2\\"), "\\test1\\test2\\");
    EXPECT_EQ(JoinPaths("test1", "test2"), "test1\\test2");
    EXPECT_EQ(JoinPaths("/test1", "test2"), "/test1\\test2");
    EXPECT_EQ(JoinPaths("/test1", "/test2"), "/test1/test2");
    EXPECT_EQ(JoinPaths("/test1", "/test2/"), "/test1/test2/");
    EXPECT_EQ(JoinPaths("/test1", "/test2/", "test3.ext"),"/test1/test2/test3.ext");
}

TEST(IsVectorContains, Nominal)
{
    EXPECT_TRUE(IsVectorContains<int>({ 1, 2, 3 }, 1));
    EXPECT_FALSE(IsVectorContains<int>({ 2, 3 }, 1));
}

TEST(IsVectorContainsDuplicateValues, Nominal)
{
    EXPECT_FALSE(IsVectorContainsDuplicateValues<int>({}));
    EXPECT_FALSE(IsVectorContainsDuplicateValues<int>({ 1 }));
    EXPECT_FALSE(IsVectorContainsDuplicateValues<int>({ 1, 2 }));
    EXPECT_FALSE(IsVectorContainsDuplicateValues<int>({ 1, 2, 3 }));
    EXPECT_TRUE(IsVectorContainsDuplicateValues<int>({ 1, 1, 2, 3 }));
    EXPECT_TRUE(IsVectorContainsDuplicateValues<int>({ 1, 1, 2, 2, 3 }));
    EXPECT_TRUE(IsVectorContainsDuplicateValues<int>({ 1, 2, 3, 3 }));
    EXPECT_FALSE(IsVectorContainsDuplicateValues<string>({ "a" }));
    EXPECT_FALSE(IsVectorContainsDuplicateValues<string>({ "a", "b" }));
    EXPECT_TRUE(IsVectorContainsDuplicateValues<string>({ "a", "a" }));
}

TEST(StringToVector, Nominal)
{
    const vector<int> list1 = StringToVector<int>("1, 2, 3 , 4,5,6 ");
    EXPECT_EQ(list1.size(), 6);
    EXPECT_EQ(list1[0], 1);
    EXPECT_EQ(list1[1], 2);
    EXPECT_EQ(list1[2], 3);
    EXPECT_EQ(list1[3], 4);
    EXPECT_EQ(list1[4], 5);
    EXPECT_EQ(list1[5], 6);
    const vector<int> list2 = StringToVector<int>("1; 2; 3 ; 4;5;6 ");
    EXPECT_EQ(list2.size(), 6);
    EXPECT_EQ(list2[0], 1);
    EXPECT_EQ(list2[1], 2);
    EXPECT_EQ(list2[2], 3);
    EXPECT_EQ(list2[3], 4);
    EXPECT_EQ(list2[4], 5);
    EXPECT_EQ(list2[5], 6);
    const vector<int> list3 = StringToVector<int>("1;, 2;; 3 ; 4;5;6 ");
    EXPECT_EQ(list3.size(), 6);
    EXPECT_EQ(list3[0], 1);
    EXPECT_EQ(list3[1], 2);
    EXPECT_EQ(list3[2], 3);
    EXPECT_EQ(list3[3], 4);
    EXPECT_EQ(list3[4], 5);
    EXPECT_EQ(list3[5], 6);
}

TEST(VectorToString, Nominal) 
{
    EXPECT_EQ(VectorToString<int>({}), "");
    EXPECT_EQ(VectorToString<int>({ 1 }), "1");
    EXPECT_EQ(VectorToString<int>({ 1, 2, 3 }), "1, 2, 3");
}

TEST(StringReplace, Nominal) 
{
    EXPECT_EQ(StringReplace("test", "-", ""), "test");
    EXPECT_EQ(StringReplace("test", "t", "a"), "aesa");
    EXPECT_EQ(StringReplace("test", "t", "---"), "---es---");
    EXPECT_EQ(StringReplace("test", "", "a"), "test");
    EXPECT_EQ(StringReplace("test", "", ""), "test");
    EXPECT_EQ(StringReplace("ttt", "ttt", "+++"), "+++");
}

TEST(StringGetAfter, Nominal) 
{
    EXPECT_EQ(StringGetAfter("test", ""), "test");
    EXPECT_EQ(StringGetAfter("test", "notinit"), "");
    EXPECT_EQ(StringGetAfter("test", "e"), "st");
    EXPECT_EQ(StringGetAfter("test, multiple tests", "test"), "s");
    EXPECT_EQ(StringGetAfter("", ""), "");
    EXPECT_EQ(StringGetAfter("path/to/dataset/sub1/image.png", "sub1/"), "image.png");
}

TEST(StringSplit, Nominal) 
{
    const vector<string> list1 = StringSplit("1,2,3,4,5 , 6", ",");
    EXPECT_EQ(list1.size(), 6);
    EXPECT_EQ(list1[0], "1");
    EXPECT_EQ(list1[1], "2");
    EXPECT_EQ(list1[2], "3");
    EXPECT_EQ(list1[3], "4");
    EXPECT_EQ(list1[4], "5 ");
    EXPECT_EQ(list1[5], " 6");
    const vector<string> list2 = StringSplit("1,2,3,4,5 , 6", "");
    EXPECT_EQ(list2.size(), 1);
    EXPECT_EQ(list2[0], "1,2,3,4,5 , 6");
    const vector<string> list3 = StringSplit("1,,2,,3,4,5 , 6", ",");
    EXPECT_EQ(list3.size(), 6);
    EXPECT_EQ(list3[0], "1");
    EXPECT_EQ(list3[1], "2");
    EXPECT_EQ(list3[2], "3");
    EXPECT_EQ(list3[3], "4");
    EXPECT_EQ(list3[4], "5 ");
    EXPECT_EQ(list3[5], " 6");
    const vector<string> list4 = StringSplit("1,,2,,3,4,5 , 6", ",,");
    EXPECT_EQ(list4.size(), 6);
    EXPECT_EQ(list4[0], "1");
    EXPECT_EQ(list4[1], "2");
    EXPECT_EQ(list4[2], "3");
    EXPECT_EQ(list4[3], "4");
    EXPECT_EQ(list4[4], "5 ");
    EXPECT_EQ(list4[5], " 6");
    const vector<string> list5 = StringSplit("1,,2,,3,4,5 , 6", ", ");
    EXPECT_EQ(list5.size(), 6);
    EXPECT_EQ(list5[0], "1");
    EXPECT_EQ(list5[1], "2");
    EXPECT_EQ(list5[2], "3");
    EXPECT_EQ(list5[3], "4");
    EXPECT_EQ(list5[4], "5");
    EXPECT_EQ(list5[5], "6");
    const vector<string> list6 = StringSplit(",1,,2,,3,4,5 , 6 ", ", ");
    EXPECT_EQ(list6.size(), 8);
    EXPECT_EQ(list6[0], "");
    EXPECT_EQ(list6[1], "1");
    EXPECT_EQ(list6[2], "2");
    EXPECT_EQ(list6[3], "3");
    EXPECT_EQ(list6[4], "4");
    EXPECT_EQ(list6[5], "5");
    EXPECT_EQ(list6[6], "6");
    EXPECT_EQ(list6[7], "");
}

TEST(IsStringStartsWith, Nominal) 
{
    EXPECT_FALSE(IsStringStartsWith("", ""));
    EXPECT_FALSE(IsStringStartsWith("a", ""));
    EXPECT_FALSE(IsStringStartsWith("", "a"));
    EXPECT_TRUE(IsStringStartsWith("a", "a"));
    EXPECT_TRUE(IsStringStartsWith("aa", "a"));
    EXPECT_TRUE(IsStringStartsWith("aa", "aa"));
    EXPECT_TRUE(IsStringStartsWith("aaaaaaaaa", "aa"));
}

TEST(IsStringContains, Nominal) {
    EXPECT_TRUE(IsStringContains("", ""));
    EXPECT_TRUE(IsStringContains("a", ""));
    EXPECT_TRUE(IsStringContains("a", "a"));
    EXPECT_TRUE(IsStringContains("ab", "a"));
    EXPECT_TRUE(IsStringContains("ab", "ab"));
    EXPECT_FALSE(IsStringContains("", "a"));
    EXPECT_FALSE(IsStringContains("ab", "c"));
}

#define CHECK_EQUAL_RESULT(find_func1, coeffs1, find_func2, coeffs2) \
{                                                                    \
    Eigen::VectorXd real1;                                           \
    Eigen::VectorXd imag1;                                           \
    const bool success1 = find_func1(coeffs1, real1, imag1);         \
    Eigen::VectorXd real2;                                           \
    Eigen::VectorXd imag2;                                           \
    const bool success2 = find_func2(coeffs2, real2, imag2);         \
    EXPECT_EQ(success1, success2);                                   \
    if (success1)                                                    \
    {                                                                \
        EXPECT_EQ(real1, real2);                                     \
        EXPECT_EQ(imag1, imag2);                                     \
    }                                                                \
}

TEST(EvaluatePolynomial, Nominal) 
{
    EXPECT_EQ(EvaluatePolynomial((Eigen::VectorXd(5) << 1, -3, 3, -5, 10).finished(), 1), 1 - 3 + 3 - 5 + 10);
    EXPECT_NEAR(EvaluatePolynomial((Eigen::VectorXd(4) << 1, -3, 3, -5).finished(), 2.0), 1 * 2 * 2 * 2 - 3 * 2 * 2 + 3 * 2 - 5, 1e-6);
}

TEST(FindLinearPolynomialRoots, Nominal)
{
    Eigen::VectorXd real;
    Eigen::VectorXd imag;
    EXPECT_TRUE(FindLinearPolynomialRoots(Eigen::Vector2d(3, -2), real, imag));
    EXPECT_EQ(real(0), 2.0 / 3.0);
    EXPECT_EQ(imag(0), 0);
    EXPECT_NEAR(EvaluatePolynomial(Eigen::Vector2d(3, -2), complex<double>(real(0), imag(0))).real(), 0.0, 1e-6);
    EXPECT_NEAR(EvaluatePolynomial(Eigen::Vector2d(3, -2), complex<double>(real(0), imag(0))).imag(), 0.0, 1e-6);
    EXPECT_FALSE(FindLinearPolynomialRoots(Eigen::Vector2d(0, 1), real, imag));
}

TEST(FindQuadraticPolynomialRootsReal, Nominal) 
{
    Eigen::VectorXd real;
    Eigen::VectorXd imag;
    Eigen::Vector3d coeffs(3, -2, -4);
    EXPECT_TRUE(FindQuadraticPolynomialRoots(coeffs, real, imag));
    EXPECT_TRUE(real.isApprox(Eigen::Vector2d(-0.868517092, 1.535183758), 1e-6));
    EXPECT_EQ(imag, Eigen::Vector2d(0, 0));
    EXPECT_NEAR(EvaluatePolynomial(coeffs, complex<double>(real(0), imag(0))).real(), 0.0, 1e-6);
    EXPECT_NEAR(EvaluatePolynomial(coeffs, complex<double>(real(1), imag(1))).imag(), 0.0, 1e-6);
}

TEST(FindQuadraticPolynomialRootsComplex, Nominal)
{
    Eigen::VectorXd real;
    Eigen::VectorXd imag;
    const Eigen::Vector3d coeffs(0.276025076998578, 0.679702676853675, 0.655098003973841);
    EXPECT_TRUE(FindQuadraticPolynomialRoots(coeffs, real, imag));
    EXPECT_TRUE(real.isApprox(Eigen::Vector2d(-1.231233560813707, -1.231233560813707), 1e-6));
    EXPECT_TRUE(imag.isApprox(Eigen::Vector2d(0.925954520440279, -0.925954520440279), 1e-6));
    EXPECT_NEAR(EvaluatePolynomial(coeffs, complex<double>(real(0), imag(0))).real(), 0.0, 1e-6);
    EXPECT_NEAR(EvaluatePolynomial(coeffs, complex<double>(real(1), imag(1))).imag(), 0.0, 1e-6);
}

TEST(FindPolynomialRootsDurandKerner, Nominal) 
{
    Eigen::VectorXd real;
    Eigen::VectorXd imag;
    Eigen::VectorXd coeffs(5);
    coeffs << 10, -5, 3, -3, 1;
    EXPECT_TRUE(FindPolynomialRootsDurandKerner(coeffs, real, imag));
    Eigen::VectorXd trueReal(4);
    trueReal << -0.201826, -0.201826, 0.451826, 0.451826;
    EXPECT_TRUE(real.isApprox(trueReal, 1e-6));
    Eigen::VectorXd trueImag(4);
    trueImag << -0.627696, 0.627696, 0.160867, -0.160867;
    EXPECT_TRUE(imag.isApprox(trueImag, 1e-6));
}

TEST(FindPolynomialRootsDurandKernerLinearQuadratic, Nominal)
{
    CHECK_EQUAL_RESULT(FindPolynomialRootsDurandKerner, Eigen::Vector2d(1, 2), FindLinearPolynomialRoots, Eigen::Vector2d(1, 2));
    CHECK_EQUAL_RESULT(FindPolynomialRootsDurandKerner, (Eigen::VectorXd(4) << 0, 0, 1, 2).finished(), FindLinearPolynomialRoots, Eigen::Vector2d(1, 2));
    CHECK_EQUAL_RESULT(FindPolynomialRootsDurandKerner, Eigen::Vector3d(1, 2, 3), FindQuadraticPolynomialRoots, Eigen::Vector3d(1, 2, 3));
    CHECK_EQUAL_RESULT(FindPolynomialRootsDurandKerner, (Eigen::VectorXd(5) << 0, 0, 1, 2, 3).finished(), FindQuadraticPolynomialRoots, Eigen::Vector3d(1, 2, 3));
}

TEST(FindPolynomialRootsCompanionMatrix, Nominal) 
{
    Eigen::VectorXd real;
    Eigen::VectorXd imag;
    Eigen::VectorXd coeffs(5);
    coeffs << 10, -5, 3, -3, 1;
    EXPECT_TRUE(FindPolynomialRootsCompanionMatrix(coeffs, real, imag));
    Eigen::VectorXd trueReal(4);
    trueReal << -0.201826, -0.201826, 0.451826, 0.451826;
    EXPECT_TRUE(real.isApprox(trueReal, 1e-6));
    Eigen::VectorXd trueImag(4);
    trueImag << 0.627696, -0.627696, 0.160867, -0.160867;
    EXPECT_TRUE(imag.isApprox(trueImag, 1e-6));
}

TEST(FindPolynomialRootsCompanionMatrixLinearQuadratic, Nominal) 
{
    CHECK_EQUAL_RESULT(FindPolynomialRootsCompanionMatrix, Eigen::Vector2d(1, 2), FindLinearPolynomialRoots, Eigen::Vector2d(1, 2));
    CHECK_EQUAL_RESULT(FindPolynomialRootsCompanionMatrix, (Eigen::VectorXd(4) << 0, 0, 1, 2).finished(), FindLinearPolynomialRoots, Eigen::Vector2d(1, 2));
    CHECK_EQUAL_RESULT(FindPolynomialRootsCompanionMatrix, Eigen::Vector3d(1, 2, 3), FindQuadraticPolynomialRoots, Eigen::Vector3d(1, 2, 3));
    CHECK_EQUAL_RESULT(FindPolynomialRootsCompanionMatrix, (Eigen::VectorXd(5) << 0, 0, 1, 2, 3).finished(), FindQuadraticPolynomialRoots, Eigen::Vector3d(1, 2, 3));
}

TEST(FindPolynomialRootsCompanionMatrixZeroSolution, Nominal) 
{
    Eigen::VectorXd real;
    Eigen::VectorXd imag;
    Eigen::VectorXd coeffs(5);
    coeffs << 10, -5, 3, -3, 0;
    EXPECT_TRUE(FindPolynomialRootsCompanionMatrix(coeffs, real, imag));
    Eigen::VectorXd trueReal(4);
    trueReal << 0.692438, -0.0962191, -0.0962191, 0;
    EXPECT_TRUE(real.isApprox(trueReal, 1e-6));
    Eigen::VectorXd trueImag(4);
    trueImag << 0, 0.651148, -0.651148, 0;
    EXPECT_TRUE(imag.isApprox(trueImag, 1e-6));
}

TEST(PRNGSeed, Nominal) 
{
    EXPECT_TRUE(PRNG == nullptr);
    SetPRNGSeed();
    EXPECT_TRUE(PRNG != nullptr);
    SetPRNGSeed(0);
    EXPECT_TRUE(PRNG != nullptr);
    thread thread([]() 
        {
        EXPECT_TRUE(PRNG == nullptr);
        SetPRNGSeed();
        EXPECT_TRUE(PRNG != nullptr);
        SetPRNGSeed(0);
        EXPECT_TRUE(PRNG != nullptr);
        });
    thread.join();
}

TEST(Repeatability, Nominal) 
{
    SetPRNGSeed(0);
    vector<int> numbers1;
    for (size_t i = 0; i < 100; i++) 
    {
        numbers1.push_back(GetRandomUniformInteger(0, 10000));
    }
    SetPRNGSeed(1);
    vector<int> numbers2;
    for (size_t i = 0; i < 100; i++)
    {
        numbers2.push_back(GetRandomUniformInteger(0, 10000));
    }
    SetPRNGSeed(0);
    vector<int> numbers3;
    for (size_t i = 0; i < 100; i++)
    {
        numbers3.push_back(GetRandomUniformInteger(0, 10000));
    }
    EXPECT_EQ(numbers1, numbers3);
    bool isAllEqual = true;
    for (size_t i = 0; i < numbers1.size(); i++)
    {
        if (numbers1[i] != numbers2[i])
        {
            isAllEqual = false;
        }
    }
    EXPECT_FALSE(isAllEqual);
}

TEST(GetRandomUniformInteger, Nominal)
{
    SetPRNGSeed();
    for (size_t i = 0; i < 1000; i++)
    {
        EXPECT_GE(GetRandomUniformInteger(-100, 100), -100);
        EXPECT_LE(GetRandomUniformInteger(-100, 100), 100);
    }
}

TEST(GetRandomUniformReal, Nominal)
{
    SetPRNGSeed();
    for (size_t i = 0; i < 1000; i++) 
    {
        EXPECT_GE(GetRandomUniformReal(-100.0, 100.0), -100.0);
        EXPECT_LE(GetRandomUniformReal(-100.0, 100.0), 100.0);
    }
}

TEST(GetRandomGaussian, Nominal)
{
    SetPRNGSeed(0);
    const double mean = 1.0;
    const double sigma = 1.0;
    const size_t numValues = 100000;
    vector<double> values;
    for (size_t i = 0; i < numValues; i++)
    {
        values.push_back(GetRandomGaussian(mean, sigma));
    }
    EXPECT_LE(abs(GetMean(values) - mean), 1e-3);
    EXPECT_LE(abs(GetStdDev(values) - sigma), 1e-3);
}

TEST(ShuffleNone, Nominal)
{
    SetPRNGSeed();
    vector<int> numbers(0);
    Shuffle(0, numbers);
    numbers = { 1, 2, 3, 4, 5 };
    vector<int> shuffledNumbers = numbers;
    Shuffle(0, shuffledNumbers);
    EXPECT_EQ(numbers, shuffledNumbers);
}

TEST(ShuffleAll, Nominal)
{
    SetPRNGSeed(0);
    vector<int> numbers(1000);
    iota(numbers.begin(), numbers.end(), 0);
    vector<int> shuffledNumbers = numbers;
    Shuffle(1000, shuffledNumbers);
    size_t numShuffled = 0;
    for (size_t i = 0; i < numbers.size(); i++) 
    {
        if (numbers[i] != shuffledNumbers[i]) 
        {
            numShuffled++;
        }
    }
    EXPECT_GT(numShuffled, 0);
}

TEST(FeatureMatchHashing, Nominal)
{
    unordered_set<pair<size_t, size_t>, MatchPairHash, MatchPairEqual> set;
    set.emplace(1520, 20892);
    EXPECT_EQ(set.size(), 1);
    set.emplace(1520, 20892);
    EXPECT_EQ(set.size(), 1);
    EXPECT_EQ(set.count(make_pair(0, 0)), 0);
    EXPECT_EQ(set.count(make_pair(1520, 20892)), 1);
    EXPECT_EQ(set.count(make_pair(1520, 20895)), 0);
    set.emplace(20892, 1519);
    EXPECT_EQ(set.size(), 2);
    EXPECT_EQ(set.count(make_pair(20892, 1518)), 0);
    EXPECT_EQ(set.count(make_pair(1519, 20892)), 1);
    EXPECT_EQ(set.count(make_pair(20892, 1520)), 1);
}

TEST(SignOfNumber, Nominal) 
{
    EXPECT_EQ(SignOfNumber(0), 0);
    EXPECT_EQ(SignOfNumber(-0.1), -1);
    EXPECT_EQ(SignOfNumber(0.1), 1);
    EXPECT_EQ(SignOfNumber(numeric_limits<float>::quiet_NaN()), 0);
    EXPECT_EQ(SignOfNumber(numeric_limits<float>::infinity()), 1);
    EXPECT_EQ(SignOfNumber(-numeric_limits<float>::infinity()), -1);
}

TEST(Clamp, Nominal)
{
    EXPECT_EQ(Clamp(0, -1, 1), 0);
    EXPECT_EQ(Clamp(0, 0, 1), 0);
    EXPECT_EQ(Clamp(0, -1, 0), 0);
    EXPECT_EQ(Clamp(0, -1, 1), 0);
    EXPECT_EQ(Clamp(0, 1, 2), 1);
    EXPECT_EQ(Clamp(0, -2, -1), -1);
    EXPECT_EQ(Clamp(0, 0, 0), 0);
}

TEST(DegToRad, Nominal)
{
    EXPECT_EQ(DegToRad(0.0f), 0.0f);
    EXPECT_EQ(DegToRad(0.0), 0.0);
    EXPECT_LT(abs(DegToRad(180.0f) - M_PI), 1e-6f);
    EXPECT_LT(abs(DegToRad(180.0) - M_PI), 1e-6);
}

TEST(RadToDeg, Nominal)
{
    EXPECT_EQ(RadToDeg(0.0f), 0.0f);
    EXPECT_EQ(RadToDeg(0.0), 0.0);
    EXPECT_LT(abs(RadToDeg(M_PI) - 180.0f), 1e-6f);
    EXPECT_LT(abs(RadToDeg(M_PI) - 180.0), 1e-6);
}

TEST(GetMedian, Nominal)
{
    EXPECT_EQ(GetMedian<int>({ 1, 2, 3, 4 }), 2.5);
    EXPECT_EQ(GetMedian<int>({ 1, 2, 3, 100 }), 2.5);
    EXPECT_EQ(GetMedian<int>({ 1, 2, 3, 4, 100 }), 3);
    EXPECT_EQ(GetMedian<int>({ -100, 1, 2, 3, 4 }), 2);
    EXPECT_EQ(GetMedian<int>({ -1, -2, -3, -4 }), -2.5);
    EXPECT_EQ(GetMedian<int>({ -1, -2, 3, 4 }), 1);
    // Test integer overflow scenario.
    EXPECT_EQ(GetMedian<int8_t>({ 100, 115, 119, 127 }), 117);
}

TEST(Percentile, Nominal)
{
    EXPECT_EQ(Percentile<int>({ 0 }, 0), 0);
    EXPECT_EQ(Percentile<int>({ 0 }, 50), 0);
    EXPECT_EQ(Percentile<int>({ 0 }, 100), 0);
    EXPECT_EQ(Percentile<int>({ 0, 1 }, 0), 0);
    EXPECT_EQ(Percentile<int>({ 0, 1 }, 50), 1);
    EXPECT_EQ(Percentile<int>({ 0, 1 }, 100), 1);
    EXPECT_EQ(Percentile<int>({ 0, 1, 2 }, 0), 0);
    EXPECT_EQ(Percentile<int>({ 0, 1, 2 }, 50), 1);
    EXPECT_EQ(Percentile<int>({ 0, 1, 2 }, 100), 2);
    EXPECT_EQ(Percentile<int>({ 0, 1, 1, 2 }, 0), 0);
    EXPECT_EQ(Percentile<int>({ 0, 1, 1, 2 }, 33), 1);
    EXPECT_EQ(Percentile<int>({ 0, 1, 1, 2 }, 50), 1);
    EXPECT_EQ(Percentile<int>({ 0, 1, 1, 2 }, 66), 1);
    EXPECT_EQ(Percentile<int>({ 0, 1, 1, 2 }, 100), 2);
}

TEST(Mean, Nominal)
{
    EXPECT_EQ(GetMean<int>({ 1, 2, 3, 4 }), 2.5);
    EXPECT_EQ(GetMean<int>({ 1, 2, 3, 100 }), 26.5);
    EXPECT_EQ(GetMean<int>({ 1, 2, 3, 4, 100 }), 22);
    EXPECT_EQ(GetMean<int>({ -100, 1, 2, 3, 4 }), -18);
    EXPECT_EQ(GetMean<int>({ -1, -2, -3, -4 }), -2.5);
    EXPECT_EQ(GetMean<int>({ -1, -2, 3, 4 }), 1);
}

TEST(Variance, Nominal)
{
    EXPECT_LE(abs(GetVariance<int>({ 1, 2, 3, 4 }) - 1.66666666), 1e-6);
    EXPECT_LE(abs(GetVariance<int>({ 1, 2, 3, 100 }) - 2401.66666666), 1e-6);
    EXPECT_LE(abs(GetVariance<int>({ 1, 2, 3, 4, 100 }) - 1902.5), 1e-6);
    EXPECT_LE(abs(GetVariance<int>({ -100, 1, 2, 3, 4 }) - 2102.5), 1e-6);
    EXPECT_LE(abs(GetVariance<int>({ -1, -2, -3, -4 }) - 1.66666666), 1e-6);
    EXPECT_LE(abs(GetVariance<int>({ -1, -2, 3, 4 }) - 8.66666666), 1e-6);
}

TEST(StdDev, Nominal)
{
    EXPECT_LE(abs(sqrt(GetVariance<int>({ 1, 2, 3, 4 })) - GetStdDev<int>({ 1, 2, 3, 4 })), 1e-6);
    EXPECT_LE(abs(sqrt(GetVariance<int>({ 1, 2, 3, 100 })) - GetStdDev<int>({ 1, 2, 3, 100 })), 1e-6);
}

TEST(NextCombination, Nominal)
{
    vector<int> list{ 0 };
    EXPECT_FALSE(NextCombination(list.begin(), list.begin() + 1, list.end()));
    list = { 0, 1 };
    EXPECT_FALSE(NextCombination(list.begin(), list.begin() + 2, list.end()));
    EXPECT_EQ(list[0], 0);
    EXPECT_TRUE(NextCombination(list.begin(), list.begin() + 1, list.end()));
    EXPECT_EQ(list[0], 1);
    EXPECT_FALSE(NextCombination(list.begin(), list.begin() + 1, list.end()));
    EXPECT_EQ(list[0], 0);
    list = { 0, 1, 2 };
    EXPECT_EQ(list[0], 0);
    EXPECT_EQ(list[1], 1);
    EXPECT_EQ(list[2], 2);
    EXPECT_TRUE(NextCombination(list.begin(), list.begin() + 2, list.end()));
    EXPECT_EQ(list[0], 0);
    EXPECT_EQ(list[1], 2);
    EXPECT_EQ(list[2], 1);
    EXPECT_TRUE(NextCombination(list.begin(), list.begin() + 2, list.end()));
    EXPECT_EQ(list[0], 1);
    EXPECT_EQ(list[1], 2);
    EXPECT_EQ(list[2], 0);
    EXPECT_FALSE(NextCombination(list.begin(), list.begin() + 2, list.end()));
    EXPECT_EQ(list[0], 0);
    EXPECT_EQ(list[1], 1);
    EXPECT_EQ(list[2], 2);
}

TEST(NChooseK, Nominal)
{
    EXPECT_EQ(NChooseK(1, 0), 1);
    EXPECT_EQ(NChooseK(2, 0), 1);
    EXPECT_EQ(NChooseK(3, 0), 1);

    EXPECT_EQ(NChooseK(1, 1), 1);
    EXPECT_EQ(NChooseK(2, 1), 2);
    EXPECT_EQ(NChooseK(3, 1), 3);

    EXPECT_EQ(NChooseK(2, 2), 1);
    EXPECT_EQ(NChooseK(2, 3), 0);

    EXPECT_EQ(NChooseK(3, 2), 3);
    EXPECT_EQ(NChooseK(4, 2), 6);
    EXPECT_EQ(NChooseK(5, 2), 10);

    EXPECT_EQ(NChooseK(500, 3), 20708500);
}

TEST(TruncateCast, Nominal)
{
    EXPECT_EQ((TruncateCast<int, int8_t>(-129)), -128);
    EXPECT_EQ((TruncateCast<int, int8_t>(128)), 127);
    EXPECT_EQ((TruncateCast<int, uint8_t>(-1)), 0);
    EXPECT_EQ((TruncateCast<int, uint8_t>(256)), 255);
    EXPECT_EQ((TruncateCast<int, uint16_t>(-1)), 0);
    EXPECT_EQ((TruncateCast<int, uint16_t>(65536)), 65535);
}

TEST(DecomposeMatrixRQ, Nominal)
{
    for (int i = 0; i < 100000; i++)
    {
        const Eigen::Matrix4d A = Eigen::Matrix4d::Random();

        Eigen::Matrix4d R, Q;
        DecomposeMatrixRQ(A, R, Q);

        EXPECT_TRUE(R.bottomRows(4).isUpperTriangular());
        EXPECT_TRUE(Q.isUnitary());
        EXPECT_NEAR(Q.determinant(), 1.0, 1e-6);
        EXPECT_TRUE(A.isApprox(R * Q, 1e-6));
    }
}

int main(int argc, char** argv)
{

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}