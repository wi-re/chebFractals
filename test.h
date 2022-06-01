#pragma once
#include <cassert>

#define assertTrue(x) {printf("assertTrue(%s)\n", #x); assert(x);}
#define assertFalse(x){ printf("assertFalse(%s)\n", #x); assert(!(x));}
#define assertLessEqual(x,comp) {printf("assertLessEqual(%s, %s)\n", #x, #comp); assert((x) <= (comp));}
#define assertGreaterEqual(x,comp) {printf("assertGreaterEqual(%s, %s)\n", #x, #comp);assert((x) >= (comp));}
#define assertEqual(x,comp) {printf("assertEqual(%s, %s)\n", #x, #comp);assert((x) == (comp));}
#define assertNotEqual(x,comp) {printf("assertNotEqual(%s, %s)\n", #x, #comp);assert(!((x) == (comp)));}
#define assertThrows(c)  {printf("assertThrows(%s)\n", #c); try{c; printf("%s\t",#c);assert(1==2);}catch(...){}}

void testInterval();
void testDomain();
void testChebTech();
void testChebTech_functions();
void testArithmetic();
void testRoots();
void testUtilities();

void testFunction_Construction();
void testFunction_Properties();
void testFunction_ClassUsage();
void testFunction_Evaluation();
void testFunction_Calculus(); // DOT MISSING
void testFunction_Roots();


void testFunction_Arithmetic();
void testFunction();
