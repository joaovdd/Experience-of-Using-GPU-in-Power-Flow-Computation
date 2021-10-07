
#define DOUBLE_MODE true

#if DOUBLE_MODE
typedef double float_type;
#else
typedef float float_type;
#endif

#if DOUBLE_MODE
typedef cuDoubleComplex complex_type;
#else
typedef cuFloatComplex complex_type;
#endif

#if DOUBLE_MODE
#define MK_COMPLEX_FUNCTION(x, y) make_cuDoubleComplex(x, y) 
#else
#define MK_COMPLEX_FUNCTION(x, y) make_cuFloatComplex(x, y)
#endif
