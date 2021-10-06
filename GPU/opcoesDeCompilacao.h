// defina como verdadeiro para que seja usado double,
// caso contr�rio ser� utilizado float 
// n�o esquecer de fazer mudan�a no outro arquivo tamb�m!!!!!!!!!!
#define DOUBLE_MODE true

// *****************************************************************************************
// N�o alterar a partir daqui***************************************************************
// *****************************************************************************************

// defini��o do tipo de ponto flutuante
#if DOUBLE_MODE
typedef double float_type;
#else
typedef float float_type;
#endif

// defini��o do tipo de n�mero complexo
#if DOUBLE_MODE
typedef cuDoubleComplex complex_type;
#else
typedef cuFloatComplex complex_type;
#endif

// defini��o do constritor complexo
#if DOUBLE_MODE
#define MK_COMPLEX_FUNCTION(x, y) make_cuDoubleComplex(x, y) 
#else
#define MK_COMPLEX_FUNCTION(x, y) make_cuFloatComplex(x, y)
#endif
