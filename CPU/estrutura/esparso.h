#pragma once
#include <vector>
#include "sistema.h"

//usado apenas pela versão GPU

struct h_sparse {
	std::vector<float_type> spJval;
	/* int* */ std::vector<subMatJ> spJsubmatType;
	/* int* */ std::vector<int> cooColIndSubMatJ;
	/* int* */ std::vector<int> cooRowIndSubMatJ;
	/* int* */ std::vector<int> cooColIndJ;
	/* int* */ std::vector<int> cooRowIndJ;
	std::vector<int> csrRowPtrJ;
	int nnzJ;

	std::vector<int> Hpos;
	std::vector<int> Lpos;
	std::vector<int> Mpos;
	std::vector<int> Npos;
};

struct d_sparse {
	float_type* spJval;
	subMatJ* spJsubmatType;
	int* cooColIndSubMatJ;
	int* cooRowIndSubMatJ;
	int* cooColIndJ;
	int* cooRowIndJ;
	int* csrRowPtrJ;
	int nnzJ;

	int* Hpos;
	int* Lpos;
	int* Mpos;
	int* Npos;
	int nnzH;
	int nnzL;
	int nnzM;
	int nnzN;
};