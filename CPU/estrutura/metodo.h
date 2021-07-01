#pragma once
enum metodoType {
	esparso,         // CPU e GPU
	esparsoSimples,  // apenas CPU
	denso,           // CPU e GPU
	denso_LAPACKE,   // apenas CPU
	nda,
	hibridoA,		 // GPU
	hibridoB		 // GPU
};

enum output_benchmarkType {
	all,
	screen,
	file
};