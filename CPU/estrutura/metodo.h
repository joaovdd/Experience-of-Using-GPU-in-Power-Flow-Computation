#pragma once
enum metodoType {
	esparso,         
	esparsoSimples,  
	denso,           
	denso_LAPACKE,   
	nda,
	hibridoA,		 
	hibridoB		 
};

enum output_benchmarkType {
	all,
	screen,
	file
};