#include "../global.h"

#include "config.h"

void loadFile() {
	std::string str;
	std::ifstream t("data.txt");

	if (t) {
		std::ostringstream ss;
		ss << t.rdbuf(); 
		str = ss.str();
	}
	else {
		printf("ERRO: A ABERTURA DO ARQUIVO FALHOU!");
		std::cin.get();
		return;
	}

	std::string aux = "laconic_mode";
	int inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	std::string laconic_mode_aux = aux;
	if (aux == "OFF") {
		global::laconic_mode = 0;
	}
	else {
		global::laconic_mode = 1;
	}

	aux = "metodo";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	if (!global::laconic_mode) { std::cout << "metodo         = " << aux << std::endl; }
	if (aux == "denso") {
		global::metodo = denso;
	}
	else if (aux == "denso_LAPACKE") {
		global::metodo = denso_LAPACKE;
	}
	else if (aux == "esparso") {
		global::metodo = esparso;
	}
	else if (aux == "esparsoSimples") {
		global::metodo = esparsoSimples;
	}

	else {
		global::metodo = nda;
	}

	aux = "lim_inj_reat";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	if (!global::laconic_mode && str.find("lim_inj_reat") != std::string::npos) { std::cout << "lim_inj_reat   = " << aux << std::endl; }
	if (aux == "OFF") {
		global::lim_inj_reat = 0;
	}
	else if (aux == "ON") {
		global::lim_inj_reat = 1;
	}

	aux = "openmp";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	if (!global::laconic_mode) { std::cout << "openmp         = " << aux << std::endl; }
	if (aux == "OFF") {
		global::openmp = 0;
	}
	else {
		global::openmp = 1;
	}

	aux = "v_inicial";
	inicio = str.find('=', str.find(aux) + aux.size()) + 1, 
		fim = str.find(';', inicio);
	global::v_inicial = atof(str.substr(inicio, fim - inicio).c_str());

	if (!global::laconic_mode) { printf("v_inicial      = %f\n", global::v_inicial); }

	aux = "theta_inicial";
	inicio = str.find('=', str.find(aux) + aux.size()) + 1, 
		fim = str.find(';', inicio);
	global::theta_inicial = atof(str.substr(inicio, fim - inicio).c_str());

	if (!global::laconic_mode) { printf("theta_inicial  = %f\n", global::theta_inicial); }

	aux = "tolerancia";
	inicio = str.find('=', str.find(aux) + aux.size()) + 1, 
		fim = str.find(';', inicio);
	global::tol = atof(str.substr(inicio, fim - inicio).c_str());
	if (!global::laconic_mode) { std::cout << "tolerancia     = " << global::tol << std::endl; }

	aux = "no_max_iter";
	inicio = str.find('=', str.find(aux) + aux.size()) + 1, 
		fim = str.find(';', inicio);
	global::no_max_iter = std::stoi(str.substr(inicio, fim - inicio));
	if (!global::laconic_mode) { std::cout << "no_max_iter    = " << global::no_max_iter << std::endl; }

	aux = "arq_entrada";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);

	global::arq_entrada = "sistemas/";
	global::arq_entrada.append(str.substr(inicio, fim - inicio));
	if (!global::laconic_mode) { std::cout << "arq_entrada    = " << str.substr(inicio, fim - inicio) << std::endl; }

	aux = "verbose_mode";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	if (!global::laconic_mode && str.find("verbose_mode") != std::string::npos) { std::cout << "verbose_mode   = " << aux << std::endl; }
	if (aux == "OFF") {
		global::verbose_mode = 0;
	}
	else if (aux == "ON") {
		global::verbose_mode = 1;
	}

	if (!global::laconic_mode) { std::cout << "laconic_mode   = " << laconic_mode_aux << std::endl; }

	aux = "output_bench";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	if (!global::laconic_mode) { std::cout << "output_bench   = " << aux << std::endl; }
	if (aux == "ALL") {
		global::output_benchmark = output_benchmarkType::all;
	}
	else if (aux == "FILE") {
		global::output_benchmark = output_benchmarkType::file;
	}
	else {
		global::output_benchmark = output_benchmarkType::screen;
	}

	aux = "output_ans";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	if (!global::laconic_mode) { std::cout << "output_ans     = " << aux << std::endl; }
	if (aux == "OFF") {
		global::output_ans = 0;
	}
	else {
		global::output_ans = 1;
	}

	aux = "output_iter";
	inicio = str.find('=', str.find(aux) + aux.size()) + 2, 
		fim = str.find(';', inicio);
	aux = str.substr(inicio, fim - inicio);
	if (!global::laconic_mode && str.find("output_iter") != std::string::npos) { std::cout << "output_iter    = " << aux << std::endl; }
	if (aux == "OFF") {
		global::output_processo_iterativo = 0;
	}
	else if (aux == "ON") {
		global::output_processo_iterativo = 1;
	}

	t.close();
}