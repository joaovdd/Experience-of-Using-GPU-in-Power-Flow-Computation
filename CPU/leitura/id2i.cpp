#include "id2i.h"

// retorna numero da barra com identificador id
int id2i(int id, sistemaType& sistema, barraType& barra) {
	// TODO: pesquisar por id em barra.id[] eretornar posição i
	for (int i = 0; i < sistema.nB; i++) {
		if (barra.id[i] == id) {
			return i + 1;
		}
	}
	throw "Barra inexistente!";
}