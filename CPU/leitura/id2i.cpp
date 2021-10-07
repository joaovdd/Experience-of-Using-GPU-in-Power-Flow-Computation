#include "id2i.h"

int id2i(int id, sistemaType& sistema, barraType& barra) {
	for (int i = 0; i < sistema.nB; i++) {
		if (barra.id[i] == id) {
			return i + 1;
		}
	}
	throw "Barra inexistente!";
}