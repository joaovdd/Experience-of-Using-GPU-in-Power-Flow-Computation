#include "benchmarks.h"

const char* benchmarkTypeString[] = {
	"processoIterativo",
	"jacobiano",
	"jacobianoStencil",
	"calcPQ"
};

namespace global {
	trackerType tracker;
}