#include "graph/Graph.h"
#include "community/DynComDetector.h"

int main() {
    DynComDetector detector("../data/test_snapshots/graph_snapshots/", 9);
    detector.initialize(0.2);
    detector.perform_dynamic_detection();
    detector.evaluate_modularity();
    detector.write_result("../data/test_snapshots/communities/");
    return 0;
}
