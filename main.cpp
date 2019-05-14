#include "LinearAlgebra/LinearAlgebra.h"

int main(int argc, char* argv[]) {

  LinearAlgebra init(&argc, &argv);
  std::cout << "Hello from process " << init.rank() << " of " << init.size() << std::endl;

  return 0;
}
