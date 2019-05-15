
#include "LinearAlgebra/Vector.h"
#include "math.h"

int main(int argc, char* argv[]) {

  LinearAlgebra init(argc, argv);
  
  Vector v(15, init);
  Vector w(15, init);
  w.scale(0.01);

  v.add(w);
  v.print();

  double norm = v.length();

  if (init.rank() == 0){
    linalg_assert( fabs(norm - 32.1776552906) < 1.0e-6 )
    std::string msg;
    msg += "[Rank " + std::to_string(init.rank()) + "] ";
    msg += std::to_string(norm);
    msg += "  Error: " + std::to_string( fabs(norm - 32.1776552906) );
    info(msg.c_str())
  }

  return 0;
}
