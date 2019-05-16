
#include "Nucleus.h"

Nucleus::Nucleus(){}

void Nucleus::OnTestStart(const ::testing::TestInfo& info) {

  // Set up Rank info string
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int nDigits = floor(log10(size) + 1);
  std::string rankString = std::to_string(rank);
  rankString.insert(rankString.begin(), nDigits-rankString.length(), ' ');
  rankString = "[Rank " + rankString + " ]";

  // Setup start info
  std::string message = "\e[32m";
  message += rankString + "[ RUN      ] \e[0m";
  message += info.test_case_name();
  message += ".";
  message += info.name();
  message += "\n";

/*
  if (rank == 0){
    std::cout << message.c_str() << std::flush;
    for (unsigned int proc = 1; proc < size; proc++){
      MPI_Recv(&message[0], message.length(), MPI_CHAR, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::cout << message.c_str() << std::flush;
    }
  } else {
    MPI_Send(message.data(), message.length(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }
*/
  this->OutputMessage(message, rank, size);
}

void Nucleus::OutputMessage(std::string message, const unsigned int rank, const unsigned int size) const {
  if (rank == 0){
    std::cout << message.c_str() << std::flush;
    for (unsigned int proc = 1; proc < size; proc++){
      MPI_Recv(&message[0], message.length(), MPI_CHAR, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::cout << message.c_str() << std::flush;
    }
  } else {
    MPI_Send(message.data(), message.length(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }
}

std::string Nucleus::PrintFullTestCommentIfPresent(const ::testing::TestInfo& test_info) {
  const char* const type_param = test_info.type_param();
  const char* const value_param = test_info.value_param();

  std::string result;
  if (type_param != NULL || value_param != NULL) {
    result += ", where ";
    // printf(", where ");
    if (type_param != NULL) {
      result += "TypeParam = ";
      result += type_param;
      // printf("%s = %s", "TypeParam", type_param);
      if (value_param != NULL)
        result += " and ";
        // printf(" and ");
    }
    if (value_param != NULL) {
      result += "GetParam() = ";
      result += value_param;
      // printf("%s = %s", "GetParam()", value_param);
    }
  }

  return result;
}

void Nucleus::OnTestEnd(const ::testing::TestInfo& info) {
  // Set up Rank info string
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int nDigits = floor(log10(size) + 1);
  std::string rankString = std::to_string(rank);
  rankString.insert(rankString.begin(), nDigits-rankString.length(), ' ');
  rankString = "[Rank " + rankString + " ]";

  std::string message;
  if (info.result()->Passed()){
    message = "\e[32m" + rankString + "[       OK ] \e[0m";
  } else {
    message = "\e[31m" + rankString + "[  FAILED  ] \e[0m";
  }
  message += info.test_case_name();
  message += ".";
  message += info.name();

  if (info.result()->Failed())
    message += this->PrintFullTestCommentIfPresent(info);

  message += " (" + std::to_string(info.result()->elapsed_time()) + " ms)\n";

  this->OutputMessage(message, rank, size);
}

void Nucleus::OnTestPartResult(const ::testing::TestPartResult& result) {
  if (result.type() == ::testing::TestPartResult::kSuccess)
    return;
  std::cout << "Part Result" << std::endl;

}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  delete listeners.Release(listeners.default_result_printer());
  listeners.Append(new Nucleus);

  int result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
