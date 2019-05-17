
#include "Nucleus.h"

Nucleus::Nucleus(){}

std::tuple<std::string, unsigned int, unsigned int> Nucleus::getMPIprefix(const color& c, const std::string& step, const align& alignment) const {

  // Figure out coloring
  std::string colorStart, colorEnd = "";
  switch (c){
    case color::GREEN:
      colorStart = "\e[32m";
      colorEnd = "\e[0m";
      break;
    case color::RED:
      colorStart = "\e[31m";
      colorEnd = "\e[0m";
      break;
    case color::YELLOW:
      colorStart = "\e[33m";
      colorEnd = "\e[0m";
      break;
    default:
      break;
  }

  // Set up Rank info string
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int nDigits = floor(log10(size) + 1);
  std::string procString = std::to_string(size);
  procString.insert(procString.begin(), nDigits-procString.length(), ' ');
  procString = "[ " + procString + " ";
  if (size == 1)
    procString += "Proc";
  else
    procString += "Procs";
  procString += " ]";

  // Format step info
  std::string stepString;
  if (step.length() == 1){
    stepString.reserve(10);
    for (unsigned int i = 0; i < 10; i++)
      stepString += step;
    stepString = "[" + stepString + "] ";
  } else {
    // stepString = std::toupper(step);
    stepString = step;
    std::transform(stepString.begin(), stepString.end(), stepString.begin(), ::toupper);
    unsigned int length = stepString.length();
    switch (alignment){
      case align::LEFT:
        stepString.append(8-length, ' ');
        break;
      case align::RIGHT:
        stepString.insert(stepString.begin(), 8-length, ' ');
        break;
      case align::CENTER:
        stepString.insert(stepString.begin(), (8-length)/2, ' ');
        stepString.append((8-length)/2, ' ');
        break;
      default:
        break;
    }
    stepString = "[ " + stepString + " ] ";
  }

  std::string message;
  message = colorStart + procString + stepString + colorEnd;

  return std::make_tuple(message, rank, size);

}

void Nucleus::OnTestStart(const ::testing::TestInfo& info) {

  auto prefix = this->getMPIprefix(color::GREEN, "run", align::LEFT);

  std::string message = std::get<0>(prefix);
  message += info.test_case_name();
  message += ".";
  message += info.name();
  message += "\n";

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));
}

void Nucleus::OutputMessage(std::string message,
                            const unsigned int rank,
                            const unsigned int size,
                            const bool outputAllProcs) const {
  if (rank == 0){
    std::cout << message.c_str() << std::flush;
    if (outputAllProcs)
      for (unsigned int proc = 1; proc < size; proc++){
        MPI_Recv(&message[0], message.length(), MPI_CHAR, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << message.c_str() << std::flush;
      }
  } else {
    if (outputAllProcs)
      MPI_Send(message.data(), message.length(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }
}

std::string Nucleus::PrintFullTestCommentIfPresent(const ::testing::TestInfo& test_info) {
  const char* const type_param = test_info.type_param();
  const char* const value_param = test_info.value_param();

  std::string result;
  if (type_param != NULL || value_param != NULL) {
    result += ", where ";
    if (type_param != NULL) {
      result += "TypeParam = ";
      result += type_param;
      if (value_param != NULL)
        result += " and ";
    }
    if (value_param != NULL) {
      result += "GetParam() = ";
      result += value_param;
    }
  }

  return result;
}

void Nucleus::OnTestEnd(const ::testing::TestInfo& info) {

  std::tuple<std::string, unsigned int, unsigned int> prefix;
  if (info.result()->Passed())
    prefix = this->getMPIprefix(color::GREEN, "ok", align::RIGHT);
  else
    prefix = this->getMPIprefix(color::RED, "failed", align::CENTER);

  std::string message = std::get<0>(prefix);
  message += info.test_case_name();
  message += ".";
  message += info.name();

  if (info.result()->Failed())
    message += this->PrintFullTestCommentIfPresent(info);

  message += " (" + std::to_string(info.result()->elapsed_time()) + " ms)\n";

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));
}

void Nucleus::OnTestPartResult(const ::testing::TestPartResult& result) {
  if (result.type() == ::testing::TestPartResult::kSuccess)
    return;


  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int nDigits = floor(log10(size) + 1);
  std::string rankString = std::to_string(rank);
  rankString.insert(rankString.begin(), nDigits-rankString.length(), ' ');
  rankString = "[ Rank " + rankString + " ]";

  std::string message;
  message += "\e[31m" + rankString + "\e[0m\n";

  std::string file_name(result.file_name() == NULL ? "unknown file" : result.file_name());
  message += file_name + ":";
  if (result.line_number() > 0)
    message += std::to_string(result.line_number()) + ":";
  message += " ";
  if (result.type() == ::testing::TestPartResult::kFatalFailure)
    message += "Failure\n";
  else if (result.type() == ::testing::TestPartResult::kSuccess)
    message += "Success";
  else
    message += "Unknown result type";
  
  message += result.message();
  message += "\n";

  this->OutputMessage(message, rank, size, true);

}

void Nucleus::OnTestCaseStart(const ::testing::TestCase& test_case) {
  auto prefix = this->getMPIprefix(color::GREEN, "-", align::FULL);
  std::string message = std::get<0>(prefix);
  message += formatTestCount(test_case.test_to_run_count()) + " from ";
  message += test_case.name();
  if (test_case.type_param() != NULL){
    message += ", where TypeParam = ";
    message += test_case.type_param();
  }
  message += "\n";

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));

}

void Nucleus::OnTestCaseEnd(const ::testing::TestCase& test_case) {
  auto prefix = this->getMPIprefix(color::GREEN, "-", align::FULL);
  std::string message = std::get<0>(prefix);
  message += formatTestCount(test_case.test_to_run_count());
  message += " from ";
  message += test_case.name();
  message += " (" + std::to_string(test_case.elapsed_time()) + " ms total)\n\n";

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));
}

void Nucleus::OnTestIterationStart(const ::testing::UnitTest& unit_test, int iteration){
  auto prefix = this->getMPIprefix(color::GREEN, "=", align::FULL);

  std::string message = std::get<0>(prefix);
  message += "Running ";
  message += formatTestCount(unit_test.test_to_run_count());
  message += " from ";
  message += formatTestCaseCount(unit_test.test_case_to_run_count());
  message += ".\n";

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));
}

std::string Nucleus::PrintFailedTests(const ::testing::UnitTest& unit_test) {
  std::string result;
  if (unit_test.failed_test_count() == 0)
    return result;

  auto prefix = this->getMPIprefix(color::RED, "failed", align::CENTER);

  for (int i = 0; i < unit_test.total_test_case_count(); i++){
    const ::testing::TestCase& test_case = *unit_test.GetTestCase(i);
    if (!test_case.should_run() || (test_case.failed_test_count() == 0))
      continue;
    for (int j = 0; j < test_case.total_test_count(); j++){
      const ::testing::TestInfo& test_info = *test_case.GetTestInfo(j);
      if (!test_info.should_run() || test_info.result()->Passed())
        continue;
      result += std::get<0>(prefix);
      result += test_case.name();
      result += ".";
      result += test_info.name();
      result += this->PrintFullTestCommentIfPresent(test_info);
      result += "\n";
    }
  }

  return result;
}

void Nucleus::OnTestIterationEnd(const ::testing::UnitTest& unit_test, int iteration){
  auto prefix = this->getMPIprefix(color::GREEN, "=", align::FULL);

  // Print Passed tests
  std::string message = std::get<0>(prefix);
  message += formatTestCount(unit_test.test_to_run_count());
  message += " from ";
  message += formatTestCaseCount(unit_test.test_case_to_run_count());
  message += " ran.";
  message += " (" + std::to_string(unit_test.elapsed_time()) + " ms total)\n";
  prefix = this->getMPIprefix(color::GREEN, "passed", align::CENTER);
  message += std::get<0>(prefix);
  message += formatTestCount(unit_test.successful_test_count());
  message += ".\n";

  // Print Failed tests
  if(!unit_test.Passed()){
    prefix = this->getMPIprefix(color::RED, "failed", align::CENTER);
    message += std::get<0>(prefix);
    int test_count = unit_test.failed_test_count();
    message += formatTestCount(test_count);
    message += ", listed below:\n";
    message += this->PrintFailedTests(unit_test);
    message += "\n  " + std::to_string(test_count) + " FAILED " + (test_count == 1 ? "TEST" : "TESTS") + "\n";
  }
  

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));
}

void Nucleus::OnEnvironmentsSetUpStart(const ::testing::UnitTest& unit_test){
  auto prefix = this->getMPIprefix(color::GREEN, "-", align::FULL);
  std::string message = std::get<0>(prefix);
  message += "Global test environment set-up.\n";

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));
}

void Nucleus::OnEnvironmentsTearDownStart(const ::testing::UnitTest& unit_test){
  auto prefix = this->getMPIprefix(color::GREEN, "-", align::FULL);
  std::string message = std::get<0>(prefix);
  message += "Global test environment tear-down\n";

  this->OutputMessage(message, std::get<1>(prefix), std::get<2>(prefix));
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
