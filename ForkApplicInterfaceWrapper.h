#ifndef _FORK_APPLIC_INTERFACE_WRAPPER_H_
#define _FORK_APPLIC_INTERFACE_WRAPPER_H_

// Dakota includes
#include<dakota_data_types.hpp>
#include<ProblemDescDB.hpp>
#include<ForkApplicInterface.hpp>

class ForkApplicInterfaceWrapper : public Dakota::ForkApplicInterface {

public:

  // constructor
  ForkApplicInterfaceWrapper(const Dakota::ProblemDescDB& problem_db);
  // destructor
  ~ForkApplicInterfaceWrapper() override;

};

inline ForkApplicInterfaceWrapper::~ForkApplicInterfaceWrapper()
{

}

#endif
