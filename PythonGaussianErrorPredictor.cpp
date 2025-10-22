/*  _______________________________________________________________________

    Dakota: Explore and predict with confidence.
    Copyright 2014-2024
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/************************************************************************
 * MODIFICATIONS IN THIS DERIVED FILE
 *
 * Derived from:
 *   - Dakota's dakota::surrogates::Python
 *   - Upstream commit: 0daaaa2237bd79445349c47d157a4a2b73db7452
 *
 * Modifications by Aditya Narkhede (AN), 2025:
 *   - Changed name from PythonGaussianErrorPredictor to PythonGaussianErrorPredictor
 *   - Removed dakota::surrogates namespace
 *   - Modified initialize_python()
 *   - Added a "variance" function to call the corresponding
 *     "predict_variance" from the python module.
 *
 * Additional Notes:
 *   - Keeping the original comments.
 *
 * Licensing:
 *   This file continues to be licensed under the GNU Lesser General Public
 *   License v2.1 or (at your option) any later version. See the bundled
 *   license file `LICENSE.LGPL-2.1` or
 *   <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 ************************************************************************/

#include <PythonGaussianErrorPredictor.h>
#include <surrogates_tools.hpp>

using namespace Eigen;
using namespace dakota::util;
using namespace dakota::surrogates;

PythonGaussianErrorPredictor::
PythonGaussianErrorPredictor(const std::string& module_and_class_name) :
  moduleAndClassName(module_and_class_name),
  ownPython(false),
  pyModuleActive(false)
{
  initialize_python();
}


PythonGaussianErrorPredictor::
PythonGaussianErrorPredictor(const MatrixXd& samples,
                             const MatrixXd& response,
                             const std::string& module_and_class_name) :
  moduleAndClassName(module_and_class_name),
  ownPython(false),
  pyModuleActive(false)
{
  initialize_python();
  build(samples, response);
}


void
PythonGaussianErrorPredictor::initialize_python()
{
  // Consider adding meaningful parameters... RWH
  configOptions.set("verbosity", 1, "console output verbosity");

  ownPython = false;
  pyModuleActive = false;
  if (!Py_IsInitialized()) {
    py::initialize_interpreter();
    ownPython = true;
    if (!Py_IsInitialized())
      throw(std::runtime_error(
            "Error: Could not initialize PythonGaussianErrorPredictor for surrogates use."));
  }
  if (!pyModuleActive) {
    try {
      size_t p = moduleAndClassName.find_last_of(".");
      if( std::string::npos == p )
        throw(std::runtime_error(
              "Invalid surrogate python module_and_class_name.\n\tUse \"module.classname\""));
      auto module_name = moduleAndClassName.substr(0, p);
      auto class_name  = moduleAndClassName.substr(p+1);
      py::object pyModule = py::module_::import(module_name.c_str());
      // AN: This calls the __init__ function of the class.
      pySurrogate = pyModule.attr(class_name.c_str())();
    }
    catch(py::error_already_set &e) {
      if (e.matches(PyExc_ModuleNotFoundError)) {
        std::cerr << "Could not load the required module '"
                  << moduleAndClassName << "'" << std::endl;
        throw;
      }
      else {
        std::cerr << "Caught a python exception:\n"
                  << e.what() << std::endl;
      }
    }
    pyModuleActive = true;
  }

  // Check that needed methods exist in the module
  // ... We could allow the user to register these... RWH
  // predict_variance option added by AN
  std::vector<std::string> req_attrs = { "construct", "predict", "predict_variance" };
  bool is_module_valid = true;
  for( auto const & req_at : req_attrs ) {
    try {
      py::object py_fn = pySurrogate.attr(req_at.c_str());
    }
    catch(py::error_already_set &e) {
      if (e.matches(PyExc_AttributeError)) {
        std::cerr << "Module '" << moduleAndClassName << "' does not "
          << "contain required method '" << req_at << "'"
          << std::endl;
      }
      is_module_valid = false;;
    }
  }
  if( !is_module_valid )
    throw(std::runtime_error("Invalid python module for surrogates"));
}

void
PythonGaussianErrorPredictor::build(const MatrixXd& samples,
                                    const MatrixXd& response)
{
  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  verbosity = configOptions.get<int>("verbosity");

  if (verbosity > 0) {
    if (verbosity == 1) {
      std::cout << "\nBuilding PythonGaussianErrorPredictor surrogate\n\n";
    } else if (verbosity == 2) {
      std::cout << "\nBuilding PythonGaussianErrorPredictor surrogate with module.method\n"
                << moduleAndClassName << "." << "construct" << "\n";
    } else
      throw(
          std::runtime_error("Invalid verbosity int for PythonGaussianErrorPredictor surrogate"));
  }
  // Hard-coded method for now; could expose to user - RWH
  const std::string fn_name("construct");
  py::object py_surr_builder = pySurrogate.attr(fn_name.c_str());
  py_surr_builder(samples, response);

  isField = (response.cols() > 1);
}

bool
PythonGaussianErrorPredictor::diagnostics_available()
{
  return !isField;
}


VectorXd
PythonGaussianErrorPredictor::value(const MatrixXd& eval_points)
{

  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  // Hard-coded method for now; could expose to user - RWH
  const std::string fn_name("predict");
  py::object py_surr_eval = pySurrogate.attr(fn_name.c_str());

  auto vals = py_surr_eval(eval_points).cast<VectorXd>();

  return vals;//.col(0);
  //return py_surr_eval(eval_points).cast<VectorXd>();
}


VectorXd
PythonGaussianErrorPredictor::values(const MatrixXd& eval_points)
{

  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  // Hard-coded method for now; could expose to user - RWH
  const std::string fn_name("predict");
  py::object py_surr_eval = pySurrogate.attr(fn_name.c_str());

  auto vals = py_surr_eval(eval_points).cast<MatrixXd>();

  return vals.row(0);
  //return py_surr_eval(eval_points).cast<VectorXd>();
}


MatrixXd
PythonGaussianErrorPredictor::gradient(const MatrixXd& eval_points)
{
  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  // Hard-coded method for now; could expose to user - RWH
  // We could add a check for this method (attribute) above in the
  // req_attrs if we knew it was needed at the time of our construction.
  const std::string fn_name("gradient");
  py::object py_surr_grad;
  try {
    py_surr_grad = pySurrogate.attr(fn_name.c_str());
  }
  catch(py::error_already_set &e) {
    if (e.matches(PyExc_AttributeError)) {
      std::cerr << "Module '" << moduleAndClassName << "' does not "
        << "contain required method '" << fn_name << "'"
        << std::endl;
      throw;
    }
  }

  return py_surr_grad(eval_points).cast<MatrixXd>();
}


MatrixXd
PythonGaussianErrorPredictor::hessian(const MatrixXd& eval_point)
{

  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  // Hard-coded method for now; could expose to user - RWH
  // We could add a check for this method (attribute) above in the
  // req_attrs if we knew it was needed at the time of our construction.
  const std::string fn_name("hessian");
  py::object py_surr_hess;
  try {
    py_surr_hess = pySurrogate.attr(fn_name.c_str());
  }
  catch(py::error_already_set &e) {
    if (e.matches(PyExc_AttributeError)) {
      std::cerr << "Module '" << moduleAndClassName << "' does not "
        << "contain required method '" << fn_name << "'"
        << std::endl;
      throw;
    }
  }

  return py_surr_hess(eval_point).cast<MatrixXd>();
}

VectorXd
PythonGaussianErrorPredictor::variance(const MatrixXd& eval_points)
{

  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  const std::string fn_name("predict_variance");
  py::object py_surr_variance;
  try {
    py_surr_variance = pySurrogate.attr(fn_name.c_str());
  }
  catch(py::error_already_set &e) {
    if (e.matches(PyExc_AttributeError)) {
      std::cerr << "Module '" << moduleAndClassName << "' does not "
        << "contain required method '" << fn_name << "'"
        << std::endl;
      throw;
    }
  }

  auto vals = py_surr_variance(eval_points).cast<VectorXd>();

  return vals;

}

double
PythonGaussianErrorPredictor::loss(const MatrixXd& test_samples,
                                   const MatrixXd& test_response)
{

  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  const std::string fn_name("loss");
  py::object py_surr_loss;
  try {
    py_surr_loss = pySurrogate.attr(fn_name.c_str());
  }
  catch(py::error_already_set &e) {
    if (e.matches(PyExc_AttributeError)) {
      std::cerr << "Module '" << moduleAndClassName << "' does not "
        << "contain required method '" << fn_name << "'"
        << std::endl;
      throw;
    }
  }

  double val = py_surr_loss(test_samples, test_response).cast<double>();

  return val;

}

void
PythonGaussianErrorPredictor::save_model()
{
  assert( pyModuleActive );
  assert( Py_IsInitialized() );

  const std::string fn_name("save");
  py::object py_surr_save;
  try {
    py_surr_save = pySurrogate.attr(fn_name.c_str());
  }
  catch(py::error_already_set &e) {
    if (e.matches(PyExc_AttributeError)) {
      std::cerr << "Module '" << moduleAndClassName << "' does not "
        << "contain required method '" << fn_name << "'"
        << std::endl;
      throw;
    }
  }

  py_surr_save();
}

